/*
    Copyright (C) 2009,2011 Stefan Enroth

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "utils.h"

void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters)
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}


/*
 * returns the 'next' position to start from.
 *
 */ 
string::size_type readAndTokenize(const string& str,
				  vector<string>& tokens,
				  const string& delimiters,
				  string::size_type firstPos)
{
  
  // Skip delimiters at beginning.
  firstPos = str.find_first_not_of(delimiters + "\n",firstPos);
  // Set the end.
  string::size_type lastPos  = str.find_first_of("\n", firstPos);
  if(lastPos == string::npos)
    return(str.length() + 1);
  string::size_type pos     = str.find_first_of(delimiters + "\n", firstPos+1);
  if(pos == string::npos)
    return(str.length() + 1);
  //cerr<<"RAT: "<<firstPos<<" "<<lastPos<<" "<<pos<<endl;
  while (pos < lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(firstPos, pos - firstPos));
      //cerr<<str.substr(firstPos, pos - firstPos)<<" ";
      // Skip delimiters.  Note the "not_of"
      firstPos = str.find_first_not_of(delimiters + "\n", pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters + "\n", firstPos);
    } 
  // Push the last one.  
  tokens.push_back(str.substr(firstPos, pos - firstPos));
  //cerr<<str.substr(firstPos, pos - firstPos)<<endl;
  //cout<<"RAT: "<<tokens[1]<<" "<<tokens[2]<<" "<<tokens[3]<<" "<<tokens[4]<<endl;
  return(lastPos);
}

/************************************************************
 * 
 * Random stuff. this should be uniform in the given span...
 * 
 *
 ************************************************************/

int uniform_range(int low,int high)
{
  int range = high - low + 1;
  int ret = low + (int)((double)range * rand()/(RAND_MAX + 1.0));
  return (ret);
}

int parseBEDline(string line,inputLine * res,int seqCol,int startCol,
				 int endCol,int strandCol, int nameCol)
{
  vector<string> data;
  float tmpVal = 1.0,tmp; 
  data.clear();
  Tokenize (line,data," \t"); // need to split on both " " and "\t" to be able to parse headers.
  if(data.size() == 0)
    return(data.size());
  res->header = 0;    
  if(strcmp((data[0]).c_str(),"track") == 0) // header line, should start with "track name=..."
    {
      res->header = 1;
      return(data.size());
    }
  res->mapped = 1;     // only mapped fragments
  res->seq    = data[seqCol];
  res->strand = (data[strandCol]== "+" ? 1:-1);
  res->pos    = atoi(data[startCol].c_str()) +1; // BED format is 0-based. change to 1-based.
  res->len    = atoi(data[endCol].c_str()) - atoi(data[startCol].c_str()); // crds are [x,y) so no need for +1 in length.
  res->name   = data[nameCol];
  res->val = (storageType)tmpVal;
  return(data.size());
}


/*
 * returns the number of lines with mapped reads 
 */
int initControlBEDlite(ifstream * inf,map<string,seqStats> * seqmap,int seqCol,int startCol,
					   int endCol,int strandCol,bool verbose)
{
  int nlines = 0;
  int *tmpVals;
  seqStats tmpStat;
  string line;
  inputLine tmp;
  pair<map<string,seqStats>::iterator,bool> ret;
  pair<map<int,int*>::iterator,bool> retVal;
  
  // initialize tmpStat
  tmpStat.truncR = 0;
  tmpStat.truncF = 0;
  tmpStat.truncC = 0;
   
  while(!(inf->eof()))
    { 
	  //cerr<<nlines<<endl;
      if((nlines % 1000) == 0) 
		cerr<<"-\r";
      if((nlines % 2000) == 0) 
		cerr<<"/\r";
      if((nlines % 3000) == 0) 
		cerr<<"|\r";
      if((nlines % 4000) == 0) 
		cerr<<"\\\r";
      if((nlines % 5000) == 0) 
		cerr<<"-\r";
      
      getline(*inf,line); 
	  // skip possible empty lines 
      if(parseBEDline(line,&tmp,seqCol,startCol,endCol,strandCol) < 1) 
		continue;
	  if(tmp.header)
		continue;
	  
      nlines++;
      // try and insert into seqmap. one entry per choromosome.
      tmpStat.minPos = max(tmp.pos,0); // allow 0-based coordinates. 
      tmpStat.maxPos = tmp.pos + tmp.len-1; 
      tmpStat.countF = 0;
      tmpStat.countR = 0;

	  if(tmp.strand == 1) // forward strand
		{
		  tmpStat.countF++;
		}
	  else
		{
		  tmpStat.countR++;
		}
	  //cerr<<"trying to add data"<<endl;
	  ret = seqmap->insert(pair<string,seqStats> (tmp.seq,tmpStat));
      
      if(!ret.second) // already a key with this value.
		{
		  if(tmpStat.minPos < (*seqmap)[tmp.seq].minPos)
			(*seqmap)[tmp.seq].minPos = tmpStat.minPos;
		  
		  if(tmpStat.maxPos > (*seqmap)[tmp.seq].maxPos)
			(*seqmap)[tmp.seq].maxPos = tmpStat.maxPos;
		  
		  if(tmp.strand == 1) // forward strand
			{
			  //cerr<<"trying to add data on forward to existing"<<endl;
			  (*seqmap)[tmp.seq].countF++;
			}
		  else
			{
			  //cerr<<"trying to add data on reverse to existing"<<endl;
			  (*seqmap)[tmp.seq].countR++;
			}
		}else{
		//cerr<<"new sequence added."<<endl;
	  }
    }
  // reset the file-stream  
  inf->clear(); 
  inf->seekg(0,ios::beg);
  return(nlines);
}

