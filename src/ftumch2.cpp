/*
    Copyright (C) 2011 Stefan Enroth

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
#include <gsl/gsl_multifit.h>

/************************************************************
 *
 * main 
 *
 ************************************************************/

int main(int argc, char* argv[]) 
{
  // parameters that you can set. 
  string delim      = "\t ";
  string chipFile   = "";
  vector<string> ctrlFiles;
  string outFile    = "";
  int readLen       = 50;
  int chunkSize     = 100000;
  int windowSize    = 5;
  int interval      = 5;
  bool talk         = false;

  string errorLine = "usage " + 
    string(argv[0]) + 
    " [Parameters]\n" +  
    "\t-i  <infile, BED-formated file containing the ChIP-reads, sorted on chromosome and position.>\n" +
	"\t-c  <space\\tab separataed list of infile(s), BED-formated file(s) \n" + 
    "\t    containing the control-reads (e.g. Input/IgG et cetera), sorted as the file given in '-i' \n" +  
    "\t-o  <outfile, BED-formated file of resulting reads after normalization, \n" + 
	"\t    with read lengths as defined by -l>\n" + 
    "\t-rl <read length, defaults to 50 >\n" + 
    "\t-cs <chunk size, number of bp considered at a time when building the model>\n" + 
	"\t-ws <window size, at every point used to build the model a window of +/- \n" + 
	"\t    this size is averaged to create an observed data point.>\n" +
    "\t-iv <interval, the step size determining the distance between points \n" + 
	"\t     used as observations in the regression model.>\n" +
    "\t-v  <set verbose>\n" + 
 	"example: \n" +  string(argv[0]) + " -i myreads.bed -c input.bed igg.bed noise.bed -o normalized.bed -rl 50 -cs 100000 -ws 5 -iv 5 \n"
   ;


  bool fail = false;
  bool ctrlfiles = false;
  string failmessage = "";
  
  for (int i=1;i<argc;i++)
    {
      if(strcmp(argv[i],"-i") == 0)
		{
		  chipFile.assign(argv[++i]);
		  ctrlfiles = false;
		}
	  else if(strcmp(argv[i],"-o") == 0)
		{
		  outFile.assign(argv[++i]);
    	  ctrlfiles = false;
		}
	  else if(strcmp(argv[i],"-c") == 0)
		{
		  ctrlfiles = true;
		}
      else if(strcmp(argv[i],"-rl") == 0)
		{
		  readLen = atoi(argv[++i]);
		  ctrlfiles = false;
		}
	  else if(strcmp(argv[i],"-cs") == 0)
		{
		  chunkSize = atoi(argv[++i]);
    	  ctrlfiles = false;
		}
	  else if(strcmp(argv[i],"-ws") == 0)
		{
		  windowSize = atoi(argv[++i]);
    	  ctrlfiles = false;
		}
	  else if(strcmp(argv[i],"-iv") == 0)
		{
		  interval = atoi(argv[++i]);
		  ctrlfiles = false;
		}
	  else if(strcmp(argv[i],"-v") == 0)
		{
		  talk = true;
		  ctrlfiles = false;
		}
	  else
		{
    
		  if(ctrlfiles) // assume that all things not parsable after -c are control files. Check for existance/readability below.
			{
			  ctrlFiles.push_back(argv[i]);
			}
		  else
			{
			  failmessage.assign("Unknown argument: ");
			  failmessage.append(argv[i]);
			  failmessage.append("\n");
			  fail = true;
			}
		}
    }
  
  // Check infile and readability. 

  if(chipFile == "")
    {
      failmessage.append("infile (-i) must be specified.\n");
      fail = true;
    }
  
  ifstream inf;
  inf.open(chipFile.c_str());
  
  if(!inf)
    {
      failmessage.append("Could not open infile '" + chipFile + "' (does the file exist?)\n");
      fail = true;
    }
  
  // Check control files. 
  if(ctrlFiles.size() < 1)
    {
      failmessage.append("at least one control file (-c) must be specified.\n");
      fail = true;
    }

  
  ifstream infc[ctrlFiles.size()];
  if(!fail)
	for (int i = 0;i<ctrlFiles.size();i++)
	  {
		infc[i].open(ctrlFiles[i].c_str());
		if(!infc[i])
		  {
			failmessage.append("Could not open ctrlfile '" + ctrlFiles[i]  + "' (does the file exist?)\n");
			fail = true;
		  }
	  }
	  
  // Check outfile and readability. 

  if(outFile == "")
    {
      failmessage.append("outfile (-o) must be specified.\n");
      fail = true;
    }
  ofstream outf;  
  if(!fail)
	outf.open(outFile.c_str(),ios::trunc);
  if(!outf)
    {
      failmessage.append("Could not open outfile '" + outFile + "' (do we have permission ?)\n");
      fail = true;
    }

  // are we ok so far? 
  if (fail)
    {
      cerr << endl << failmessage.c_str() << endl << errorLine << endl;
	  //try and close any opened files
	  inf.close();
	  for (int i = 0;i<ctrlFiles.size();i++)  
		infc[i].close();
	  outf.close(); 
	  return(1);
    }
  
  /*
   * Get some initial parameters. 
   */
  map <string,seqStats> seqMapChip; 
  map <string,seqStats> *seqMapCtrls;
  seqMapCtrls = new map<string,seqStats>[ctrlFiles.size()];
  
  map <string,seqStats>::iterator it;
  map <int,int*>::iterator valIt;
  
  // Read the reference sequences and the range of each file
  cout<<"Reading BED-files."<<endl;
  int nlinesChIP, nlinesCtrl=0;
  cout<<"ChIP file.."<<endl;
  nlinesChIP = initControlBEDlite(&inf,&seqMapChip,0,1,2,5,true);
  cout<<"Control file(s) .."<<endl;
  for (int i = 0;i<ctrlFiles.size();i++)  
	nlinesCtrl += initControlBEDlite(&infc[i],&seqMapCtrls[i],0,1,2,5,true);
  
  cout<<"ChIP-data consists of "<<nlinesChIP<<" mapped fragments."<<endl;
  cout<<"Control-data consists of "<<nlinesCtrl<<" mapped fragments."<<endl;
  
  // print some stats. 
  cout <<"ChIP Read Statistics::"<<endl;
  cout <<setw(10)<<"Name\t"<<setw(10)<<"minCrd\t"<<setw(10)<<"maxCrd\t"<<setw(10)<<"F_counts\t"<<setw(10)<<"R_counts\t"<<endl;
  for ( it=seqMapChip.begin() ; it != seqMapChip.end(); it++ )
    {
      cout <<setw(10)<< (*it).first << "\t" <<setw(10)<< (*it).second.minPos << "\t" << setw(10)<<(*it).second.maxPos<<"\t";
      cout <<setw(10)<< (*it).second.countF << "\t" <<setw(10)<< (*it).second.countR<<endl;
    }
  cout <<"Control Statistics::"<<endl;
  for (int i = 0;i<ctrlFiles.size();i++)
	{
	  cout<<ctrlFiles[i]<<endl;
	  cout <<setw(10)<<"Name\t"<<setw(10)<<"minCrd\t"<<setw(10)<<"maxCrd\t"<<setw(10)<<"F_counts\t"<<setw(10)<<"R_counts\t"<<endl;
	  for ( it=seqMapCtrls[i].begin() ; it != seqMapCtrls[i].end(); it++ )
		{
		  cout <<setw(10)<< (*it).first << "\t" <<setw(10)<< (*it).second.minPos << "\t" << setw(10)<<(*it).second.maxPos<<"\t";
		  cout <<setw(10)<< (*it).second.countF << "\t" <<setw(10)<< (*it).second.countR<<endl;
		}
	}
  
  cout<<"Processing reads in chunks of "<<chunkSize<<" bp."<<endl;
  int lowPos,highPos;
  int chunkRange;
  int chunkObs;
  int obsCount;
  int *winSumF = new int[ctrlFiles.size()+1];
  int *winSumR = new int[ctrlFiles.size()+1];
  double chisqF,chisqR;
  gsl_matrix *XF, *covF,*XR, *covR;
  gsl_vector *yF,*cF,*rF,*yR,*cR,*rR;
  double *resiF,*resiR;
  int *cntF,*cntR;
  int posOff;
  

  // initialize. 
  string line;
  inputLine chipLine;
  inputLine *ctrlLines;
  ctrlLines = new inputLine[ctrlFiles.size()];
  
  // read first line from each file. check position and chr. assume that the files are ordered inside chr. i.e don't loop
  // over chr by seqMap but over info in the files. retrieve min/max position from the seqMap depending on file contents. 
  // also assume that the chr ordering is the same in chip & control files.
  getline(inf,line);
  parseBEDline(line,&chipLine,0,1,2,5);
  for(int i=0;i<ctrlFiles.size();i++)
	{
	  getline(infc[i],line);
	  parseBEDline(line,&ctrlLines[i],0,1,2,5);
	}

  // initialize with the "first" chromosome and its min/max pos.
  string currChr = chipLine.seq;
  int chrMinPos,chrMaxPos;
  int chrMinPosCtrl,chrMaxPosCtrl;
  int currLine = 1;
  int memNeeded;
  int chrPosChip,chrPosCtrl;
  int ctrlIndex;

  // tmp. storage for the chip/control signals.
  unsigned short *chipF,*chipR,*ctrlF,*ctrlR;

  // introduce curr pos, curr Chr etc. and a loop on !EOF in the chip file. 
  // no point in normalizing where there are no signals in chip...

  while(currLine <= nlinesChIP)
    {
	  chrMinPos = seqMapChip.find(currChr)->second.minPos;
	  chrMaxPos = seqMapChip.find(currChr)->second.maxPos;;
	  if(talk)
		cout<<"ChIP: "<<chipLine.seq<<" "<<chrMinPos<<" "<<chrMaxPos<<endl;
	  chrPosChip = chrMaxPos - chrMinPos +1;
	  // check the min/max for this chr in ctrl-data
	  chrMinPosCtrl = INT_MAX;
	  chrMaxPosCtrl = -1;
	  for(int i = 0;i<ctrlFiles.size();i++)
		{
		  if(seqMapCtrls[i].count(currChr))
			{
			  chrMinPosCtrl = min(chrMinPosCtrl,seqMapCtrls[i][currChr].minPos);
			  chrMaxPosCtrl = max(chrMaxPosCtrl,seqMapCtrls[i][currChr].maxPos);
			} 
		}
	  if(talk)
		cout<<"Control: "<<chipLine.seq<<" "<<chrMinPosCtrl<<" "<<chrMaxPosCtrl<<endl;
	  chrPosCtrl = chrMaxPosCtrl - chrMinPosCtrl +1;
	  memNeeded = sizeof(unsigned short)*(chrPosChip + ctrlFiles.size()*chrPosCtrl);
	  // allocate memory to hold the entire chromosome, do the regression in chunks. 
	  try{
		cout<<"Trying to allocate: ";	
		if(memNeeded > 1000000000)
		  cout<<memNeeded/1000000000<<" Gb for "<<currChr<<".";
		else if (memNeeded > 1000000)
		  cout<<memNeeded/1000000<<" Mb for "<<currChr<<".";
		else if (memNeeded > 1000)
		  cout<<memNeeded/1000<<" kb for "<<currChr<<".";
		else
		  cout<<memNeeded<<" bytes for raw signals"<<currChr<<".";
		
		chipF = new unsigned short[chrPosChip];
		chipR = new unsigned short[chrPosChip];
		ctrlF = new unsigned short[chrPosCtrl*ctrlFiles.size()]; // these will need to be accessed in a "[i + chrPosCtrl*j]"-type of fashion.
		ctrlR = new unsigned short[chrPosCtrl*ctrlFiles.size()];
		cout<<" Done."<<endl;
	  }catch (std::bad_alloc  &f){
		cerr<<string(argv[0])<<" couldn't allocate as much memory as it wanted. Failure: '"<<f.what()<<endl;
		// close files. 
		inf.close();
		for (int i = 0;i<ctrlFiles.size();i++)  
		  infc[i].close();
		outf.close();	
		
		delete[] chipF;
		delete[] chipR;
		delete[] ctrlF;
		delete[] ctrlR;
		delete[] resiF;
		delete[] resiR;
		delete[] cntF;
		delete[] cntR;
		return(-1);
	  }

	  // make sure it's all zeroes.
	  for(int i=0;i<chrPosChip;i++)
		{
		  chipF[i] = 0;
		  chipR[i] = 0;
		}
	  for(int i=0;i<chrPosCtrl;i++)
		for(int j = 0;j<ctrlFiles.size();j++)
		  {
			ctrlF[i + j*chrPosCtrl] = 0;
			ctrlR[i + j*chrPosCtrl] = 0;
		  }
			
	  
	  // read in the sought chip-data
	  while((chipLine.seq == currChr)  && !(inf.eof())) // chip-file
		{
		  //cout<<chipLine.seq<<"\t"<<line<<endl;
		  // update previous line's data.
		  if(chipLine.strand == 1)
			chipF[chipLine.pos-chrMinPos]++;
		  else
			chipR[chipLine.pos-chrMinPos+chipLine.len]++;
		  // read in the nextline.
		  getline(inf,line);
		  parseBEDline(line,&chipLine,0,1,2,5);
		  currLine++;
		}
	   if((chipLine.seq == currChr)  && (inf.eof())) // chip-file, last read, ok chr, use.
		{
		  //cout<<"Last line of the ChipFIle"<<endl;
		  //cout<<chipLine.seq<<"\t"<<line<<endl;
		  if(chipLine.strand == 1)
			chipF[chipLine.pos-chrMinPos]++;
		  else
			chipR[chipLine.pos-chrMinPos+chipLine.len-1]++;
		  
	}
	  
	  // read in the sought ctrl-data
	  for(int i = 0;i<ctrlFiles.size();i++)
		{
		  // is there data at all for this chr in this control file?
		  if(seqMapCtrls[i].count(currChr) == 1)
			{
			  // we're assuming that the chromosomes are in the same order in the chip & ctrl files. 
			  // cases:
			  // chr on current line is not the same as in chip
              //    => we know that we should have chr data on this chr & that chrs comes in the same order. this can prob. only 
			  //       happen for a chr-specific chromosome, e.g. its safe to read past and check again. 
			  // chr on current line is the same as in chip
			  //   => this is good. last time (either preFirst or not) should have read prev. chr completely. so just start reading until we hit 
			  //      another chr. 
			  // 
			 
			  if(ctrlLines[i].seq != currChr)
				{
				  // read past th "wrong" chromosome(s). 
				  while(!infc[i].eof() && ctrlLines[i].seq != currChr)
					{
					  getline(infc[i],line);
					  parseBEDline(line,&ctrlLines[i],0,1,2,5);
					}
				}
			  
			  // now we have the first line of the the correct chr in 'ctrlLines[i]'
			  // Read in the complete chr and store the data accordingly.
			  while(!infc[i].eof() && ctrlLines[i].seq == currChr)
				{
				  if(ctrlLines[i].strand == 1)
					ctrlF[ctrlLines[i].pos-chrMinPosCtrl + i*chrPosCtrl]++;
				  else
					ctrlR[ctrlLines[i].pos-chrMinPosCtrl+ctrlLines[i].len + i*chrPosCtrl-1]++;
				  getline(infc[i],line);
				  parseBEDline(line,&ctrlLines[i],0,1,2,5);
				}
			}
		  
		}
	  cout<<"Analysing "<<currChr<<endl;
	  currChr = chipLine.seq; // store "next" chromosome
	  // now all data for this chr is read. Start analysing in chunks. 
	  lowPos = chrMinPos;
	  while(lowPos < chrMaxPos) // loop over this chromosome data in chunks.
		{ 
		  if(!talk)
			{
			  cout<<lowPos<<" of "<<chrMaxPos<<"\r";
			}
		  highPos = lowPos + chunkSize-1;
		  if(highPos >= (chrMaxPos - 0.5*chunkSize)) // less than 0.8 of a chunk left. merge.
			highPos =  chrMaxPos;
		  chunkRange = highPos - lowPos + 1;
		  if(talk)
			cout<<"["<<lowPos<<","<<highPos<<"]\tsize: "<<chunkRange<<endl;
		  
		  resiF = new double[chunkRange];
		  resiR = new double[chunkRange];
		  cntF = new int[chunkRange];
		  cntR = new int[chunkRange];
		  for(int i=0;i<chunkRange;i++)
			{
			  resiF[i] = 0.0;
			  resiR[i] = 0.0;
			  cntF[i] = 0;
			  cntR[i] = 0;
			}
		  
		  // for each chunk, step forward in 'interval' steps and average signals in that window.
		  chunkObs = (chunkRange-2*windowSize)/interval + 1;
		  if (talk) 
			cout<<"\tsampling this chunk at "<<chunkObs<<" positions."<<endl;  
		  // Storage for the signals on '+'
		  XF = gsl_matrix_alloc (chunkObs, ctrlFiles.size());
		  yF = gsl_vector_alloc (chunkObs);
		  rF = gsl_vector_alloc (chunkObs);
		  cF = gsl_vector_alloc (ctrlFiles.size());
		  covF = gsl_matrix_alloc (ctrlFiles.size(), ctrlFiles.size());
		  // Storage for the signals on '-'
		  XR = gsl_matrix_alloc (chunkObs, ctrlFiles.size());
		  yR = gsl_vector_alloc (chunkObs);
		  rR = gsl_vector_alloc (chunkObs);
		  cR = gsl_vector_alloc (ctrlFiles.size());
		  covR = gsl_matrix_alloc (ctrlFiles.size(), ctrlFiles.size());
		  
		  // loop over the signals in interval steps and average in +/- windowSize. fill in the matrices. 
		  obsCount = 0;
		  for (int i = lowPos+windowSize;i<(highPos-windowSize);)
			{
			  // collect sums over each signal in the sough window
			  for (int j = 0;j < ctrlFiles.size() + 1;j++)
				{
				  winSumF[j] = 0;
				  winSumR[j] = 0;
				}
			  for (int j = -windowSize;j<=windowSize;j++)
				{				
				  winSumF[0] += chipF[i - chrMinPos + j];
				  winSumR[0] += chipR[i - chrMinPos + j];
				  for(int k = 0;k<ctrlFiles.size();k++)
					{
					  ctrlIndex = i - chrMinPosCtrl + j + k*chrPosCtrl;
					  if(ctrlIndex >= 0 && ctrlIndex <chrPosCtrl) // is there ctrl data for this position?
						{
						  winSumF[1+k] += ctrlF[ctrlIndex];
						  winSumR[1+k] += ctrlR[ctrlIndex];
						}
					}
				}
			  // the chip signal
			  gsl_vector_set (yF, obsCount, (double)winSumF[0]/(double)(2*windowSize+1));
			  gsl_vector_set (yR, obsCount, (double)winSumR[0]/(double)(2*windowSize+1));
			  
			  // the control signals
			  for (int j = 0;j < ctrlFiles.size();j++)
				{
				  gsl_matrix_set (XF, obsCount, j, (double)winSumF[j+1]/(double)(2*windowSize+1));
				  gsl_matrix_set (XR, obsCount, j, (double)winSumR[j+1]/(double)(2*windowSize+1));
				}
			  obsCount++;
			  i+=interval; 
			}
		  // fit the models.
		  gsl_multifit_linear_workspace * work  = gsl_multifit_linear_alloc (chunkObs, ctrlFiles.size());
		  /* 
		   * '+' Strand
		   */
		  gsl_multifit_linear (XF, yF, cF, covF,&chisqF, work);
		  if(talk)
			{
			  cout<<"\t'+' chisq: "<<chisqF<<"\t"<<"c's:";
			  for (int j = 0;j < ctrlFiles.size();j++)
				cout<<gsl_vector_get(cF,j)<<" ";
			  cout<<endl;
			}
		  /* 
		   * '-' Strand
		   */
		  gsl_multifit_linear (XR, yR, cR, covR,&chisqR, work);
		  if(talk)
			{
			  cout<<"\t'-' chisq: "<<chisqR<<"\t"<<"c's:";
			  for (int j = 0;j < ctrlFiles.size();j++)
				cout<<gsl_vector_get(cR,j)<<" ";
			  cout<<endl;
			}
		  
		  gsl_multifit_linear_free (work);
		  // calculate residuals.
		  if(talk)
			cout<<"\tCaclulating residuals..";
		  gsl_multifit_linear_residuals (XF,yF,cF,rF);
		  gsl_multifit_linear_residuals (XR,yR,cR,rR);
		  if(talk)
			cout<<"done."<<endl<<"\tRebuilding signal..";
		  // rebuild a per-bp-signal 
		  for (int i=0;i<chunkObs;i++)
			{
			  // center of this observation. 
			  posOff = 1+(i+1)*interval;
			  if(posOff > chunkRange) // outside of our chunk, should never happen.
				continue;
			  if(gsl_vector_get(rF,i) > 0.5)  // original R-code used 'round' on the residuals, ceiling(x-0.5) does the same thing
				for (int j = -windowSize;j<=windowSize;j++)
				  {
					resiF[j+posOff] += ceil(gsl_vector_get(rF,i)-0.5);
					cntF[j+posOff] += 1;
				  }
			  if(gsl_vector_get(rR,i) > 0.5)
				for (int j = -windowSize;j<=windowSize;j++)
				  {
					resiR[j+posOff] += ceil(gsl_vector_get(rR,i)-0.5);
					cntR[j+posOff] += 1;
				  }
			}
		  
		  if(talk)
			cout<<"done."<<endl<<"\tWriting output..";
		  for (int i=0;i<chunkRange;i++)
			{
			  if(cntF[i] > 0)
				{
				  resiF[i] = resiF[i]/(double)cntF[i];
				  if(resiF[i] > 0) 
					{
					  for(int j=0;j<ceil(resiF[i]);j++)
						{
						  outf<<currChr<<"\t"<<lowPos + i-1<<"\t"<<lowPos+i+readLen-2<<"\tDUMMY\t"; // bed is zero-based, halfopen (ie.-1/-2)
						  outf<<resiF[i]<<"\t+"<<endl;
						}
					}
				}
			  
			  if(cntR[i] > 0)
				{
				  resiR[i] = resiR[i]/(double)cntR[i];
				  if(resiR[i] > 0)
					{
					  for(int j=0;j<ceil(resiR[i]);j++)
						{
						  outf<<currChr<<"\t"<<lowPos + i-readLen-2<<"\t"<<lowPos+i-1<<"\tDUMMY\t";
						  outf<<resiR[i]<<"\t-"<<endl;
						}
					}

				}
			}
		  if(talk)
			cout<<"done."<<endl;
		  lowPos = highPos+1;
		  delete[] resiF;
		  delete[] resiR;
		  delete[] cntF;
		  delete[] cntR;
		  gsl_matrix_free(XF);
		  gsl_vector_free(yF);
		  gsl_vector_free(rF);
		  gsl_vector_free(cF);
		  gsl_matrix_free(covF);

		  gsl_matrix_free(XR);
		  gsl_vector_free(yR);
		  gsl_vector_free(rR);
		  gsl_vector_free(cR);
		  gsl_matrix_free(covR);
		}
	  if(!talk)
		cout<<endl;
	}


  // close files. 
  inf.close();
  for (int i = 0;i<ctrlFiles.size();i++)  
	infc[i].close();
  outf.close();
  string statFname = "readStats.txt";
  bool writeStats = true;
  ofstream ofc;
  ofc.open(statFname.c_str(),ios::trunc);
  if (ofc.fail())
    {
      failmessage.clear();
      failmessage.append("ERROR: Output file \"");
      failmessage.append(statFname.c_str());
      failmessage.append("\" could not be created, skipping.\n");
      writeStats = false;
    }
  
  if(writeStats)
    {
	  ofc <<"Chip reads"<<endl<<"Name\t"<<"minCrd\t"<<"maxCrd\t"<<"F_counts\t"<<"R_counts\t"<<endl;
      for ( it=seqMapChip.begin() ; it != seqMapChip.end(); it++ )
		{
		  ofc << (*it).first << "\t" << (*it).second.minPos << "\t" << (*it).second.maxPos<<"\t";
		  ofc << (*it).second.countF << "\t" << (*it).second.countR;
		  ofc <<endl;
		}
	  for (int i = 0;i<ctrlFiles.size();i++)
		{
		  ofc<<ctrlFiles[i]<<endl;
		  ofc <<"Control reads"<<endl<<"Name\t"<<"minCrd\t"<<"maxCrd\t"<<"F_counts\t"<<"R_counts\t"<<endl;
		  for ( it=seqMapCtrls[i].begin() ; it != seqMapCtrls[i].end(); it++ )
			{
			  ofc << (*it).first << "\t" << (*it).second.minPos << "\t" << (*it).second.maxPos<<"\t";
			  ofc << (*it).second.countF << "\t" << (*it).second.countR;
			  ofc <<endl;
			}
		}
    }else{
    cerr<<failmessage.c_str()<<endl;
  }
  ofc.close();
  
  return(0);
}



  
