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



#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <iterator>
#include <map>
// itoa & random stuff
#include <ctime>
#include <sys/time.h>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cfloat>
#include <cstring>

#include <math.h>

typedef unsigned short storageType; 

using namespace std;

/*
 * spits the current string on delimiters
 *
 */ 
void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters = " ");

/*
 * as Tokenize but returns the 'next' position to start from.
 *
 */ 
string::size_type readAndTokenize(const string& str,
				  vector<string>& tokens,
				  const string& delimiters = " ",
				  string::size_type firstPos = 0);

/*
 * used to keep track of some statistics of each sequence.
 *  
 */
struct seqStats{
  int minPos;
  int maxPos;
  int truncR;
  int truncF;
  int truncC;
  int countF;
  int countR;
  map <int,int*> histF;
  map <int,int*> histR;
};

/*
 *
 * used to store data from a single parsed input-line. 
 *
 */
struct inputLine{
  int strand;      // 1/-1
  int pos;         // 1-based leftmost position.
  int len;         // length of mapped fragment, determined by seq. 
  string seq;      // name of reference sequence.  
  string name;     // name of the BED line.
  int mapped;      // was this read mapped? (1/0)
  int header;      // is this line a header line (1/0), obsoletes the above. 
  storageType val; // if input specified values are used. 
};

/*
 *  
 *  Extract the necessary information from a single line of the BED-file
 *
 * Column-numbers are 0-based.
 */
int parseBEDline(string str,inputLine*,int seqCol,int startCol,
				 int endCol,int strandCol,int nameCol = 3);     


/*
 *
 * returns number of lines in infile + high/low coordinates of all seqs. 
 */

int initControlBEDlite(ifstream * inf,map<string,seqStats> * seqmap,int seqCol,int startCol,
					   int endCol,int strandCol,bool verbose);

#endif
