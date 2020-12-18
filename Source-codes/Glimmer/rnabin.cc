////////////////////////////////////////////////////////////////
//  Author:  Olga Troyanskaya
//
//    File:  Glimmer/rnabin.cc
//
//  Version:  2.03  9 Dec 2002
//
//    Copyright (c) 1999, 2000, 2002 by Olga Troyanskaya, Arthur Delcher,
//    Steven Salzberg, Simon Kasif, and Owen White.  All rights reserved.
//    Redistribution is not permitted without the express written
//    permission of the authors.
//
//    Desc:  a function for Glimmer that finds RNA binding sites prior 
//           to start sites passed to it.
//
//
//  to compile:  use g++ or other C++ compiler
//               the program requires  rnabin.h 
////////////////////////////////////////////////////////////////


#include "rnabin.h"


void RNAbin (char * Data, long int genomeLen, long int coords[][2], long int newStartsArr [], long int N) {
  
  //variables
 
  char highSeq [N] [MAX_SEQ_LENGTH+1];
  int complement[N];  //1 if complement, 0 if not
  long int coordsL[N]; //storer index of where in coords array the coordinates of start and end position of the low scoring sequences are
  char sequences[N][MAX_SEQ_LENGTH+1];
  double probMatrix[MAX_SEQ_LENGTH][4]; //matrix[col][row], row = P(A), P(C), P(G), P(T), col = position of nucleotide, matrix contains P(nucleotide) at each position
  double oldMatrix[BINDING_SIZE][4]; //used to store old copy of probMatrix
  double reducedMatrix[BINDING_SIZE][4]; //matrix reduced from probMatrix by reduceMatrix function
  int alignedStartInd[N]; //contains start coords for best aligned window
  double alignedSeqScores[N]; //array stores scores for each aligned sequence with respect to the final matrix  
  //int seqNum=0;
  int numGenesHigh, numGenesLow;
  double mc;  //largest single probability change in the matrix
  double treshold = DEFAULT_TRESHOLD;
  long int numSequences = N-1;

  retrieveSeqs(Data, genomeLen, numSequences, coords, complement, sequences);

  fillProbMatrix(MAX_SEQ_LENGTH, numSequences, sequences, probMatrix);
  reduceMatrix(MAX_SEQ_LENGTH, probMatrix, reducedMatrix);
  
  do {
    align(numSequences, MAX_SEQ_LENGTH, alignedStartInd, sequences, reducedMatrix);
    fillReducedMatrix(BINDING_SIZE, numSequences, alignedStartInd, sequences, reducedMatrix, oldMatrix);
    mc = maxChange(BINDING_SIZE, reducedMatrix, oldMatrix);
  } while (mc > treshold);
 
  computeScores(BINDING_SIZE, numSequences, alignedStartInd, alignedSeqScores, sequences, reducedMatrix);
  sortSeq(numSequences, numGenesHigh, numGenesLow, alignedSeqScores, highSeq, coordsL, sequences);

  //repeat process but for high scoring sequences only
  fillProbMatrix(MAX_SEQ_LENGTH, numGenesHigh, highSeq, probMatrix);
  reduceMatrix(MAX_SEQ_LENGTH, probMatrix, reducedMatrix);
  
  do {
    align(numGenesHigh, MAX_SEQ_LENGTH, alignedStartInd, highSeq, reducedMatrix);
    fillReducedMatrix(BINDING_SIZE, numGenesHigh, alignedStartInd, highSeq, reducedMatrix, oldMatrix);
    mc = maxChange(BINDING_SIZE, reducedMatrix, oldMatrix);
  } while (mc > treshold);

#ifdef DEBUG  
  cout << endl << endl << "for high scoring sequences: " << endl;
  cout << "num: " << numGenesHigh << endl;
  matrixPrint (reducedMatrix, BINDING_SIZE);
  printConsensus(alignedStartInd, BINDING_SIZE, reducedMatrix);

  cout << endl << "for low scoring sequences: " << endl;
  cout << "num: " << numGenesLow << endl;
#endif

  findBetterStart(Data, numGenesLow, coordsL, coords, BINDING_SIZE, DEFAULT_MIN_GENE_LEN, reducedMatrix, sequences, complement, numSequences, newStartsArr);
}



void retrieveSeqs(char * Data, long int Len, int genesNum, long int coords[][2], int complement[], char sequences [][MAX_SEQ_LENGTH+1]) {
  int Index = 0, i = 0, j = 0, End = 1;
    
  for (Index = 0; Index <= genesNum; Index++) {
    complement[Index] = (coords[Index][S] > coords[Index][E]);
    if (!complement[Index]) { //gene on non-complement strand
    	i= coords[Index][S] - 21;
	End = coords[Index][S]-1;
        if  (End < 1)
            End += Len;
	j=0;
	do
	  {
	    i ++;
	    if  (i < 1)  //to handle wraparounds
	        i += Len;
            if  (i > Len)
                i -= Len;
	    sequences[Index][j] = Data[i];
	    j++;
	  }  while  (i != End);
	sequences[Index][j] = '\0';
    } else { //gene on complementary strand
    	i = coords[Index][S] + 21;
	End = coords[Index][S] + 1;
        if  (End > Len)
            End -= Len;
	j=0;
	do
	  {
	    i --;
            if  (i > Len)
                i -= Len;
            if  (i < 1)
                i += Len;
	    sequences[Index][j] = Complement(Data[i]);
	    j++;
	  }  while  (i != End); 
	sequences[Index][j] = '\0';
    }
  }
}

void fillProbMatrix(int strLen, int numSequences, char sequences[][MAX_SEQ_LENGTH+1], double probMatrix[][4]) {
  
  int ind = 0; 
  int position = 0; //position in sequence being examined
  double numA=0, numC=0, numG=0, numT=0; //number of each nucleotide in a specified position in sequence
  double numberOfSeqs = numSequences+1;
  
  for (position = 0; position < strLen; position++) {
    for (ind = 0; ind <= numSequences; ind++) {
      switch (sequences[ind][position])  { 
         case 'A':
         case 'a':
	   numA++;
	   break;
         case 'C':
         case 'c':
	   numC++;
	   break;
         case 'G':
         case 'g':
	   numG++;
	   break;
         case 'T':
         case 't':
	   numT++;
	   break;
      }//switch
    }//for(index...
    probMatrix[position][nA] = numA/numberOfSeqs;
    probMatrix[position][nC] = numC/numberOfSeqs;
    probMatrix[position][nG] = numG/numberOfSeqs;
    probMatrix[position][nT] = numT/numberOfSeqs;
    numA = 0;
    numC = 0;
    numG = 0;
    numT = 0;
  }//for(position...
 
}

void reduceMatrix(int strLen, double probMatrix[][4], double reducedMatrix[][4]) {
  int row = 0;
  int col = 0; //column
  double max = 0;
  double scores[strLen];  //array to store scores in
  int numGroups = strLen - BINDING_SIZE + 1; //the number of groups of 8 cols to examine 
  int groupStartCol = 0;  //first column in group with highest skew
  
  //calculate the multiple of max P(nucl) at each column for all possible groups of 8 adjacent columns to find the group with highest bias (highest score)
  for (col = 0; col < strLen; col++) {
    for (row = 0; row < 4; row++) { //find max P(nucl) in column
      if (probMatrix[col][row] >= max) {
	max = probMatrix[col][row];
      }
    }
 
    scores[col] = max;
    max = 0;
  }

  double groupScore = 1.0, maxGroupScore = 0.0;
  int lcv = 0;
  
  for (col = 0; col < numGroups; col++) {
    for (lcv = col; lcv <= (col + BINDING_SIZE - 1); lcv++ ) {
      groupScore = scores[lcv] * groupScore;
    }  
 
    if (groupScore >= maxGroupScore) {
      //if this group is more skewed then any previously examined group
      maxGroupScore = groupScore;
      groupStartCol = col;
    }
    groupScore = 1;
  }
   
  for (col = 0; col < BINDING_SIZE; col++) {
    for (row = 0; row < 4; row++) {
      reducedMatrix[col][row] = probMatrix[groupStartCol + col][row];
    }
  }
}


void align(int geneNum, int strLen, int * winStartInd,  char sequences[][MAX_SEQ_LENGTH+1], double reducedMatrix[][4]) {
  //functions recomputes the matrix by finding better alignment of sequences
  int numWindows = strLen - BINDING_SIZE + 1; //the number of windows to examine 
  double maxWindScore = 0;
  //score is a multiple of P(N) for all nucleotides in a window
  //N is the nucleotide in the position examined and P(N) is the P of N being in that position 
  //P(N) is from reducedMatrix
  
  double score = 1;
  int win=0, position=0, lcv=0, ind=2000;
  
  for (ind = 0; ind <= geneNum; ind++) {
    for (win = 0; win < numWindows; win++) {
      for (lcv = win; lcv <= (win + BINDING_SIZE - 1); lcv++ ) {
	switch (sequences[ind][lcv])  { 
	 case 'A':
         case 'a':
	   score *= reducedMatrix[position][nA];
	   break;
         case 'C':
         case 'c':
	   score *= reducedMatrix[position][nC];
	   break;
         case 'G':
         case 'g':
	   score *= reducedMatrix[position][nG];
	   break;
         case 'T':
         case 't':
	   score *= reducedMatrix[position][nT];
	   break;
	};  //switch
	position++;
      }
      if (score > maxWindScore) {
	maxWindScore = score;
	winStartInd[ind] = win;
      }
      score = 1;
      position = 0;
     
    }    
    maxWindScore = 0;
  }
}


void fillReducedMatrix(int strLen, int numSeq, int * winStartInd, char sequences[][MAX_SEQ_LENGTH+1], double reducedMatrix[][4], double oldMatrix [][4]) {
  //recomputes reducedMatrix according to new alignment and copies the old matrix into oldMatrix
  
  
  int position = 0; //position in sequence being examined
  int ind=0; //index of sequence in sequences array
  double numA=0, numC=0, numG=0, numT=0; //number of each nucleotide in a specified position in sequence
  double numberOfSeqs = numSeq + 1;
  int lcv = 0;

  matrixCopy(reducedMatrix, oldMatrix, BINDING_SIZE); //store old probability matrix in oldMatrix
  
  for (lcv = 0; lcv < strLen; lcv++) {
    for (ind = 0; ind < numberOfSeqs; ind++) {
      position = winStartInd[ind] + lcv; //find position in aligned sequence
      switch (sequences[ind][position])  { 
         case 'A':
         case 'a':
	   numA++;
	   break;
         case 'C':
         case 'c':
	   numC++;
	   break;
         case 'G':
         case 'g':
	   numG++;
	   break;
         case 'T':
         case 't':
	   numT++;
	   break;
      }//switch
    }//for(index...
    reducedMatrix[lcv][nA] = numA/numberOfSeqs;
    reducedMatrix[lcv][nC] = numC/numberOfSeqs;
    reducedMatrix[lcv][nG] = numG/numberOfSeqs;
    reducedMatrix[lcv][nT] = numT/numberOfSeqs;
    numA = 0;
    numC = 0;
    numG = 0;
    numT = 0;
  }//for(position...
}


double maxChange(int maxCol, double reducedMatrix[][4], double oldMatrix[][4]) {
  double change = 0, mChange = 0;
  int row =0;
  int col = 0;
  
  for (col = 0; col < maxCol; col++) {
    for (row = 0; row < 4; row++) {
      change = fabs(reducedMatrix[col][row] - oldMatrix[col][row]); 
      if (change >= mChange) {
	mChange = change;
      }
    }
  }
  return mChange;
}

void printConsensus(int * winStartInd, int maxCol, double reducedMatrix[][4]) {
  int row = 0, col = 0;
  int bestMatchNuclRow = 0;
 
  for (col = 0; col < maxCol; col++) {
    for (row = 0; row < 4; row++) {
      if (reducedMatrix[col][row] > reducedMatrix[col][bestMatchNuclRow]) {
	bestMatchNuclRow = row;
      }
    }
    switch (bestMatchNuclRow) {
    case nA:
      cout << "A";
      break;
    case nC:
      cout << "C";
      break;
    case nG:
      cout << "G";
      break;
    case nT:
      cout << "T";
      break;
    }
    bestMatchNuclRow = 0;
  }
  cout << endl;
}


void computeScores(int strLen, int numGenes, int* winStartInd, double* secScores, char sequences[][MAX_SEQ_LENGTH+1], double reducedMatrix[][4]) {
  int ind = 0;
  int lcv = 0;
  int position = 0;
  double score = 1;
  for (ind = 0; ind <= numGenes; ind++) {
    for (lcv = 0; lcv < strLen; lcv++) {
       position = winStartInd[ind] + lcv; //find position in aligned sequence
        switch (sequences[ind][position])  { 
         case 'A':
         case 'a':
	   score *= reducedMatrix[lcv][nA];
	   break;
         case 'C':
         case 'c':
	   score *= reducedMatrix[lcv][nC];
	   break;
         case 'G':
         case 'g':
	   score *= reducedMatrix[lcv][nG];
	   break;
         case 'T':
         case 't':
	   score *= reducedMatrix[lcv][nT];
	   break;
      }//switch
    }
    secScores[ind] = score;
    score = 1;
  }
} 

int Cmp (const void * A, const void * B)
  //comparison function for qsort to sort indexedScores by descending score
{
   double A_sc = (* ((Scores *) A)).score;
   double B_sc = (* ((Scores *) B)).score;

   if  (A_sc > B_sc) { 
     return  -1;
   }
   else if  (A_sc < B_sc) {
     return  1;
   }
   else {
       return  0;
   }
}




void sortSeq
    (int numGenes, int &numGenesH, int &numGenesL, double* secScores,
     char highS [] [MAX_SEQ_LENGTH+1],  long int coordsL[],
     char sequences[][MAX_SEQ_LENGTH+1]) {
  int ind, posH=0, posL=0, cutOffH = 0, cutOffL = 0;
  Scores * indexedScores;

  indexedScores = (Scores *) malloc ((numGenes+1) * sizeof (Scores));
  assert (indexedScores != NULL);

  for (ind = 0; ind <= numGenes; ind++) {
    indexedScores[ind].score = secScores[ind];
    indexedScores[ind].index = ind;
  }

  qsort(indexedScores, (numGenes+1), sizeof(Scores), Cmp);

  cutOffH = (numGenes+1)*DEFAULT_PERCENT_HIGH/100;
  cutOffL = (numGenes+1)*(100 - DEFAULT_PERCENT_LOW)/100;

  for (ind = 0; ind <= cutOffH; ind++) {
    memmove (highS[posH], sequences[(indexedScores[ind].index)], MAX_SEQ_LENGTH + 1);
    posH++;
  }

  for (ind = numGenes; ind >= cutOffL; ind--) {
    coordsL[posL] = (indexedScores[ind].index);
    posL++;
  }

  numGenesH = posH-1;
  numGenesL = posL-1;

  free (indexedScores);
}



double scoreRegion(char *region, int strLen, int &winStartInd, double reducedMatrix[][4]) {
  int numWindows = strLen - BINDING_SIZE + 1; //the number of windows to examine 
  
  double maxWindScore = 0;
  //score is a multiple of P(N) for all nucleotides in a window
  //N is the nucleotide in the position examined and P(N) is the P of N being in that position 
  //P(N) is from reducedMatrix

  double score = 1;
  int win=0, position=0, lcv=0;
  int skippedZero = FALSE;

	  for (win = 0; win < numWindows; win++) {
	    for (lcv = win; lcv <= (win + BINDING_SIZE - 1); lcv++ ) {
	      switch (region[lcv])  { 
	      case 'A':
	      case 'a':
		score *= reducedMatrix[position][nA];
		break;
	      case 'C':
	      case 'c':
		score *= reducedMatrix[position][nC];
		break;
	      case 'G':
	      case 'g':
		score *= reducedMatrix[position][nG];
		break;
	      case 'T':
	      case 't':
		score *= reducedMatrix[position][nT];
		break;
	      };  //switch
	      position++;
	    }
	    
	    if (score > maxWindScore) {
	      maxWindScore = score;
	      winStartInd = win;
	    }
	    skippedZero = FALSE;
	    position = 0;
	    score = 1;
 	  }
	  return maxWindScore;
}


void findBetterStart(char * data, int numL, long int * coordsL, long int coords[][2], int strLen, int minGeneLength, double reducedMatrix[][4], char sequences[][MAX_SEQ_LENGTH+1], int *complement, long int numSequences, long int *newStartsArr)

{
  int indexOfSeq = 0, coordInGenome = 1, i=0, lcv = 0, start = FALSE;
  double maxScore = 0, score =0, origScore=0;
  long int newStartCoord = 1, s=1, e=1, c = 1, newStart = 1;
  int winStart = 0;
  char seq[24];
#ifdef DEBUG
  int count = 0;
  char pattern[9];
  char newstart[4];
  char oldstart[4];
#endif 
  long int maxDistFromStart = DEFAULT_MAX_DIST_START;
  long int dist = 0, quit = FALSE; 
 

  for (lcv = 0; lcv <= numSequences; lcv++) {
    newStartsArr[lcv] = coords[lcv][S]; //init newStarts to original coords
  }
  lcv = 0;

  for (indexOfSeq = 0; indexOfSeq <= numL; indexOfSeq++) {
    
    maxScore = 0;
    s = coords[(coordsL[indexOfSeq])][S];
    e = coords[(coordsL[indexOfSeq])][E];
    origScore = scoreRegion(sequences[(coordsL[indexOfSeq])], MAX_SEQ_LENGTH, winStart, reducedMatrix);

#ifdef DEBUG    
    //for printout purposes - can remove
    origSeq = sequences[(coordsL[indexOfSeq])];
    int y;
    for (y = 0; y < BINDING_SIZE; y++) {
      origSeqMatch[y] = origSeq[(winStart + y)];
    }
    origSeqMatch[y] = '\0';
    //end printout region
#endif
    
    newStart = s;

    if (!complement[(coordsL[indexOfSeq])]) {  //process sequences that are not complements
      for (c = s+3; c < (e - minGeneLength); ) {
	start = FALSE;

#if  0   // Olga's old code
         if ((data[c] == 't') || (data[c] == 'a') || (data[c] == 'g')){ 
	  if (((data[c+1] == 't') && (data[c+2] == 'g'))) {
	   start = TRUE;
	 }
       }
#endif

// Debra's new version
       if ((data[s]=='a') && (data[s+1]=='t') && (data[s+2]=='g')) {
          if ((data[c]=='a') && (data[c+1]=='t') && (data[c+2]=='g'))
          start = TRUE;
       } // for atg old start only new atg start is valid

       if ((data[s]=='g') && (data[s+1]=='t') && (data[s+2]=='g')) {
          if ((data[c]=='a') || (data[c]=='g')) {
             if ((data[c+1]=='t') && (data[c+2]=='g'))
             start = TRUE;
          }
       } // for gtg old start, either atg or gtg new start is valid

      if ((data[s]=='t') && (data[s+1]=='t') && (data[s+2]=='g')) {
         if ((data[c]=='a') || (data[c]=='g') ||(data[c]=='t')) {
            if ((data[c+1]=='t') && (data[c+2]=='g'))
            start = TRUE;
         }
      } // for ttg old start, all new starts, atg, ttg and gtg, are valid


       if (start) {
	 
	 for (coordInGenome = c; coordInGenome < (c + MAX_SEQ_LENGTH); coordInGenome++) {
	   seq[i] = data[(coordInGenome-MAX_SEQ_LENGTH)];
	   i++;
	 }
	 seq[i] = '\0';
	 score = scoreRegion(seq, MAX_SEQ_LENGTH, winStart, reducedMatrix);

	 if (score > maxScore) { 
	   dist = abs((s-c));
	   if (dist <= maxDistFromStart) { 

	   	   maxScore = score;
		   newStart = s;
		   newStartCoord = c;
#ifdef DEBUG
		   //for test printout only
		   oldstart[0] = data[s];
		   oldstart[1] = data[s+1];
		   oldstart[2] = data[s+2];
		   oldstart[3] = '\0';
		   newstart[0] = data[c];
		   newstart[1] = data[c+1];
		   newstart[2] = data[c+2];
		   newstart[3] = '\0';
		   //end for test printout
		   
		   for (count = 0; count < BINDING_SIZE; count++) {
		     pattern[count] = seq[(winStart + count)];
		   }
		   pattern[count]='\0';
#endif
	   }
	   else quit = TRUE;
	 }
	 i = 0;
       }
       if (quit) break; //break out when went over the maxDistFromStart
       c+=3; //increment c by 3 to keep in the same reading frame
     }
     quit = FALSE;
   } else { //process sequences from complementary strand
     
	 for (c = s-3; c > (e + minGeneLength); ) {
	   start = FALSE;

#if  0  // Olga's old code
	   if ((data[c] == 't') || (data[c] == 'a') || (data[c] == 'c')){ 
	     if (((data[c-2] == 'c') && (data[c-1] == 'a'))) {
	       start = TRUE;
	     }
	   } 
#endif
	   
// Debraj's new version
        if ((data[s]=='t') && (data[s-1]=='a') && (data[s-2]=='c')) {
          if ((data[c]=='t') && (data[c-1]=='a') && (data[c-2]=='c'))
          start = TRUE;
        } // for atg old start only new atg start is valid

        if ((data[s]=='c') && (data[s-1]=='a') && (data[s-2]=='c')) {
          if ((data[c]=='t') || (data[c]=='c')) {
             if ((data[c-1]=='a') && (data[c-2]=='c'))
             start=TRUE;
          }
       } // for gtg old start, either atg or gtg new start is valid

       if ((data[s]=='a') && (data[s-1]=='a') && (data[s-2]=='c')) {
         if ((data[c]=='t') || (data[c]=='c') ||(data[c]=='a')) {
            if ((data[c-1]=='a') && (data[c-2]=='c'))
            start=TRUE;
         }
       }


	   if (start) {
	     i=0;
	     for (coordInGenome = c+MAX_SEQ_LENGTH; coordInGenome > c; coordInGenome--) {
	       seq[i] = Complement(data[coordInGenome]);
	       i++;
    	     }
	     seq[MAX_SEQ_LENGTH] = '\0';
	     score = scoreRegion(seq, MAX_SEQ_LENGTH, winStart, reducedMatrix);

	     if (score > maxScore) { 
	       dist = abs((s-c));
	       if (dist <= maxDistFromStart) {
		 maxScore = score;
		 newStart = s;
		 newStartCoord = c;
#ifdef DEBUG		 
		 //for test printout only
		 oldstart[0] = Complement(data[s]);
		 oldstart[1] = Complement(data[s-1]);
		 oldstart[2] = Complement(data[s-2]);
		 oldstart[3] = '\0';
		 newstart[0] = Complement(data[c]);
		 newstart[1] = Complement(data[c-1]);
		 newstart[2] = Complement(data[c-2]);
		 newstart[3] = '\0';
		 //end for test printout
 
		 for (count = 0; count < BINDING_SIZE; count++) {
		   pattern[count] = seq[(winStart + count)];
		 }
		 pattern[count]='\0';
#endif
	       }
	       else quit = TRUE;
	     }
	     i = 0;
	   }
	   if (quit) break;
 	   c-=3; //decrement c by 3 to keep in the same reading frame
 	 }
	 quit = FALSE;
   }


   if (maxScore > origScore)
     {
       newStartsArr[(coordsL[indexOfSeq])] = newStartCoord;
#ifdef DEBUG
       dist = abs((s-newStartCoord));
       cout << "old: " << origSeqMatch << "\t" << s << "\t" << oldstart << "\t" << dist << "\tnew: " << pattern << "\t" << newStartCoord << "\t" << newstart << "\t" << e  << endl;
       t++;
#endif
     }
  }
#ifdef DEBUG
  cout << "TOTAL NUMBER OF NEW STARTS FOUND: " << t << endl;
#endif
}



