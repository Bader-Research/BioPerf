////////////////////////////////////////////////////////////////
//  Author:  Olga Troyanskaya
//
//    File:  Glimmer/rnabin.h
//
//  Version:  2.03  9 Dec 2002
//
//    Copyright (c) 1999, 2000, 2002 by Olga Troyanskaya, Arthur Delcher,
//    Steven Salzberg, Simon Kasif, and Owen White.  All rights reserved.
//    Redistribution is not permitted without the express written
//    permission of the authors.
//
//  Desc: Function and variable declarations for the RNAbin function 
//////////////////////////////////////////////////////////////////

#ifndef _RNABIN_H_INCLUDED
#define _RNABIN_H_INCLUDED


#include  "delcher.h"
#include  "gene.h"


#define nA 0
#define nC 1
#define nG 2
#define nT 3
#define S 0
#define E 1
#define  DEFAULT_MIN_GENE_LEN  120
#define DEFAULT_MAX_DIST_START 150
#define DEFAULT_PERCENT_HIGH 20   //percent of scores used to recalculate the matrix
#define DEFAULT_PERCENT_LOW 10
#define DEFAULT_TRESHOLD 0.00001

// *******************LOOK HERE FOR DEBUGGING****************** 
//#define DEBUG //uncomment this if you want the program to print results out to the terminal window

  typedef struct {
    double score; //score
    int index; //index of sequence to which score corresponds
  } Scores;

//constants
const int MAX_GENE_NUMBER = 5000;
const int MAX_SEQ_LENGTH = 20; 
const int BINDING_SIZE = 7; //length of RNA binding site

//!!!!---FUNCTION DECLARATIONS---!!!!
//helper functions (implementations in RNAbin.h)
//void fileRead (const char *); //read the file 
void matrixCopy(double [][4], double [][4], int); //copies first argument into second argument, second argument is the number of cols in matrix
void matrixPrint(double [][4], int); //prints out the matrix, second argument is the number of cols

//other functions
void fillProbMatrix(int, int, char [][MAX_SEQ_LENGTH+1], double [][4]); //calculates probabilities of each nucleotide being at each position, and stores those probabilites in the matrix
void reduceMatrix(int, double [][4], double [][4]); //reduces the matrix to BINDING_SIZE
void align(int, int, int*, char [][MAX_SEQ_LENGTH+1], double [][4]);  //aligns sequences and puts aligned coors into alignedStartInd
void fillReducedMatrix(int, int, int*, char [][MAX_SEQ_LENGTH+1], double [][4], double [][4]); //recomputes reducedMatrix and stores the old reducedMatrix in oldMatrix
double maxChange(int, double [][4], double [][4]); //returns the largest change between new and old matrices, argument is num of cols
void printConsensus(int*, int, double [][4]); //prints out consensus sequence, first argument is alignedStartInd[], second argument is size of consensus (BINDING_SIZE)
void computeScores(int, int, int*, double*, char [][MAX_SEQ_LENGTH+1], double [][4]); //computes scores for each aligned sequence using reducedMatrix and stores scores in array passed to it (alignedSeqScores), first argument is string length, second argument is gene number, 3rd argument - alignedStartInd, 4th argument - scores array.
void sortSeq(int, int&, int&, double*, char [][MAX_SEQ_LENGTH+1], long int*, char [][MAX_SEQ_LENGTH+1]); //passed numGenes, scores array, finds low cutoff score and copies all sequences above that score into a separate array, which it is passed as a second argument, and all sequences below that score into a separates array, which is passed as the third argument
double scoreRegion(char *, int, int&, double [][4]); //scores the region that's passed as an argument agains reducedMatrix, the second argument is the length, third argument is winStartInd
void findBetterStart(char *, int, long int *, long int [][2], int, int, double [][4], char [][MAX_SEQ_LENGTH+1], int *, long int, long int *);
void retrieveSeqs(char *, long int, int, long int [][2], int *, char [][MAX_SEQ_LENGTH+1]);

//!!!!---IMPLEMENTATIONS OF SOME HELPER FUNCTIONS---!!!!

void matrixCopy(double fromMtrx[][4], double  toMtrx[][4], int maxCol) {
  int row =0;
  int col = 0;
  for (col = 0; col < maxCol; col++) {
    for (row = 0; row < 4; row++) {
      toMtrx[col][row] = fromMtrx[col][row]; 
    }
  }
}

void matrixPrint(double  matrix[][4], int maxCol) {
  for (int j=0; j<4; j++) {
    for (int i=0; i < maxCol; i++) {
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }
}

#endif //!_RNAbin_h



