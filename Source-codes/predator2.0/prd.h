#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>

#define BUFSZ         1024
#define NAcd            20
#define NPROP            7
#define NTHRESH    2*NPROP
#define MAX_LEN      10000
#define MAX_FIELD       50
#define MBRWIN          13
#define MBRWIN22 (MBRWIN/2)
#define MBRNBEST        25
#define MAXDB         1000
#define WINBETA          4
#define WINALPHA         7
#define PARDIST          5 
#define ANTIDIST         4
#define MINHEL           6
#define MINBET           4
#define MINCOIL          3
#define MAX_ASSIGN    1000
#define MAX_FIELD       50
#define MAX_SEQ       5000
#define OUTWID          50
#define INFTHR         0.1
#define NBIN             5
#define MAXREL        NBIN*NBIN*NBIN*NBIN*NBIN*NBIN*NBIN

#define Minimum(x,y)              ((x)<(y) ? x : y)
#define Maximum(x,y)              ((x)<(y) ? y : x)
#define Eps                       0.000001
#define Inf(x)                    ((-x)*log((x)))
#define INSIDE(x,left,right)      ( (x >= (left)) && (x <= (right)) )


#define TRUE                      1
#define FALSE                     0
#define ERR                      -1

#define NMARKS 50

typedef char BUFFER[BUFSZ+1];
#ifndef ICMH__
typedef char bool;
#endif
#define NFILETYPES 5
enum FILETYPE {Fasta, Clustal, MSF, Stride, Dssp, GCG, EMBL, PIR, GDE, Hssp, Extern, Unknown};

typedef struct { 
                 BUFFER InputFile, OutFile, OutSeqFile, PredId;
                 BUFFER StrideFile, DsspFile, DatabaseFile;
		 bool Dssp, Info, Long, All, OutSeq;
                 bool Verbose, OrigAli, LookUp, Single;
                 char ChainId;
                 float Safety, NonRed;
                 int GapVal, DelVal;
		 int NAli, MinAliLength;
	       } COMMAND;

typedef struct {
                 char *Id, *De, *Se, *AliSe;
		 char *Pred, *Asn, *KnownAsn, *Hml[MAXDB], PdbChain;
                 int Len, AliLen, *AliCnt, NHml;
                 float **Prop, PercCorr, *Rel, *PHel, *PBeta, *PCoil;
		 bool ToBePredicted, Selected;
                 enum FILETYPE FileType;
	       } SEQ;

typedef struct {
                 float *Percent, *HCri;
		 int NAli, *Nc, *Score, *Accepted;
		 int *Len, *Beg1, *Beg2, *End1, *End2;
		 char **A1, **A2;
	       } ALN;

typedef struct {
                 int NAli;
                 float *Score, *HCri, *Percent;
		 int *ResN, *AliLength;
		 SEQ **Seq;
	       } PROJ;

#define DIRDELIM '/'

/* #define MSDOS 1 */

#if MSDOS
#define DIRDELIM '\\'
#endif
