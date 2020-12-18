#include "prd.h"
#include "prdprot.h"

#define MINFLOAT -99999.0
#define _Min(a, b)       (((a) < (b)) ?  (a) : (b))
#define _MAX2(a, b)       ((a) > (b)) ?  (a) : (b)
#define _MAX3(a, b, c)   ((a) > (_MAX2((b), (c)))) ? (a) : (_MAX2((b), (c)))
#define m(x,y) Mat[x][y]
#define VerboseOutput 0
#define Penalty(x) (GapI+((x)-1)*GapE)
#define SMITHWATERMAN 1

static char acid[] = "ARNDCQEGHILKMFPSTWYVXZB";

float AlignGlobal(char *SeqX, char *SeqY, int LenX, int LenY, float CompMatrix[NAcd][NAcd],
                  float GapI, float GapE, char *AliX, char *AliY)
{

  float **Matrix;
  register int i;
  char *SeqXX, *SeqYY;
  
  Matrix = FloatMatrix(LenY,LenX);
  SeqXX = (char *)ckalloc(LenX*sizeof(char));
  SeqYY = (char *)ckalloc(LenY*sizeof(char));

  for( i=0; i<LenX; i++ )
    SeqXX[i] = AAChar2Bin(SeqX[i]);
  for( i=0; i<LenY; i++ )
    SeqYY[i] = AAChar2Bin(SeqY[i]);
  
  SumUpFast (SeqXX, SeqYY, LenX, LenY, GapI, GapE, CompMatrix, Matrix);
  traceback (SeqXX, SeqYY, LenX, LenY, GapI, GapE, Matrix, AliX, AliY);

  FreeFloatMatrix(Matrix,LenY);
  free(SeqXX);
  free(SeqYY);

  return(PercentId(AliX,AliY));
}


void 
SumUpFast (char *SeqX, char *SeqY, int LenX, int LenY, float GPI, float GPE, float CompMatrix[NAcd][NAcd], float **Matrix)
{
  int y, x;
  float maxy, Dia;
/*  float max_val; */
/*  int max_x, max_y;*/
  float *maxx;

  maxx = (float *)ckalloc((LenX+1)*sizeof(float));
  
  for (y = 0; y < LenY; y++)
    for (x = 0; x < LenX; x++)
      Matrix[y][x] = CompMatrix[SeqY[y]][SeqX[x]];
  /* DumpMatrix2(LenY-10, LenY-1, LenX-10, LenX-1); */

  for (x = LenX - 1; x >= 0; x--) {
    maxx[x] = MINFLOAT;
  }

  for (y = LenY - 2; y >= 0; y--) {
    maxy = MINFLOAT;
    for (x = LenX - 2; x >= 0; x--) {
      if (y < LenY - 2) {
	maxx[x + 1] = _MAX2 (Matrix[y + 2][x + 1] - GPI,
			     maxx[x + 1] - GPE);
      }
      Dia = Matrix[y + 1][x + 1];
      Matrix[y][x] += _MAX3 (Dia, maxy, maxx[x + 1]);
#if SMITHWATERMAN
/*
  if (Matrix[y][x] < 0) {
	Matrix[y][x] = 0;
      }
      if (Matrix[y][x] > max_val)
	max_val = Matrix[y][x], max_x = x, max_y = y;
        */
#endif /* SMITHWATERMAN */
      maxy = _MAX2 (Matrix[y + 1][x + 1] - GPI,
		    maxy - GPE);
    }
  }
  free(maxx);
  
}



void 
traceback (char *SeqX, char *SeqY, int LenX, int LenY, float GapI, float GapE, float **Matrix, char *AliX, char *AliY )
{
  int x, y, x1, y1, x2, y2;
  int MaxX, MaxY, MaxPosX, MaxPosY;
  float Max;
  char *pX, *pY;
/*  int Gap, i;*/

/*  int NGap, NIdentical, NAligned; */

  /* printf("traceback\n"); */
  MaxX = LenX - 1;
  MaxY = LenY - 1;
  pX = AliX;
  pY = AliY;
  x2 = y2 = 0;
  x = MaxX;
  y = MaxY;
  Max = -999;
  for (x = 0; x <= MaxX; x++) {
    if (Matrix[0][x] > Max) {
      Max = Matrix[0][x];
      MaxPosX = x;
      MaxPosY = 0;
    }
  }

  for (y = 0; y <= MaxY; y++) {
    if (Matrix[y][0] > Max) {
      Max = Matrix[y][0];
      MaxPosX = 0;
      MaxPosY = y;
    }
  }

  while (x2 < MaxPosX) {
    *pY = '.';
    pY++;
    *pX = AABin2Char (SeqX[x2]);
    x2++;
    pX++;
  }
  while (y2 < MaxPosY) {
    *pX = '.';
    pX++;
    *pY = AABin2Char (SeqY[y2]);
    y2++;
    pY++;
  }

  *pY = AABin2Char (SeqY[y2]);
  y2++;
  *pX = AABin2Char (SeqX[x2]);
  x2++;
  pX++;
  pY++;

  /* printf("%3d %3d\n", MaxPosX, MaxPosY); */
  x = MaxPosX + 1;
  y = MaxPosY + 1;
  while ((x <= MaxX) && (y <= MaxY)) {
    Max = Matrix[y][x];
    MaxPosX = x;
    MaxPosY = y;
/*    Gap = 0;*/

    for (y1 = y + 1; y1 <= MaxY; y1++) {
      if (Matrix[y1][x] - Penalty (y1 - y) > Max) {
	Max = Matrix[y1][x] - Penalty (y1 - y);
	MaxPosX = x;
	MaxPosY = y1;
/*	Gap = 1;*/
      }
    }

    for (x1 = x + 1; x1 <= MaxX; x1++) {
      if (Matrix[y][x1] - Penalty (x1 - x) > Max) {
	Max = Matrix[y][x1] - Penalty (x1 - x);
	MaxPosX = x1;
	MaxPosY = y;
/*	Gap = 1;*/
      }
    }

/*
    if (Gap)
      NGap++;
      if (Gap) printf("Gap  "); */

    while (x2 < MaxPosX) {
      *pY = '.';
      pY++;
      *pX = AABin2Char (SeqX[x2]);
      x2++;
      pX++;
    }
    while (y2 < MaxPosY) {
      *pX = '.';
      pX++;
      *pY = AABin2Char (SeqY[y2]);
      y2++;
      pY++;
    }

    *pY = AABin2Char (SeqY[y2]);
    y2++;
    *pX = AABin2Char (SeqX[x2]);
    x2++;
    pX++;
    pY++;

    /* printf("%3d %3d\n", MaxPosX, MaxPosY); */

    x = MaxPosX + 1;
    y = MaxPosY + 1;
  }
  *pY = 0;
  *pX = 0;

/*
  for (i = 0; i < strlen (AliX); i++) {
    if (isalpha (AliX[i]) && isalpha (AliY[i])) {
      NAligned++;
      if (AliX[i] == AliY[i])
	NIdentical++;
    }
  } */
}

char 
AAChar2Bin (char AA)
{
  char *p;
  char BB;

  if( AA == 'X' )
    BB = 'A';
  else
  if( AA == 'B' )
    BB = 'N';
  else
  if( AA == 'Z' )
    BB = 'Q';
  else
    BB = AA;
    
  p = strchr (acid, toupper (BB));
  if (p) {
    return (char) (p - acid);
  } else {
    printf ("Error in Sequence: Wrong Character: %d=|%c|\n", AA, AA);
    exit (-1);
  }
  return 0;
}

int 
AABin2Char (char AA)
{
  return acid[AA];
}


