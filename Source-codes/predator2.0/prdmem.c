#include "prd.h"
#include "prdprot.h"

extern float **Bridge1, **Bridge2, **Bridge3;
extern float **Alpha, ***MbrMat, ThresholdsM[NTHRESH], ThresholdsS[NTHRESH];
extern int    *SeqDb[MAXDB], *StrDb[MAXDB], NStr;
extern char  *NmDb[MAXDB];
extern float **Rel;
extern int NRel;

void *ckalloc(size_t bytes)
{
  register void *ret;
  void die(char *format, ... );
  
  if( (ret = malloc(bytes)) == 0) die("Out of  memory\n");

  return ret;	
}

void *myrealloc(void *ptr, size_t bytes)
{
  register void *ret;
  
  if( (ret = realloc(ptr, bytes)) == 0 ) die("Out of  memory\n");

  return ret;	
}

float **FloatMatrix(int M, int N)
{
  int m;
  float **Matrix;

  Matrix = (float **)ckalloc(M*sizeof(float *));

  for( m=0; m<M; m++ ) Matrix[m] = (float *)ckalloc(N*sizeof(float));

  return(Matrix);
}

void FreeFloatMatrix(float **Matrix, int M)
{
  int m;

  for( m=0; m<M; m++ ) free(Matrix[m]);

  free(Matrix);

}

float ***FloatCube(int M, int N, int K)
{
  int m, n, k;
  float ***Cube;

  Cube = (float ***)ckalloc(M*sizeof(float **));

  for( m=0; m<M; m++ ) {
    Cube[m] = (float **)ckalloc(N*sizeof(float *));
    for( n=0; n<N; n++ )
      Cube[m][n] = (float *)ckalloc(K*sizeof(float));
  }

  for( m=0; m<M; m++ )
    for( n=0; n<N; n++ )
      for( k=0; k<K; k++ )
	Cube[m][n][k] = 0.0;

  return(Cube);
}

void FreeFloatCube(float ***Cube, int M, int N)
{
  int m, n;

  for( m=0; m<M; m++ ) {
    for( n=0; n<N; n++ )
      free(Cube[m][n]);
    free(Cube[m]);
  }

  free(Cube);
}

char **CharMatrix(int M, int N)
{
  int m;
  char **Matrix;

  Matrix = (char **)ckalloc(M*sizeof(char *));

  for( m=0; m<M; m++ ) Matrix[m] = (char *)ckalloc(N*sizeof(char));

  return(Matrix);
}

void FreeCharMatrix(char **Matrix, int M)
{
  int m;

  for( m=0; m<M; m++ ) free(Matrix[m]);

  free(Matrix);

}


int **IntMatrix(int M, int N)
{
  int m;
  int **Matrix;

  Matrix = (int **)ckalloc(M*sizeof(int *));

  for( m=0; m<M; m++ ) Matrix[m] = (int *)ckalloc(N*sizeof(int));

  return(Matrix);
}

void FreeIntMatrix(int **Matrix, int M)
{
  int m;

  for( m=0; m<M; m++ ) free(Matrix[m]);

  free(Matrix);

}

void MainAlloc(COMMAND **Cmd, SEQ ***Seq)
{
  Bridge1 = FloatMatrix(NAcd,NAcd);
  Bridge2 = FloatMatrix(NAcd,NAcd);
  Bridge3 = FloatMatrix(NAcd,NAcd);
  Alpha   = FloatMatrix(NAcd,NAcd);
  MbrMat  = FloatCube(MBRWIN,NAcd,NAcd);
  (*Cmd)  = (COMMAND *)ckalloc(sizeof(COMMAND));
  (*Seq)  = (SEQ **)ckalloc(MAX_SEQ*sizeof(SEQ *));
  Rel = (float **)ckalloc(MAXREL*sizeof(float *));
}

void MainDealloc(COMMAND **Cmd, SEQ ***Seq)
{
  register int i;
  
  free(*Cmd);

  free((*Seq));

  FreeFloatMatrix(Bridge1,NAcd);
  FreeFloatMatrix(Bridge2,NAcd);
  FreeFloatMatrix(Bridge3,NAcd);
  FreeFloatMatrix(Alpha,NAcd);
  FreeFloatCube(MbrMat,MBRWIN,NAcd);

  for( i=0; i<NStr; i++ ) {
    free(SeqDb[i]);
    free(StrDb[i]);
    free(NmDb[i]);
  }
  for( i=0; i<NRel; i++ )
    free(Rel[i]);
  free(Rel);
}

void FreeSeq(SEQ **Seq)
{
  register int i;
  
  if( (*Seq)->Id != NULL )
    free((*Seq)->Id);
  if( (*Seq)->De != NULL )
    free((*Seq)->De);
  if( (*Seq)->Se != NULL )
    free((*Seq)->Se);
  if( (*Seq)->AliSe != NULL )
    free((*Seq)->AliSe);
  if( (*Seq)->Pred != NULL )
    free((*Seq)->Pred);
  if( (*Seq)->Asn != NULL )
    free((*Seq)->Asn);
  if( (*Seq)->KnownAsn != NULL )
    free((*Seq)->KnownAsn);
  if( (*Seq)->AliCnt != NULL )
    free((*Seq)->AliCnt);
  if( (*Seq)->Rel != NULL )
    free((*Seq)->Rel);
  if( (*Seq)->PHel != NULL )
    free((*Seq)->PHel);
  if( (*Seq)->PBeta != NULL )
    free((*Seq)->PBeta);
  if( (*Seq)->PCoil != NULL )
    free((*Seq)->PCoil);
  FreeFloatMatrix((*Seq)->Prop,NPROP);
  for( i=0; i<(*Seq)->NHml; i++ )
    free((*Seq)->Hml[i]);
  free((*Seq));
}
