#define YES                       1
#define NO                        0
#include "prd.h"
#include "prdprot.h"

extern float **Bridge1, **Bridge2, **Bridge3;
extern float **Alpha, ***MbrMat, ThresholdsM[NTHRESH], ThresholdsS[NTHRESH];
extern int    *SeqDb[MAXDB], *StrDb[MAXDB], NStr;
extern char  *NmDb[MAXDB];

extern float **Rel;
extern int NRel;

extern char   *AMINOACIDS;
extern char   *ALNSYMBOLS;
extern char   *STRUCTURES;

extern bool FirstCall;
extern int MaxNameLength;

bool Hml[MAXDB];

void PredictSS(SEQ *Seq, int From, int To, char *Name, COMMAND *Cmd)
{

  int *ISeq, Length, NameLength;
  register int i, j;

  if( Cmd->Info ) {

    FirstCall = 0;
    
    if( FirstCall ) {

      printf("\n");
      for( i=0; i<MaxNameLength+3; i++ )
        printf(" ");
      printf("0%% ");

      for( i=4; i<NMARKS-4; i++ )
        printf(".");
      printf(" 100%%\n");
      fflush(stdout);
      FirstCall = FALSE;
    }

    NameLength = (int)strlen(Name);

    printf(" Info>  ");
    for( i=0; i<MaxNameLength; i++ ) {
      if( i<NameLength )
        printf("%c",Name[i]);
      else
        printf(" ");
    }
    printf(" : ");
    fflush(stdout);
  }
  
  Length = To-From;

  ISeq   = (int *)ckalloc(Seq->Len*sizeof(int));
  Char2Int(Seq->Se,ISeq,Seq->Len);

  for( i=0; i<Seq->Len; i++ ) {
    Seq->KnownAsn[i] = '?';
    for( j=0; j<NPROP; j++ )
      Seq->Prop[j][i] = 0.0;
  }

  PredictAnti(ISeq,Seq->Len,Seq->Prop[0]);
  PredictPar(ISeq,Seq->Len,Seq->Prop[1]);
  PredictAlp(ISeq,Seq->Len,Seq->Prop[2]);
  PredictNN(ISeq+From,Length,Seq->Prop[3]+From,Seq->Prop[4]+From,Seq->Prop[5]+From,Seq->KnownAsn+From,Cmd);
  PredictTurn(ISeq+From,Length,Seq->Prop[6]+From);
  free(ISeq);

}

void StoreHomol(SEQ **Seq)
{

  register int i;

  for( i=0; i<NStr; i++ )
    if( Hml[i] ) {
      (*Seq)->Hml[(*Seq)->NHml] = (char *)ckalloc(((int)strlen(NmDb[i])+1)*sizeof(char));
      strcpy((*Seq)->Hml[(*Seq)->NHml],NmDb[i]);
      (*Seq)->NHml++;
    }
}

void PredictAnti(int *ISeq, int SeqLen, float *Anti)
{
  float Sum1, Sum2, Mx, Max;
  register int k, i, ii, j, ind1, ind2;

  for( i=0; i<SeqLen; i++ )
    Anti[i] = 0.0;

  for( i=WINBETA-1; i<SeqLen; i++ ) {

    Max = 0.0;

    for( j=i+ANTIDIST; j<SeqLen-WINBETA; j++ ) {

      for( k=0; k<WINBETA; k++ ) {
        
        ind1 = ISeq[i-k];
	ind2 = ISeq[j+k];
        
	if( (k%2) == 0 ) {
	  if( k == 0 ) {
	    Sum1 = Bridge1[ind1][ind2];
	    Sum2 = Bridge2[ind1][ind2];
	  }
	  else {
	    Sum1 += Bridge1[ind1][ind2];
	    Sum2 += Bridge2[ind1][ind2];
	  }
	}
	else {
	  Sum1 += Bridge2[ind1][ind2];
	  Sum2 += Bridge1[ind1][ind2];
	}
	
      }
      Mx = Maximum(Sum1,Sum2)/(float)WINBETA;
      if( Max < Mx )
        Max = Mx;
      
    }

    for( ii=WINBETA-1; ii<=i-WINBETA-ANTIDIST; ii++ ) {
      j = i-ANTIDIST;
      for( k=0; k<WINBETA; k++ ) {
        
        ind1 = ISeq[ii-k];
	ind2 = ISeq[j+k];
        
	if( (k%2) == 0 ) {
	  if( k == 0 ) {
	    Sum1 = Bridge1[ind1][ind2];
	    Sum2 = Bridge2[ind1][ind2];
	  }
	  else {
	    Sum1 += Bridge1[ind1][ind2];
	    Sum2 += Bridge2[ind1][ind2];
	  }
	}
	else {
	  Sum1 += Bridge2[ind1][ind2];
	  Sum2 += Bridge1[ind1][ind2];
	}
	
      }
      Mx = Maximum(Sum1,Sum2)/(float)WINBETA;
      if( Max < Mx )
        Max = Mx;
    }

    if( i >= WINBETA )
      Anti[i-(WINBETA)] = Max;
    
  }
  
  TakeBestBeta(Anti,SeqLen);
}


void PredictPar(int *ISeq, int SeqLen, float *Par)
{
  float Sum1, Sum2, Mx;
  register int k, i, ii, j, ind1, ind2;
  int OffSet;

  for( i=0; i<SeqLen; i++ )
    Par[i] = 0.0;

  OffSet = WINBETA + PARDIST;

  for( i=0; i<SeqLen-WINBETA; i++ ) {

    for( j=i+OffSet; j<SeqLen-WINBETA; j++ ) {
      
      for( k=0; k<WINBETA; k++ ) {
        
	ind1 = ISeq[i+k];
	ind2 = ISeq[j+k];

	if( (k%2) == 0 ) {
	  if( k == 0 ) {
	    Sum1 = Bridge3[ind1][ind2];
	    Sum2 = Bridge3[ind2][ind1];
	  }
	  else {
	    Sum1 += Bridge3[ind1][ind2];
	    Sum2 += Bridge3[ind2][ind1];
	  }
	}
	else {
	  Sum1 += Bridge3[ind2][ind1];
	  Sum2 += Bridge3[ind1][ind2];
	}
	
      }
      Mx = Maximum(Sum1,Sum2)/(float)WINBETA;
      if( Par[i] < Mx )
        Par[i] = Mx;
    }

    for( ii=0; ii<=i-OffSet; ii++ ) {
      j = i;
      for( k=0; k<WINBETA; k++ ) {
        
	ind1 = ISeq[ii+k];
	ind2 = ISeq[j+k];

	if( (k%2) == 0 ) {
	  if( k == 0 ) {
	    Sum1 = Bridge3[ind1][ind2];
	    Sum2 = Bridge3[ind2][ind1];
	  }
	  else {
	    Sum1 += Bridge3[ind1][ind2];
	    Sum2 += Bridge3[ind2][ind1];
	  }
	}
	else {
	  Sum1 += Bridge3[ind2][ind1];
	  Sum2 += Bridge3[ind1][ind2];
	}
	
      }
      Mx = Maximum(Sum1,Sum2)/(float)WINBETA;
      if( Par[i] < Mx )
        Par[i] = Mx;
    }
  }
  
  TakeBestBeta(Par,SeqLen);

}

void PredictAlp(int *ISeq, int SeqLen, float *Alp)
{

  register int i, j;

  for( i=0; i<SeqLen-WINALPHA-4; i++ ) {
    Alp[i] = 0.0;
    for( j=0; j<WINALPHA; j++ )
      Alp[i] += Alpha[ISeq[i+j]][ISeq[i+j+4]];
    
    Alp[i] /= (float)WINALPHA;
    
  }
  TakeBestAlpha(Alp,SeqLen);
}

void Derive(float **Mat, int N, float *Pred)
{
  
  register int i, j;

  for( j=0; j<N; j++ ) {
    Pred[j] = -100.0;
    for( i=1; i<N; i++ )
      if( Pred[j] < Mat[i][j] )
 	Pred[j] = Mat[i][j];
  }
}

void TakeBestBeta(float *Pred, int N)
{

  register int i, k;
  float *Out;

  Out = (float *)ckalloc(N*sizeof(float));

  for( i=0; i<N; i++ )
    Out[i] = 0.0;

  for( i=0; i<N-WINBETA; i++ ) {

    for( k=i+1; k<i+WINBETA-1 && k<N; k++ )
      if( Pred[i] > Out[k] )
	Out[k] = Pred[i];
  }

  for( i=0; i<N; i++ )
    Pred[i] = Out[i];

  free(Out);
}

void TakeBestAlpha(float *Pred, int N)
{

  register int i, k;
  float *Out;

  Out = (float *)ckalloc(N*sizeof(float));

  for( i=0; i<N; i++ )
    Out[i] = 0.0;

  for( i=0; i<N-WINALPHA-4; i++ ) {

    for( k=i+2; k<i+WINALPHA+2 && k<N; k++ )
      if( Pred[i] > Out[k] )
	Out[k] = Pred[i];
  }

  for( i=0; i<N; i++ )
    Pred[i] = Out[i];

  free(Out);
}
	

void MakeSymAnti(float **Matrix, int N)
{

  register int i, j;
  
  for( i=WINBETA; i<N-WINBETA; i++ )
    for( j=i+ANTIDIST; j<N-WINBETA; j++ )
      Matrix[j+WINBETA][i-WINBETA] = Matrix[i][j];
}
  
void MakeSymPar(float **Matrix, int N)
{

  register int i, j;
  int OffSet;

  OffSet = WINBETA+PARDIST;
  
  for( i=0; i<N-OffSet-WINBETA; i++ )
    for( j=i+OffSet; j<N-WINBETA; j++ )
      Matrix[j][i] = Matrix[i][j];
}


void PredictNN(int *ISeq, int SeqLen, float *AlpNN, float *BetaNN, float *CoilNN, char *KnownAsn, COMMAND *Cmd)
{
  register int i, ii, j, k, PdbCount;
  float Score[MBRNBEST], CurrScore, CurrMax;
  int Structure[MBRNBEST], Position[MBRNBEST], Homol[MAXDB];
  int CurrPred, IndMax, CurrLen, *p1, *p2, IPr, NDots = 1, NPrinted = 0;
  float PrintStep;
  float *f0, *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9, *f10, *f11, *f12;
      
  for( PdbCount=0; PdbCount<NStr; PdbCount++ ) {
    Homol[PdbCount] = Homologous(ISeq,SeqLen,SeqDb[PdbCount]+1,SeqDb[PdbCount][0]);
    Hml[PdbCount] = FALSE;
  }

  PrintStep = (float)(SeqLen - MBRWIN)/(float)NMARKS;
  if( PrintStep < 0.5 ) {
    NDots = 1.0/PrintStep;
    PrintStep = 1.0;
  }
  IPr = 0;
  
  for( i=0; i<SeqLen-MBRWIN; i++ ) {

    if( Cmd->Info && NPrinted < NMARKS-1 && (i == 0 || ( i!= 0 ) && i - IPr >= (int)PrintStep ) ) {
      for( j=0; j<NDots; j++ ) {
        printf(".");
        NPrinted++;
      }        
      fflush(stdout);
      IPr = i;
    }
    
    p1 = ISeq+i;
    
    f0 = MbrMat[0][p1[0]];
    f1 = MbrMat[1][p1[1]];
    f2 = MbrMat[2][p1[2]];
    f3 = MbrMat[3][p1[3]];
    f4 = MbrMat[4][p1[4]];
    f5 = MbrMat[5][p1[5]];
    f6 = MbrMat[6][p1[6]];
    f7 = MbrMat[7][p1[7]];
    f8 = MbrMat[8][p1[8]];
    f9 = MbrMat[9][p1[9]];
    f10 = MbrMat[10][p1[10]];
    f11 = MbrMat[11][p1[11]];
    f12 = MbrMat[12][p1[12]];
    
    CurrPred = i+MBRWIN22;
    
    for( j=0; j<MBRNBEST; j++ ) {
      Score[j] = 1.5;
      Structure[j] = 0;
      Position[j] = j;
    }
    
    CurrMax = 1.5;
    IndMax = 0;
    
    for( PdbCount=0; PdbCount<NStr; PdbCount++ ) {
      
      if( Homol[PdbCount] )
	continue;

      p2 = SeqDb[PdbCount]+1;
      
      CurrLen = SeqDb[PdbCount][0]-MBRWIN+1;

      for( j=1; j<CurrLen; j++ ) {
        
        if( (CurrScore = f0[p2[0]] + f1[p2[1]] + f2[p2[2]] + f3[p2[3]] + f4[p2[4]] +
             f5[p2[5]] + f6[p2[6]] + f7[p2[7]] + f8[p2[8]]) < CurrMax &&
            (CurrScore += f9[p2[9]]) < CurrMax &&
            (CurrScore += f10[p2[10]]) < CurrMax &&
            (CurrScore += f11[p2[11]]) < CurrMax &&
            (CurrScore += f12[p2[12]]) < CurrMax ) {
         
          Score[IndMax] = CurrScore;
          Structure[IndMax] = PdbCount;
          Position[IndMax] = j;
          
          for( IndMax = 0, CurrMax = Score[0], ii=1; ii<MBRNBEST; ii++ )
            if( Score[ii] > CurrMax ) {
              CurrMax = Score[ii];
              IndMax = ii;
            }
        }
       p2++;
      }
      
    }

    Reason(Score,Structure,Position,&AlpNN[CurrPred],&BetaNN[CurrPred],&CoilNN[CurrPred]);
    
  }

  if( Cmd->Info ) {
    for( i=NPrinted; i<NMARKS-1; i++ )
      printf(".");
    printf("!\n");
    fflush(stdout);
  }

  for( PdbCount=0; PdbCount<NStr; PdbCount++ )
    if( Homol[PdbCount] )
      Hml[PdbCount] = LookUp(KnownAsn,ISeq,SeqLen,SeqDb[PdbCount]+1,StrDb[PdbCount]+1,SeqDb[PdbCount][0]);

}

void Reason(float *Score, int *Structure, int *Position,float *Score_H, float *Score_E, float *Score_C)
{

  register int i;
  int Count_H, Count_E, Count_C;
  int CurrStr;

  Count_H = Count_E = Count_C = 0;
  *Score_H = *Score_E = *Score_C = 0.0;
  
  for( i=0; i<MBRNBEST; i++ ) {

    if( Structure[i] == ERR )
      die("\n\nThe database you are using is too small. Please contact Dmitrij Frishman (FRISHMAN@EMBL-heidelberg.de)\n");
    
    CurrStr = StrDb[Structure[i]][Position[i]+MBRWIN22];

    switch(CurrStr) {
    case 0: 
      Count_H++;
      *Score_H += Score[i];
      break;
      
    case 1: 
      Count_E++;
      *Score_E += Score[i];
      break;
      
    case 2: 
      Count_C++;
      *Score_C += Score[i];
      break;
    }
  }

  if( *Score_H > Eps ) {
    *Score_H /= (float)Count_H;
    *Score_H = ((float)Count_H) / (((float)MBRNBEST) * (*Score_H)) ;
  }

  if( *Score_E > Eps ) {
    *Score_E /= (float)Count_E;
    *Score_E = ((float)Count_E) / (((float)MBRNBEST) * (*Score_E)) ;
  }

  if( *Score_C > Eps ) {
    *Score_C /= (float)Count_C;
    *Score_C = ((float)Count_C) / (((float)MBRNBEST) * (*Score_C)) ;
  }
}
      

bool Homologous(int *ISeq, int SeqLen, int *IDbSeq, int DbSeqLen)
{

  register int i, j, k;
  int win = 6, TotalCount = 0, Ln;

  Ln = DbSeqLen-win;

  for( i=0; i<SeqLen-win; i++ ) {
    for( j=0; j<Ln; j++ ) {
      for( k=0; k<win; k++ )
	if( ISeq[i+k] != IDbSeq[j+k] )
	  break;
      if( k == win ) {
	TotalCount++;
	break;
      }
    }
  }

  if( TotalCount > 12 )
    return(TRUE);

  return(FALSE);
}

void LnMat(float **Matrix)
{
  
  register int i, j;
  
  for( i=0; i<NAcd; i++ )
    for( j=0; j<NAcd; j++ ) {
      if( Matrix[i][j] < 0.0005 )
	Matrix[i][j] = 0.0005;
      Matrix[i][j] = log(Matrix[i][j]);
    }
}

void NormalizeMatrixArtif(float **Matrix, float Min, float Max)
{

  register int i, j;

  for( i=0; i<NAcd; i++ )
    for( j=0; j<NAcd; j++ ) {
      if( Matrix[i][j] > Max )
	Matrix[i][j] = Max;
		if( Matrix[i][j] < Min )
	Matrix[i][j] = Min;
      Matrix[i][j] = (Matrix[i][j]-Min)/(Max-Min);
    }
}


void Decide(SEQ *Seq, COMMAND *Cmd)
{
  register int i, j;
  float Sum, *Thresholds;

  if( Cmd->Single )
    Thresholds = ThresholdsS;
  else
    Thresholds = ThresholdsM;

  for( i=0; i<Seq->Len; i++ ) {
    
    if( (  Seq->AliCnt[i] && Seq->Prop[2][i] >= Thresholds[2]) ||
        ( !Seq->AliCnt[i] && Seq->Prop[2][i] >= Thresholds[9]) ) 
      Seq->Pred[i] = 'h';
    else
    if( (  Seq->AliCnt[i] &&
          ( Seq->Prop[0][i] >= Thresholds[0] || Seq->Prop[1][i] >= Thresholds[1]) ) ||
        ( !Seq->AliCnt[i] &&
          ( Seq->Prop[0][i] >= Thresholds[7] || Seq->Prop[1][i] >= Thresholds[8]) ) )
      Seq->Pred[i] = 'e';
    else
      Seq->Pred[i] = 'c';
  }

  for( i=7; i<Seq->Len-8-MINHEL+1; i++ ) {
      Sum = 0.0;
      for( j=0; j<MINHEL; j++ )
	Sum += Seq->Prop[3][i+j];
      Sum /= (float)MINHEL;
      if( Seq->AliCnt[i] && Sum > Thresholds[4] || !Seq->AliCnt[i] && Sum > Thresholds[11] ) {
	for( j=1; j<MINHEL-1; j++ )
	  Seq->Pred[i+j] = 'h';
      }
    }
  
    for( i=7; i<Seq->Len-8-MINBET+1; i++ ) {
      Sum = 0.0;
      for( j=0; j<MINBET; j++ )
	Sum += Seq->Prop[4][i+j];
      Sum /= (float)MINBET;
      if( Seq->AliCnt[i] && Sum > Thresholds[5] || !Seq->AliCnt[i] && Sum > Thresholds[12] ) {
	for( j=1; j<MINBET-1; j++ )
	  Seq->Pred[i+j] = 'e';
      }
    }
	  
    for( i=7; i<Seq->Len-8-MINCOIL+1; i++ ) {
      Sum = 0.0;
      for( j=0; j<MINCOIL; j++ )
	Sum += Seq->Prop[5][i+j];
      Sum /= (float)MINCOIL;
      if( Seq->AliCnt[i] && Sum > Thresholds[6] || !Seq->AliCnt[i] && Sum > Thresholds[13] ) {
	for( j=1; j<MINCOIL-1; j++ )
	  Seq->Pred[i+j] = 'c';
      }
    }
  
  for( i=0; i<Seq->Len; i++ )
    if( Seq->AliCnt[i] && 
        Seq->Prop[6][i+0] > Thresholds[3] && 
        Seq->Prop[6][i+1] > Thresholds[3] && 
        Seq->Prop[6][i+2] > Thresholds[3] && 
        Seq->Prop[6][i+3] > Thresholds[3]   ||
        !Seq->AliCnt[i] && 
        Seq->Prop[6][i+0] > Thresholds[10] && 
        Seq->Prop[6][i+1] > Thresholds[10] && 
        Seq->Prop[6][i+2] > Thresholds[10] && 
        Seq->Prop[6][i+3] > Thresholds[10] ) {
      Seq->Pred[i+1] = 'c';
      Seq->Pred[i+2] = 'c';
    }
  

  Seq->Pred[0] = Seq->Pred[1] = Seq->Pred[2] = Seq->Pred[Seq->Len-1] =
    Seq->Pred[Seq->Len-2] = Seq->Pred[Seq->Len-3] = 'c';
  
  CorrectAsn(Seq->Pred,Seq->Len,'e','c',2);
  CorrectAsn(Seq->Pred,Seq->Len,'h','c',4);

  if( Cmd->LookUp ) {
    for( i=0; i<Seq->Len; i++ )
      if( Seq->KnownAsn[i] != '?')
        Seq->Pred[i] = Seq->KnownAsn[i];
  }
}

void PredictTurn(int *ISeq, int SeqLen, float *OutTurn)
{

  static float P[20][4] = {
    { 0.81, 0.96, 0.66, 0.89 },
    { 0.69, 0.93, 0.75, 0.93 },
    { 1.54, 1.02, 2.14, 1.06 },
    { 1.56, 1.24, 1.86, 0.99 },
    { 1.42, 0.73, 0.98, 1.20 },
    { 0.89, 0.94, 0.93, 1.01 },
    { 0.87, 1.35, 0.92, 0.89 },
    { 1.09, 1.04, 2.14, 1.64 },
    { 1.25, 0.95, 1.16, 0.93 },
    { 0.66, 0.61, 0.42, 0.68 },
    { 0.73, 0.67, 0.47, 0.78 },
    { 0.80, 1.22, 0.94, 1.10 },
    { 0.70, 0.48, 0.41, 0.68 },
    { 0.98, 0.66, 0.96, 0.95 },
    { 1.48, 2.45, 0.63, 0.96 },
    { 1.29, 1.23, 1.06, 1.03 },
    { 1.08, 0.79, 0.94, 1.20 },
    { 0.62, 0.65, 0.76, 0.79 },
    { 1.04, 0.75, 0.83, 1.07 },
    { 0.72, 0.70, 0.54, 0.84 } };

  register int i, k;
  int AACode1, AACode2, AACode3, AACode4;
  float *Pred;

  Pred = (float *)ckalloc(SeqLen*sizeof(float));

  for( i=0; i<SeqLen-4; i++ ) {
    if( (AACode1 = ISeq[i])   < 20 &&
        (AACode2 = ISeq[i+1]) < 20 &&
        (AACode3 = ISeq[i+2]) < 20 &&
        (AACode4 = ISeq[i+3]) < 20 )
      Pred[i] = ( P[AACode1][0] + P[AACode2][1] + P[AACode3][2] + P[AACode4][3] )/4.0;
    else
      Pred[i] = 0.0;
  }
  
  for( i=0; i<SeqLen; i++ )
    OutTurn[i] = 0.0;

  for( i=0; i<SeqLen-4; i++ )
    for( k=0; k<4; k++ )
      if( Pred[i] > OutTurn[i+k] )
	OutTurn[i+k] = Pred[i];

  free(Pred);
}

void Normalize(void)
{

  LnMat(Bridge1);
  LnMat(Bridge2);
  LnMat(Bridge3);
  LnMat(Alpha);
  
  NormalizeMatrixArtif(Bridge1,-8.0,-2.5);
  NormalizeMatrixArtif(Bridge2,-8.0,-2.5);
  NormalizeMatrixArtif(Bridge3,-8.0,-2.5);
  NormalizeMatrixArtif(Alpha,-8.0,0.0);
}

void PrintMatrix(float **Matrix, int N, int M, char *Text)
{
  register int i, j;

  for( i=0; i<N; i++ ) {
    printf("%s ",Text);
    for( j=0; j<M; j++ )
      printf("%7.3f ",Matrix[i][j]) ;
    printf("\n");
  }
}
  

bool LookUp(char *KnownAsn, int *ISeq, int SeqLen, int *IDbSeq, int *IDbStr, int DbSeqLen)
{

  register int i, j, k;
  int win = 7;
  bool Used = FALSE;

  for( i=0; i<SeqLen-win+1; i++ ) {
    for( j=0; j<DbSeqLen-win; j++ ) {
      for( k=0; k<win; k++ )
	if( ISeq[i+k] != IDbSeq[j+k] )
	  break;
      if( k == win ) {
        Used = TRUE;
        for( k=0; k<win; k++ )
          KnownAsn[i+k] = STRUCTURES[IDbStr[j+k]];
	break;
      }
    }
  }

  return(Used);
}

void Assess(SEQ *Seq, COMMAND *Cmd)
{
  char *KnownAsn, *Pred;
  register int i;

  if( strcmp(Cmd->StrideFile,"") || strcmp(Cmd->DsspFile,"") ) {
    
    KnownAsn = (char *)ckalloc(Seq->Len*sizeof(char));
    Pred     = (char *)ckalloc(Seq->Len*sizeof(char));
    
    for( i=0; i<Seq->Len; i++ ) {
      KnownAsn[i] = toupper(Seq->Asn[i]);
      if( KnownAsn[i] != 'H' && KnownAsn[i] != 'E' )
        KnownAsn[i] = 'C';
      Pred[i]     = toupper(Seq->Pred[i]);
    }    
    
    CorrectAsn(KnownAsn,Seq->Len,'E','C',2);
    CorrectAsn(KnownAsn,Seq->Len,'H','C',4);
    
    Seq->PercCorr = PercentCorrect(Pred,KnownAsn,Seq->Len);
    
    free(KnownAsn);
    free(Pred);
  }
  else
    Seq->PercCorr = -1.0;
}

    

float PercentCorrect(char *TestAsn, char *KnownAsn, int Length)
{
  int Res, Count=0;;

  for( Res=0; Res<Length; Res++ )
    if( KnownAsn[Res] == TestAsn[Res] )
      Count++;

  return( (float)((float)Count/(float)Length) );
}


void ProbAndRel(SEQ *Seq)
{
  register int i, j;

  for( i=0; i<Seq->Len; i++ ) {
    for( j=0; j<NRel; j++ ) {
      if( Seq->Prop[0][i] >= Rel[j][4]  && Seq->Prop[0][i] <= Rel[j][5]  &&
          Seq->Prop[1][i] >= Rel[j][6]  && Seq->Prop[1][i] <= Rel[j][7]  &&
          Seq->Prop[2][i] >= Rel[j][8]  && Seq->Prop[2][i] <= Rel[j][9]  &&
          Seq->Prop[6][i] >= Rel[j][10] && Seq->Prop[6][i] <= Rel[j][11] &&
          Seq->Prop[3][i] >= Rel[j][12] && Seq->Prop[3][i] <= Rel[j][13] &&
          Seq->Prop[4][i] >= Rel[j][14] && Seq->Prop[4][i] <= Rel[j][15] &&
          Seq->Prop[5][i] >= Rel[j][16] && Seq->Prop[5][i] <= Rel[j][17] ) {
        Seq->Rel[i]   = Rel[j][0];
        Seq->PHel[i]  = Rel[j][1];
        Seq->PBeta[i] = Rel[j][2];
        Seq->PCoil[i] = Rel[j][3];
        break;
      }
    }
  }
}

    

      
