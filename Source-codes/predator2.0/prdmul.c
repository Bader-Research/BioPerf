#include "prd.h"
#include "prdprot.h"

void Multiple(SEQ **Seq, int NSeq, COMMAND *Cmd)
{

  ALN   ***Aln;
  PROJ  **Proj;
  register int i, j, k;

  Aln = (ALN ***)ckalloc(MAX_SEQ*sizeof(ALN **));
    
  for( i=0; i<NSeq; i++ ) {
    
    if( !Seq[i]->ToBePredicted )
      continue;
    
    if( Cmd->Info && Cmd->All ) {
      printf(" Info>  Generating the final prediction for sequence %s",Seq[i]->Id);
      fflush(stdout);
    }

    Aln[i] = (ALN **)ckalloc(NSeq*sizeof(ALN *));
    AllocProj(&Proj,Seq[i]->Len,NSeq);
      
    if( !Cmd->Single ) {
      
      for( j=0; j<NSeq; j++ )
        
        if( i != j ) {
          
          AllocAln(&Aln[i][j],Seq[i]->Len,Seq[j]->Len,Cmd->NAli);
          
          if( Cmd->OrigAli == TRUE &&
              (Seq[0]->FileType == MSF ||
               Seq[0]->FileType == Clustal) ) 
            Borrow(Seq[i],Seq[j],Aln[i][j]); /* Not implemented! */
          else
            Align(Seq[i]->Se,Seq[i]->Len,Seq[j]->Se,Seq[j]->Len,Aln[i][j],
                  Cmd->DelVal,Cmd->GapVal,Cmd->NAli,Cmd->MinAliLength,Cmd->Safety);
          
          for( k=0; k<Aln[i][j]->NAli; k++ ) {
            if( Aln[i][j]->Accepted[k] ) {
              
/*
          if( Cmd->Verbose )
            PrintAln(Aln[i][j],i,j,k);
                  */
              
              Project(Proj,Aln[i][j],Seq[i],Seq[j],k,Cmd);
            }
          }
          DeallocAln(&Aln[i][j],Cmd->NAli);
        }
/* if( Cmd->Verbose )
   PrintProp(Seq,i,Proj,Cmd); */
    }
    
    if( Cmd->Info && Cmd->All )
      printf("\n");
    
    Combine(Seq,i,Proj,Cmd);
    
    free(Aln[i]);
    DeallocProj(&Proj,Seq[i]->Len);

  }

  free(Aln);
  
}

void SelectSequences(SEQ **Seq, int NSeq, COMMAND *Cmd)
{
  register int i;

  if( NSeq == 1 ) {
    if( Seq[0]->Len < 30 ) 
      die("Sequence %s is too short\n",Seq[0]->Id);
    Seq[0]->ToBePredicted = TRUE;
    return;
  }
    
  if( Cmd->All ) {
    for( i=0; i<NSeq; i++ ) {
      if( Seq[i]->Len < 30 ) 
	die("Sequence %s is too short\n",Seq[i]->Id);
      Seq[i]->ToBePredicted = TRUE;
    }
    return;
  }

  if( strlen(Cmd->PredId) != 0 ) {
    for( i=0; i<NSeq; i++ )
      if( !strcmp(Cmd->PredId,Seq[i]->Id) ) {
	if( Seq[i]->Len < 30 ) 
	  die("Sequence %s is too short\n",Seq[i]->Id);
	Seq[i]->ToBePredicted = TRUE;
	break;
      }
	 if( i == NSeq )
      die("Sequence %s not found in file %s\n",Cmd->PredId,Cmd->InputFile);
    else
      return;
  }

  if( Seq[0]->Len < 30 ) 
    die("Sequence %s is too short\n",Seq[0]->Id);
  Seq[0]->ToBePredicted = TRUE;
  
  return;
}


void AllocAln(ALN **Aln, int Len1, int Len2, int NAli)
{

  int MaxAliLen;

  MaxAliLen = 10*Maximum(Len1,Len2);

  *Aln             = (ALN *)    ckalloc(sizeof(ALN));
  (*Aln)->Nc       = (int *)    ckalloc(NAli*sizeof(int));
  (*Aln)->Score    = (int *)    ckalloc(NAli*sizeof(int));
  (*Aln)->HCri     = (float *)  ckalloc(NAli*sizeof(float));
  (*Aln)->Len      = (int *)    ckalloc(NAli*sizeof(int));
  (*Aln)->Beg1     = (int *)    ckalloc(NAli*sizeof(int));
  (*Aln)->Beg2     = (int *)    ckalloc(NAli*sizeof(int));
  (*Aln)->End1     = (int *)    ckalloc(NAli*sizeof(int));
  (*Aln)->End2     = (int *)    ckalloc(NAli*sizeof(int));
  (*Aln)->Percent  = (float *)  ckalloc(NAli*sizeof(float));
  (*Aln)->Accepted = (int *)    ckalloc(NAli*sizeof(int));

  (*Aln)->A1       = CharMatrix(NAli,MaxAliLen);
  (*Aln)->A2       = CharMatrix(NAli,MaxAliLen);
  (*Aln)->NAli     = 0;
}

void DeallocAln(ALN **Aln, int NAli)
{

  free((*Aln)->Nc);
  free((*Aln)->Score);
  free((*Aln)->HCri);
  free((*Aln)->Len);
  free((*Aln)->Beg1);
  free((*Aln)->Beg2);
  free((*Aln)->End1);
  free((*Aln)->End2);
  free((*Aln)->Percent);
  free((*Aln)->Accepted);

  FreeCharMatrix((*Aln)->A1,NAli);
  FreeCharMatrix((*Aln)->A2,NAli);
  
  free(*Aln);
}

void PrintAln(ALN *Aln, int Seq1, int Seq2, int AlnN)
{
  register int m;

  printf("ALN %3d <-> %3d %3d %7.5f %7.5f %d %d\n",
	 Seq1,Seq2,AlnN,Aln->Percent[AlnN],Inf(Aln->Percent[AlnN]/100.0),Aln->Nc[AlnN],Aln->Score[AlnN]);
  printf("ALN %3d <-> %3d %3d %-5d ",Seq1,Seq2,AlnN,Aln->Beg1[AlnN]+1);
  for( m=0; m<Aln->Len[AlnN]; m++ )
    printf("%c",Aln->A1[AlnN][m]);
  printf("%5d\n",Aln->End1[AlnN]+1);
  printf("ALN %3d <-> %3d %3d %-5d ",Seq1,Seq2,AlnN,Aln->Beg2[AlnN]+1);
  for( m=0; m<Aln->Len[AlnN]; m++ )
    printf("%c",Aln->A2[AlnN][m]);
  printf("%5d\n\n",Aln->End2[AlnN]+1);
}

void Project(PROJ **Proj, ALN *Aln, SEQ *Seq1, SEQ *Seq2, int AlnN, COMMAND *Cmd)
{
  int Pos1, Pos2, From, To;
  register int m;


  for( m=0; m<Seq1->Len; m++ )
    Proj[m]->ResN[Proj[m]->NAli] = ERR;
  
  Pos1 = Aln->Beg1[AlnN];
  Pos2 = Aln->Beg2[AlnN];
  
  for( m=0; m<Aln->Len[AlnN]; m++ ) {

    if( Aln->A1[AlnN][m] == '-' )
      Pos2++;
    else
    if( Aln->A2[AlnN][m] == '-' )
      Pos1++;
    else {
      Proj[Pos1]->Score[Proj[Pos1]->NAli] = (float)Aln->Score[AlnN];
      Proj[Pos1]->HCri[Proj[Pos1]->NAli] = Aln->HCri[AlnN];
      Proj[Pos1]->AliLength[Proj[Pos1]->NAli] = Aln->Len[AlnN];
      Proj[Pos1]->Percent[Proj[Pos1]->NAli] = Aln->Percent[AlnN];
      Proj[Pos1]->Seq[Proj[Pos1]->NAli] = Seq2;
      Proj[Pos1]->ResN[Proj[Pos1]->NAli] = Pos2;
      Proj[Pos1]->NAli++;
      Pos1++;
      Pos2++;
    }
  }

  if( Seq2->ToBePredicted == FALSE ) {
    From = Aln->Beg2[AlnN]-MBRWIN;
    if( From < 0 )
      From = 0;
    To   = Aln->End2[AlnN]+MBRWIN;
    if( To > Seq2->Len )
      To = Seq2->Len;
    
    PredictSS(Seq2,From,To,Seq2->Id,Cmd);
  }  
}

void AllocProj(PROJ ***Proj, int Len, int NSeq)
{
  register int j;

  (*Proj) = (PROJ **)ckalloc(Len*sizeof(PROJ *));
  
  for( j=0; j<Len; j++ ) {
    (*Proj)[j]            = (PROJ *)   ckalloc(sizeof(PROJ));
    (*Proj)[j]->NAli      = 0;
    (*Proj)[j]->Seq       = (SEQ **)   ckalloc(NSeq*sizeof(SEQ *));
    (*Proj)[j]->Score     = (float *)  ckalloc(NSeq*sizeof(float));
    (*Proj)[j]->HCri      = (float *)  ckalloc(NSeq*sizeof(float));
    (*Proj)[j]->Percent   = (float *)  ckalloc(NSeq*sizeof(float));
    (*Proj)[j]->ResN      = (int *)    ckalloc(NSeq*sizeof(int));
    (*Proj)[j]->AliLength = (int *)    ckalloc(NSeq*sizeof(int));
  }
}


void DeallocProj(PROJ ***Proj, int Len)
{
  register int j;

  for( j=0; j<Len; j++ ) {
    free((*Proj)[j]->Seq);
    free((*Proj)[j]->Score);
    free((*Proj)[j]->HCri);
    free((*Proj)[j]->Percent);
    free((*Proj)[j]->ResN);
    free((*Proj)[j]->AliLength);
    free((*Proj)[j]);
  }
  free(*Proj);
}


void Get3D(SEQ **Seq, int NSeq, COMMAND *Cmd)
{

  SEQ **Known3D;
  char *AsnFile;
  int NSeq3D = 0, ActiveChain = ERR;
  register int i, j;
  float NonRed = -1.0; /* Update 2.1 */

  for( i=0; i<NSeq; i++ )
    if( Seq[i]->Asn == NULL )
      Seq[i]->Asn = (char *)ckalloc(Seq[i]->Len*sizeof(char));
      
  if( strcmp(Cmd->StrideFile,"") || strcmp(Cmd->DsspFile,"") ) {

    if( strcmp(Cmd->StrideFile,"") )
      AsnFile = Cmd->StrideFile;
    else
      AsnFile = Cmd->DsspFile;
    
    Known3D = (SEQ **)ckalloc(MAX_SEQ*sizeof(SEQ *));

    NSeq3D = ReadSeq(Known3D,AsnFile,NonRed);
    
    for( i=0; i<NSeq3D; i++ )
      if( Known3D[i]->PdbChain == Cmd->ChainId ) {
        ActiveChain = i;
        break;
      }
    if( ActiveChain == ERR )
      die("PDB chain %c not found in %s\n",Cmd->ChainId,AsnFile);

    for( i=0; i<NSeq; i++ ) {
      
      if( !Seq[i]->ToBePredicted )
        continue;

      Concord(Seq[i],Known3D[ActiveChain]);

      for( j=0; j<Seq[i]->Len; j++ )
        if( Seq[i]->Se[j] != Known3D[ActiveChain]->Se[j] )
          die("Sequences %s and %s differ\n",Cmd->InputFile,AsnFile);
      else
        Seq[i]->Asn[j] = Known3D[ActiveChain]->Asn[j];
      
    }

    for( i=0; i<NSeq3D; i++ )
      FreeSeq(&Known3D[i]);
    free(Known3D);
    
  }
  else
    for( i=0; i<NSeq; i++ )
      for( j=0; j<Seq[i]->Len; j++ )
        Seq[i]->Asn[j] = '?';

}

void PrintProp(SEQ **Seq, int SeqN, PROJ **Proj, COMMAND *Cmd)
{
  register int j, k, m;

  printf("NAME %s\n",Cmd->InputFile);
  for( j=0; j<Seq[SeqN]->Len; j++ ) {
    printf("PRED %4d %c %c  %s > ",j,Seq[SeqN]->Se[j],Seq[SeqN]->Asn[j],Cmd->InputFile);
    for( m=0; m<NPROP; m++ )
      printf("%7.5f ",Seq[SeqN]->Prop[m][j]);
    for( k=0; k<Proj[j]->NAli; k++ )
      if( Proj[j]->ResN[k] != ERR ) {
        printf(" | %7.5f %7.5f %7.5f %7.5f %c %2d ",
               Proj[j]->Score[k],Proj[j]->Percent[k],Inf((Proj[j]->Percent[k]/100.0)),
               Proj[j]->HCri[k],Proj[j]->Seq[k]->Se[Proj[j]->ResN[k]],Proj[j]->AliLength[k]);
        for( m=0; m<NPROP; m++ )
          printf("%7.5f ",Proj[j]->Seq[k]->Prop[m][Proj[j]->ResN[k]]);
      }
      else
        printf("000.000 - | ");
    printf("\n");
    
  }
  printf("\n");
}

void Combine(SEQ **Seq, int SeqN, PROJ **Proj, COMMAND *Cmd)
{
  register int j, k, m;
  float WeightedSum, SumOfWeights, INF;
  float C = 0.8;

  if( Cmd->Verbose )
    printf("NEXT %s\n",Cmd->InputFile);
  
  for( j=0; j<Seq[SeqN]->Len; j++ ) {

    if( Cmd->Verbose )
      printf("COMB %4d %c %c  %s > ",j,Seq[SeqN]->Se[j],Seq[SeqN]->Asn[j],Cmd->InputFile);
    
    for( m=0; m<NPROP; m++ ) {
       
      if( !Cmd->Single && Seq[SeqN]->Prop[m][j] > 0.00001 ) {
        
        WeightedSum  = 0.0;
        SumOfWeights = 0.0;
        
        for( k=0; k<Proj[j]->NAli; k++ ) {
          
          INF = Inf((Proj[j]->Percent[k]/100.0));
          
          if(  Proj[j]->ResN[k] != ERR && INF > INFTHR &&
               Proj[j]->Seq[k]->Prop[m][Proj[j]->ResN[k]] > 0.0001 ) {
            
            WeightedSum  += INF*Proj[j]->Seq[k]->Prop[m][Proj[j]->ResN[k]];
            SumOfWeights += INF;
            
            if( m == 0 )
              Seq[SeqN]->AliCnt[j]++;
        
          }
        }

        Seq[SeqN]->Prop[m][j] = (C*Seq[SeqN]->Prop[m][j]+WeightedSum)/(SumOfWeights+C);
      }
      if( Cmd->Verbose )
        printf("%7.5f ",Seq[SeqN]->Prop[m][j]);
    }
    
    if( Cmd->Verbose )
      printf("%d \n",Seq[SeqN]->AliCnt[j]);
  }
}


float HsspCrit(char *Seq1, char *Seq2, int Len)
{

  register int i;
  int Length = 0;
  float Crit;
  

  for( i=0; i<Len; i++ )
    if( Seq1[i] != '-' && Seq2[i] != '-' )
      Length++;
  
  Crit = (float)(290.15*pow((float)Length,(-0.562)));
  return(Crit);
}

void Unlink(int *Beg1, int *Beg2, int *End1, int *End2, int *AliLen, int N, int *Accepted)
{
  int **Overlap;
  register int i;
 
  for( i=0; i<N; i++ )
    Accepted[i] = 1;

  Overlap = IntMatrix(N,N);

  if( FindOverlappedAlignments(Beg1,Beg2,End1,End2,N,Overlap) != 0 )
    DoUnlink(Overlap,AliLen,N,Accepted);
  
  FreeIntMatrix(Overlap,N);

}

int FindOverlappedAlignments(int *Beg1, int *Beg2, int *End1, int *End2, int N, int **Overlap)
{

  register int i, j;
  int NOvl = 0;

  for( i=0; i<N-1; i++ )
    for( j=i+1; j<N; j++ )
      if( INSIDE(Beg1[i],Beg1[j],End1[j]) || INSIDE(End1[i],Beg1[j],End1[j]) ||
	  INSIDE(Beg1[j],Beg1[i],End1[i]) || INSIDE(End1[j],Beg1[i],End1[i]) ||
          INSIDE(Beg2[i],Beg2[j],End2[j]) || INSIDE(End2[i],Beg2[j],End2[j]) ||
	  INSIDE(Beg2[j],Beg2[i],End2[i]) || INSIDE(End2[j],Beg2[i],End2[i]) ) {
	Overlap[i][j] = 1;
	NOvl++;
      }
      else
	Overlap[i][j] = 0;
  
  return(NOvl);
}

void DoUnlink(int **Overlap, int *AliLen, int N, int *Accepted)
{

  register int i, j;
  int MaxLen = 0, MaxAli = -1;

  for( i=0; i<N-1; i++ )
    for( j=i+1; j<N; j++ )
		if( Overlap[i][j] ) {
	if( AliLen[i] > MaxLen ) {
	  MaxLen = AliLen[i];
	  MaxAli = i;
	}
	if( AliLen[j] > MaxLen ) {
	  MaxLen = AliLen[j];
	  MaxAli = j;
	}
      }
	  
  if( MaxAli == -1 )
    return;

  for( i=0; i<MaxAli; i++ )
    if( Overlap[i][MaxAli] ) {
      Overlap[i][MaxAli] = 0;
      Accepted[i] = 0;
    }
  for( i=MaxAli+1; i<N; i++ )
    if( Overlap[MaxAli][i] ) {
      Overlap[MaxAli][i] = 0;
      Accepted[i] = 0;
    }
 
    DoUnlink(Overlap,AliLen,N,Accepted);
}


void NoDash(SEQ **Seq, int NSeq)
{
  register int i, j;
  int NewLen;
  char *Tmp;

  if( Seq[0]->FileType != MSF && Seq[0]->FileType != GDE && Seq[0]->FileType != Clustal )
    return;

  for( i=0; i<NSeq; i++ ) {
    strncpy(Seq[i]->AliSe,Seq[i]->Se,Seq[i]->Len);
    Seq[i]->AliLen = Seq[i]->Len;
    NewLen = 0;
    Tmp = (char *)ckalloc((Seq[i]->Len+1)*sizeof(char));
    for( j=0; j<Seq[i]->Len; j++ )
      if( Seq[i]->Se[j] != '-' && Seq[i]->Se[j] != '.' )
        Tmp[NewLen++] = Seq[i]->Se[j];
      else
      if( Seq[i]->Se[j] == '.' )
        Seq[i]->Se[j] = '-';
    strncpy(Seq[i]->Se,Tmp,NewLen);
    Seq[i]->Len = NewLen;
    free(Tmp);
    
  }
}

void Borrow(SEQ *Seq1, SEQ *Seq2, ALN *Aln)
{
  register int i;
  int IdentCount = 0;
  
  if( Seq1->AliLen != Seq2->AliLen )
    die("Aligned sequences have different length\n");
  
  Aln->Len[0] = 0;
  
  Aln->NAli        = 1;
  Aln->Beg1[0]     = 0;
  Aln->Beg2[0]     = 0;
  Aln->End1[0]     = Seq1->Len;
  Aln->End2[0]     = Seq2->Len;
  Aln->Accepted[0] = TRUE;

  for( i=0; i<Seq1->AliLen; i++ )
    if( Seq1->AliSe[i] != '-' || Seq2->AliSe[i] != '-' ) {
      Aln->A1[0][Aln->Len[0]] = Seq1->AliSe[i];
      Aln->A2[0][Aln->Len[0]] = Seq2->AliSe[i];
      if( Aln->A1[0][Aln->Len[0]] == Aln->A2[0][Aln->Len[0]] )
        IdentCount++;
      Aln->Len[0]++;
    }
      
  Aln->Percent[0] = (float)IdentCount/(float)Aln->Len[0];
  Aln->Score[0]   = 0.0; /* Not used! */
  Aln->Nc[0]      = 0; /* Not used! */
  
  Aln->HCri[0] = HsspCrit(Aln->A1[0],Aln->A2[0],Aln->Len[0]);
}

  
void NonRed(SEQ **Seq, int NSeq, COMMAND *Cmd)
{

  float CompMatrix[NAcd][NAcd] = { /* Gonnet matrix */
    { 7.6,  4.6,  4.9,  4.9,  5.7,  5.0,  5.2,  5.7,  4.4,  4.4,
      4.0,  4.8,  4.5,  2.9,  5.5,  6.3,  5.8,  1.6,  3.0,  5.3 },
    { 4.6,  9.9,  5.5,  4.9,  3.0,  6.7,  5.6,  4.2,  5.8,  2.8,
      3.0,  7.9,  3.5,  2.0,  4.3,  5.0,  5.0,  3.6,  3.4,  3.2 },
    { 4.9,  5.5,  9.0,  7.4,  3.4,  5.9,  6.1,  5.6,  6.4,  2.4,
      2.2,  6.0,  3.0,  2.1,  4.3,  6.1,  5.7,  1.6,  3.8,  3.0 },
    { 4.9,  4.9,  7.4,  9.9,  2.0,  6.1,  7.9,  5.3,  5.6,  1.4,
      1.2,  5.7,  2.2,  0.7,  4.5,  5.7,  5.2,  0.0,  2.4,  2.3 },
    { 5.7,  3.0,  3.4,  2.0, 16.7,  2.8,  2.2,  3.2,  3.9,  4.1,
      3.7,  2.4,  4.3,  4.4,  2.1,  5.3,  4.7,  4.2,  4.7,  5.2 },
    { 5.0,  6.7,  5.9,  6.1,  2.8,  7.9,  6.9,  4.2,  6.4,  3.3,
      3.6,  6.7,  4.2,  2.6,  5.0,  5.4,  5.2,  2.5,  3.5,  3.7 },
    { 5.2,  5.6,  6.1,  7.9,  2.2,  6.9,  8.8,  4.4,  5.6,  2.5,
      2.4,  6.4,  3.2,  1.3,  4.7,  5.4,  5.1,  0.9,  2.5,  3.3 },
    { 5.7,  4.2,  5.6,  5.3,  3.2,  4.2,  4.4, 11.8,  3.8,  0.7,
      0.8,  4.1,  1.7,  0.0,  3.6,  5.6,  4.1,  1.2,  1.2,  1.9 },
    { 4.4,  5.8,  6.4,  5.6,  3.9,  6.4,  5.6,  3.8, 11.2,  3.0,
      3.3,  5.8,  3.9,  5.1,  4.1,  5.0,  4.9,  4.4,  7.4,  3.2 },
    { 4.4,  2.8,  2.4,  1.4,  4.1,  3.3,  2.5,  0.7,  3.0,  9.2,
      8.0,  3.1,  7.7,  6.2,  2.6,  3.4,  4.6,  3.4,  4.5,  8.3 },
    { 4.0,  3.0,  2.2,  1.2,  3.7,  3.6,  2.4,  0.8,  3.3,  8.0,
      9.2,  3.1,  8.0,  7.2,  2.9,  3.1,  3.9,  4.5,  5.2,  7.0 },
    { 4.8,  7.9,  6.0,  5.7,  2.4,  6.7,  6.4,  4.1,  5.8,  3.1,
      3.1,  8.4,  3.8,  1.9,  4.6,  5.3,  5.3,  1.7,  3.1,  3.5 },
    { 4.5,  3.5,  3.0,  2.2,  4.3,  4.2,  3.2,  1.7,  3.9,  7.7,
      8.0,  3.8,  9.5,  6.8,  2.8,  3.8,  4.6,  4.2,  5.0,  6.8 },
    { 2.9,  2.0,  2.1,  0.7,  4.4,  2.6,  1.3,  0.0,  5.1,  6.2,
      7.2,  1.9,  6.8, 12.2,  1.4,  2.4,  3.0,  8.8, 10.3,  5.3 },
    { 5.5,  4.3,  4.3,  4.5,  2.1,  5.0,  4.7,  3.6,  4.1,  2.6,
      2.9,  4.6,  2.8,  1.4, 12.8,  5.6,  5.3,  0.2,  2.1,  3.4 },
    { 6.3,  5.0,  6.1,  5.7,  5.3,  5.4,  5.4,  5.6,  5.0,  3.4,
      3.1,  5.3,  3.8,  2.4,  5.6,  7.4,  6.7,  1.9,  3.3,  4.2 },
    { 5.8,  5.0,  5.7,  5.2,  4.7,  5.2,  5.1,  4.1,  4.9,  4.6,
      3.9,  5.3,  4.6,  3.0,  5.3,  6.7,  7.7,  1.7,  3.3,  5.2 },
    { 1.6,  3.6,  1.6,  0.0,  4.2,  2.5,  0.9,  1.2,  4.4,  3.4,
      4.5,  1.7,  4.2,  8.8,  0.2,  1.9,  1.7, 19.4,  9.3,  2.6 },
    { 3.0,  3.4,  3.8,  2.4,  4.7,  3.5,  2.5,  1.2,  7.4,  4.5,
      5.2,  3.1,  5.0, 10.3,  2.1,  3.3,  3.3,  9.3, 13.0,  4.1 },
    { 5.3,  3.2,  3.0,  2.3,  5.2,  3.7,  3.3,  1.9,  3.2,  8.3,
      7.0,  3.5,  6.8,  5.3,  3.4,  4.2,  5.2,  2.6,  4.1,  8.6 }
  };

  register int i, j;
  char *Ali1, *Ali2;
  float GapI = 6.0, GapE = 0.8, PerId;

  if( Cmd->NonRed < 0.0 )
    return;

  Ali1 = (char *)ckalloc(MAX_LEN*10*sizeof(char));
  Ali2 = (char *)ckalloc(MAX_LEN*10*sizeof(char));

  for( i=0; i<NSeq; i++ ) {

    if( !Seq[i]->Selected )
      continue;

    for( j=i+1; j<NSeq; j++ ) {

      if( !Seq[j]->Selected )
        continue;

      PerId = 100.0;
      
      if( Cmd->NonRed < 40.0 && VeryClose(Seq[i]->Se,Seq[j]->Se,Seq[i]->Len,Seq[j]->Len) ||
          (PerId = AlignGlobal(Seq[i]->Se,Seq[j]->Se,Seq[i]->Len,
                               Seq[j]->Len,CompMatrix,GapI,GapE,Ali1,Ali2) ) >= Cmd->NonRed )
        Seq[j]->Selected = FALSE;

      if( Cmd->Info )
        printf("%-19.19s <-> %-19.19s : %7.3f\n",Seq[i]->Id,Seq[j]->Id,PerId);
      
    }
    
  }

  Cmd->OutSeq = TRUE;
  OutSeq(Seq,NSeq,Cmd);

  free(Ali1);
  free(Ali2);

}

float PercentId(char *Ali1, char *Ali2)
{
  int Len, Count = 0;
  register int i;
  
  if( (Len = (int)strlen(Ali1)) != (int)strlen(Ali2) )
    die("Alignments differ in length\n");
  
  for( i=0; i<Len; i++ )
    if( Ali1[i] != '-' && Ali1[i] == Ali2[i] )
      Count++;
    
  return(100.0*(float)Count/(float)Len);
    
}

  
bool VeryClose(char *Seq1, char *Seq2, int Len1, int Len2)
{
  char *sq1, *sq2;
  int ln1, ln2, window = 6, Required, Count = 0;
  float perc = 0.6;
  register int i, j, k;
  

  if( Len1 <= Len2 ) {
    sq1 = Seq1;
    sq2 = Seq2;
    ln1 = Len1;
    ln2 = Len2;
  }
  else {
    sq1 = Seq2;
    sq2 = Seq1;
    ln1 = Len2;
    ln2 = Len1;
  }

  Required = (int)(perc*ln1);

  for( i=0; i<ln1-window; i++ )
    for( j=0; j<ln2-window; j++ ) {
      for( k=0; k<window; k++ )
        if( sq1[i+k] != sq2[j+k] )
          break;
      if( k == window ) {
        Count++;
        break;
      }
    }

  if( Count > Required )
    return(TRUE);
  else
    return(FALSE);
  
}


