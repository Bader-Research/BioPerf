#include "prd.h"
#include "prdprot.h"

void Output(SEQ *Seq, COMMAND *Cmd)
{

  FILE *Out;
  register int i, j;
  int NBlocks, Tail, CurrentWidth, ind;
  static bool FirstCall = TRUE;

  if( !strlen(Cmd->OutFile) )
    Out = stdout;
  else 
  if(  FirstCall == TRUE  && (Out = fopen(Cmd->OutFile,"w")) == 0 ||
       FirstCall == FALSE && (Out = fopen(Cmd->OutFile,"a")) == 0  )
    die("Can not open output file %s\n",Cmd->OutFile);

  FirstCall = FALSE;
  
  if( Seq->NHml ) {

    fprintf(Out," Info>  Identical 7-residue fragments found in:\n");
    for( i=0; i<Seq->NHml; i++ )
      fprintf(Out," Info>  %s",Seq->Hml[i]);
  }

  fprintf(Out,"\n> %s\n",Seq->Id);

  NBlocks = Seq->Len/OUTWID;
  Tail = Seq->Len % OUTWID;
  if( Tail ) 
    NBlocks++;
  else
    Tail = OUTWID;

  for( i=0; i<NBlocks; i++ ) {

    CurrentWidth = ( i<(NBlocks-1) ? OUTWID : Tail );
    ind = i*OUTWID;

    fprintf(Out,"     ");
    for( j=0; j<CurrentWidth; j++ )
      if( j != 0 && ((j+1)%10) == 0 )
	fprintf(Out,".");
      else
	fprintf(Out," ");
    fprintf(Out,"\n");

    fprintf(Out,"%-4d ",ind+1);
    for( j=0; j<OUTWID; j++ )
      if( j<CurrentWidth )
	fprintf(Out,"%c",Seq->Se[ind+j]);
      else
	fprintf(Out," ");
    fprintf(Out," %4d\n",ind+CurrentWidth);

    fprintf(Out,"     ");
    for( j=0; j<CurrentWidth; j++ )
	fprintf(Out,"%c",C2_(Seq->Pred[ind+j]));
    fprintf(Out,"\n");

    if( strcmp(Cmd->StrideFile,"") || strcmp(Cmd->DsspFile,"") ) {
      fprintf(Out,"     ");
      for( j=0; j<CurrentWidth; j++ )
	fprintf(Out,"%c",Seq->Asn[ind+j]);
      fprintf(Out,"\n");
    }

    fprintf(Out,"\n");
  }

  if( Seq->PercCorr >= 0.0 )
    fprintf(Out,"Percent correct in three states: %7.3f in chain %s%s\n",Seq->PercCorr,Seq->Id,Seq->De);
  
  if( Out != stdout )
    fclose(Out);
}

void OutputLong(SEQ *Seq, COMMAND *Cmd)
{

  FILE *Out;
  register int i;
  static bool FirstCall = TRUE;

  if( !strlen(Cmd->OutFile) )
    Out = stdout;
  else 
  if(  FirstCall == TRUE  && (Out = fopen(Cmd->OutFile,"w")) == 0 ||
       FirstCall == FALSE && (Out = fopen(Cmd->OutFile,"a")) == 0 )
    die("Can not open output file %s\n",Cmd->OutFile);

  FirstCall = FALSE;
  
  if( Seq->NHml ) {

    fprintf(Out," Info>  Identical 7-residue fragments found in:\n");
    for( i=0; i<Seq->NHml; i++ )
      fprintf(Out," Info>  %s",Seq->Hml[i]);
  }

  fprintf(Out,"\nNAME %s\n",Seq->Id);
  fprintf(Out,"HEADER  |- Residue -|  Pred  Rel      NAli   Asn\n");

  for( i=0; i<Seq->Len; i++ ) {
    fprintf(Out,"PRED %4d    %3s    %c  %c   %7.3f ",
            i+1,OneToThree(Seq->Se[i]),Seq->Se[i],Seq->Pred[i],Seq->Rel[i]);
    if( strcmp(Cmd->StrideFile,"") || strcmp(Cmd->DsspFile,"") )
      fprintf(Out,"   %-4d   %c\n",Seq->AliCnt[i],Seq->Asn[i]);
    else
      fprintf(Out,"   %-4d   ?\n",Seq->AliCnt[i]);
  }
    
  if( Seq->PercCorr >= 0.0 )
    fprintf(Out,"Percent correct in three states: %7.3f in chain %s%s\n",Seq->PercCorr,Seq->Id,Seq->De);

  if( Out != stdout )
    fclose(Out);
}

void OutSeq(SEQ **Seq, int NSeq, COMMAND *Cmd)
{

  register int i, j;
  FILE *fi;

  if( !Cmd->OutSeq )
    return;
  
  if( (int)strlen(Cmd->OutSeqFile) == 0 )
    fi = stdout;
  else
    if( (fi = fopen(Cmd->OutSeqFile,"a")) == NULL )
      die("Can not open output sequence file\n");
  

  for( i=0; i<NSeq; i++ ){
    if( Seq[i]->Selected && (Cmd->ChainId == '$' || Seq[i]->PdbChain == Cmd->ChainId ) ) {
      fprintf(fi,"> %s ",Seq[i]->Id);
      if( Seq[i]->PdbChain != '$')
        fprintf(fi,"%c ",Seq[i]->PdbChain);
      else
        fprintf(fi,"%s ",Seq[i]->De);
      fprintf(fi,"\n");
      for( j=0; j<Seq[i]->Len; j++ ) {
        if( j%60 == 0 && j != 0 )
          fprintf(fi,"\n");
        fprintf(fi,"%c",Seq[i]->Se[j]);
      }
/*      if( j%61 != 0 )*/
        fprintf(fi,"\n");
    }
  }
  if( fi != stdout )
    fclose(fi);
  
  exit(1);
  
}

void OutputEXT(SEQ *Seq, int N, char **PredEXT, float **RelEXT, float **ProbHEXT, float **ProbEEXT, float **ProbCEXT)
{

  register int i;

  for( i=0; i<Seq->Len; i++ ) {
    PredEXT[N][i] = C2_(Seq->Pred[i]);
    RelEXT[N][i] = Seq->Rel[i];
    ProbHEXT[N][i] = Seq->PHel[i];
    ProbEEXT[N][i] = Seq->PBeta[i];
    ProbCEXT[N][i] = Seq->PCoil[i];
  }

}

char C2_(char SStr)
{
  if( toupper(SStr) == 'C' )
    return('_');
  else
    return( toupper(SStr) );
}

char _2C(char SStr)
{
  if( toupper(SStr) == '_' )
    return('C');
  else
    return( toupper(SStr) );
}


