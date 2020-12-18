#include "prd.h"
#include "prdprot.h"
/* #include <console.h> */  /* For Macintosh only, see readme.mac */

/* Propensities of amino acids to form beta-bridges of type I, II, and III */
float **Bridge1, **Bridge2, **Bridge3;

/* Propensities of amino acids to form hydrogen bonds i,i+4 */
float **Alpha;

/* Memory-based reasoning (or nearest neighbour) matrix.
   See Zhang et al., 1992 for details */
float ***MbrMat;

/* Recognition thresholds for single and multiple sequence prediction */
float ThresholdsM[NTHRESH], ThresholdsS[NTHRESH];

/* Sequence/structure database */
int   *SeqDb[MAXDB], *StrDb[MAXDB], NStr;
char  *NmDb[MAXDB];

/* Reliability and probability data */
float **Rel;
int NRel;

char     *AMINOACIDS = "ARNDCQEGHILKMFPSTWYVXZB";
char     *ALNSYMBOLS = "ARNDCQEGHILKMFPSTWYVXZB-.";
char     *STRUCTURES = "HEC";

bool FirstCall = TRUE;
int MaxNameLength = 0;

int main(int argc, char **argv)
{
  int NSeq;
  register int i;
  COMMAND  *Cmd;
  SEQ  **Seq;

  /* argc = ccommand(&argv); */ /* For Macintosh only, see readme.mac */ 

  MainAlloc(&Cmd,&Seq);
  
  ProcessOptions(argv,argc,Cmd);

  NSeq = ReadSeq(Seq,Cmd->InputFile,Cmd->NonRed); /* Update 2.1 */

  if( NSeq == 1 )
    Cmd->Single = TRUE;

  NoDash(Seq,NSeq); /* Unalign if a multiple alignment is supplied */
  
  NonRed(Seq,NSeq,Cmd); /* Make sequence set non-redundant and die */

  OutSeq(Seq,NSeq,Cmd); /* Output sequence(s) and die */

  SelectSequences(Seq,NSeq,Cmd);

  Get3D(Seq,NSeq,Cmd); /* Get known assignment (STRIDE or DSSP) */

  if( !ReadTrainData(Cmd) )
    die("");

  for( i=0; i<NSeq; i++ )
    if( Seq[i]->ToBePredicted == TRUE ) {
      PredictSS(Seq[i],0,Seq[i]->Len,Seq[i]->Id,Cmd);
      StoreHomol(&Seq[i]);
    }

  Multiple(Seq,NSeq,Cmd);

  for( i=0; i<NSeq; i++ ) {

    if( Seq[i]->ToBePredicted == TRUE ) {
      
      Decide(Seq[i],Cmd);
      
      ProbAndRel(Seq[i]);
      
      Assess(Seq[i],Cmd);
      
      if( Cmd->Long )
        OutputLong(Seq[i],Cmd);
      else
        Output(Seq[i],Cmd);
    }
  }

  for( i=0; i<NSeq; i++ )
    FreeSeq(&Seq[i]);

  MainDealloc(&Cmd,&Seq);

  return 0;
}

void ProcessOptions(char **List, int ListLength, COMMAND *Cmd)
{

  int i, InpFile = 0;
  char OPTION;

  if( Uniq(List,ListLength) == 0 )
    PrintHelp("All options must be unique\n");

  DefaultCmd(Cmd);

  for( i=1; i<ListLength; i++ ) {
    
    if( *List[i] == '-' ) {
      
      OPTION = toupper((*(List[i]+1)));
      
      if(OPTION == 'D') Cmd->Dssp = TRUE;
      else if( OPTION == 'L' ) Cmd->Long    = TRUE;
      else if( OPTION == 'A' ) Cmd->All     = TRUE;
      else if( OPTION == 'V' ) Cmd->Verbose = TRUE;
      else if( OPTION == 'R' ) Cmd->OrigAli = TRUE;
      else if( OPTION == 'U' ) Cmd->LookUp  = FALSE;
      else if( OPTION == 'S' ) Cmd->Single  = TRUE;
      else if( OPTION == 'H' ) Cmd->Info    = TRUE;
      else if( OPTION == 'N' ) Cmd->NonRed  = atof(List[i]+2);
      else
        if( OPTION == 'O' ) {
	strcpy(Cmd->OutSeqFile,List[i]+2);
	Cmd->OutSeq = TRUE;
      }
      else if( OPTION == 'F' ) strcpy(Cmd->OutFile,List[i]+2);
      else if( OPTION == 'I' ) strcpy(Cmd->PredId,List[i]+2);
      else if( OPTION == 'X' ) strcpy(Cmd->StrideFile,List[i]+2);
      else if( OPTION == 'Y' ) strcpy(Cmd->DsspFile,List[i]+2);
      else if( OPTION == 'Z' ) Cmd->ChainId = toupper( *(List[i]+2) );
      else if( OPTION == 'B' ) strcpy(Cmd->DatabaseFile,List[i]+2);
      else
	PrintHelp("");
    }
    else {
      strcpy(Cmd->InputFile,List[i]);
      InpFile++;
    }
  }

  if( InpFile > 1 )
    PrintHelp("Only one input file is allowed\n");
  else
  if( !InpFile )
    PrintHelp("You must specify input file\n"); 
  
  if( strcmp(Cmd->StrideFile,"") && strcmp(Cmd->DsspFile,"") )
    PrintHelp("Please use only SSTRIDE or DDSSP assignment \n"); 
  else
  if( (strcmp(Cmd->StrideFile,"") || strcmp(Cmd->DsspFile,"")) && Cmd->ChainId == '$' )
    PrintHelp("You must specify PDB chain\n");
  else
  if( Cmd->All && Cmd->Single )
    PrintHelp("Options -a and -s are not compatible\n");

}

void PrintHelp(char *Msg)
{

  if( strcmp(Msg,"") )
    fprintf(stderr,"\nERROR: %s",Msg);
      
  fprintf(stderr,"\nProtein secondary structure prediction from single or multiple sequences\n");
  fprintf(stderr,"Version 2.1, February 1997\n");
  fprintf(stderr,"Please cite:\n");
  fprintf(stderr,"    Frishman, D. and Argos, P. (1997) Proteins, 27, 329-335.\n");
  fprintf(stderr,"    Frishman, D. and Argos, P. (1996) Prot. Eng., 9, 133-142.\n");
  fprintf(stderr,"Questions/bug reports to frishman@mips.biochem.mpg.de\n");
  fprintf(stderr,"\nUsage        : predator [Options] InputFile [ > file ]\n");
  fprintf(stderr,"Input formats: FASTA, CLUSTAL, MSF, DSSP, STRIDE\n\n");
  fprintf(stderr,"General options:  \n");
  fprintf(stderr,"  -fFileName  Output file (default is terminal)\n");
  fprintf(stderr,"  -l          Long output (default is short output)\n");
  fprintf(stderr,"  -o          Output sequence(s) and die\n");
  fprintf(stderr,"  -h          Indicate progress by dots\n");
  fprintf(stderr,"\nSelection of sequences to predict:  \n");
  fprintf(stderr,"  -a          Make prediction for All sequences in the input file\n");
  fprintf(stderr,"  -iSeqId     Make prediction for the sequence SeqId\n");
  fprintf(stderr,"By default prediction is made for the first sequence in the set\n");
  fprintf(stderr,"\nPrediction options:\n");
  fprintf(stderr,"  -s          Perform single sequence prediction. Ignore other sequences in the set.\n");
  fprintf(stderr,"  -r          Preserve the original alignment in the CLUSTAL or MSF file (do not unalign)\n");
  fprintf(stderr,"  -u          Do not copy assignment directly from the PDB database if found\n");
  fprintf(stderr,"  -d          Use DSSP target assignment (default is STRIDE, see documentation for details)\n");
  fprintf(stderr,"  -bFileName  Use database file FileName\n");
  fprintf(stderr,"\nComparison with the known assignment (for testing):\n");
  fprintf(stderr,"  -xFileName  Read STRIDE file\n");
  fprintf(stderr,"  -yFileName  Read DSSP file\n");
  fprintf(stderr,"  -zChain     PDB Chain (must be specified if option -x or -y are used)\n");
  fprintf(stderr,"\nAdditional functions\n");
  fprintf(stderr,"  -nPercentId Find a subset of sequences with no more than PercentId identity\n");
  fprintf(stderr,"              between any pair of sequences (quick and dirty algorithm)\n");
  fprintf(stderr,"\nOptions are position and case insensitive\n");
  fprintf(stderr,"\nExamples:\n\n");
  fprintf(stderr,"  predator -l globins.aln -fglobins.pred\n");
  fprintf(stderr,"  Predict secondary structure of the first sequence in the multiple alignment\n");
  fprintf(stderr,"  Create long output and write it in the file globins.pred\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"  predator 5ruba.seq -x5rub.str -za\n");
  fprintf(stderr,"  Read sequence from the file 5ruba.seq, make prediction and compare it\n");
  fprintf(stderr,"  with the known assignment for the chain A from the file 5rub.str\n");
  fprintf(stderr,"  NOTE: chain \" \" must be specified as \"-\"; e.g. -z-\n");

  exit(0);

}

void DefaultCmd(COMMAND *Cmd)
{

  Cmd->Dssp         = FALSE;
  Cmd->Info         = FALSE;
  Cmd->Long         = FALSE;
  Cmd->All          = FALSE;
  Cmd->OutSeq       = FALSE;
  Cmd->Verbose      = FALSE;
  Cmd->OrigAli      = FALSE;
  Cmd->LookUp       = TRUE;
  Cmd->Single       = FALSE;

  Cmd->Safety       =  -5;
  Cmd->DelVal       =  -14;
  Cmd->GapVal       =  -4;
  
  Cmd->NAli         = 10;
  Cmd->MinAliLength = 10;
  
  strcpy(Cmd->OutFile,"");
  strcpy(Cmd->StrideFile,"");
  strcpy(Cmd->DsspFile,"");
  strcpy(Cmd->DatabaseFile,"");
  strcpy(Cmd->OutSeqFile,"");
  strcpy(Cmd->PredId,"");

  Cmd->ChainId = '$';
  Cmd->NonRed  = -1.0;
}

int Uniq(char **List, int ListLength)
{
  int i, j;
  
  for( i=1; i<ListLength-1; i++ ) {
    if( *List[i] != '-' ) continue;
    for( j=i+1; j<ListLength; j++ ) {
      if( *List[j] != '-' ) continue;
      if( !strcmp(List[i],List[j] ) ) return(FALSE);
    }
  }
  
  return(TRUE);
}





