#include "prd.h"
#include "prdprot.h"
#include "unistd.h"

extern float **Bridge1, **Bridge2, **Bridge3;
extern float **Alpha, ***MbrMat, ThresholdsM[NTHRESH], ThresholdsS[NTHRESH];
extern int   *SeqDb[MAXDB], *StrDb[MAXDB], NStr;
extern char  *NmDb[MAXDB];

extern float **Rel;
extern int NRel;

extern char  *AMINOACIDS;
extern char  *ALNSYMBOLS;
extern char  *STRUCTURES;

extern int MaxNameLength;

bool ReadTrainData(COMMAND *Cmd)
{

  FILE *fi;
  BUFFER FileName, Tmp;
  char *PRE_DIR;

  if( strcmp(Cmd->DatabaseFile,"") )
    strcpy(FileName,Cmd->DatabaseFile);
  else {
    if( (PRE_DIR = getenv("PRE_DIR")) != 0 )
      strcpy(FileName,PRE_DIR);
    else
      strcpy(FileName,".");
    
    if( Cmd->Dssp )
      sprintf(Tmp,"%cdssp.dat",DIRDELIM);
    else
      sprintf(Tmp,"%cstride.dat",DIRDELIM);
    
    strcat(FileName,Tmp);
  }

  if( (fi = fopen(FileName,"r")) == 0 )
    return(wrer("Can not open %s\n",FileName));
  
  if( !ReadDatabase(fi)     ||
      !ReadProp(fi,Bridge1) ||
      !ReadProp(fi,Bridge2) ||
      !ReadProp(fi,Bridge3) ||
      !ReadProp(fi,Alpha)   ||
      !ReadMbr(fi)          ||
      !ReadThresholds(fi)   ||
      !ReadRel(fi,Cmd) )
    return(wrer("Error reading %s\n",FileName));
  
  Normalize();

  fclose(fi);

  return(TRUE);
}


bool ReadProp(FILE *fi, float **Mat)
{

  BUFFER Buffer;
  char *Field[MAX_FIELD];
  int i, j, NFields;

  i = 0;
  while( i<NAcd ) {
    
    if( !fgets(Buffer,BUFSZ,fi) )
      return(FALSE);
    
    if( Buffer[0] == '#' || !(NFields = SplitString(Buffer,Field,NAcd) ) )
      continue;
    
    if( NFields != NAcd )
      return(FALSE);
    
    for( j=0; j<NAcd; j++ )
      Mat[i][j] = atof(Field[j]);
    
    i++;
  }

  return(TRUE);
}

bool ReadMbr(FILE *fi)
{

  BUFFER Buffer;
  char *Field[MAX_FIELD];
  int i, j, k, NLines, LineCount=0, NFields;
  
  
  NLines = MBRWIN*NAcd;
  i = j = 0;
  while( LineCount < NLines ) {

    if( !fgets(Buffer,BUFSZ,fi) )
      return(FALSE);

    if( Buffer[0] == '#' || !(NFields = SplitString(Buffer,Field,NAcd) ) )
      continue;

    if( NFields != NAcd )
      return(FALSE);

    for( k=0; k<NAcd; k++ )
      MbrMat[i][j][k] = atof(Field[k]);

    j++;
    if( j == NAcd ) {
      j = 0;
      i++;
    }
    LineCount++;
  }

  return(TRUE);
}


bool ReadDatabase(FILE *fi)
{

  char Buffer[MAX_LEN], *p;
  int CurrLen, Flag = 0, RecordType, RecordTypePrev;
  register int i;
  
  for( i=0; i<MAXDB; i++ ) {
    NmDb[i] = NULL;
    SeqDb[i] = NULL;
    StrDb[i] = NULL;
  }
  
  NStr = 0;
  RecordTypePrev = 3;
  
  while( fgets(Buffer,MAX_LEN,fi) != 0 ) {

    if( Buffer[0] == 'N' && Buffer[1] == 'M' )
      RecordType = 0;
    else
    if( Buffer[0] == 'S' && Buffer[1] == 'Q' )
      RecordType = 1;
    else
    if( Buffer[0] == 'S' && Buffer[1] == 'S' )
      RecordType = 2;
    else
      RecordType = 3;

    if( RecordType == 0 && RecordTypePrev != 2 &&  RecordTypePrev != 2 ||
        RecordType == 1 && RecordTypePrev != 0 ||
        RecordType == 2 && RecordTypePrev != 1 ) {
      for( i=0; i<NStr; i++ ) {
        if( NmDb[i] != NULL )
          free(NmDb[i]);
        if( SeqDb[i] != NULL )
          free(SeqDb[i]);
        if( StrDb[i] != NULL )
          free(StrDb[i]);
        return(FALSE);
      }
    }
        
    if( RecordType == 3 ) {
      if( Flag == 1 )
        break; /* Done reading sequence/structure pairs */
      else
        continue; /* They are still ahead! */
    }

    Flag = 1; /* Began reading sequence/structure pairs */
    
    p = Buffer+3;
    CurrLen = (int)strlen(p);
    switch(RecordType) {
      case 0:
        NmDb[NStr] = (char *)ckalloc((CurrLen+2)*sizeof(char));
        strcpy(NmDb[NStr],p);
        break;
      case 1:
        SeqDb[NStr] = (int *)ckalloc((CurrLen+2)*sizeof(int));
        SeqDb[NStr][0] = CurrLen;
        Char2Int(p,SeqDb[NStr]+1,CurrLen);
        break;
      case 2:
        StrDb[NStr] = (int *)ckalloc(strlen(Buffer+2)*sizeof(int));
        StrDb[NStr][0] = CurrLen;
        Str2Int(p,StrDb[NStr++]+1,CurrLen);
    }
    RecordTypePrev = RecordType;
  }
  
  return(TRUE);
}
      
bool ReadThresholds(FILE *fi)
{
  BUFFER Buffer;
  char *Field[MAX_FIELD];
  int NFields, i, Flag = 0;

  while( fgets(Buffer,BUFSZ,fi ) ) {
    
    if( Buffer[0] == '#' || !(NFields = SplitString(Buffer,Field,NAcd) ) ) 
        continue;
    
    if( NFields != NTHRESH )
      return(FALSE);

    if( Flag == 0 ) {
      for( i=0; i<NTHRESH; i++ )
        ThresholdsM[i] = atof(Field[i]);
      Flag = 1;
    } else {
      for( i=0; i<NTHRESH; i++ )
        ThresholdsS[i] = atof(Field[i]);
    break;
    }
  }
  return(TRUE);
}


bool ReadRel(FILE *fi, COMMAND *Cmd)
{

  char Buffer[MAX_LEN], *p;
  char *Field[MAX_FIELD];
  register int i;
  int N = 4+2*NPROP;
  
  NRel = 0;
  
  while( fgets(Buffer,MAX_LEN,fi) != 0 ) {

    if( !Cmd->Single && (Buffer[0] != 'R' || Buffer[1] != 'M') ||
         Cmd->Single && (Buffer[0] != 'R' || Buffer[1] != 'S') )
      continue;
      
    p = Buffer+3;

    if( SplitString(p,Field,N) != N ) 
      return(FALSE);
    
    Rel[NRel] = (float *)ckalloc(N*sizeof(float));

    for( i=0; i<N; i++ )
      Rel[NRel][i] =  atof(Field[i]);
    
    NRel++;

    if( NRel == MAXREL )
      return(TRUE);
  }

  return(TRUE);
} 
      
int ReadSeq(SEQ **Seq, char *InputFile, float NonRed) /*Update 2.1 */
{
  FILE *fi;
  int NSeq;
  int (*(Rd)[NFILETYPES])(SEQ **Seq, FILE *fi) =
  { ReadFasta, ReadClustal, ReadMSF, ReadStride, ReadDssp };
  register int i;
  enum FILETYPE FileType;
  
  if( !(fi = fopen(InputFile,"r")) )
    die("Can not open input file %s\n",InputFile);

  FileType = GetFileType(fi,InputFile);

  fseek(fi,0,SEEK_SET);

  if( (NSeq = Rd[(int)(FileType)](Seq,fi)) < 1 )
    die("No sequences in %s\n",InputFile);
  
  for( i=0; i<NSeq; i++ ) {
    Seq[i]->FileType = FileType;
    Init2Seq(Seq[i],NonRed);
  }

/*  printf("Input file is %d with %d sequences\n",FileType,NSeq);*/

  fclose(fi);

  for( i=0; i<NSeq; i++ )
    if( (int)strlen(Seq[i]->Id) > MaxNameLength )
      MaxNameLength = (int)strlen(Seq[i]->Id);

  return(NSeq);
}

#define NFIRSTLINES 10

enum FILETYPE GetFileType(FILE *fi, char *FileName)
{
  char Lines[NFIRSTLINES][BUFSZ];
  int LineCount = 0;
  register int i, j;

  while(  LineCount < NFIRSTLINES && fgets(Lines[LineCount],BUFSZ,fi) ) {
    for( j=0; j<(int)strlen(Lines[LineCount]); j++ )
      Lines[LineCount][j] = tolower(Lines[LineCount][j]);
    LineCount++;
  }

  if( LineCount == 0 )
    die("Empty file %s\n",FileName);

  if( Lines[0][0] == '>' || Lines[1][0] == '>' || Lines[2][0] == '>' ||
      Lines[3][0] == '>' || Lines[4][0] == '>' || Lines[5][0] == '>' ||
      Lines[6][0] == '>' || Lines[7][0] == '>' || Lines[8][0] == '>' )
    return(Fasta);

/*
  if( Lines[0][0] == '%' )
    return(GDE);
    */

  if( Lines[0][0] == 'r' && Lines[0][1] == 'e' && Lines[0][2] == 'm' )
    return(Stride);

  if(strstr(Lines[0],"sander") || ( LineCount > 1 && strstr(Lines[1],"sander")) )
    return(Dssp);

  if( strstr(Lines[0],"clustal") )
    return(Clustal);

  for( i=0; i<LineCount; i++ )
    if( strstr(Lines[i],"msf") != NULL && strstr(Lines[i],"..") != NULL )
      return(MSF);

  die("Unknown file type %s\n",FileName);
  return(Unknown);
}

int ReadMSF(SEQ **Seq, FILE *fi)
{
  BUFFER Buffer;
  char *Field[MAX_FIELD];
  int NSeq = 0, SeqCount = 0;
  bool Start = FALSE;
  register int i;

  while( fgets(Buffer,BUFSZ,fi) ) {
    
    
    if( Buffer[0] == '/' && Buffer[1] == '/' ) {
      Start = TRUE;
      if( NSeq == 0 )
        die("No sequence names found in MSF file\n");
      continue;
    }
    
    if( Start == FALSE ) {
      /* Update 2.1 */
      if( SplitString(Buffer,Field,2) < 2 )
        continue;
      /* End update 2.1 */
      if( !strcmp(Field[0],"Name:") )
        Init1Seq(&Seq[NSeq++],Field[1],Field[1]);
      continue;
    }
    else {

      if( (int)strspn(Buffer, " 0123456789\n") == (int)strlen(Buffer) )
        continue;
        
      SplitString(Buffer,Field,1);

      if( strcmp(Field[0],Seq[SeqCount]->Id) )
        die("Error reading MSF file\n");

      /* Update 2.1 */
      for( i=1; i<(int)strlen(Buffer); i++ )
        if( Buffer[i] == ' ' && Buffer[i-1] != ' ' )
          break;
      /* End of update 2.1 */

      Add(&Seq[SeqCount++]->Se,Buffer+i+1,ALNSYMBOLS);

      if( SeqCount == NSeq )
        SeqCount = 0;
    }
  }

  return(NSeq);

}

int ReadClustal(SEQ **Seq, FILE *fi)
{
  BUFFER Buffer;
  char *Field[MAX_FIELD];
  int BlockCount = -1, LineCount, NSeq;
  bool Block = FALSE;

  if( !fgets(Buffer,BUFSZ,fi) )  /* Skip header line */
      return(FALSE);                 

  while( fgets(Buffer,BUFSZ,fi) ) {

    if( Buffer[0] == ' ' )       /* Skip marker lines */
      continue;

    if( (int)strlen(Buffer) < 2 ) {   /* End of sequence block */
      Block = FALSE;
      continue;
    }

    if( Block == FALSE ) {          /* Beginning of sequence block */
      if( BlockCount == 0 )
        NSeq = LineCount;
      else
      if( BlockCount != -1 && NSeq != LineCount )
        die("Inconsistent number of sequences in alignment blocks\n");
      BlockCount++;
      LineCount = -1;
      Block = TRUE;
    }
    
    if( SplitString(Buffer,Field,10) != 2 )
      die("Strange line %s in input file\n",Buffer);
    
    LineCount++;
    
    if( BlockCount == 0 )
      Init1Seq(&Seq[LineCount],Field[0],Field[0]);
    
    Add(&Seq[LineCount]->Se,Field[1],ALNSYMBOLS);
  }
  
  return(++LineCount);
}


int ReadFasta(SEQ **Seq, FILE *fi)
{
  BUFFER Buffer;
  char *Field[MAX_FIELD];
  int NSeq = -1, NFields;
  register int i;

  while( fgets(Buffer,BUFSZ,fi) ) {

    if( Buffer[0] == '#' || (int)strlen(Buffer) < 2 )
      continue;
    
    if( Buffer[0] == '>' ) {
      
      NSeq++;
      NFields = SplitString(Buffer+1,Field,100);
      
      if( NFields == 0 )
        Init1Seq(&Seq[NSeq],"","");
      else
      if( NFields == 1 )
        Init1Seq(&Seq[NSeq],Field[0],"");
      else {
        Init1Seq(&Seq[NSeq],Field[0],Field[1]);
        for( i=2; i<NFields; i++ ) {
          Add(&Seq[NSeq]->De," ","");
          Add(&Seq[NSeq]->De,Field[i],"");
        }
      }
      
    }
    else
      Add(&Seq[NSeq]->Se,Buffer,ALNSYMBOLS);
  }
      
  return(++NSeq);
}

int ReadStride(SEQ **Seq, FILE *fi)
{
  
  BUFFER Buffer, Buffer1;
  char *Field[MAX_FIELD];
  int NSeq = -1;
  int StrideWidth = 50;
  register int i;


  while( fgets(Buffer,BUFSZ,fi) ) {
    
    if( Buffer[0] == 'L' && Buffer[1] == 'O' && Buffer[2] == 'C' )
      break;

    if( Buffer[0] == 'C' && Buffer[1] == 'H' && Buffer[2] == 'N' ) {
      
      NSeq++;
      if( SplitString(Buffer+1,Field,3) < 3 )
        die("Error in stride file\n");
      Init1Seq(&Seq[NSeq],Field[1],Field[1]);
      Seq[NSeq]->PdbChain = *Field[2];
    }
    else
    if( Buffer[0] == 'S' && Buffer[1] == 'E' && Buffer[2] == 'Q' ) {
	  
      if( !fgets(Buffer1,BUFSZ,fi) || Buffer1[0] != 'S' || Buffer1[1] != 'T' || Buffer1[2] != 'R' )
        die("Error in stride file\n");

      for( i=10; i<10+StrideWidth; i++ ) {
        if( Buffer[i] == ' ' )
          break;
        if( Buffer1[i] == ' ' || Buffer1[i] == 'T' || Buffer1[i] == 'B' ||
            Buffer1[i] == 'I' || Buffer1[i] == 'G' )
          Buffer1[i] = 'C';
        
      }   

      Buffer[i] = '\0';
      Buffer1[i] = '\0';
            
      Add(&Seq[NSeq]->Se,Buffer+10,AMINOACIDS);
      Add(&Seq[NSeq]->Asn,Buffer1+10,STRUCTURES);
    }
  }

  return(++NSeq);
}


int ReadDssp(SEQ **Seq, FILE *fi)
{
  int NSeq = -1, Start = 0, i;
  bool DuplicateChain = FALSE;
  BUFFER Buffer, Header;
  char *Tmp, *Tmp1, CurrentChainId;
  
  Tmp  = (char *)ckalloc(MAX_LEN*sizeof(char));
  Tmp1 = (char *)ckalloc(MAX_LEN*sizeof(char));

  while( fgets(Buffer,BUFSZ,fi) != NULL && !DuplicateChain ) {

    if( Buffer[2] == '#' ) Start = 1;
    else
    if( Buffer[0] == 'H' && Buffer[1] == 'E' && Buffer[2] == 'A' && Buffer[3] == 'D' ) {
      strncpy(Header,Buffer+62,4);
      Header[4] = '\0';
    }
    else
    if( Start == 1 ) {

      if( Buffer[13] == '!' ) continue;

      if( NSeq > 0 ) {
        for( i=0; i<NSeq; i++ ) 
          if( Seq[i]->PdbChain == CurrentChainId || 
	      (Seq[i]->PdbChain == '-' && CurrentChainId == '-') ) {
            DuplicateChain = TRUE;
            fprintf(stderr,"# Duplicate chain(s) in DSSP file\n");
	    break;
	  }
      }
      
      if( DuplicateChain )
        break;
      
      if( NSeq == -1 || Buffer[11] != CurrentChainId ) {

        if( NSeq != -1 ) {
          if( Seq[NSeq]->Len ) {
            Seq[NSeq]->Se = (char *)ckalloc(Seq[NSeq]->Len*sizeof(char));
            Seq[NSeq]->Asn = (char *)ckalloc(Seq[NSeq]->Len*sizeof(char));
            strncpy((Seq[NSeq]->Se),Tmp,Seq[NSeq]->Len);
            strncpy((Seq[NSeq]->Asn),Tmp1,Seq[NSeq]->Len);
          }
          else
            die("No sequence in input file\n");
        }
          
        NSeq++;
	Init1Seq(&Seq[NSeq],Header,Header);
        Seq[NSeq]->Len = 0;
        CurrentChainId = Buffer[11];
        Seq[NSeq]->PdbChain = ((CurrentChainId == ' ') ? ('-') : CurrentChainId);
      }
      
      if( islower(Buffer[13]) ) Buffer[13] = 'C';
      
      Tmp[Seq[NSeq]->Len] = Buffer[13];
      
      if( Buffer[16] == ' ' )
        Tmp1[Seq[NSeq]->Len] = 'C';
      else
        Tmp1[Seq[NSeq]->Len] = Buffer[16];
      
      Seq[NSeq]->Len++;
    }
  }

  if( NSeq != -1 && Seq[NSeq]->Len ) {
    Seq[NSeq]->Se = (char *)ckalloc(Seq[NSeq]->Len*sizeof(char));
    Seq[NSeq]->Asn = (char *)ckalloc(Seq[NSeq]->Len*sizeof(char));
    strncpy((Seq[NSeq]->Se),Tmp,Seq[NSeq]->Len);
    strncpy((Seq[NSeq]->Asn),Tmp1,Seq[NSeq]->Len);
  }
  else
    die("No sequence in input file\n");

  free(Tmp);
  free(Tmp1);
  
  return(++NSeq);
}

void Init1Seq(SEQ **Seq, char *Id, char *De)
{
  (*Seq) = (SEQ *)ckalloc(sizeof(SEQ));
  (*Seq)->Id = (char *)ckalloc((strlen(Id)+1)*sizeof(char));
  (*Seq)->De = (char *)ckalloc((strlen(De)+1)*sizeof(char));
  strcpy((*Seq)->Id,Id);
  strcpy((*Seq)->De,De);
  (*Seq)->Se       = NULL;
  (*Seq)->AliSe    = NULL;
  (*Seq)->Asn      = NULL;
  (*Seq)->KnownAsn = NULL;
  (*Seq)->PdbChain = '$';
  (*Seq)->Len      = 0;
  (*Seq)->AliLen   = 0;
  (*Seq)->NHml     = 0;
}

void Init2Seq(SEQ *Seq, float NonRed) /* Update 2.1 */
{
  register int i;
  
  if( Seq->Len == 0 )
    Seq->Len = (int)strlen(Seq->Se);

  Seq->AliSe    = (char *)ckalloc(Seq->Len*sizeof(char));

  if( NonRed < 0.0 ) {
    Seq->Pred     = (char *)ckalloc(Seq->Len*sizeof(char));
    Seq->KnownAsn = (char *)ckalloc(Seq->Len*sizeof(char));
    Seq->AliCnt   = (int *)ckalloc(Seq->Len*sizeof(int));
    Seq->Rel      = (float *)ckalloc(Seq->Len*sizeof(float));
    Seq->PHel     = (float *)ckalloc(Seq->Len*sizeof(float));
    Seq->PBeta    = (float *)ckalloc(Seq->Len*sizeof(float));
    Seq->PCoil    = (float *)ckalloc(Seq->Len*sizeof(float));
    Seq->Prop     = FloatMatrix(NPROP,Seq->Len);

    for( i=0; i<Seq->Len; i++ ) {
      Seq->AliCnt[i] = 0;
      Seq->Rel[i] = 0.0;
      Seq->PHel[i] = 0.0;
      Seq->PBeta[i] = 0.0;
      Seq->PCoil[i] = 0.0;
    }
  }

  Seq->Selected = TRUE;
  Seq->ToBePredicted = FALSE;
}


void Add(char **Dest, char *Source, char *Alphabet)
{
  int SourceLength, OldDestLength, NewDestLength, AlpLength, Good = 1;
  register int i, j = 0;

  SourceLength = (int)strlen(Source);
  AlpLength = (int)strlen(Alphabet);

  for( i=0; i<SourceLength; i++ ) 
    if( isalpha(Source[i]) ) /* Update 2.1 to fix idiotic bug in gcc for Linux */
    Source[i] = toupper(Source[i]);

  if( !AlpLength )
    Good += SourceLength;
  else {
    for( i=0; i<SourceLength; i++ ) 
      if( strchr(Alphabet,Source[i]) )
	Good++;
  }

  if( (*Dest) == NULL ) {
    OldDestLength = 0;
    NewDestLength = Good;
  }
  else {
    OldDestLength = (int)strlen(*Dest);
    NewDestLength = OldDestLength+Good;
  }

  if( (*Dest) != NULL )
    (*Dest) = (char *)myrealloc((char *)(*Dest),NewDestLength*sizeof(char));
  else
    (*Dest) = (char *)ckalloc(NewDestLength*sizeof(char));


  for( i=0; i<SourceLength; i++ )
    if( !AlpLength || strchr(Alphabet,Source[i]) )
      (*Dest)[OldDestLength + (j++) ] = Source[i];
  (*Dest)[OldDestLength +j] = '\0';
}


void Concord(SEQ *Seq, SEQ *Known3D)
{
  register int i;
  
  if( Seq->Se[1] == Known3D->Se[0] &&
      Seq->Se[2] == Known3D->Se[1] &&
      Seq->Se[3] == Known3D->Se[2] &&
      Seq->Se[4] == Known3D->Se[3] &&
      Seq->Se[5] == Known3D->Se[4] &&
      Seq->Se[6] == Known3D->Se[5] ) {
    for( i=0; i<Seq->Len-1; i++ )
      Seq->Se[i] = Seq->Se[i+1];
    Seq->Len--;
  }
  
  if( Seq->Len > Known3D->Len )
    Seq->Len = Known3D->Len;
  else
  if( Seq->Len < Known3D->Len )
    Known3D->Len = Seq->Len;
}


int ReadSeqEXT(char **SeqEXT, char **NamesEXT, int *LenEXT, int NSeqEXT, int SeqToPredictEXT, SEQ **Seq)
{

  register int i, j;
  float Dummy = -1.0; /* Update 2.1 */

  for( i=0; i<NSeqEXT; i++ ) {
    Init1Seq(&Seq[i],NamesEXT[i],"");
    Seq[i]->Se = (char *)ckalloc(LenEXT[i]*sizeof(char));
    for( j=0; j<LenEXT[i]; j++ )
      Seq[i]->Se[j] = toupper(SeqEXT[i][j]);
    Seq[i]->Len = LenEXT[i];
    Seq[i]->FileType = Extern;
    Init2Seq(Seq[i],Dummy);
    if( i == SeqToPredictEXT || SeqToPredictEXT == -1 )
      Seq[i]->ToBePredicted = TRUE;
  }    
  
  for( i=0; i<NSeqEXT; i++ )
    if( (int)strlen(Seq[i]->Id) > MaxNameLength )
      MaxNameLength = (int)strlen(Seq[i]->Id);

  return(NSeqEXT);
}

