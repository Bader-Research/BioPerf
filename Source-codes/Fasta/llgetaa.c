
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa34t25d2 $ - $Id: llgetaa.c,v 1.14 2003/11/02 16:41:22 wrp Exp $ */

/*
   Feb, 1998 - version for prss 

   March, 2001 - modifications to support comp_thr.c: use libpos to indicate
   whether the score is shuffled==1 or unshuffled==0.  This simplifies
   complib.c and makes comp_thr.c possible

   modified version of nxgetaa.c that generates random sequences
   for a library
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "mm_file.h"

#include "uascii.h"
#include "structs.h"

#define XTERNAL
#include "upam.h"
#undef XTERNAL

#define YES 1
#define NO 0
#define MAXLINE 512

#ifndef min
#define min(x,y) ((x) > (y) ? (y) : (x))
#endif

#if !defined(_CARBON) && defined(__MWERKS__)

extern char prompt[];
#include <StandardFile.h>
extern StandardFileReply freply;
extern FSSpec q1Spec;
extern void STFileDlog(char *prompt,StandardFileReply *freply,
						SFTypeList plist,int nl);
#endif

int nsfnum;	/* number of superfamily numbers */
int sfnum[10];	/* superfamily number from types 0 and 5 */
int nsfnum_n;
int sfnum_n[10];

static int use_stdin=0;
static char llibstr0[256];
static char llibstr1[256];
static char o_line[256];

#define NO_FORMAT 0
#define FASTA_FORMAT 1
#define GCG_FORMAT 2
static int seq_format=NO_FORMAT;
static char seq_title[200];

getseq(char *filen, int *qascii,
       unsigned char *seq, int maxs, char *libstr,
       long *sq0off)
{
  FILE *fptr;
  char line[512],*bp;
  int i, j, n;
  int ic;
  int sstart, sstop, sset=0;
  int have_desc = 0;
  int llen, l_offset;

  seq_title[0]='\0';

  sstart = sstop = -1;
#ifndef DOS
  if ((bp=strchr(filen,':'))!=NULL) {
#else
  if ((bp=strchr(filen+3,':'))!=NULL) {
#endif
    *bp='\0';
    if (*(bp+1)=='-') sscanf(bp+2,"%d",&sstop);
    else sscanf(bp+1,"%d-%d",&sstart,&sstop);
    sset=1;
  }

  if (!use_stdin) {
    if (strcmp(filen,"-") && strcmp(filen,"@")) {
      if ((fptr=fopen(filen,"r"))==NULL) {
	fprintf(stderr," could not open %s\n",filen);
	return 0;
      }
    }
    else {
      fptr = stdin;
      use_stdin = 1;
    }
  }
  else {
    fptr = stdin;
    have_desc = 1;
    strncpy(llibstr1,o_line,sizeof(llibstr1));
    l_offset = 0;
  }

  if (sset==1) {
    filen[strlen(filen)]=':';
    if (*sq0off==1 || sstart>1) *sq0off = sstart;
  }

  n=0;
  while(fgets(line,sizeof(line),fptr)!=NULL) {
    if (line[0]=='>') {
      if (have_desc) {
	strncpy(o_line,line,sizeof(o_line));
	goto last;
      }
      l_offset = 0;
      seq_format = FASTA_FORMAT;
#ifdef STAR_X
      qascii['*'] = qascii['X'];
#endif
      sfnum[0] = nsfnum = 0;

      if ((bp=(char *)strchr(line,'\n'))!=NULL) *bp='\0';
      strncpy(seq_title,line+1,sizeof(seq_title));
      strncpy(llibstr0,line+1,sizeof(llibstr0));

      if ((bp=(char *)strchr(line,' '))!=NULL) *bp='\0';
      strncpy(libstr,line+1,12);
      libstr[12]='\0';
    }
    else if (seq_format==NO_FORMAT) {
      seq_format = GCG_FORMAT;
      qascii['*'] = qascii['X'];
      l_offset = 10;
      llen = strlen(line);
      while (strncmp(&line[llen-3],"..\n",(size_t)3) != 0) {
	if (fgets(line,sizeof(line),fptr)==NULL) return 0;
	llen = strlen(line);
      }
      if ((bp=(char *)strchr(line,' '))!=NULL) *bp='\0';
      else if ((bp=(char *)strchr(line,'\n'))!=NULL) *bp='\0';
      strncpy(libstr,line,12);
      libstr[12]='\0';
      if (fgets(line,sizeof(line),fptr)==NULL) return 0;
    }

    if (seq_format==GCG_FORMAT && strlen(line)<l_offset) continue;

    if (line[0]!='>'&& line[0]!=';') {
      for (i=l_offset; (n<maxs)&&
	     ((ic=qascii[line[i]&AAMASK])<EL); i++)
	if (ic<NA) seq[n++]= ic;
      if (ic == ES) break;
    }
    else {
      if (have_desc) {
	strncpy(o_line,line,sizeof(o_line));
	goto last;
      }
      else {
	have_desc = 1;
      }
    }
  }

 last:
  if (n==maxs) {
    fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
    fflush(stderr);
  }
  if ((bp=strchr(libstr,'\n'))!=NULL) *bp = '\0';
  if ((bp=strchr(libstr,'\r'))!=NULL) *bp = '\0';
  seq[n]= EOSEQ;

  if (fptr!=stdin) fclose(fptr);

  if (sset) {
    if (sstart <= 0) sstart = 1;
    if (sstop <= 0) sstop = n;
    sstart--;
    sstop--;
    for (i=0, j=sstart; j<=sstop; i++,j++)
      seq[i] = seq[j];
    n = sstop - sstart +1;
    seq[n]=EOSEQ;
  }

  return n;
}

gettitle(filen,title,len)
  char *filen, *title; int len;
{
  FILE *fptr;
  char line[512];
  char *bp;
  int ll,sset;
#ifdef MSDOS
  char *strpbrk();
#endif
  sset = 0;

  if (use_stdin) {
    if (use_stdin == 1) {
      use_stdin++;
      strncpy(title,llibstr0,len);
    }
    else
      strncpy(title,llibstr1,len);
    return strlen(title);
  }

  if ((bp=strchr(filen,':'))!=NULL) { *bp='\0'; sset=1;}
	  
  if ((fptr=fopen(filen,"r"))==NULL) {
    fprintf(stderr," file %s was not found\n",filen);
    fflush(stderr);
    return 0;
  }

  if (sset==1) filen[strlen(filen)]=':';

  while(fgets(line,sizeof(line),fptr)!=0) {
    if (line[0]=='>'|| line[0]==';') goto found;
  }
  fclose(fptr);
  title[0]='\0';
  return 0;

 found:
#ifdef MSDOS
  bp = strpbrk(line,"\n\r");
#else
  bp = strchr(line,'\n');
#endif
  if (bp!=NULL) *bp = 0;
  strncpy(title,line,len);
  title[len-1]='\0';
  fclose(fptr);
  return strlen(title);
}	

FILE *libf=NULL;

long lpos;
char lline[MAXLINE];
int lfflag=0;	/* flag for CRLF in EMBL CDROM files */
#define LFCHAR '\015'  /* for MWC 5.5 */


int agetlib(); void aranlib();	/* pearson fasta format */

/*	the following is from fgetgb.c */

#ifdef __MWERKS__
SFTypeList llist={'TEXT',0L,0L,0L};
#define getenv mgetenv
#define LLN 1
#endif

/* a file name for openlib may now include a library type suffix */
/* only opens fasta format files */

static char libn_save[MAX_FN];
static int ldna_save=0;
static int do_shuffle;
static int shuff_cnt=10;
static int w_flag = 0;
#ifdef DEBUG
static FILE *dfile=NULL;
#endif
static unsigned char *aa_save;
static int n1_save;
static int i_even;

/* lmf_str * is used here for compatibility with the "normal" openlib,
   but is largely unnecessary */

struct lmf_str *
openlib(char *lname, int ldnaseq, int *sascii, struct mngmsg m_msg)
{
  char rline[10],libn[MAX_FN], *bp, *getenv ();
  int wcnt, ll, opnflg;
  char dfname[MAX_FN];
  int libtype;
  struct lmf_str *m_fptr;

   if (m_msg.shuff_wid > 0) w_flag = m_msg.shuff_wid;
   if (m_msg.shuff_max > shuff_cnt) shuff_cnt = m_msg.shuff_max;

#ifdef DEBUG
  if (m_msg.dfile[0]!='\0') {
    strncpy(dfname,m_msg.dfile,sizeof(dfname));
    strncat(dfname,"_rlib",sizeof(dfname));
    dfile = fopen(dfname,"w");
  }
#endif

  wcnt = 0;
  libtype = 0;

#if !defined(_CARBON) && defined(__MWERKS__)
  HSetVol(NULL,q1Spec.vRefNum,q1Spec.parID);
#endif

  strncpy(libn_save,lname,sizeof(libn_save));

  /* now allocate a buffer for the opened text file */
  if ((m_fptr = calloc(1,sizeof(struct lmf_str)))==NULL) {
    fprintf(stderr," cannot allocate lmf_str (%ld) for %s\n",
	    sizeof(struct lmf_str),lname);
    return NULL;
  }

  strncpy(m_fptr->lb_name,lname,MAX_FN);
  m_fptr->lb_name[MAX_FN-1]='\0';

  m_fptr->sascii = sascii;
  m_fptr->getlib = agetlib;
  m_fptr->ranlib = aranlib;
  m_fptr->mm_flg = 0;

  do_shuffle = 0;
  irand(0);		/* initialize the random number generator */

  return m_fptr;
}

void
closelib()
{
  if (libf!=NULL) {
    fclose(libf);
    libf = NULL;
  }
#ifdef DEBUG
  if (dfile) fclose(dfile);
#endif
}

static int ieven=0;

int
agetlib(unsigned char *seq, 
	int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont, 
	struct lmf_str *lf_fd)
{
  long sq1_off;
  int i;

  if (!do_shuffle) {
    do_shuffle = 1;
    
    if ((n1_save = getseq(libn_save,lf_fd->sascii,
			  seq,maxs,libstr,&sq1_off)) < 1)
      return n1_save;

    if ((aa_save = (unsigned char *)calloc(n1_save+1,sizeof(unsigned char)))==
	NULL) fprintf(stderr," cannot allocate %d for saved sequence\n",
		       n1_save);

    memcpy((void *)aa_save,(void *)seq,n1_save);
    *libpos = 0;
    return n1_save;
  }
  else {	/* return a shuffled sequence - here we need a window size; */
    if (shuff_cnt-- <= 0 ) return -1;
    if (w_flag > 0) wshuffle(aa_save,seq,n1_save,w_flag,&ieven);
    else shuffle(aa_save,seq,n1_save);
    seq[n1_save] = EOSEQ;
#ifdef DEBUG
    if (dfile!=NULL) {
      fprintf(dfile,">%d\n",shuff_cnt);
      for (i=0; i<n1_save; i++) {
	if (aa[seq[i]]>0) fputc(aa[seq[i]],dfile);
	else {fprintf(stderr,"error aa0[%d]: %d %d\n",
		      i,seq[i],aa[seq[i]]);}
	if (i%60 == 59) fputc('\n',dfile);
      }
      fputc('\n',dfile);
    }
#endif
    *libpos = 1;
    return n1_save;
  }
}

void
aranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp;
  int ll;

  if (seq_title[0]=='>' || seq_title[0]==';') {
    strncpy(str,&seq_title[1],cnt);
  }
  else {
    strncpy(str,seq_title,cnt);
  }
  str[cnt-1]='\0';
  if ((bp = strchr(str,'\001'))!=NULL) *bp='\0';
  else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
  else str[cnt-1]='\0';
}

/*
void
revcomp(unsigned char *seq, int n, int *c_nt)
{
  unsigned char tmp;
  int i, ni;


  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = c_nt[seq[i]];
    seq[i] = c_nt[seq[ni]];
    seq[ni] = tmp;
  }
  if ((n%2)==1) {
    i = n/2;
    seq[i] = c_nt[seq[i]];
  }
}
*/

struct lmf_str *
re_openlib(struct lmf_str *om_fptr, int outtty)
{
  return om_fptr;
}

int re_getlib(unsigned char *aa1, int n1, int maxt3, int loff, int cont,
	      int term_code, long *loffset, long *l_off,
	      struct lmf_str *m_file_p)
{
  *loffset = 1;
  *l_off = 0;
  return n1;
}

