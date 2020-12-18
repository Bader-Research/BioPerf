/*      align.c
	protein driver for linear sequence comparison method

	March, 1995 - increased space for initseq() for long long gaps,
	changed default gap penalty to -14, -4 for protein.

*/

#include <stdio.h>
#include <stdlib.h>

/*#include <ctype.h>*/

char *refstr="\nPlease cite:\n X. Huang and W. Miller (1991) Adv. Appl. Math. 12:373-381\n";

#define TRUE 1
#define FALSE 0

#ifndef BIGMEM
#define MAXTST 2000	/* longest test sequence */
#define MAXLIB 2000
#define MAXDIAG (MAXTST+MAXLIB)
#else
#define MAXTST 10000
#define MAXLIB 50000
#define MAXDIAG (MAXTST+MAXLIB)
#endif

FILE *outfd;		/* fd for output file */

/* globals for matching */

long lmark;		/* position in library file from ftell() */

char libstr[21];	/* partial title from library sequence */
char name0[11], name1[11];	/* for labeling output */

char *aa0=NULL, *aa1=NULL;	/* amino acid sequence data */
char *seqc0, *seqc1;	/* aligned sequences */
long sq0off=1, sq1off=1;

int dnaseq, lcont;
int bktup, bkfact, scfact, bestoff, bestscale, histint, bestmax;

int nseq=2;		/* same sequence twice */
int maxn;		/* max space for lib sequence */
int n0, n1;	/* length of aa0, length of aa1, n0+n1,	diagonal offset */
long loffset = 0l;		/* offset into sequence */

/*  the following are defaults for values that are read by
    pam.c from *.mat if SMATRIX is defined */

int nshow; char rline[20],sline[20];

/* output options */
int showall,markx, llen;

char ttitle[60], ltitle[60];

int smark[4] = {-10000,-10000,-10000,-10000};

int min0,min1,max0,max1,mins;
#ifdef TPLOT
char lvstr[40];
#endif

extern int optind;
char *libenv, *aaenv, *smptr;
char smstr[40];

#include "upam.gbl"		/* includes pam array */
int reformat_seq_file(char *in, char *out1,char *out2);

main(argc, argv)
        int argc; char **argv;
{
	int a;
        char tname[40], lname[40], qline[40];
	int K;
	int itemp, iln, nln;
	char *getenv(), *cptr, *bp, *strchr();
	float percent;
	char ced_seq0[10000];
	char ced_seq1[10000];
	int ced_use_id;
	char seq1_file[1000];
	char seq2_file[1000];

	initenv(argc,argv);

	if ((aa0=calloc(MAXTST+MAXLIB,sizeof(char)))==0) {
		fprintf(stderr," cannot allocate sequence array\n");
		exit(1);
		}
	maxn = MAXTST+MAXLIB;

        sprintf ( seq1_file, "%s.1", argv[1]);
	sprintf ( seq2_file, "%s.2", argv[1]);
	
	reformat_seq_file( argv[1], seq1_file, seq2_file);

	
	       
	/*NOTE: 
	  lalign seq1 seq2 V1 V2
	  V1: Positive=Number of local aln=  V1
	  Negative=Min Score Local aln= -V1
	  NULL    =>Number of local aln=10
	  
	  V2:
	  1: use_id as weight
	  0: use raw_score as weight (default)
	*/
	strncpy(tname,seq1_file,40);
	if ((n0=getseq(tname,aa0,maxn,&dnaseq, ced_seq0))==0) 
	    {
	    fprintf(stderr," %s : %s sequence not found\n",tname,sqtype);
	    exit(1);
	    }
	resetp(dnaseq);
	dnaseq = -1;
	strncpy(lname,seq2_file,40);
	K=10;
	ced_use_id=0;	
	if (argc>=4) 
	    {
	    sscanf(argv[3],"%d",&K);
	    if ( K==0)K=10;
	    }
	if ( argc>=5)
	    sscanf(argv[4],"%d",&ced_use_id);
		    

	if (strcmp(tname,lname)==0) nseq=1;
	else nseq=2;

 	strncpy(name0,tname,6);
 	gettitle(tname,ttitle,50);
 	if (strlen(ttitle)>0)
		if (*ttitle=='>') strncpy(name0,&ttitle[1],6);
		else strncpy(name0,ttitle,6);
	else
		strncpy(name0,tname,6);
	name0[6]='\0';
	if ((bp=strchr(name0,' '))!=NULL) *bp='\0';


	if ( strcmp(argv[2], "stdout")==0)outfd=stdout;
	else outfd = fopen ( argv[2],"w");

	aa1 = aa0 + n0 + 2;
	maxn -= n0 + 3;

	n1=getseq(lname,aa1,maxn,&dnaseq, ced_seq1);
	

	gettitle(lname,ltitle,50);
	if (strlen(ltitle)>0)
	    if (*ltitle=='>') strncpy(name1,&ltitle[1],6);
	    else strncpy(name1,ltitle,6);
	else strncpy(name1,lname,6);
	name1[6]='\0';
	if ((bp=strchr(name1,' '))!=NULL) *bp='\0';
	




	initseq(min(n0,n1)*2);
	if ( strcmp ( ttitle+1, ltitle+1)!=0)
	   {
	   fprintf ( outfd, "2\n");
	   fprintf ( outfd, "%s %d %s\n",ttitle+1, n0, ced_seq0);
	   fprintf ( outfd, "%s %d %s\n",ltitle+1, n1, ced_seq1);
	   fprintf ( outfd, "#0 1\n");
	   }
	else
	   {
	   fprintf ( outfd, "1\n");
	   fprintf ( outfd, "%s %d %s\n",ttitle+1, n0, ced_seq0);
	   fprintf ( outfd, "#0 0\n");
	   }
	initpam2();	/* convert 1-d pam to 2-d pam2 */

	SIM(aa0-1,aa1-1,n0,n1,K,pam2,-(gdelval-ggapval),-ggapval,nseq, ced_use_id, &outfd);
	fprintf ( outfd, "CPU 0\n");	
	remove ( seq1_file);
	remove ( seq2_file);
	fclose ( outfd);
	return;
}

extern int *sascii, nascii[], aascii[];

initenv(argc,argv)
	int argc;
	char **argv;
{
	char *cptr, *getenv();
	int copt, getopt();
	extern char *optarg;

	libenv="\0";
	aaenv="\0";

	sascii = aascii;
	pam = abl50;
	gdelval = -14;
	ggapval = -4;
	strncpy(smstr,"BLOSUM50",sizeof(smstr));
	smptr=smstr;
	sq = aa;
	hsq = haa;
	nsq = naa;
	dnaseq = 0;

	showall = 1;

	if ((cptr=getenv("LINLEN"))!=NULL) sscanf(cptr,"%d",&llen);
	else llen = 60;
	if (llen>=200) llen=200-1;
	markx=0;
	if ((cptr=getenv("MARKX"))==NULL) markx=0;
	else sscanf(cptr,"%d",&markx);

	if (dnaseq>=0) {
		if ((smptr=getenv("SMATRIX"))!=NULL && initpam(smptr)) {
			dnaseq = -1;
			}
		else smptr=smstr;
		}
	}

resetp(dnaseq)
     int dnaseq;
{
  if (dnaseq==1) {
    fprintf(stderr," resetting to DNA matrix\n");
    pam = npam;
    strncpy(smstr,"DNA",sizeof(smstr));
    smptr = smstr;
    gdelval = -16;
  }
}

int match, mismh;

initpam2()
{
	int i, j, k, tmp;

	match = -1000; mismh = 1000;
	k=0;
	for (i=0; i<nsq; i++)
		for (j=0; j<=i; j++) {
			tmp=pam2[j][i] = pam2[i][j] = pam[k++];
			if (tmp>match) match=tmp;
			if (tmp<mismh) mismh=tmp;
		      }
	}

int smin0, smin1, smins;	/* set bounds for discons */

calcons(aa0,n0,aa1,n1,res,nc,nident)
     char *aa0, *aa1;
     int n0, n1;
     int *res;
     int *nc;
     int *nident;
{
  int i0, i1;
  int op, nid, lenc, nd, ns, itmp;
  char *sp0, *sp1;
  int *rp;
  
  /* first fill in the ends */
  min0--; min1--;

  smin0 = min0;
  smin1 = min1;
  smins = mins = 0;

/* now get the middle */

  sp0 = seqc0+mins;
  sp1 = seqc1+mins;
  rp = res;
  lenc = nid = op = 0;
  i0 = min0;
  i1 = min1;
  
  while (i0 < max0 || i1 < max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      *sp0 = sq[aa0[i0++]];
      *sp1 = sq[aa1[i1++]];
      lenc++;
      if (*sp0 == *sp1) nid++;
      else if ((dnaseq==1) && (*sp0=='T' && *sp1=='U') ||
	       (*sp0=='U' && *sp1=='T')) nid++;
      sp0++; sp1++;
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {
	*sp0++ = '-';
	*sp1++ = sq[aa1[i1++]];
	op--;
	lenc++;
      }
      else {
	*sp0++ = sq[aa0[i0++]];
	*sp1++ = '-';
	op++;
	lenc++;
      }
    }
  }

  *nident = nid;
  *nc = lenc;
/*	now we have the middle, get the right end */
  nd = 0;
  return mins+lenc+nd;
}

initseq(seqsiz)		/* initialize arrays */
	int seqsiz;
{
	seqc0=calloc(seqsiz,sizeof(char));
	seqc1=calloc(seqsiz,sizeof(char));
	if (seqc0==NULL || seqc1==NULL)
		{fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
		 exit(1);}
	}

freeseq()
{
	free(seqc0); free(seqc1);
	}

int reformat_seq_file(char *in, char *out1,char *out2)
    {
    int a, b, c, d;
    FILE *fp;
    FILE *fp1;
    FILE *fp2;

    fp=fopen ( in, "r");
    fp1=fopen ( out1, "w");
    fp2=fopen ( out2, "w");

    
    while ( (c=fgetc(fp))!='>');
    
    
    fprintf ( fp1, ">");
    while ( (c=fgetc(fp))!='>' && c!=EOF){fprintf ( fp1, "%c",c);}
    fclose ( fp1);
    
    fprintf ( fp2, ">");
    while ( (c=fgetc(fp))!='>' && c!=EOF){fprintf ( fp2, "%c",c);}
    fclose ( fp2);
    
    return;
    }
