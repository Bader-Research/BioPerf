/*	May, June 1987	- modified for rapid read of database

	June 2, 1987 - added TFASTA
	March 30, 1988 - combined ffgetaa, fgetgb;
	April 8, 1988 - added PIRLIB format for unix
	Feb 4, 1989 - added universal subroutines for libraries
	November, 1992 - added ncbi search cd support

	copyright (c) 1987,1988,1989,1992 William R. Pearson

	getnt.c	associated subroutines for matching sequences */

/*
8-April-88
	The compile time #define PIRLIB allows this routine to be used
	to read protein and DNA sequence libraries in the NBRF/PIR
	VAX/VMS library format.  That is:

	>P1;LCBO
	This is a line of description
	GTYH ... the sequence starts on this line

	This may ease conversion from UWGCG format libraries. It
	has not been extensively tested.

	In addition, sequence libraries with a '>' in the 4th position
	are recognized as NBRF format libraries for consistency with
	UWGCG

	February 4, 1988 - this starts a major revision of the getaa
	routines.  The goal is to be able to seach the following format
	libraries:

		0 - normal FASTA format
		1 - full Genbank tape format
		2 - NBRF/PIR CODATA format
		3 - EMBL/Swiss-prot format
		4 - Intelligentics format
		5 - NBRF/PIR VMS format
		9 - compressed genbank
		10 - NCBI SEARCH format
		11 - NCBI setdb/blastp (1.3.2) AA

	see file altlib.h to confirm numbers

	This is done with a new global variable and a requirement for the
	FASTLIBS file.  The FASTLIBS file will now indicate both the sequence
	type (protein = 0, DNA = 1) and the file format (the numbers shown
	above, although intelligenetics may become an alternative to Pearson).
	This will be done by always using a function pointer for getlib and
	ranlib(), and setting up a bunch of different getlib() and ranlib()
	functions.  Openlib() will be substantially simplified.
*/

/* 	Nov 12, 1987	- this version checks to see if the sequence
	is DNA or protein by asking whether > 85% is A, C, G, T

	May 5, 1988 - modify the DNA/PROTEIN checker by re-reading
	DNA sequences in order to check for 'U'.
*/

#include <stdio.h>
#include <strings.h>
#include "uascii.gbl"

#ifdef VMS
#define PIRLIB
#endif

#define XTERNAL
#include "upam.gbl"
#undef XTERNAL

#define TRUE 1
#define FALSE 0
#define MAXLINE 512

#define MAXR 15
int lascii[] = {ES, 0, 1, 7,
		 2, 5, 9,13,
		 3, 8, 6,12,
		10,11,14,15};

#define LAMASK 15

getseq(filen,seq,maxs,dnaseq, ced_seq)
	char *filen, *seq;
	int maxs, *dnaseq;
	char *ced_seq;
{
	
	int ced_len=0;
        FILE *fptr;
	char line[512];
	int i, j, n;
	int ic;

	if ((fptr=fopen(filen,"r"))==NULL) {
		fprintf(stderr," could not open %s\n",filen);
		return 0;
		}
	n=0;
	while(fgets(line,512,fptr)!=NULL) {
#ifdef PIRLIB
		if (line[0]=='>'&& (line[3]==';'||line[3]=='>'))
		    fgets(line,512,fptr);
		    
		else
#endif
		if (line[0]!='>'&& line[0]!=';') {
		    for (i=0; (n<maxs)&&
			((ic=sascii[line[i]&AAMASK])<EL); i++)
				if (ic<NA)
				    {
				    seq[n++]= ic;
				    ced_seq[ced_len++]= line[i];
				    }
		    if (ic == ES) break;
		    }
		
		}
	if (n==maxs) {
		fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
		fflush(stderr);
		}
	seq[n]= EOSEQ;
	ced_seq[ced_len]='\0';
	
	if (*dnaseq ==0 && (float)scanseq(seq,n,"ACGT")/(float)n > 0.85) {
		*dnaseq = 1;
				/* convert from protein to DNA sequence */
		sascii = nascii;
		fseek(fptr,0l,0);
		n=0;
		while(fgets(line,512,fptr)!=NULL) {
#ifdef PIRLIB
			if (line[0]=='>'&& (line[3]==';'||line[3]=='>'))
				fgets(line,512,fptr);
			else
#endif
			if (line[0]!='>'&& line[0]!=';') {
			    for (i=0; (n<maxs)&&
				((ic=sascii[line[i]&AAMASK])<EL); i++)
					if (ic<NA)
					    {seq[n++]= ic;
					    
					    }
			    if (ic == ES) break;
			    }
			}
		if (n==maxs) {
			fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
			fflush(stderr);
			}
		seq[n]= EOSEQ;

		sq = nt;
		nsq = nnt;
		hsq = hnt;
		pam = npam;
		strcpy(sqnam,"nt");
		strcpy(sqtype,"DNA");
		}

	fclose(fptr);

	return n;
	}

gettitle(filen,title,len)
	char *filen, *title; int len;
{
	FILE *fptr;
	char line[512];
	char *bp;
	int ll;
#ifdef MSDOS
	char *strpbrk();
#else
	char *strchr();
#endif
	if ((fptr=fopen(filen,"r"))==NULL) {
		fprintf(stderr," file %s was not found\n",filen);
		fflush(stderr);
		return 0;
		}

	while(fgets(line,512,fptr)!=0) {
		if (line[0]=='>'|| line[0]==';') goto found;
		}
	fclose(fptr);
	title[0]='\0';
	return 0;

found:
#ifdef PIRLIB
	if (line[0]=='>'&&(line[3]==';'||line[3]=='>')) {
		if ((bp = strchr(line,'\n'))!=NULL) *bp='\0';
		ll=strlen(line); line[ll++]=' '; line[ll]='\0';
		fgets(&line[ll],512-ll,fptr);
	}
#endif
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

#ifndef VMS
FILE *libf=NULL;
#else
int libf = -1;
#endif
#ifdef NOLIB
int leof = 0;
#endif

long lpos;
char lline[MAXLINE];


#ifndef NOLIB
#include "genbank.h"
#include "altlib.h"
extern int ldnaseq;
int (*getlib)();
int (*ranlib)();
#define GETLIB agetlib
#define RANLIB aranlib
struct slibhdr libhdr;
struct seqhdr namrec;
char seqrec[SQRLEN];
int jsave, recsav;
#else
#define LASTLIB 11
#define BINARYGB 9
#define GETLIB getlib
#define RANLIB ranlib
#endif

/*	the following is from fgetgb.c */

#include <fcntl.h>
#ifndef O_RAW
#ifdef O_BINARY
#define O_RAW O_BINARY
#else
#define O_RAW 0
#endif		/* O_BINARY */
#endif		/* O_RAW */
int libfd= -1;
#ifndef NOLIB
extern int deftype;	/* default library type */
extern int outtty;	/* flag for no interaction */
#ifndef UNIX
#define RBSTR "rb"	/* read file in binary mode */
#else
#define RBSTR "r"
#endif
#else
int deftype=0;
int outtty=1;
#endif
int libtype;		/* current open library type */
int sfnum;		/* superfamily number from types 0 and 5 */

/* a file name for openlib may now include a library type suffix */

openlib(lname,libenv)
	char *lname, *libenv;
{
	char rline[10],libn[120], *strchr(), *bp;
	long ftell();
	int wcnt, ll, opnflg;

	if (lname[0]=='#') return -9;
	wcnt = 0;

#ifndef NOLIB
	if (strlen(libenv)!=0) {
		strncpy(libn,libenv,120);
#ifdef UNIX
		strncat(libn,"/",120);
#endif
		strncat(libn,lname,120-strlen(libn));
		}
	else strncpy(libn,lname,120);
#else
	strncpy(libn,lname,120);
#endif

	/* check for library type */
	if ((bp=strchr(libn,' '))!=NULL) {
	    *bp='\0';
	    sscanf(bp+1,"%d",&libtype);
	    if (libtype<0 || libtype >= LASTLIB) {
		fprintf(stderr," invalid library type: %d (>%d)- resetting\n%s\n",
			libtype,LASTLIB,lname);
		libtype=deftype;
		}
	    }
	else libtype=deftype;

#ifndef NOLIB
	getlib=getliba[libtype];
	ranlib=ranliba[libtype];

l1:	if (libtype<=LASTTXT) opnflg=((libf=fopen(libn,"r"))!=NULL);
	else if (libtype==NCBISRC) opnflg=((libf=fopen(libn,RBSTR))!=NULL);
	else if (libtype==NCBIBL13) opnflg=(ncbl_openlib(libn)!= -1);

	if (!opnflg) {
#else
l1:	if ((libf=fopen(libn,"r"))==NULL) {
#endif
  	   if (outtty) {
		fprintf(stderr," cannot open %s library\n",libn);
		fprintf(stderr," enter new file name or <RET> to quit ");
		fflush(stderr);
		if (fgets(libn,120,stdin)==NULL) return -1;
		if ((bp=strchr(libn,'\n'))!=0) *bp='\0';
		if (strlen(libn)==0) return 0;
		if (++wcnt > 10) return -1;
		goto l1;
	      }
	   else return 0;
	 }
#ifndef NOLIB
	if (libtype<=LASTTXT) {
		lpos = ftell(libf);
		if (fgets(lline,MAXLINE,libf)==NULL) return -1;
		}
	else if (libtype==NCBISRC) {
	  if (ncbi_openlib(deftype)== -1) return -1;
	  ncbi_opendef(libn);
	}
#else		/* NOLIB */
	lpos = ftell(libf);
	if (fgets(lline,MAXLINE,libf)==NULL) return -1;
	leof = 0;
#endif		/* NOLIB */
	return 1;
	}

closelib()
{
  if (libf!=NULL) {
    fclose(libf);
    libf = NULL;
  }
}

GETLIB(seq,maxs,libstr,libpos,lcont)
	char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int ll;
	int ic;
	register char *cp;
	register char *seqp;
	register int *ap;
	char *seqm, *seqm1, *linep, *strchr(), *bp;

	seqp = seq;
	seqm = &seq[maxs-9];
	seqm1 = seqm-1;
#ifndef TFASTA
	ap = sascii;
#else
	ap = nascii;
#endif
	if (*lcont==0) {
#ifndef NOLIB
		while (lline[0]!='>' && lline[0]!=';') {
			lpos = ftell(libf);
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
		}
#ifdef SUPERFAMNUM
		if ((bp=strchr(lline,SFCHAR))!=NULL) {
			*bp='\0';
			sscanf(bp+1,"%d",&sfnum);
		      }
		else sfnum=0;
#else
		sfnum = 0;
#endif
		
		strncpy(libstr,lline+1,20);

		libstr[10]='\0';
		*libpos = lpos;
#else	/* NOLIB */
		if (leof) return 0;
		*libpos = lpos;
		if (lline[0]=='>' || lline[0]==';') {
			strncpy(libstr,lline+1,20);
			libstr[10]='\0';
			}
		else {
			libstr[0]='\0';
			strncpy(seqp,lline,(int)(seqm-seqp));
			for (cp=seqp; seqp<seqm1; ) {
				if ((*seqp++=ap[*cp++])<NA) continue;
				if (*--seqp>NA) break;
				}
			if (*seqp==ES) goto done;
			}
#endif
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets(seqp,(int)(seqm-seqp),libf)!=NULL) {
		if (*seqp=='>') goto new;
		if (*seqp==';') {
			if (strchr(seqp,'\n')==NULL) goto cont;
			continue;
			}
		for (cp=seqp; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		lpos = ftell(libf);
		}
	goto done;
new:	strncpy(lline,seqp,MAXLINE);
	lline[MAXLINE-1]='\0';
	if (strchr(seqp,'\n')==NULL) fgets(lline,MAXLINE-strlen(lline),libf);
	goto done;

cont:
	fgets(lline,MAXLINE,libf);
	seqm1 = seqp;

done:	if (seqp>=seqm1) {
		(*lcont)++;
		}
	else {
#ifdef NOLIB
	leof = 1;
#endif
	*lcont=0;
		}


	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

RANLIB(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	int ll;
	char *strchr();

	fseek(libf, seek, 0);
	fgets(lline,MAXLINE,libf);

	if (lline[0]=='>' || lline[0]==';') {
		strncpy(str,lline+1,cnt);
		str[cnt-1]='\0';
#ifdef SUPERFAMNUM
		if ((bp = strchr(str,SFCHAR))!=NULL) *bp='\0';
		else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
		else str[cnt-1]='\0';
#else
		if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
		else str[cnt-1]='\0';
#endif
		}
	else {
		str[0]='\0';
		}
#ifdef NOLIB
	leof=0;
#endif
	}

#ifndef NOLIB
ggetlib(seq,maxs,libstr,libpos,lcont)
	char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	int i, j, m, n;
	int arec, nrec, rrec, trec;
	register int *lap;
	register char *sptr, *tptr;
	char *ttptr;
	int mc2, tmp;
	char sname[11];
	long lseek();

	sptr = seq;
	lap = lascii;

	if (*lcont==0) {
		*libpos = lseek(libfd,0L,1);
		if (read(libfd,(char *)&namrec,sizeof(namrec))==0) return 0;
		recsav = namrec.rcnt[0] + (namrec.rcnt[1]<<8) - 1;
		for (j=0; j<11; j++) 
			if (!(*sptr++ = namrec.seqstart[j])) goto done;
		maxs -= 22;
		strncpy(sname,namrec.seqnam,10);
		sname[10]='\0';
		strncpy(libstr,sname,20);
		libstr[10]='\0';
		}

	arec = (maxs-2)/(2*SQRLEN);
	rrec = min(arec,recsav);
	if (rrec>0
	    && (trec=read(libfd,(char *)sptr,rrec*SQRLEN))!=rrec*SQRLEN)
		goto error;
	if (rrec == recsav) *lcont=0;
	else {(*lcont)++; recsav -= rrec;}

	sptr += rrec*SQRLEN;

done:
	tptr = ttptr = seq + 2*(int)(sptr - seq);
	while (sptr>seq) {
		*--tptr = lap[*--sptr&LAMASK];
		*--tptr = lap[(*sptr>>4)&LAMASK];
		}

	tptr = ttptr;
	while (*--tptr>MAXR);
	n = (int)(tptr-seq)+1;
	seq[n]= EOSEQ;
	return (n);

error:	fprintf(stderr," error reading %10s %4d %4d %4d\n",sname,trec,rrec*SQRLEN,SQRLEN);
	fflush(stderr);
	return (-1);
	}

extern int ixstat;

granlib(str,cnt,seek)
	char *str; int cnt;
	long seek;
{
	int ctmp;
	long lseek();

	lseek(libfd,seek,0);
	if (read(libfd,(char *)&namrec,sizeof(namrec))==0) return 0;
	strncpy(str,namrec.seqnam,10);
	str[10]='\0';
	if (ixstat>0) idxann(str,&str[10],cnt-10);
	lseek(libfd,seek,0);
	}

char *cpsave;

lgetlib(seq,maxs,libstr,libpos,lcont)
	char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int i, n, ll;
	int ic;
	register char *cp;
	register char *seqp;
	register int *ap;
	char *seqm, *seqm1, *linep, *strchr();

	seqp = seq;
	seqm = &seq[maxs-11];
	seqm1 = seqm-1;
#ifndef TFASTA
	ap = sascii;
#else
	ap = nascii;
#endif
	if (*lcont==0) {
		while (lline[0]!='L' || lline[1]!='O' || 
				strncmp(lline,"LOCUS",5)) { /* find LOCUS */
			lpos = ftell(libf);
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
			}
		strncpy(libstr,&lline[11],20);
		libstr[10]='\0';
		*libpos=lpos;
		while (lline[0]!='O' || lline[1]!='R' ||
				strncmp(lline,"ORIGIN",6)) { /* find ORIGIN */
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
			}
		}
	else {
		for (cp= cpsave; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA) continue;
			if (*(--seqp)>NA) break;
			}
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets(lline,sizeof(lline),libf)!=NULL) {
		if (lline[0]=='/') goto new;
		for (cp= &lline[10]; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		}
	goto done;
new:	lpos = ftell(libf);
	fgets(lline,sizeof(lline),libf);

done:	if (seqp>=seqm1) {
		cpsave = cp;
		(*lcont)++;
		}
	else *lcont=0;

	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

lranlib(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	int ll;
	char *strchr();

	fseek(libf, seek, 0);
	fgets(lline,MAXLINE,libf);

	strncpy(str,&lline[12],10);
	str[10]='\0';
	fgets(lline,sizeof(lline),libf);
	while (lline[0]!='D' || lline[1]!='E' || strncmp(lline,"DEFINITION",10))
		fgets(lline,sizeof(lline),libf);
	strncpy(&str[10],&lline[11],cnt-10);
	str[cnt-1]='\0';
	if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

	fseek(libf,seek,0);
	fgets(lline,MAXLINE,libf);
	}

pgetlib(seq,maxs,libstr,libpos,lcont)
	char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int i, n, ll;
	int ic;
	register char *cp;
	register char *seqp;
	register int *ap;
	char *seqm, *seqm1, *linep, *strchr();

	seqp = seq;
	seqm = &seq[maxs-11];
	seqm1 = seqm-1;
#ifndef TFASTA
	ap = sascii;
#else
	ap = nascii;
#endif
	if (*lcont==0) {
		while (lline[0]!='E' || lline[1]!='N' || strncmp(lline,"ENTRY",5))
		{ /* find ENTRY */
			lpos = ftell(libf);
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
			}
		strncpy(libstr,&lline[16],8);
		libstr[8]='\0';
		*libpos = lpos;
		while (lline[0]!='S' || lline[2]!='Q' || strncmp(lline,"SEQUENCE",8))
		{ /* find SEQUENCE */
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
			}
		fgets(lline,sizeof(lline),libf); /* get the extra line */
		}
	else {
		for (cp= cpsave; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA) continue;
			if (*(--seqp)>NA) break;
			}
		if (*seqp==ES) goto done;
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets(lline,sizeof(lline),libf)!=NULL) {
		if (lline[0]=='/') goto new;
		for (cp= &lline[8]; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA) continue;
			if (*(--seqp)>NA) break;
		      };
		if (*seqp==ES) goto done;
		}
	goto done;
new:	lpos = ftell(libf);
	fgets(lline,sizeof(lline),libf);

done:	if (seqp>=seqm1) {
		cpsave = cp;
		(*lcont)++;
		}
	else *lcont=0;

	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

pranlib(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	int ll;
	char *strchr();

	fseek(libf, seek, 0);
	fgets(lline,MAXLINE,libf);

	strncpy(str,&lline[16],8);
	str[8]='\0';
	fgets(lline,sizeof(lline),libf);
	while (lline[0]!='T' || lline[1]!='I' || strncmp(lline,"TITLE",5))
		fgets(lline,sizeof(lline),libf);
	strncpy(&str[8],&lline[16],cnt-9);
	str[cnt-1]='\0';
	if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

	fseek(libf,seek,0);
	fgets(lline,MAXLINE,libf);
	}

long seqsiz;

egetlib(seq,maxs,libstr,libpos,lcont)
	char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int ll;
	int ic;
	register char *cp;
	register char *seqp;
	register int *ap;
	char *seqm, *seqm1, *linep, *strchr();
	char id[11];  /* Holds Identifier */

	seqp = seq;
	seqm = &seq[maxs-11];
	seqm1 = seqm-1;
#ifndef TFASTA
	ap = sascii;
#else
	ap = nascii;
#endif
	if (*lcont==0) {
		while (lline[0]!='I' || lline[1]!='D') { /* find ID */
			lpos = ftell(libf);
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
			}
		sscanf(&lline[5],"%s",id);
		sprintf(libstr,"%-10.10s",id);
		*libpos = lpos;
		while (lline[0]!='S' || lline[1]!='Q') { /* find ORIGIN */
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
			}
		sscanf(&lline[14],"%ld",&seqsiz);
		}
	else {
		for (cp= cpsave; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA) continue;
			if (*(--seqp)>NA) break;
			}
		if (*seqp==ES) goto done;
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets(lline,sizeof(lline),libf)!=NULL) {
		if (lline[0]=='/') goto new;
		lline[70]='\0';
		for (cp= &lline[5]; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		}
	goto done;
new:	lpos = ftell(libf);
	fgets(lline,sizeof(lline),libf);
	goto done;

done:	if (seqp>=seqm1) {
		cpsave = cp;
		(*lcont)++;
		seqsiz -= (long)(seqp-seq);
		}
	else *lcont=0;

	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
/*	if (*lcont==0 && (long)(seqp-seq)!=seqsiz)
		printf("%s read %d of %d\n",libstr,(int)(seqp-seq),seqsiz);
*/
	return (int)(seqp-seq);
	}

eranlib(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	char id[11];  /* Holds Identifier */
	int ll;
	char *strchr();

	fseek(libf, seek, 0);
	fgets(lline,MAXLINE,libf);

	sscanf(&lline[5],"%s",id);
	sprintf(str,"%-10.10s ",id);
	fgets(lline,sizeof(lline),libf);
	while (lline[0]!='D' || lline[1]!='E') fgets(lline,sizeof(lline),libf);
	strncpy(&str[11],&lline[5],cnt-11);
	str[cnt-1]='\0';
	if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
	if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

	fseek(libf,seek,0);
	fgets(lline,MAXLINE,libf);
	}

igetlib(seq,maxs,libstr,libpos,lcont)
	char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int i, n, ll;
	int ic;
	register char *cp;
	register char *seqp;
	register int *ap;
	char *seqm, *seqm1, *linep, *bp, *strchr();

	seqp = seq;
	seqm = &seq[maxs-9];
	seqm1 = seqm-1;
#ifndef TFASTA
	ap = sascii;
#else
	ap = nascii;
#endif
	if (*lcont==0) {
		while (lline[0]!=';') {
			lpos = ftell(libf);
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
			}
		*libpos = lpos;
		while (lline[0]==';') fgets(lline,sizeof(lline),libf);
		strncpy(libstr,lline+1,10);
		libstr[9]='\0';
		if((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets(seqp,(int)(seqm-seqp),libf)!=NULL) {
		if (*seqp=='>') goto new;
		if (*seqp==';') {
			if (strchr(seqp,'\n')==NULL) goto cont;
			continue;
			}
		for (cp=seqp; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		lpos = ftell(libf);
		}
	goto done;
new:	strncpy(lline,seqp,MAXLINE);
	lline[MAXLINE-1]='\0';
	if (strchr(seqp,'\n')==NULL) fgets(lline,MAXLINE-strlen(lline),libf);
	goto done;

cont:
	fgets(lline,MAXLINE,libf);
	seqm1 = seqp;

done:	if (seqp>=seqm1) {
		(*lcont)++;
		}
	else {
	*lcont=0;
		}


	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

iranlib(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	int ll;
	char *strchr();
	char tline[120];

	fseek(libf, seek, 0);
	fgets(lline,MAXLINE,libf);

	if (lline[0]=='>' || lline[0]==';') {
		strncpy(tline,lline+1,sizeof(tline));
		str[cnt-1]='\0';
		if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
		else str[cnt-1]='\0';
		}
	else {
		tline[0]='\0';
		}

	while (lline[0]==';') fgets(lline,sizeof(lline),libf);
	if ((bp=strchr(lline,'\n'))!=NULL) *bp=0;
	if ((bp=strchr(lline,' '))!=NULL) *bp=0;
	strncpy(str,lline,cnt);
	strncat(str,"  ",cnt-strlen(str));
	strncat(str,tline,cnt-strlen(str));
	str[cnt-1]='\0';
	
	fseek(libf,seek,0);
	fgets(lline,MAXLINE,libf);
	}

vgetlib(seq,maxs,libstr,libpos,lcont)
	char *seq;
	int maxs;
	char *libstr;
	long *libpos;
	int *lcont;
{
	long ftell();
	int i, n, ll;
	int ic;
	register char *cp;
	register char *seqp;
	register int *ap;
	char *seqm, *seqm1, *linep, *strchr(), *bp;

	seqp = seq;
	seqm = &seq[maxs-9];
	seqm1 = seqm-1;
#ifndef TFASTA
	ap = sascii;
#else
	ap = nascii;
#endif
	if (*lcont==0) {
		while (lline[0]!='>' && lline[0]!=';') {
			lpos = ftell(libf);
			if (fgets(lline,MAXLINE,libf)==NULL) return 0;
		}
		if ((bp=strchr(lline,SFCHAR))!=NULL) {
			*bp='\0';
			sscanf(bp+1,"%d",&sfnum);
		      }
		else sfnum=0;
		if ((bp=strchr(lline,'\n'))!=NULL) *bp='\0';
		strncpy(libstr,&lline[4],20);
		fgets(lline,MAXLINE,libf);
		libstr[10]='\0';
		*libpos = lpos;
		}

	lline[0]='\0';
	while (seqp<seqm1 && fgets(seqp,(int)(seqm-seqp),libf)!=NULL) {
		if (*seqp=='>') goto new;
		if (*seqp==';') {
			if (strchr(seqp,'\n')==NULL) goto cont;
			continue;
			}
		for (cp=seqp; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		lpos = ftell(libf);
		}
	goto done;
new:	strncpy(lline,seqp,MAXLINE);
	lline[MAXLINE-1]='\0';
	if (strchr(seqp,'\n')==NULL) fgets(lline,MAXLINE-strlen(lline),libf);
	goto done;

cont:
	fgets(lline,MAXLINE,libf);
	seqm1 = seqp;

done:	if (seqp>=seqm1) {
		(*lcont)++;
		}
	else {

	*lcont=0;
		}


	*seqp = EOSEQ;
	if ((int)(seqp-seq)==0) return 1;
	return (int)(seqp-seq);
	}

vranlib(str,cnt,seek)
	char *str; int cnt; long seek;
{
	char *bp;
	int ll;
	char *strchr();

	fseek(libf, seek, 0);
	fgets(lline,MAXLINE,libf);

	if (lline[0]=='>'&&(lline[3]==';'||lline[3]=='>')) {
		strncpy(str,&lline[4],cnt);

		if ((bp = strchr(str,':'))!=NULL) *bp='\0';
		if ((bp=strchr(str,'\r'))!=NULL) *bp='\0';
		else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
		else str[cnt-1]='\0';

		fgets(lline,MAXLINE,libf);
		if ((bp=strchr(lline,'\r'))!=NULL) *bp='\0';
		if ((bp=strchr(lline,'\n'))!=NULL) *bp='\0';
		strncat(str," ",cnt);
		strncat(str,lline,cnt-strlen(str));
		}
	else {
		str[0]='\0';
		}

	fseek(libf,seek,0);
	fgets(lline,MAXLINE,libf);
	}

#endif	/* NOLIB */

scanseq(seq,n,str)
	char *seq, *str;
	int n;
{
	int tot,i;
	char aaray[MAXSQ];		/* this must be set > nsq */
	
	for (i=0; i<MAXSQ; i++)  aaray[i]=0;
	for (i=0; i<strlen(str); i++) aaray[sascii[str[i]]]=1;
	for (i=tot=0; i<n; i++) tot += aaray[seq[i]];
	return tot;
	}

revcomp(seq,n)
	char *seq; int n;
{
  char tmp;
  int i, ni;

  for (i=0; i< n; i++)
    if (sq[seq[i]]=='A') seq[i] = nascii['T'];
    else if (sq[seq[i]]=='C') seq[i] = sascii['G'];
    else if (sq[seq[i]]=='G') seq[i] = sascii['C'];
    else if (sq[seq[i]]=='T') seq[i] = sascii['A'];
    else if (sq[seq[i]]=='R') seq[i] = sascii['Y'];
    else if (sq[seq[i]]=='Y') seq[i] = sascii['R'];
    else if (sq[seq[i]]=='M') seq[i] = sascii['K'];
    else if (sq[seq[i]]=='K') seq[i] = sascii['M'];
    else if (sq[seq[i]]=='D') seq[i] = sascii['H'];
    else if (sq[seq[i]]=='H') seq[i] = sascii['D'];
    else if (sq[seq[i]]=='V') seq[i] = sascii['B'];
    else if (sq[seq[i]]=='B') seq[i] = sascii['V'];

  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = seq[i];
    seq[i] = seq[ni];
    seq[ni] = tmp;
  }
}

max(arg1,arg2)
	int arg1, arg2;
{
	return (arg1>arg2) ? arg1 : arg2;
	}

#ifdef VMS
memcpy(ar0, ar1, n)
	char *ar0, *ar1; unsigned n;
{
	while (n--) *ar0++ = *ar1++;
	}

openidx() {}
/* newname generates a new filename with prefix oname and suffix suff */

newname(nname,oname,suff,maxn)
	char *nname, *oname, *suff;
	int maxn;
{
	char *tptr;
	if (*oname!='@') strncpy(nname,oname,maxn);
	else strncpy(nname,oname+1,maxn);

	for (tptr=nname; *tptr!='.'&& *tptr; tptr++); /* get to '.' or EOS */
	*tptr++='.'; *tptr='\0';
	strncat(nname,suff,maxn);
	}
#endif

#ifdef MACLSC
memcpy(ar0, ar1, n)
	char *ar0, *ar1; unsigned n;
{
	while (n--) *ar0++ = *ar1++;
	}
#endif
