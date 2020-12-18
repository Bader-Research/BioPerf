/* copyright (c) 1998, 1999 William R. Pearson and the U. of Virginia */

/* $Name: fa34t25d2 $ - $Id: dropfs.c,v 1.20 2004/11/10 21:06:03 wrp Exp $ */

/* this code implements the "fasts" algorithm, which compares a set of
   protein fragments to a protein sequence.  Comma's are used to separate
   the sequence fragments, which need not be the same length.

   The expected input is:

   >mgstm1
   MGDAPDFD,
   MILGYW,
   MLLEYTDS

   The fragments do not need to be in the correct order (which is
   presumably unknown from the peptide sequencing.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"
/*
#define XTERNAL
#include "upam.h"
*/
#define EOSEQ 0

/* globals for fasta */
#define MAXWINDOW 64

#define ESS 49
#define NMAP_X 23 /* for 'X' */
#define NMAP_Z 24 /* for '*' */
#define MAXHASH 32
#define NMAP MAXHASH+1

#ifndef MAXSAV
#define MAXSAV 25
#endif
#define MAX_SEG 25

static char *verstr="3.36 June 2000";
extern void w_abort(char *, char *);
int shscore(const unsigned char *aa0, const int n0, int **pam2, int nsq);

#ifdef TFASTS
extern int aatran(const unsigned char *ntseq, unsigned char *aaseq, const int maxs, const int frame);
#endif

struct dstruct	{	/* diagonal structure for saving current run */
   int     score;	/* hash score of current match */
   int     start;	/* start of current match */
   int     stop;	/* end of current match */
   struct savestr *dmax;   /* location in vmax[] where best score data saved */
};

struct savestr {
   int     score;		/* pam score with segment optimization */
   int     score0;		/* pam score of best single segment */
   int     gscore;		/* score from global match */
   int     dp;			/* diagonal of match */
   int     start;		/* start of match in lib seq */
   int     stop;		/* end of match in lib seq */
};

struct f_struct {
  struct dstruct *diag;
  struct savestr vmax[MAXSAV];	/* best matches saved for one sequence */
  struct savestr *vptr[MAXSAV];
  struct savestr *lowmax;
  int nsave;
  int ndo;
  int noff;
  int nm0, nmoff[MAX_SEG];		/* offset number, start */
  int *aa0b, *aa0e, *aa0i;
  unsigned char *aa0t;
  int hmask;			/* hash constants */
  int *pamh1;			/* pam based array */
  int *pamh2;			/* pam based kfact array */
  int *link, *harr;		/* hash arrays */
  int kshft;			/* shift width */
  int nsav, lowscor;		/* number of saved runs, worst saved run */
#ifdef TFASTS
  unsigned char *aa1x;
  int n10;
#endif
  struct bdstr *bss;
  struct swstr *ss;
  struct swstr *r_ss;
  int *waa;
  int *res;
  int max_res;
};

#ifdef ALLOCN0
void savemax (struct dstruct *, int, struct f_struct *);
#else
void savemax (struct dstruct *, struct f_struct *);
#endif

int spam(const unsigned char *, const unsigned char *, int, struct savestr *, int **, struct f_struct *);
int sconn(struct savestr **, int nsave, int cgap, int pgap, struct f_struct *);
void kpsort(struct savestr **, int);
void kssort(struct savestr **, int);
int sconn_a(unsigned char *, int, struct f_struct *, int cgap, int pgap, struct a_struct *);
void kpsort(struct savestr **, int);


/* initialize for fasta */

void
init_work (unsigned char *aa0, const int n0, 
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
  int mhv;
   int hmax;
   int i0, ib, hv;
   int pamfact;
   int btemp;
   struct f_struct *f_str;
   /* these used to be globals, but do not need to be */
   int ktup, fact, kt1;

   int maxn0;
   int *pwaa;
   int i, j, q;
   struct bdstr *bss;
   struct swstr *ss, *r_ss;
   int *waa;
   int *res;

   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   btemp = shscore(aa0,n0,ppst->pam2[0],ppst->nsq) / 3;
   ppst->param_u.fa.cgap = btemp;

   ppst->param_u.fa.pgap = ppst->gdelval + ppst->ggapval;
   pamfact = ppst->param_u.fa.pamfact;
   ppst->param_u.fa.ktup = ktup = 1;

   if (pamfact == -1) pamfact = 0;
   else if (pamfact == -2) pamfact = 1;

   /* fasts3 cannot work with lowercase symbols as low complexity;
      thus, NMAP must be disabled; this depends on aascii['X']  */
   if (ppst->hsq[NMAP_X] == NMAP ) {ppst->hsq[NMAP_X]=1;}
   if (ppst->hsq[NMAP_Z] == NMAP ) {ppst->hsq[NMAP_Z]=1;}
   /* this does not work in a threaded environment */
   /*    else {fprintf(stderr," cannot find 'X'==NMAP\n");} */

   for (i0 = 1, mhv = -1; i0 <= ppst->nsq; i0++)
      if (ppst->hsq[i0] < NMAP && ppst->hsq[i0] > mhv)  mhv = ppst->hsq[i0];

   if (mhv <= 0) {
      fprintf (stderr, " maximum hsq <=0 %d\n", mhv);
      exit (1);
   }

   for (f_str->kshft = 0; mhv > 0; mhv /= 2)
      f_str->kshft++;

/*      kshft = 2;	*/
   kt1 = 0;
   hv = 1;
   for (i0 = 0; i0 < ktup; i0++) hv = hv << f_str->kshft;
   hmax = hv;
   f_str->hmask = (hmax >> f_str->kshft) - 1;

   if ((f_str->aa0t = (unsigned char *) calloc(n0+1, sizeof(char))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0t array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0b = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0b array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0e = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0e array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0i = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0i array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->harr = (int *) calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate hash array: hmax: %d hmask: %d\n",
	      hmax, f_str->hmask);
     exit (1);
   }
   if ((f_str->pamh1 = (int *) calloc (ppst->nsq+1, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate pamh1 array\n");
     exit (1);
   }
   if ((f_str->pamh2 = (int *) calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate pamh2 array\n");
     exit (1);
   }
   if ((f_str->link = (int *) calloc (n0, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate hash link array");
     exit (1);
   }

   for (i0 = 0; i0 < hmax; i0++) f_str->harr[i0] = -1;
   for (i0 = 0; i0 < n0; i0++) f_str->link[i0] = -1;

   /* encode the aa0 array */
   /* no initial loop because kt1 = 0 */

   hv = 0;
   f_str->nmoff[0] = 0;
   f_str->nm0 = 1;
   for (i0 = 0; i0 < n0; i0++)
   {
     if (aa0[i0] == ESS || aa0[i0] == 0) {
       aa0[i0] = EOSEQ;	/* set ESS to 0 */
       /*       fprintf(stderr," converted ',' to 0\n");*/
       f_str->nmoff[f_str->nm0++] = i0+1; 
       if (f_str->nm0 >= MAX_SEG) {
	 fprintf(stderr," too many fragments (> %d)\n",MAX_SEG);
	 exit(1);
       }
       hv = 0;
       continue;
     }

      hv = ((hv & f_str->hmask) << f_str->kshft) + ppst->hsq[aa0[i0]];
      f_str->link[i0] = f_str->harr[hv];
      f_str->harr[hv] = i0;
      f_str->pamh2[hv] = ppst->pam2[0][aa0[i0]][aa0[i0]];
   }
   f_str->nmoff[f_str->nm0] = n0+1;
   ppst->param_u.fa.cgap /= f_str->nm0;

   /* setup aa0b[], aa0e[], which specify the begining and end of each
      segment */

   for (ib = i0 = 0; i0 < n0; i0++) {
     if (aa0[i0]==EOSEQ) ib++;
     f_str->aa0b[i0] =  f_str->nmoff[ib];
     f_str->aa0e[i0] =  f_str->nmoff[ib+1]-2;
     f_str->aa0i[i0] =  ib;
     /*
     fprintf(stderr,"%2d %c: %2d %2d %2d\n",i0,ppst->sq[aa0[i0]],
	     f_str->aa0b[i0],f_str->aa0e[i0],f_str->aa0i[i0]);
     */
   }

   /*
   for (i0=1; i0<=ppst->nsq; i0++) {
     fprintf(stderr," %c: %2d ",ppst->sq[i0],f_str->harr[i0]);
     hv = f_str->harr[i0];
     while (hv >= 0) {
       fprintf(stderr," %2d",f_str->link[hv]);
       hv = f_str->link[hv];
     }
     fprintf(stderr,"\n");
   }
   */

/* this has been modified from 0..<ppst->nsq to 1..<=ppst->nsq because the
   pam2[0][0] is now undefined for consistency with blast
*/
   for (i0 = 1; i0 <= ppst->nsq; i0++)
     f_str->pamh1[i0] = ppst->pam2[0][i0][i0];

   f_str->ndo = 0;
   f_str->noff = n0-1;
#ifndef ALLOCN0
   if (f_str->diag==NULL) 
     f_str->diag = (struct dstruct *) calloc ((size_t)MAXDIAG,
					      sizeof (struct dstruct));
#else
   if (f_str->diag==NULL) 
     f_str->diag = (struct dstruct *) calloc ((size_t)n0,
					      sizeof (struct dstruct));
#endif

   if (f_str->diag == NULL) {
      printf (" cannot allocate diagonal arrays: %ld\n",
	      (long) MAXDIAG * (long) (sizeof (struct dstruct)));
      exit (1);
   }

#ifdef TFASTS
   if ((f_str->aa1x =(unsigned char *)calloc((size_t)ppst->maxlen+2,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "cannot allocate aa1x array %d\n", ppst->maxlen+2);
     exit (1);
   }
   f_str->aa1x++;
#endif

   if ((waa= (int *)calloc ((size_t)(ppst->nsq+1)*n0,sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate waa struct %3d\n",ppst->nsq*n0);
     exit(1);
   }

   pwaa = waa;
   for (i=0; i<ppst->nsq+1; i++) {
     for (j=0;j<n0; j++) {
       *pwaa = ppst->pam2[0][i][aa0[j]];
       pwaa++;
     }
   }
   f_str->waa = waa;

   maxn0 = max(3*n0/2,MIN_RES);
   if ((res = (int *)calloc((size_t)maxn0,sizeof(int)))==NULL) {
     fprintf(stderr,"cannot allocate alignment results array %d\n",maxn0);
     exit(1);
   }
   f_str->res = res;
   f_str->max_res = maxn0;

   *f_arg = f_str;
}


/* pstring1 is a message to the manager, currently 512 */
/* pstring2 is the same information, but in a markx==10 format */
void
get_param (struct pstruct *pstr, char *pstring1, char *pstring2)
{
#ifndef TFASTS
  char *pg_str="FASTS";
#else
  char *pg_str="TFASTS";
#endif

  sprintf (pstring1, "%s (%s) function [%s matrix (%d:%d)] ktup: %d\n join: %d, gap-pen: %d/%d, width: %3d",pg_str,verstr,
	       pstr->pamfile, pstr->pam_h,pstr->pam_l, pstr->param_u.fa.ktup, pstr->param_u.fa.cgap,
	       pstr->gdelval, pstr->ggapval, pstr->param_u.fa.optwid);
  if (pstr->param_u.fa.iniflag) strcat(pstring1," init1");
  /*
  if (pstr->zsflag==0) strcat(pstring1," not-scaled");
  else if (pstr->zsflag==1) strcat(pstring1," reg.-scaled");
  */
  if (pstring2 != NULL) {
    sprintf (pstring2, "; pg_name: %s\n; pg_ver: %s\n; pg_matrix: %s (%d:%d)\n\
; pg_gap-pen: %d %d\n; pg_ktup: %d\n; pg_cgap: %d\n",
	     pg_str,verstr,pstr->pamfile, pstr->pam_h,pstr->pam_l, pstr->gdelval,
	     pstr->ggapval,pstr->param_u.fa.ktup,pstr->param_u.fa.cgap);
   }
}

void
close_work (const unsigned char *aa0, const int n0,
	    struct pstruct *ppst,
	    struct f_struct **f_arg)
{
  struct f_struct *f_str;

  f_str = *f_arg;

  if (f_str != NULL) {

    free(f_str->res);
    free(f_str->waa);
    free(f_str->diag);
    free(f_str->link);
    free(f_str->pamh2); 
    free(f_str->pamh1);
    free(f_str->harr);
    free(f_str->aa0i);
    free(f_str->aa0e);
    free(f_str->aa0b);
    free(f_str->aa0t);
    free(f_str);
    *f_arg = NULL;
  }
}

void do_fasta (const unsigned char *aa0, const int n0,
	      const unsigned char *aa1, const int n1,
	      struct pstruct *ppst, struct f_struct *f_str,
	      struct rstruct *rst, int *hoff)
{
   int     nd;		/* diagonal array size */
   int     lhval;
   int     kfact;
   register struct dstruct *dptr;
   register int tscor;
#ifndef ALLOCN0
   register struct dstruct *diagp;
#else
   register int dpos;
   int     lposn0;
#endif
   struct dstruct *dpmax;
   register int lpos;
   int     tpos;
   struct savestr *vmptr;
   int     scor, tmp;
   int     im, ib, nsave;
   int     cmps ();		/* comparison routine for ksort */
   int ktup;
   int doffset;
   int nm_u[MAX_SEG];

   ktup = 1;

   if (n1 < ktup) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     return;
   }

   if (n0+n1+1 >= MAXDIAG) {
     fprintf(stderr,"n0,n1 too large: %d, %d\n",n0,n1);
     rst->score[0] = rst->score[1] = rst->score[2] = -1;
     return;
   }

#ifdef ALLOCN0
   nd = n0;
#else
   nd = n0 + n1;
#endif

   dpmax = &f_str->diag[nd];
   for (dptr = &f_str->diag[f_str->ndo]; dptr < dpmax;)
   {
      dptr->stop = -1;
      dptr->dmax = NULL;
      dptr++->score = 0;
   }

   for (vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++)
      vmptr->score = 0;
   f_str->lowmax = f_str->vmax;
   f_str->lowscor = 0;


   /* start hashing */
   lhval = lpos = 0;

#ifndef ALLOCN0
   diagp = &f_str->diag[f_str->noff];
   for (; lpos < n1; lpos++, diagp++) {
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval]; tpos >= 0; tpos = f_str->link[tpos]) {
       if ((tscor = (dptr = &diagp[-tpos])->stop) >= 0) {
#else
   lposn0 = f_str->noff + lpos;
   for (; lpos < n1; lpos++, lposn0++) {
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval]; tpos >= 0; tpos = f_str->link[tpos]) {
       dpos = lposn0 - tpos;
       if ((tscor = (dptr = &f_str->diag[dpos % nd])->stop) >= 0) {
#endif
	 tscor++;
	 if ((tscor -= lpos) <= 0) {
	   scor = dptr->score;
	   if ((tscor += (kfact = f_str->pamh2[lhval])) < 0 && f_str->lowscor < scor)
#ifdef ALLOCN0
	     savemax (dptr, dpos, f_str);
#else
	     savemax (dptr, f_str);
#endif
	     if ((tscor += scor) >= kfact) {
	       dptr->score = tscor;
	       dptr->stop = lpos;
	     }
	     else {
	       dptr->score = kfact;
	       dptr->start = dptr->stop = lpos;
	     }
	 }
	 else {
	   dptr->score += f_str->pamh1[aa0[tpos]];
	   dptr->stop = lpos;
	 }
       }
       else {
	 dptr->score = f_str->pamh2[lhval];
	 dptr->start = dptr->stop = lpos;
       }
     }				/* end tpos */

#ifdef ALLOCN0
      /* reinitialize diag structure */

     if ((dptr = &f_str->diag[lpos % nd])->score > f_str->lowscor)
       savemax (dptr, lpos, f_str);
     dptr->stop = -1;
     dptr->dmax = NULL;
     dptr->score = 0;
#endif
   }				/* end lpos */

#ifdef ALLOCN0
   for (tpos = 0, dpos = f_str->noff + n1 - 1; tpos < n0; tpos++, dpos--) {
     if ((dptr = &f_str->diag[dpos % nd])->score > f_str->lowscor)
       savemax (dptr, dpos, f_str);
   }
#else
   for (dptr = f_str->diag; dptr < dpmax;) {
     if (dptr->score > f_str->lowscor) savemax (dptr, f_str);
     dptr->stop = -1;
     dptr->dmax = NULL;
     dptr++->score = 0;
   }
   f_str->ndo = nd;
#endif

/*
        at this point all of the elements of aa1[lpos]
        have been searched for elements of aa0[tpos]
        with the results in diag[dpos]
*/

   for (nsave=0, vmptr=f_str->vmax; vmptr< &f_str->vmax[MAXSAV]; vmptr++) {
      if (vmptr->score > 0) {
	/*
	fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d",
	     f_str->noff+vmptr->start-vmptr->dp,
	     f_str->noff+vmptr->stop-vmptr->dp,
	     vmptr->start,vmptr->stop,
	     vmptr->dp,vmptr->score);
	*/
	 vmptr->score = spam (aa0, aa1, n1, vmptr, ppst->pam2[0], f_str);
	 /*
	 fprintf(stderr," spam: %d (%d-%d)\n",vmptr->score,vmptr->start,vmptr->stop);
	 */
	 if (vmptr->score > 0) f_str->vptr[nsave++] = vmptr;
      }
   }

   if (nsave <= 0) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     f_str->nsave = 0;
     return;
   }

   /*
   fprintf(stderr,"n0: %d; n1: %d; noff: %d\n",n0,n1,f_str->noff);
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	     f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
	     f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
	     f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
	     f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
   }

   fprintf(stderr,"---\n");
   */

   kssort(f_str->vptr,nsave);

   /* make certain each seg is used only once */

   for (ib=0; ib<f_str->nm0; ib++) nm_u[ib]=0;
   for (ib=0; ib < nsave; ib++) {
     doffset = f_str->vptr[ib]->dp - f_str->noff;
     tpos=f_str->aa0i[f_str->vptr[ib]->start - doffset];
     if (nm_u[tpos] == 0) nm_u[tpos]=1;
     else f_str->vptr[ib]->score = -1;
   }

   scor = sconn (f_str->vptr, nsave, ppst->param_u.fa.cgap, 
		 ppst->param_u.fa.pgap, f_str);

   kssort(f_str->vptr,nsave);

   /* here we should use an nsave that is consistent with sconn and nm0 */

   f_str->nsave = nsave;
   if (nsave > f_str->nm0) f_str->nsave = f_str->nm0;

   rst->score[1] = f_str->vptr[0]->score;
   rst->score[0] = max (scor, f_str->vptr[0]->score);
   rst->score[2] = rst->score[0];
}

void do_work (const unsigned char *aa0, const int n0,
	      const unsigned char *aa1, const int n1,
	      int frame,
	      struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, struct rstruct *rst)
{
  int hoff, n10, i;

  rst->escore = 1.0;
  rst->segnum = rst->seglen = 1;

  if (n1 < ppst->param_u.fa.ktup) {
    rst->score[0] = rst->score[1] = rst->score[2] = -1;
    return;
  }
#ifdef TFASTS
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  if (ppst->debug_lib)
    for (i=0; i<n10; i++)
      if (f_str->aa1x[i]>ppst->nsq) {
	fprintf(stderr,
		"residue[%d/%d] %d range (%d)\n",i,n1,
		f_str->aa1x[i],ppst->nsq);
	f_str->aa1x[i]=0;
	n10=i-1;
      }

  do_fasta (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff);
#else	/* FASTA */
  do_fasta (aa0, n0, aa1, n1, ppst, f_str, rst, &hoff);
#endif

}

void do_opt (const unsigned char *aa0, const int n0,
	     const unsigned char *aa1, const int n1,
	     int frame,
	     struct pstruct *ppst,
	     struct f_struct *f_str,
	     struct rstruct *rst)
{
  int optflag, tscore, hoff, n10;

  optflag = ppst->param_u.fa.optflag;
  ppst->param_u.fa.optflag = 1;

#ifdef TFASTS  
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  do_fasta (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff);
#else	/* FASTA */
  do_fasta(aa0,n0,aa1,n1,ppst,f_str,rst, &hoff);
#endif
  ppst->param_u.fa.optflag = optflag;
}

#ifdef ALLOCN0
void
savemax (dptr, dpos, f_str)
  register struct dstruct *dptr;
  int  dpos;
  struct f_struct *f_str;
{
   register struct savestr *vmptr;
   register int i;

#else
void 
savemax (dptr, f_str)
  register struct dstruct *dptr;
  struct f_struct *f_str;
{
   register int dpos;
   register struct savestr *vmptr;
   register int i;

   dpos = (int) (dptr - f_str->diag);

#endif

/* check to see if this is the continuation of a run that is already saved */

   if ((vmptr = dptr->dmax) != NULL && vmptr->dp == dpos &&
	 vmptr->start == dptr->start)
   {
      vmptr->stop = dptr->stop;
      if ((i = dptr->score) <= vmptr->score)
	 return;
      vmptr->score = i;
      if (vmptr != f_str->lowmax)
	 return;
   }
   else
   {
      i = f_str->lowmax->score = dptr->score;
      f_str->lowmax->dp = dpos;
      f_str->lowmax->start = dptr->start;
      f_str->lowmax->stop = dptr->stop;
      dptr->dmax = f_str->lowmax;
   }

   for (vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++)
      if (vmptr->score < i)
      {
	 i = vmptr->score;
	 f_str->lowmax = vmptr;
      }
   f_str->lowscor = i;
}

/* this version of spam scans the diagonal to find the best local score,
   then resets the boundaries for a global alignment and re-scans */

int spam (const unsigned char *aa0, const unsigned char *aa1,int n1,
	  struct savestr *dmax, int **pam2,
	  struct f_struct *f_str)
{
   int     lpos, doffset;
   int     tot, mtot;
   struct {
     int     start, stop, score;
   } curv, maxv;
   register const unsigned char *aa0p, *aa1p;

   doffset = dmax->dp - f_str->noff;
   aa1p = &aa1[lpos = dmax->start];
   aa0p = &aa0[lpos - doffset];
   curv.start = lpos;

   tot = curv.score = maxv.score = 0;
   for (; lpos <= dmax->stop; lpos++) {
     tot += pam2[*aa0p++][*aa1p++];
     if (tot > curv.score) {
       curv.stop = lpos;	/* here, curv.stop is actually curv.max */
       curv.score = tot;
      }
      else if (tot < 0) {
	if (curv.score > maxv.score) {
	  maxv.start = curv.start;
	  maxv.stop = curv.stop;
	  maxv.score = curv.score;
	}
	tot = curv.score = 0;
	curv.start = lpos+1;
      }
   }

   if (curv.score > maxv.score) {
     maxv.start = curv.start;
     maxv.stop = curv.stop;
     maxv.score = curv.score;
   }

   if (maxv.score <= 0) return 0;

   /* now, reset the boundaries of the alignment using aa0b[]
      and aa0e[], which specify the residues that start and end
      the segment */
      
   maxv.start = f_str->aa0b[maxv.stop-doffset] + doffset;
   if (maxv.start < 0) maxv.start = 0;
   maxv.stop = f_str->aa0e[maxv.stop-doffset] + doffset;
   if (maxv.stop > n1) maxv.stop = n1-1;
   aa1p = &aa1[lpos = maxv.start];
   aa0p = &aa0[lpos - doffset];
   for (tot=0; lpos <= maxv.stop; lpos++) tot += pam2[*aa0p++][*aa1p++];
   maxv.score = tot;

/*	if (maxv.start != dmax->start || maxv.stop != dmax->stop)
		printf(" new region: %3d %3d %3d %3d\n",maxv.start,
			dmax->start,maxv.stop,dmax->stop);
*/
   dmax->start = maxv.start;
   dmax->stop = maxv.stop;

   return maxv.score;
}

int sconn (struct savestr **v, int n, 
       int cgap, int pgap, struct f_struct *f_str)
{
   int     i, si, cmpp ();
   struct slink
   {
      int     score;
      struct savestr *vp;
      struct slink *next;
   }      *start, *sl, *sj, *so, sarr[MAXSAV];
   int     lstart, ltmp, tstart, plstop, ptstop, ptstart, tstop;

   pgap = 0;	/* the data should have gaps */

   /*  sort the score left to right in lib pos */
   kpsort (v, n);

   start = NULL;

   /*  for the remaining runs, see if they fit */

   for (i = 0, si = 0; i < n; i++) {
/*	if the score is less than the gap penalty, it never helps */
      if (v[i]->score < cgap)  continue;

      /* the segment is worth adding; find out where? */
      lstart = v[i]->start;
      ltmp = v[i]->stop;
      tstart = lstart - v[i]->dp + f_str->noff;
      tstop = ltmp - v[i]->dp + f_str->noff;

/*	put the run in the group */
      sarr[si].vp = v[i];
      sarr[si].score = v[i]->score;
      sarr[si].next = NULL;

/*  if it fits, then increase the score

    start points to the highest scoring run
   -> next is the second highest, etc.
   put the segment into the highest scoring run that it fits into
*/
      for (sl = start; sl != NULL; sl = sl->next) {
	 ltmp = sl->vp->start;
	 plstop = sl->vp->stop;
	 ptstart = ltmp - sl->vp->dp + f_str->noff;
	 ptstop = plstop - sl->vp->dp + f_str->noff;
	 if (plstop < lstart && ( ptstop < tstart || ptstart > tstop)) {
	    sarr[si].score = sl->score + v[i]->score + pgap;
	    /*
	    fprintf(stderr,"sconn %d added %d/%d getting %d; si: %d\n",
		    i,v[i]->start, v[i]->score,sarr[si].score,si);
	    */
	    break;
	 }
      }

/*	now recalculate where the score fits */
      if (start == NULL)
	start = &sarr[si];
      else
	 for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	    if (sarr[si].score > sj->score) {
	       sarr[si].next = sj;
	       if (so != NULL)
		  so->next = &sarr[si];
	       else
		  start = &sarr[si];
	       break;
	    }
	    so = sj;
	 }
      si++;
   }

   if (start != NULL) return (start->score);
   else return (0);
}

void
kssort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->score >= v[j + gap]->score)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

void
kpsort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->start <= v[j + gap]->start)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

/* calculate the 100% identical score */
shscore(const unsigned char *aa0, const int n0, int **pam2, int nsq)
{
  int i, sum;
  for (i=0,sum=0; i<n0; i++)
    if (aa0[i] != EOSEQ && aa0[i]<=nsq) sum += pam2[aa0[i]][aa0[i]];
  return sum;
}

/* sorts alignments from right to left (back to front) based on stop */

void
krsort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->stop > v[j + gap]->stop)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

int  do_walign (unsigned char *aa0, int n0,
		unsigned char *aa1, int n1,
		int frame,
		struct pstruct *ppst, 
		struct f_struct *f_str, 
		int **ares, int *nres, struct a_struct *aln)
{
  int hoff, n10;
  struct rstruct rst;
  int ib, i;
  unsigned char *aa0t, *aa1p;
  struct savestr *vmptr;

#ifdef TFASTS
  f_str->n10 = n10 = aatran(aa1,f_str->aa1x,n1,frame);
  aa1p = f_str->aa1x;
  aln->qlrev = 0;
  aln->qlfact = 1;
  aln->llfact = aln->llmult = 3;
  if (frame > 3) aln->llrev = 1;
  aln->frame = 0;
#else
  aln->llfact = aln->llmult = aln->qlfact = 1;
  aln->llrev = aln->qlrev = 0;
  aln->frame = 0;
  n10 = n1;
  aa1p = aa1;
#endif

  /* this should be modified to do a global/local alignment for each
     of the query fragments, and then optimally align those alignments
     - currently, very long fragments can mask shorter identical
     alignments */

  do_fasta(aa0, n0, aa1p, n10, ppst, f_str, &rst, &hoff);

  /* the alignment portion takes advantage of the information left
     over in f_str after do_fasta is done.  in particular, it is
     easy to run a modified sconn() to produce the alignments.

     unfortunately, the alignment display routine wants to have
     things encoded as with bd_align and sw_align, so we need to do that.  */

  if ((aa0t = (unsigned char *)calloc(n0+1,sizeof(unsigned char)))==NULL) {
    fprintf(stderr," cannot allocate aa0t %d\n",n0+1);
    exit(1);
  }

   kssort(f_str->vptr,f_str->nsave);

   /* at some point, we want one best score for each of the segments */

   for ( ; f_str->nsave > 0; f_str->nsave--) 
     if (f_str->vptr[f_str->nsave-1]->score >0) break;

   /* copy aa0[] into f_str->aa0t[] */
   for (i=0; i<n0; i++) f_str->aa0t[i] = aa0[i];

  *nres = sconn_a (aa0t,n0,f_str,ppst->param_u.fa.cgap,
		   ppst->param_u.fa.pgap, aln);
  free(aa0t);

  *ares = f_str->res;
  return rst.score[0];
}

/* this version of sconn is modified to provide alignment information */
/* in addition, it needs to know whether a segment has been used before */

int sconn_a (unsigned char *aa0, int n0, 
	     struct f_struct *f_str,int cgap, int pgap,
	     struct a_struct *aln)
{
   int     i, si, cmpp (), n;
   unsigned char *aa0p;
   int sx, dx, doff;

   struct savestr **v;
   struct slink {
     int     score;
     struct savestr *vp;
     struct slink *snext;
     struct slink *aprev;
   } *start, *sl, *sj, *so, sarr[MAXSAV];
   int     lstart, lstop, ltmp, plstart, tstart, plstop, ptstop, ptstart, tstop;

   int *res, nres, tres;

/*	sort the score left to right in lib pos */

   v = f_str->vptr;
   n = f_str->nsave;

   pgap = 0;

   /* set things up in case nothing fits */
   if (v[0]->score <= 0) return 0;

   if (v[0]->score < cgap) {
     sarr[0].vp = v[0];
     sarr[0].score = v[0]->score;
     sarr[0].snext = NULL;
     sarr[0].aprev = NULL;
     start = &sarr[0];
   }
   else {

     krsort (v, n);	/* sort from left to right in library */

     start = NULL;

     /*	for each alignment, see if it fits */

     for (i = 0, si = 0; i < n; i++) {
       /*	if the score is less than the join threshold, skip it */

       if (v[i]->score < cgap) continue;

       lstart = v[i]->start;
       lstop = v[i]->stop;
       tstart = lstart - v[i]->dp + f_str->noff;
       tstop = lstop - v[i]->dp + f_str->noff;

       /*	put the alignment in the group */

       sarr[si].vp = v[i];
       sarr[si].score = v[i]->score;
       sarr[si].snext = NULL;
       sarr[si].aprev = NULL;

       /* 	if it fits, then increase the score */
       /* start points to a sorted (by total score) list of candidate
	  overlaps */

       for (sl = start; sl != NULL; sl = sl->snext) { 
	 plstart = sl->vp->start;
	 plstop = sl->vp->stop;
	 ptstart = plstart - sl->vp->dp + f_str->noff;
	 ptstop = plstop - sl->vp->dp + f_str->noff;
	 if (plstart > lstop && (ptstop < tstart || ptstart > tstop)) {
	   sarr[si].score = sl->score + v[i]->score + pgap;
	   sarr[si].aprev = sl;
	   break;		/* quit as soon as the alignment has been added */
	 }
       }

       /* now recalculate the list of best scores */
       if (start == NULL)
	 start = &sarr[si];	/* put the first one in the list */
       else
	 for (sj = start, so = NULL; sj != NULL; sj = sj->snext) {
	   if (sarr[si].score > sj->score) { /* new score better than old */
	     sarr[si].snext = sj;		/* snext best after new score */
	     if (so != NULL)
	       so->snext = &sarr[si];	/* prev_best->snext points to best */
	     else  start = &sarr[si];	/* start points to best */
	     break;			/* stop looking */
	   }
	   so = sj;		/* previous candidate best */
	 }
       si++;				/* increment to snext alignment */
     }
   }

   res = f_str->res;
   tres = nres = 0;
   aa0p = aa0;
   aln->min1 = start->vp->start+1;
   aln->min0 = 1;
   for (sj = start; sj != NULL; sj = sj->aprev ) {
     doff = (int)(aa0p-aa0) - (sj->vp->start-sj->vp->dp+f_str->noff);
     /*
       fprintf(stderr,"doff: %3d\n",doff);
     */
     for (dx=sj->vp->start,sx=sj->vp->start-sj->vp->dp+f_str->noff;
	  dx <= sj->vp->stop; dx++) {
       *aa0p++ = f_str->aa0t[sx++];
       tres++;
       res[nres++] = 0;
     }
     sj->vp->dp -= doff;
     if (sj->aprev != NULL) {
       if (sj->aprev->vp->start - sj->vp->stop - 1 > 0 )
	 tres += res[nres++] = (sj->aprev->vp->start - sj->vp->stop - 1);
     }

     /*
       fprintf(stderr,"t0: %3d, tx: %3d, l0: %3d, lx: %3d, dp: %3d noff: %3d, score: %3d\n",
       sj->vp->start - sj->vp->dp + f_str->noff,
       sj->vp->stop - sj->vp->dp + f_str->noff,
       sj->vp->start,sj->vp->stop,sj->vp->dp,
       f_str->noff,sj->vp->score);

       fprintf(stderr,"%3d - %3d: %3d\n",
       sj->vp->start,sj->vp->stop,sj->vp->score);
     */

     aln->max1 = sj->vp->stop+1;
     aln->max0 = aln->max1 - sj->vp->dp + f_str->noff;
   }

   /*     fprintf(stderr,"(%3d - %3d):(%3d - %3d)\n",
	    aln->min0,aln->max0,aln->min1,aln->max1);
     */
     /* now replace f_str->aa0t with aa0 */
   for (i=0; i<n0; i++) f_str->aa0t[i] = aa0[i];

   return tres;
}

int calcons(unsigned char *aa0, int n0,
	    unsigned char *aa1, int n1,
	    int *res, int nres, int *nc,
	    struct a_struct *aln, struct pstruct pst,
	    char *seqc0, char *seqc1,
	    struct f_struct *f_str)
{
  int i0, i1, nn1, n0t;
  int op, lenc, len_gap, nd, ns, itmp;
  unsigned char *aa1p;
  char *sp0, *sp1;
  int *rp;
  
#ifndef TFASTS
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* first fill in the ends */
  aln->min0--; aln->min1--;
  n0 -= (f_str->nm0-1);

  if (min(aln->min0,aln->min1)<aln->llen || aln->showall==1)
    			/* will we show all the start ?*/
    if (aln->min0>=aln->min1) {      /* aa0 extends more to left */
      aln->smins=0;
      if (aln->showall==1) aln->mins=aln->min0;
      else aln->mins = min(aln->min0,aln->llen/2);
      aancpy(seqc0,(char *)f_str->aa0t+aln->min0-aln->mins,aln->mins,pst);
      aln->smin0 = aln->min0-aln->mins;
      if ((aln->mins-aln->min1)>0) {
	memset(seqc1,' ',aln->mins-aln->min1);
	aancpy(seqc1+aln->mins-aln->min1,(char *)aa1p,aln->min1,pst);
	aln->smin1 = 0;
      }
      else {
	aancpy(seqc1,(char *)aa1p+aln->min1-aln->mins,aln->mins,pst);
	aln->smin1 = aln->min1-aln->mins;
      }
    }
    else {
      aln->smins=0;
      if (aln->showall == 1) aln->mins=aln->min1;
      else aln->mins = min(aln->min1,aln->llen/2);
      aancpy(seqc1,(char *)(aa1p+aln->min1-aln->mins),aln->mins,pst);
      aln->smin1 = aln->min1-aln->mins;
      if ((aln->mins-aln->min0)>0) {
	memset(seqc0,' ',aln->mins-aln->min0);
	aancpy(seqc0+aln->mins-aln->min0,(char *)f_str->aa0t,aln->min0,pst);
	aln->smin0 = 0;
      }
      else {
	aancpy(seqc0,(char *)f_str->aa0t+aln->min0-aln->mins,aln->mins,pst);
	aln->smin0 = aln->min0-aln->mins;
      }
    }
  else {
    aln->mins= min(aln->llen/2,min(aln->min0,aln->min1));
    aln->smins=aln->mins;
    aln->smin0=aln->min0;
    aln->smin1=aln->min1;
    aancpy(seqc0,(char *)f_str->aa0t+aln->min0-aln->mins,aln->mins,pst);
    aancpy(seqc1,(char *)aa1p+aln->min1-aln->mins,aln->mins,pst);
  }

/* now get the middle */

  sp0 = seqc0+aln->mins;
  sp1 = seqc1+aln->mins;
  rp = res;
  n0t = lenc = len_gap = aln->nident = aln->ngap_q = aln->ngap_l =aln->nfs = op = 0;
  i0 = aln->min0;
  i1 = aln->min1;
  
  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      *sp0 = pst.sq[f_str->aa0t[i0++]];
      *sp1 = pst.sq[aa1p[i1++]];
      n0t++;
      lenc++;
      if (toupper(*sp0) == toupper(*sp1)) aln->nident++;
      else if ((toupper(*sp0)=='I' && toupper(*sp1)=='L') ||
	       (toupper(*sp0)=='L' && toupper(*sp1)=='I')) aln->nident++;
      sp0++; sp1++;
    }
    else {
      if (op==0) { op = *rp++;}
      if (op>0) {
	*sp0++ = '-';
	*sp1++ = pst.sq[aa1p[i1++]];
	op--;
	len_gap++;
	lenc++;
      }
      else {
	*sp0++ = pst.sq[f_str->aa0t[i0++]];
	*sp1++ = '-';
	op++;
	n0t++;
	len_gap++;
	lenc++;
      }
    }
  }

  *nc = lenc-len_gap;
/*	now we have the middle, get the right end */

  ns = aln->mins + lenc + aln->llen;
  ns -= (itmp = ns %aln->llen);
  if (itmp>aln->llen/2) ns += aln->llen;
  nd = ns - (aln->mins+lenc);
  if (nd > max(n0t-aln->max0,nn1-aln->max1)) nd = max(n0t-aln->max0,nn1-aln->max1);
  
  if (aln->showall==1) {
    nd = max(n0t-aln->max0,nn1-aln->max1);	/* reset for showall=1 */
    /* get right end */
    aancpy(seqc0+aln->mins+lenc,(char *)f_str->aa0t+aln->max0,n0t-aln->max0,pst);
    aancpy(seqc1+aln->mins+lenc,(char *)aa1p+aln->max1,nn1-aln->max1,pst);
    /* fill with blanks - this is required to use one 'nc' */
    memset(seqc0+aln->mins+lenc+n0t-aln->max0,' ',nd-(n0t-aln->max0));
    memset(seqc1+aln->mins+lenc+nn1-aln->max1,' ',nd-(nn1-aln->max1));
  }
  else {
    if ((nd-(n0t-aln->max0))>0) {
      aancpy(seqc0+aln->mins+lenc,(char *)f_str->aa0t+aln->max0,
	     n0t-aln->max0,pst);
      memset(seqc0+aln->mins+lenc+n0t-aln->max0,' ',nd-(n0t-aln->max0));
    }
    else aancpy(seqc0+aln->mins+lenc,(char *)f_str->aa0t+aln->max0,nd,pst);
    if ((nd-(nn1-aln->max1))>0) {
      aancpy(seqc1+aln->mins+lenc,(char *)aa1p+aln->max1,nn1-aln->max1,pst);
      memset(seqc1+aln->mins+lenc+nn1-aln->max1,' ',nd-(nn1-aln->max1));
    }
    else aancpy(seqc1+aln->mins+lenc,(char *)aa1p+aln->max1,nd,pst);
  }
  
  return aln->mins+lenc+nd;
}


int calc_id(unsigned char *aa0, int n0,
	    unsigned char *aa1, int n1,
	    int *res, int nres,
	    struct a_struct *aln, struct pstruct pst,
	    struct f_struct *f_str)
{
  int i0, i1, nn1, n0t;
  int op, lenc, len_gap, nd, ns, itmp;
  unsigned char *aa1p;
  int sp0, sp1;
  int *rp;
  
#ifndef TFASTF
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* first fill in the ends */
  n0 -= (f_str->nm0-1);

  /* now get the middle */
  rp = res;
  n0t = lenc = len_gap = aln->nident = aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  i0 = aln->min0-1;
  i1 = aln->min1-1;
  
  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      sp0 = pst.sq[f_str->aa0t[i0++]];
      sp1 = pst.sq[aa1p[i1++]];
      n0t++;
      lenc++;
      if (toupper(sp0) == toupper(sp1)) aln->nident++;
    }
    else {
      if (op==0) { op = *rp++;}
      if (op>0) {
	i1++;
	op--;
	len_gap++;
	lenc++;
      }
      else {
	i0++;
	op++;
	n0t++;
	len_gap++;
	lenc++;
      }
    }
  }
  
  return lenc-len_gap;
}


#ifdef PCOMPLIB
#include "p_mw.h"
void
update_params(struct qmng_str *qm_msg, struct pstruct *ppst)
{
  ppst->n0 = qm_msg->n0;
}
#endif
