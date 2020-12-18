/* copyright (c) 1998, 1999 William R. Pearson and the U. of Virginia */

/*  - dropfs.c,v 1.1.1.1.2.1 1999/12/17 14:05:40 ajm6q Exp */

/* changed to return 2.0, rather than -1.0, for failure */

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
#include "tatstats.h"

#define EOSEQ 0
#define ESS 49
#define NMAP_X 23 /* for 'X' */
#define NMAP_Z 24 /* for '*' */
#define MAXHASH 32
#define NMAP MAXHASH+1

static char *verstr="4.2 March 2004";

int shscore(const unsigned char *aa0, const int n0, int **pam2, int nsq);

#ifdef TFAST
extern int aatran(const unsigned char *ntseq, unsigned char *aaseq, const int maxs, const int frame);
#endif

void savemax(struct dstruct *, struct f_struct *);

int spam(const unsigned char *, const unsigned char *, int, struct savestr *, int **, struct f_struct *);
int sconn(struct savestr **v,
	  int nsave,
	  int cgap,
	  int pgap,
	  struct f_struct *,
	  struct rstruct *,
	  struct pstruct *,
	  unsigned char *aa0, int n0,
	  unsigned char *aa1, int n1,
	  int opt_prob);

void kpsort(struct savestr **, int);
void kssort(struct savestr **, int);
int sconn_a(unsigned char *, int, unsigned char *, int,
	    struct f_struct *, int cgap, int pgap, struct a_struct *,
	    struct pstruct *);
void kpsort(struct savestr **, int);

/* initialize for fasta */

void
init_work (unsigned char *aa0, const int n0, 
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
   int mhv, phv;
   int hmax, nsegs;
   int i0, ib, hv, old_hv;
   int pamfact;
   int btemp;
   struct f_struct *f_str;
   /* these used to be globals, but do not need to be */
   int ktup, fact, kt1;

   int maxn0;
   int i, j, q;
   int tat_size;
   int *res;

   unsigned char *query;
   int k, l, m, n, N, length, index;

   double *tatprobptr;

   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   btemp = shscore(aa0,n0,ppst->pam2[0],ppst->nsq) / 3;
   ppst->param_u.fa.cgap = btemp;

   ppst->param_u.fa.pgap = ppst->gdelval + ppst->ggapval;
   ktup = ppst->param_u.fa.ktup;
   fact = ppst->param_u.fa.scfact;

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

   for (f_str->kshft = 0; mhv > 0; mhv /= 2) f_str->kshft++;

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

   if ((f_str->aa0ti = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0ti array; %d\n",n0+1);
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

   /* for FASTS/FASTM, we want to know when we get to the end of a peptide,
      so we can ensure that we set the end and restart */

   if ((f_str->l_end = (int *) calloc (n0, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate link end array");
     exit (1);
   }

   for (i0 = 0; i0 < hmax; i0++) f_str->harr[i0] = -1;
   for (i0 = 0; i0 < n0; i0++) f_str->link[i0] = -1;
   for (i0 = 0; i0 < n0; i0++) f_str->l_end[i0] = 0;

   /* count the number of peptides */
   nsegs = 1;
   for (i0 = 0; i0 < n0; i0++) {
     if (aa0[i0] == ESS || aa0[i0] == 0) nsegs++;
   }

   /* allocate space for peptides offsets, nm_u */
   if ((f_str->nmoff = (int *)calloc(nsegs+1, sizeof(int)))==NULL) {
     fprintf(stderr, " cannot allocat nmoff array: %d\n", nsegs);
     exit(1);
   }

   if ((f_str->nm_u = (int *)calloc(nsegs+1, sizeof(int)))==NULL) {
     fprintf(stderr, " cannot allocat nm_u array: %d\n", nsegs);
     exit(1);
   }

   /* encode the aa0 array */
   /* no initial loop because kt1 = 0 */

   phv = hv = 0;
   f_str->nmoff[0] = 0;
   f_str->nm0 = 1;
   for (i0 = 0; i0 < n0; i0++) {
     if (aa0[i0] == ESS || aa0[i0] == 0) {
       aa0[i0] = EOSEQ;	/* set ESS to 0 */
       /*       fprintf(stderr," converted ',' to 0\n");*/
       f_str->nmoff[f_str->nm0++] = i0+1; 
       f_str->l_end[i0-1] = 1;
       phv = hv = 0;
       continue;
     }

     hv = ((hv & f_str->hmask) << f_str->kshft) + ppst->hsq[aa0[i0]];
     f_str->link[i0] = f_str->harr[hv];
     f_str->harr[hv] = i0;
     f_str->pamh2[hv] = ppst->pam2[0][aa0[i0]][aa0[i0]];
   }
   f_str->l_end[n0-1] = 1;

   f_str->nmoff[f_str->nm0] = n0+1;

   /* build an integer array of the max score that can be achieved from that position -
      use in savemax to mark some segments as fixed */

   ppst->param_u.fa.cgap /= f_str->nm0;

   /* setup aa0b[], aa0e[], which specify the begining and end of each
      segment */

   for (ib = i0 = 0; i0 < n0; i0++) {
     if (aa0[i0]==EOSEQ) {
       f_str->aa0b[i0] = -1;
       f_str->aa0e[i0] = -1;
       f_str->aa0i[i0] = -1;
       ib++;
     }
     f_str->aa0b[i0] =  f_str->nmoff[ib];
     f_str->aa0e[i0] =  f_str->nmoff[ib+1]-2;
     f_str->aa0i[i0] =  ib;

     /*
     fprintf(stderr,"%2d %c: %2d %2d %2d\n",i0,ppst->sq[aa0[i0]],
	     f_str->aa0b[i0],f_str->aa0e[i0],f_str->aa0i[i0]);
     */
   }

   /* now that we have f_str->nm0 (number of peptides), we can
      allocate space for what used to be: 
      vmax[MAXSAV], *vptr[MAXSAV], and sarr[MAXSAV]
   */

   f_str->maxsav = max(MAXSAV,2*f_str->nm0);

   if ((f_str->vmax = (struct savestr *)
	calloc(f_str->maxsav,sizeof(struct savestr)))==NULL) {
     fprintf(stderr, "Couldn't allocate vmax[%d].\n",f_str->maxsav);
     exit(1);
   }

   if ((f_str->vptr = (struct savestr **)
	calloc(f_str->maxsav,sizeof(struct savestr *)))==NULL) {
     fprintf(stderr, "Couldn't allocate vptr[%d].\n",f_str->maxsav);
     exit(1);
   }

   if ((f_str->sarr = (struct slink *)
	calloc(f_str->maxsav,sizeof(struct slink)))==NULL) {
     fprintf(stderr, "Couldn't allocate sarr[%d].\n",f_str->maxsav);
     exit(1);
   }

   /* Tatusov Statistics Setup */

   /* initialize priors array. */
   if((f_str->priors = (double *)calloc(ppst->nsq+1, sizeof(double))) == NULL) {
     fprintf(stderr, "Couldn't allocate priors array.\n");
     exit(1);
   }

   calc_priors(f_str->priors, ppst, f_str, NULL, 0, ppst->pseudocts);

   /* pre-calculate the Tatusov probability array for each full segment */

   if(ppst->zsflag >= 1 && ppst->zsflag <= 3 && f_str->nm0 <= 10) {

     tat_size = (1<<f_str->nm0) -1;
     f_str->dotat = 1;
     f_str->tatprobs = (struct tat_str **) malloc((size_t)tat_size*sizeof(struct tat_str *));
     if (f_str->tatprobs == NULL) {
       fprintf (stderr, " cannot allocate tatprobs array: %ld\n",
		tat_size * sizeof(struct tat_str *));
       exit (1);
     }

     f_str->intprobs = (double **) malloc((size_t)tat_size * sizeof(double *));
     if(f_str->intprobs == NULL) {
       fprintf(stderr, "Couldn't allocate intprobs array.\n");
       exit(1);
     }

     for(k = 0, l = f_str->nm0 ; k < l ; k++) {
       query = &(aa0[f_str->nmoff[k]]);
       length = f_str->nmoff[k+1] - f_str->nmoff[k] - 1;

       /* this segment alone */
       index = (1 << k) - 1;
       generate_tatprobs(query, 0, length - 1, f_str->priors, ppst->pam2[0], ppst->nsq, &(f_str->tatprobs[index]), NULL);

       /* integrate the probabilities */
       N = f_str->tatprobs[index]->highscore - f_str->tatprobs[index]->lowscore;
       tatprobptr = (double *) calloc(N+1, sizeof(double));
       if(tatprobptr == NULL) {
	 fprintf(stderr, "Couldn't calloc tatprobptr.\n");
	 exit(1);
       }
       f_str->intprobs[index] = tatprobptr;

       for (i = 0; i <= N ; i++ ) {
	 tatprobptr[i] = f_str->tatprobs[index]->probs[i];
	 for (j = i + 1 ; j <= N ; j++ ) {
	   tatprobptr[i] += f_str->tatprobs[index]->probs[j];
	 }
       }

       /* this segment built on top of all other subcombinations */
       for(i = 0, j = (1 << k) - 1 ; i < j ; i++) {
	 index = (1 << k) + i;
	 generate_tatprobs(query, 0, length - 1, f_str->priors, ppst->pam2[0], ppst->nsq, &(f_str->tatprobs[index]), f_str->tatprobs[i]);

	 /* integrate the probabilities */
	 N = f_str->tatprobs[index]->highscore - f_str->tatprobs[index]->lowscore;
	 tatprobptr = (double *) calloc(N+1, sizeof(double));
	 if(tatprobptr == NULL) {
	   fprintf(stderr, "Couldn't calloc tatprobptr.\n");
	   exit(1);
	 }
	 f_str->intprobs[index] = tatprobptr;
       
	 for (m = 0; m <= N ; m++ ) {
	   tatprobptr[m] = f_str->tatprobs[index]->probs[m];
	   for (n = m + 1 ; n <= N ; n++ ) {
	     tatprobptr[m] += f_str->tatprobs[index]->probs[n];
	   }
	 }
       }
     }
   } else {
     f_str->dotat = 0;
     f_str->shuff_cnt = ppst->shuff_node;
   }

   /* End of Tatusov Statistics Setup */

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
   if (f_str->diag==NULL) 
     f_str->diag = (struct dstruct *) calloc ((size_t)MAXDIAG,
					      sizeof (struct dstruct));
   if (f_str->diag == NULL) {
      fprintf (stderr, " cannot allocate diagonal arrays: %ld\n",
	      (long) MAXDIAG * (long) (sizeof (struct dstruct)));
      exit (1);
   }

#ifdef TFAST
   if ((f_str->aa1x =(unsigned char *)calloc((size_t)ppst->maxlen+2,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "cannot allocate aa1x array %d\n", ppst->maxlen+2);
     exit (1);
   }
   f_str->aa1x++;
#endif

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
#ifdef FASTS
#ifndef TFAST
  char *pg_str="FASTS";
#else
  char *pg_str="TFASTS";
#endif
#endif

#ifdef FASTM
#ifndef TFAST
  char *pg_str="FASTM";
#else
  char *pg_str="TFASTM";
#endif
#endif

  sprintf (pstring1, "%s (%s) function [%s matrix (%d:%d)]",pg_str,verstr,
	       pstr->pamfile, pstr->pam_h,pstr->pam_l);
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
  int i, j;

  f_str = *f_arg;

  if (f_str != NULL) {

    free(f_str->res);
#ifdef TFAST
    free(f_str->aa1x - 1); /* because f_str->aa1x got ++'ed when allocated! */
#endif
    free(f_str->diag);
    free(f_str->l_end);
    free(f_str->link);
    free(f_str->pamh2); 
    free(f_str->pamh1);
    free(f_str->harr);
    free(f_str->vmax);
    free(f_str->vptr);
    free(f_str->sarr);
    free(f_str->aa0i);
    free(f_str->aa0e);
    free(f_str->aa0b);
    free(f_str->aa0ti);
    free(f_str->aa0t);
    free(f_str->nmoff);
    free(f_str->nm_u);

    if(f_str->dotat) {
      for(i = 0, j = (1 << f_str->nm0) - 1 ; i < j ; i++) {
	free(f_str->tatprobs[i]->probs);
	free(f_str->tatprobs[i]);
	free(f_str->intprobs[i]);
      }
      free(f_str->tatprobs);
      free(f_str->intprobs);
    }

    free(f_str->priors);
    free(f_str);
    *f_arg = NULL;
  }
}

void do_fasts (unsigned char *aa0, const int n0,
	       unsigned char *aa1, const int n1,
	       struct pstruct *ppst, struct f_struct *f_str,
	       struct rstruct *rst, int *hoff, int opt_prob)
{
   int     nd;		/* diagonal array size */
   int     lhval;
   int     kfact;
   register struct dstruct *dptr;
   register int tscor;
   register struct dstruct *diagp;
   struct dstruct *dpmax;
   register int lpos;
   int     tpos;
   struct savestr *vmptr, *vmaxmax;
   int     scor, tmp;
   int     im, ib, nsave;
   int     cmps ();		/* comparison routine for ksort */
   int ktup;
   int doffset;


   vmaxmax = &f_str->vmax[f_str->maxsav];

   ktup = 1;

   if (n1 < ktup) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     rst->escore = 1.0;
     rst->segnum = 0;
     rst->seglen = 0;
     return;
   }

   if (n0+n1+1 >= MAXDIAG) {
     fprintf(stderr,"n0,n1 too large: %d, %d\n",n0,n1);
     rst->score[0] = rst->score[1] = rst->score[2] = -1;
     rst->escore = 2.0;
     rst->segnum = 0;
     rst->seglen = 0;
     return;
   }

   nd = n0 + n1;

   dpmax = &f_str->diag[nd];
   for (dptr = &f_str->diag[f_str->ndo]; dptr < dpmax;)
   {
      dptr->stop = -1;
      dptr->dmax = NULL;
      dptr++->score = 0;
   }

   for (vmptr = f_str->vmax; vmptr < vmaxmax; vmptr++)
      vmptr->score = 0;
   f_str->lowmax = f_str->vmax;
   f_str->lowscor = 0;

   /* start hashing */
   diagp = &f_str->diag[f_str->noff];
   for (lhval=lpos=0; lpos < n1; lpos++, diagp++) {
     if (ppst->hsq[aa1[lpos]]>=NMAP) {	/* skip residue */
       lpos++ ; diagp++;
       while (lpos < n1 && ppst->hsq[aa1[lpos]]>=NMAP) {lpos++; diagp++;}
       if (lpos >= n1) break;
       lhval = 0;
     }

     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];

     for (tpos = f_str->harr[lhval]; tpos >= 0; tpos = f_str->link[tpos]) {
       if ((tscor = (dptr = &diagp[-tpos])->stop) >= 0) {
	 tscor++;	/* tscor is stop of current, increment it */
	 if ((tscor -= lpos) <= 0) {  /* tscor, the end of the current
					 match, is before lpos, so there
					 is a mismatch - this is also the
					 mismatch cost */
	   scor = dptr->score;	/* save the run score on the diag */
	   if ((tscor += (kfact = f_str->pamh2[lhval])) < 0 
	       && f_str->lowscor < scor) {
	     /* if what we will get (tscor + kfact) is < 0 and the
		score is better than the worst savemax() score, save
		it */
	     savemax (dptr, f_str);
	   }

	   /* if extending is better than starting over, extend */
	   if ((tscor += scor) >= kfact) {
	     dptr->score = tscor;
	     dptr->stop = lpos;
	     if (f_str->l_end[tpos]) {
	       savemax(dptr, f_str);
	       dptr->score = -10000;
	     }
	   }
	   else {     /* otherwise, start new */
	     dptr->score = kfact;
	     dptr->start = dptr->stop = lpos;
	   }
	 } 
	 else { /* tscor is after lpos, so extend one residue */
	   dptr->score += f_str->pamh1[aa0[tpos]];
	   dptr->stop = lpos;
	   if (f_str->l_end[tpos]) {
	     savemax(dptr, f_str);
	     dptr->score = -10000;
	   }
	 }
       }
       else {	/* start new */
	 dptr->score = f_str->pamh2[lhval];
	 dptr->start = dptr->stop = lpos;
       }
     }				/* end tpos */
   }				/* end lpos */

   for (dptr = f_str->diag; dptr < dpmax;) {
     if (dptr->score > f_str->lowscor) savemax (dptr, f_str);
     dptr->stop = -1;
     dptr->dmax = NULL;
     dptr++->score = 0;
   }
   f_str->ndo = nd;

/*
        at this point all of the elements of aa1[lpos]
        have been searched for elements of aa0[tpos]
        with the results in diag[dpos]
*/

   for (nsave=0, vmptr=f_str->vmax; vmptr< vmaxmax; vmptr++) {
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
	fprintf(stderr,"  sscore: %d %d-%d\n",vmptr->score,vmptr->start,vmptr->stop);
	*/
	if (vmptr->score > 0) f_str->vptr[nsave++] = vmptr;
      }
   }

   if (nsave <= 0) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     rst->escore = 1.0;
     rst->segnum = 0;
     rst->seglen = 0;
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

   for (ib=0; ib<f_str->nm0; ib++) f_str->nm_u[ib]=0;
   for (ib=0; ib < nsave; ib++) {
     doffset = f_str->vptr[ib]->dp - f_str->noff;
     tpos=f_str->aa0i[f_str->vptr[ib]->start - doffset];
     if (f_str->nm_u[tpos] == 0) {
       f_str->nm_u[tpos]=1;
     } else {
       f_str->vptr[ib]->score = -1;
     }
   }

   kssort(f_str->vptr,nsave);
   for (ib = nsave-1; ib >= 0; ib--)
     if (f_str->vptr[ib]->score > -1) break;
   nsave = ib+1;

   scor = sconn (f_str->vptr, nsave, ppst->param_u.fa.cgap, 
		 ppst->param_u.fa.pgap, f_str, rst, ppst, aa0, n0, aa1, n1,
		 opt_prob);

   if (rst->escore < 0.0) rst->escore = 2.0;
   kssort(f_str->vptr,nsave);

   /* here we should use an nsave that is consistent with sconn and nm0 */

   f_str->nsave = nsave;
   if (nsave > f_str->nm0) f_str->nsave = f_str->nm0;

   rst->score[1] = f_str->vptr[0]->score;
   rst->score[0] = rst->score[2] = max(scor, f_str->vptr[0]->score);

}

void do_work (unsigned char *aa0, const int n0,
	      unsigned char *aa1, const int n1,
	      int frame,
	      struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, struct rstruct *rst)
{
  int opt_prob;
  int hoff, n10, i;

  if (qr_flg==1 && f_str->shuff_cnt <= 0) {
    rst->escore = 2.0;
    rst->score[0]=rst->score[1]=rst->score[2]= -1;
    return;
  }

  if (f_str->dotat || ppst->zsflag == 4 || ppst->zsflag == 14 ) opt_prob=1;
  else opt_prob = 0;
  if (ppst->zsflag == 2 || ppst->zsflag == 12) opt_prob = 0;
  if (qr_flg) {
    opt_prob=1;
    /*    if (frame==1) */
      f_str->shuff_cnt--;
  }

  if (n1 < ppst->param_u.fa.ktup) {
    rst->score[0] = rst->score[1] = rst->score[2] = -1;
    rst->escore = 2.0;
    return;
  }
#ifdef TFAST
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

  do_fasts (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff, opt_prob);
#else	/* FASTA */
  do_fasts (aa0, n0, aa1, n1, ppst, f_str, rst, &hoff, opt_prob);
#endif

  rst->comp = rst->H = -1.0;
}

void do_opt (unsigned char *aa0, const int n0,
	     unsigned char *aa1, const int n1,
	     int frame,
	     struct pstruct *ppst, struct f_struct *f_str,
	     struct rstruct *rst)
{
  int lag, tscore, hoff, n10;

#ifdef TFAST
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  do_fasts (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff, 1);
#else	/* FASTA */
  do_fasts(aa0,n0,aa1,n1,ppst,f_str,rst, &hoff, 1);
#endif
}


/* modify savemax() so that full length 100% matches are marked
   so that they cannot be removed

   modify savemax() to split alignments that include a comma
*/

/* savemax(dptr, f_str) takes a current diagonal run (saved in dptr),
   and places it in the set of runs to be saved (in  f_str->vmax[])
*/

void 
savemax (dptr, f_str)
  register struct dstruct *dptr;
  struct f_struct *f_str;
{
  register int dpos;	/* position along the diagonal, -n0 .. n1 */
  register int i;
  register struct savestr *vmptr;
  struct savestr *vmaxmax;

  vmaxmax = &f_str->vmax[f_str->maxsav];

  dpos = (int) (dptr - f_str->diag);	/* current diagonal */

/* check to see if this is the continuation of a run that is already saved */
/* if we are at the end of the query, save it regardless */

  if ((vmptr = dptr->dmax) != NULL	/* have an active run */
      && vmptr->dp == dpos &&		/* on the correct diagonal */
      vmptr->start == dptr->start) {	/* and it starts at the same place */
    vmptr->stop = dptr->stop;	/* update the end of the match in vmax[] */

/* if the score is worse, don't update, return - if the score gets bad
   enough, it will restart in the diagonal scan */
    if ((i = dptr->score) <= vmptr->score) { return;} 

/* score is better, update */
    vmptr->score = i;
/* if the score is not the worst, return */
    if (vmptr != f_str->lowmax) { return;}
  }
  else {	/* not a continuation */
    /* save in the lowest place */
    i = f_str->lowmax->score = dptr->score;   /* 'i' is used as a bound */
    f_str->lowmax->dp = dpos;
    f_str->lowmax->start = dptr->start;
    f_str->lowmax->stop = dptr->stop;
    dptr->dmax = f_str->lowmax;
  }

  /* rescan the list for the worst score */
  for (vmptr = f_str->vmax; vmptr < vmaxmax; vmptr++)
    if (vmptr->score < i) {
      i = vmptr->score;
      f_str->lowmax = vmptr;
    }
  f_str->lowscor = i;
}

/* this version of spam scans the diagonal to find the best local score,
   then resets the boundaries for a global alignment and re-scans */

/* NOOVERHANG allows one to score any overhanging alignment as zero.
   Useful for SAGE alignments.  Normally, one allows overhangs because
   of the possibility of partial sequences.
*/

#undef NOOVERHANG

/* 
   May, 2005 - spam() has an intesting bug that occurs when two
   peptides match in order, separated by one position (the comma).  In
   this case, spam() splits the match, and only returns the better of
   the two matches.  So, if spam splits an alignment at a comma, it
   needs the ability to insert the missing match.

*/

int spam (const unsigned char *aa0, const unsigned char *aa1,int n1,
	  struct savestr *dmax, int **pam2,
	  struct f_struct *f_str)
{
   int     lpos, doffset;
   int     tot, mtot;
   struct {
     int  start, stop, score;
   } curv, maxv;
   register const unsigned char *aa0p, *aa1p;

   doffset = dmax->dp - f_str->noff;
   aa1p = &aa1[dmax->start];
   aa0p = &aa0[dmax->start - doffset];
   curv.start = lpos;

   tot = curv.score = maxv.score = 0;
   for (lpos = dmax->start; lpos <= dmax->stop; lpos++) {
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
   if (maxv.start < 0) {
     maxv.start = 0;
#ifdef NOOVERHANG
     return 0;
#endif
   }

   maxv.stop = f_str->aa0e[maxv.stop-doffset] + doffset;
   if (maxv.stop > n1) {
     maxv.stop = n1-1;
#ifdef NOOVERHANG
     return 0;
#endif
   }
   aa1p = &aa1[lpos = maxv.start];
   aa0p = &aa0[lpos - doffset];

   for (tot=0; lpos <= maxv.stop; lpos++) {
     tot += pam2[*aa0p++][*aa1p++];
   }

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
	   int cgap, int pgap, struct f_struct *f_str,
	   struct rstruct *rst, struct pstruct *ppst,
	   unsigned char *aa0, int n0,
	   unsigned char *aa1, int n1, int opt_prob)
{
   int     i, si, cmpp ();
   struct slink *start, *sl, *sj, *so, *sarr;
   int     lstart, ltmp, tstart, plstop, ptstop, ptstart, tstop;
   double  tatprob;
   int     dotat;

   sarr = f_str->sarr;

   pgap = 0;	/* the data should have gaps */

   /*  sort the score left to right in lib pos */
   kpsort (v, n);

   start = NULL;
   rst->score[0] = 0;
   rst->escore = 2.0;

/*  for the remaining runs, see if they fit */
/*  lstart/lstop -> start/stop in library sequence
    tstart/tstop -> start/stop in query sequence
    plstart/plstop ->
*/

   for (i = 0, si = 0; i < n; i++) {

     /* the segment is worth adding; find out where? */
     lstart = v[i]->start;
     ltmp = v[i]->stop;
     tstart = lstart - v[i]->dp + f_str->noff;
     tstop = ltmp - v[i]->dp + f_str->noff;

     /*	put the run in the group */
     sarr[si].vp = v[i];
     sarr[si].score = v[i]->score;
     sarr[si].next = NULL;
     sarr[si].prev = NULL;
     sarr[si].tat = NULL;

/*
  opt_prob for FASTS only has to do with using aa1 for priors,
  i.e. we always calculate tatprobs for segments in FASTS (unlike
  FASTF)
*/
     if(opt_prob) {
       sarr[si].tatprob = 
	 calc_tatusov(NULL, &sarr[si], aa0, n0, aa1, n1, 
		      ppst->pam2[0], ppst->nsq, f_str, 
		      ppst->pseudocts, opt_prob, ppst->zsflag);
       if (sarr[si].tatprob < 0.0) {
	 fprintf(stderr," negative tatprob: %lg\n",sarr[si].tatprob);
	 sarr[si].tatprob = 1.0;
       }
       sarr[si].tat = sarr[si].newtat;
     }

/*  if it fits, then increase the score

    start points to the highest scoring run
    -> next is the second highest, etc.
    put the segment into the highest scoring run that it fits into
*/
     for (sl = start; sl != NULL; sl = sl->next) {
       ltmp = sl->vp->start;
 /* plstop -> previous lstop */
       plstop = sl->vp->stop;
 /* ptstart -> previous t(query) start */
       ptstart = ltmp - sl->vp->dp + f_str->noff;
 /* ptstop -> previous t(query) stop */
       ptstop = plstop - sl->vp->dp + f_str->noff;
#ifndef FASTM
 /* if the previous library stop is before the current library start */
       if (plstop < lstart && ( ptstop < tstart || ptstart > tstop))
#else
 /* if the previous library stop is before the current library start */
       if (plstop < lstart && ptstop < tstart)
#endif
       {
	 if(!opt_prob) {
	    sarr[si].score = sl->score + v[i]->score + pgap;
	    sarr[si].prev = sl;
	    break;
	  } else {
	    tatprob = calc_tatusov(sl, &sarr[si], aa0, n0, aa1, n1, 
				   ppst->pam2[0], ppst->nsq, f_str, 
				   ppst->pseudocts, opt_prob, ppst->zsflag);
	    /* if our tatprob gets worse when we add this, forget it */
	    if(tatprob > sarr[si].tatprob) {
	      free(sarr[si].newtat->probs); /* get rid of new tat struct */
	      free(sarr[si].newtat);
	      continue; /* reuse this sarr[si] */
	    } else {
	      sarr[si].tatprob = tatprob;
	      free(sarr[si].tat->probs); /* get rid of old tat struct */
	      free(sarr[si].tat);
	      sarr[si].tat = sarr[si].newtat;
	      sarr[si].prev = sl;
	      sarr[si].score = sl->score + v[i]->score + pgap;
	      /*
		fprintf(stderr,"sconn %d added %d/%d getting %d; si: %d, tat: %g\n",
		i,v[i]->start, v[i]->score,sarr[si].score,si, tatprob);
	      */
	      break;
	    }
	  }
	}
      }
      
      /* now recalculate where the score fits */
      if (start == NULL) start = &sarr[si];
      else {
	if(!opt_prob) {
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
	} else {
	  for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	    if ( sarr[si].tatprob < sj->tatprob ||
		 ((sarr[si].tatprob == sj->tatprob) && sarr[si].score > sj->score) ) {
	      sarr[si].next = sj;
	      if (so != NULL)
		so->next = &sarr[si];
	      else
		start = &sarr[si];
	      break;
	    }
	    so = sj;
	  }
	}
      }

      si++;
   }
      
   if(opt_prob) {
     for (i = 0 ; i < si ; i++) {
       free(sarr[i].tat->probs);
       free(sarr[i].tat);
     }
   }

   if (start != NULL) {
     if(opt_prob) {
       rst->escore = start->tatprob;
     } else {
       rst->escore = 2.0;
     }

     rst->segnum = rst->seglen = 0;
     for(sj = start ; sj != NULL; sj = sj->prev) {
       rst->segnum++;
       rst->seglen += sj->vp->stop - sj->vp->start + 1;
     }
     return (start->score);
   } else {
     rst->escore = 1.0;
   }

   rst->segnum = rst->seglen = 0;
   return (0);
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

#ifdef TFAST
  f_str->n10 = n10 = aatran(aa1,f_str->aa1x,n1,frame);
  aa1p = f_str->aa1x;
  aln->qlrev = 0;
  aln->qlfact= 1;
  aln->llfact = aln->llmult = 3;
  if (frame > 3) aln->llrev = 1;
  else aln->llrev = 0;
  aln->frame = 0;
#else
  aln->llfact = aln->llmult = aln->qlfact = 1;
  aln->llrev = aln->qlrev = 0;
  aln->frame = 0;
  n10 = n1;
  aa1p = aa1;
#endif

  do_fasts(aa0, n0, aa1p, n10, ppst, f_str, &rst, &hoff, 1);

  /* the alignment portion takes advantage of the information left
     over in f_str after do_fasts is done.  in particular, it is
     easy to run a modified sconn() to produce the alignments.

     unfortunately, the alignment display routine wants to have
     things encoded as with bd_align and sw_align, so we need to do that.
     */

  /* unnecessary; do_fasts just did this */
  /*  kssort(f_str->vptr,f_str->nsave);  */

   /* at some point, we want one best score for each of the segments */

   for ( ; f_str->nsave > 0; f_str->nsave--) 
     if (f_str->vptr[f_str->nsave-1]->score >0) break;

   if ((aa0t = (unsigned char *)calloc(n0+1,sizeof(unsigned char)))==NULL) {
     fprintf(stderr," cannot allocate aa0t %d\n",n0+1);
     exit(1);
   }

   /* copy aa0[] into f_str->aa0t[] */
   for (i=0; i<n0; i++) f_str->aa0t[i] = aa0t[i] = aa0[i];
   f_str->aa0t[i] = aa0t[i] = '\0';

   *nres = sconn_a (aa0t,n0,aa1p,n10,f_str,ppst->param_u.fa.cgap,
		    ppst->param_u.fa.pgap, aln, ppst);

   free(aa0t);

   *ares = f_str->res;
   return rst.score[0];
}

/* this version of sconn is modified to provide alignment information */
/* in addition, it needs to know whether a segment has been used before */

/* sconn_a fills in the res[nres] array, but this is passed implicitly
   through f_str->res[f_str->nres] */

int sconn_a (unsigned char *aa0, int n0, unsigned char *aa1, int n1,
	     struct f_struct *f_str,int cgap, int pgap,
	     struct a_struct *aln, struct pstruct *ppst)
{
   int     i, si, cmpp (), n;
   unsigned char *aa0p;
   int sx, dx, doff, *aa0tip;

   struct savestr **v;
   struct slink *start, *sl, *sj, *so, *sarr;
   int     lstart, lstop, ltmp, plstart, tstart, plstop, ptstop, ptstart, tstop;

   int *res, nres, tres;

   double tatprob;

/*	sort the score left to right in lib pos */

   v = f_str->vptr;
   n = f_str->nsave;
   sarr = f_str->sarr;

   pgap = 0; cgap = 0;

   /* set things up in case nothing fits */
   if (v[0]->score <= 0) return 0;

   if (v[0]->score < cgap) {
     sarr[0].vp = v[0];
     sarr[0].score = v[0]->score;
     sarr[0].next = NULL;
     sarr[0].prev = NULL;
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
       sarr[si].next = NULL;
       sarr[si].prev = NULL;
       sarr[si].tat = NULL;

       sarr[si].tatprob = 
	 calc_tatusov(NULL, &sarr[si], aa0, n0, aa1, n1, 
		      ppst->pam2[0], ppst->nsq, f_str, 
		      ppst->pseudocts, 1, ppst->zsflag);
       sarr[si].tat = sarr[si].newtat;


       /* 	if it fits, then increase the score */
       /* start points to a sorted (by total score) list of candidate
	  overlaps */

       for (sl = start; sl != NULL; sl = sl->next) { 
	 plstart = sl->vp->start;
	 plstop = sl->vp->stop;
	 ptstart = plstart - sl->vp->dp + f_str->noff;
	 ptstop = plstop - sl->vp->dp + f_str->noff;
#ifndef FASTM
	 if (plstart > lstop && (ptstop < tstart || ptstart > tstop)) {
#else
         if (plstop > lstart && ptstart > tstop) {
#endif
	   /* alignment always uses probabilistic scoring ... */
	   /*   sarr[si].score = sl->score + v[i]->score + pgap;
		sarr[si].prev = sl;
		break; */		/* quit as soon as the alignment has been added */

	   tatprob = calc_tatusov(sl, &sarr[si], aa0, n0, aa1, n1, 
				  ppst->pam2[0], ppst->nsq, f_str, 
				  ppst->pseudocts, 1, ppst->zsflag);
	   /* if our tatprob gets worse when we add this, forget it */
	   if(tatprob > sarr[si].tatprob) {
	     free(sarr[si].newtat->probs); /* get rid of new tat struct */
	     free(sarr[si].newtat);
	     continue; /* reuse this sarr[si] */
	   } else {
	     sarr[si].tatprob = tatprob;
	     free(sarr[si].tat->probs); /* get rid of old tat struct */
	     free(sarr[si].tat);
	     sarr[si].tat = sarr[si].newtat;
	     sarr[si].prev = sl;
	     sarr[si].score = sl->score + v[i]->score + pgap;
	     /*
	       fprintf(stderr,"sconn %d added %d/%d getting %d; si: %d, tat: %g\n",
	       i,v[i]->start, v[i]->score,sarr[si].score,si, tatprob);
	     */
	     break;
	   }
	 }
       }

       /* now recalculate the list of best scores */
       if (start == NULL)
	 start = &sarr[si];	/* put the first one in the list */
       else
	 for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	   /* if (sarr[si].score > sj->score) { */ /* new score better than old */
	   if ( sarr[si].tatprob < sj->tatprob ||
		((sarr[si].tatprob == sj->tatprob) && sarr[si].score > sj->score) ) {
	     sarr[si].next = sj;		/* next best after new score */
	     if (so != NULL)
	       so->next = &sarr[si];	/* prev_best->next points to best */
	     else  start = &sarr[si];	/* start points to best */
	     break;			/* stop looking */
	   }
	   so = sj;		/* previous candidate best */
	 }
       si++;				/* increment to next alignment */
     }
   }

   for (i = 0 ; i < si ; i++) {
     free(sarr[i].tat->probs);
     free(sarr[i].tat);
   }

   res = f_str->res;
   tres = nres = 0;
   aa0p = aa0;
   aa0tip = f_str->aa0ti;	/* point to temporary index */
   aln->min1 = start->vp->start;
   aln->min0 = 0;

   for (sj = start; sj != NULL; sj = sj->prev ) {
     doff = (int)(aa0p-aa0) - (sj->vp->start-sj->vp->dp+f_str->noff);
     
     /* fprintf(stderr,"doff: %3d\n",doff); */
     
     for (dx=sj->vp->start,sx=sj->vp->start-sj->vp->dp+f_str->noff;
	  dx <= sj->vp->stop; dx++) {
       *aa0tip++ = f_str->aa0i[sx];	/* save index */
       *aa0p++ = f_str->aa0t[sx++];	/* save sequence at index */
       tres++;
       res[nres++] = 0;
     }
     sj->vp->dp -= doff;
     if (sj->prev != NULL) {
       if (sj->prev->vp->start - sj->vp->stop - 1 > 0 )
	 tres += res[nres++] = (sj->prev->vp->start - sj->vp->stop - 1);
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

   /*
   fprintf(stderr,"(%3d - %3d):(%3d - %3d)\n",
	   aln->min0,aln->max0,aln->min1,aln->max1);
   */
   
   /* now replace f_str->aa0t with aa0
      (f_str->aa0t is permanent, aa0 is not)*/
   for (i=0; i<n0; i++) f_str->aa0t[i] = aa0[i];

   return tres;
}

#include "a_mark.h"

int calcons(unsigned char *aa0, int n0,
	    unsigned char *aa1, int n1,
	    int *res, int nres, int *nc,
	    struct a_struct *aln, struct pstruct pst,
	    char *seqc0, char *seqc1, char *seqca,
	    struct f_struct *f_str)
{
  int i0, i1, nn1, n0t;
  int op, lenc, len_gap, nd, ns, itmp;
  unsigned char *aa1p;
  char *sp0, *sp1, *spa;
  int *rp;
  
#ifndef TFAST
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* first fill in the ends */
  aln->min0; aln->min1;
  n0 -= (f_str->nm0-1);

  if (min(aln->min0,aln->min1)<aln->llen || aln->showall==1)
    			/* will we show all the start ?*/
    if (aln->min0>=aln->min1) {                        /* aa0 extends more to left */
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

  memset(seqca,M_BLANK,aln->mins);

/* now get the middle */

  spa = seqca+aln->mins;
  sp0 = seqc0+aln->mins;
  sp1 = seqc1+aln->mins;
  rp = res;
  n0t = lenc = len_gap = aln->nident = aln->nsim = aln->ngap_q = aln->ngap_l = op = 0;
  i0 = aln->min0;
  i1 = aln->min1;
  
  /* op is the previous "match/insert" operator; *rp is the current
     operator or repeat count */

  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {	/* previous was match (or start), current is match */
      op = *rp++;		/* get the next match/insert operator */

      /* get the alignment symbol */
      if ((itmp=pst.pam2[0][f_str->aa0t[i0]][aa1p[i1]])<0) { *spa = M_NEG; }
      else if (itmp == 0) { *spa = M_ZERO;}
      else {*spa = M_POS;}
      if (*spa == M_ZERO || *spa == M_POS) { aln->nsim++;}

      *sp0 = pst.sq[f_str->aa0t[i0++]];	/* get the residues for the consensus */
      *sp1 = pst.sq[aa1p[i1++]];
      n0t++;
      lenc++;
      if (toupper(*sp0) == toupper(*sp1)) {aln->nident++; *spa = M_IDENT;}
      sp0++; sp1++; spa++;
    }
    else {	/* either op != 0 (previous was insert) or *rp != 0
		   (current is insert) */
      if (op==0) { op = *rp++;}	/* previous was match, start insert */
      				/* previous was insert - count through gap */
      *sp0++ = '-';
      *sp1++ = pst.sq[aa1p[i1++]];
      *spa++ = M_DEL;
      op--;
      len_gap++;
      lenc++;
    }
  }

  *spa = '\0';
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

int calcons_a(unsigned char *aa0, unsigned char *aa0a, int n0,
	      unsigned char *aa1, int n1,
	      int *res, int nres, int *nc,
	      struct a_struct *aln, struct pstruct pst,
	      char *seqc0, char *seqc0a, char *seqc1, char *seqca,
	      char *ann_arr, struct f_struct *f_str)
{
  int i0, i1, nn1, n0t;
  int op, lenc, len_gap, nd, ns, itmp, p_ac, fnum, o_fnum;
  unsigned char *aa1p, *aa0ap;
  char *sp0, *sp0a, *sp1, *spa;
  int *rp;
  
#ifndef TFAST
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* first fill in the ends */
  aln->min0; aln->min1;
  n0 -= (f_str->nm0-1);

  if (min(aln->min0,aln->min1)<aln->llen || aln->showall==1)
    			/* will we show all the start ?*/
    if (aln->min0>=aln->min1) {                        /* aa0 extends more to left */
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

  memset(seqca,M_BLANK,aln->mins);
  memset(seqc0a,' ', aln->mins);

/* now get the middle */

  spa = seqca+aln->mins;
  sp0 = seqc0+aln->mins;
  sp0a = seqc0a+aln->mins;
  sp1 = seqc1+aln->mins;
  rp = res;
  n0t=lenc=len_gap=aln->nident=aln->nsim=aln->ngap_q=aln->ngap_l=op=p_ac= 0;
  i0 = aln->min0;
  i1 = aln->min1;
  
  /* op is the previous "match/insert" operator; *rp is the current
     operator or repeat count */

  o_fnum = f_str->aa0ti[i0];
  aa0ap = &aa0a[f_str->nmoff[o_fnum]+i0];

  while (i0 < aln->max0 || i1 < aln->max1) {
    fnum = f_str->aa0ti[i0];
    if (op == 0 && *rp == 0) {	/* previous was match (or start), current is match */
      if (p_ac == 0) { /* previous code was a match */
	if (fnum != o_fnum) { /* continuing a match, but with a different fragment */
	  aa0ap = &aa0a[f_str->nmoff[fnum]];
	  o_fnum = fnum;
	}
      }
      else {
	p_ac = 0; o_fnum = fnum = f_str->aa0ti[i0];
	aa0ap = &aa0a[f_str->nmoff[fnum]];
      }
      op = *rp++;		/* get the next match/insert operator */

      /* get the alignment symbol */
      if ((itmp=pst.pam2[0][f_str->aa0t[i0]][aa1p[i1]])<0) { *spa = M_NEG; }
      else if (itmp == 0) { *spa = M_ZERO;}
      else {*spa = M_POS;}
      if (*spa == M_ZERO || *spa == M_POS) { aln->nsim++;}

      *sp0 = pst.sq[f_str->aa0t[i0++]];	/* get the residues for the consensus */
      *sp0a++ = ann_arr[*aa0ap++];
      *sp1 = pst.sq[aa1p[i1++]];
      n0t++;
      lenc++;
      if (toupper(*sp0) == toupper(*sp1)) {aln->nident++; *spa = M_IDENT;}
      sp0++; sp1++; spa++;
    }
    else {	/* either op != 0 (previous was insert) or *rp != 0
		   (current is insert) */
      if (op==0) { op = *rp++;}	/* previous was match, start insert */
      				/* previous was insert - count through gap */
      if (p_ac != 1) {
	p_ac = 1; fnum = f_str->aa0ti[i0];
      }

      *sp0++ = '-';
      *sp1++ = pst.sq[aa1p[i1++]];
      *spa++ = M_DEL;
      *sp0a++ = ' ';
      op--;
      len_gap++;
      lenc++;
    }
  }

  *sp0a = *spa = '\0';
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

void aaptrshuffle(unsigned char *res, int n) {

  int i, j;
  unsigned char tmp;

  for( i = n; --i; ) {

    /* j = nrand(i); if (i == j) continue; */ /* shuffle */
    j = (n - 1) - i; if (i <= j ) break; /* reverse */

    tmp = res[i];
    res[i] = res[j];
    res[j] = tmp;
  }
}

void aa0shuffle(unsigned char *aa0, int n0, struct f_struct *f_str) {

  int i;
  int j;

  for(i = 0 ; i < f_str->nm0 ; i++) { /* for each fragment */

    aaptrshuffle(&(aa0[f_str->nmoff[i]]), 
		 f_str->nmoff[i+1] - f_str->nmoff[i] - 1 );

  }

}

/* build an array of match/ins/del - length strings */
int calc_code(const unsigned char *aa0, const int n0,
	      const unsigned char *aa1, const int n1,
	      int *res, int nres,
	      struct a_struct *aln, struct pstruct pst,
	      char *al_str, int al_str_n, struct f_struct *f_str)
{
  int i0, i1, nn1;
  int op, lenc, len_gap, nd, ns, itmp;
  int p_ac, op_cnt;
  const unsigned char *aa1p;
  char tmp_cnt[20];
  char sp0, sp1, *sq;
  int *rp;
  int o_fnum,fnum = 0;

  if (pst.ext_sq_set) {
    sq = pst.sqx;
  }
  else {
    sq = pst.sq;
  }

#ifndef TFAST
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* first fill in the ends */

  rp = res;
  lenc = len_gap =aln->nident=aln->nsim=aln->ngap_q=aln->ngap_l=aln->nfs=op=p_ac = 0;
  op_cnt = 0;

  i0 = aln->min0;	/* start in aa0 (f_str->aa0t) */
  i1 = aln->min1;	/* start in aa1 */
  tmp_cnt[0]='\0';
  
  o_fnum = f_str->aa0ti[i0] + 1;	/* fragment number */
  while (i0 < aln->max0 || i1 < aln->max1) {
    fnum = f_str->aa0ti[i0]+1;
    if (op == 0 && *rp == 0) {	/* previous was match, this is match */
      if (p_ac == 0) {	/* previous code was a match */
	if (fnum == o_fnum) { op_cnt++;}
	else {		/* continuing a match, but with a different fragment */
	  update_code(al_str,al_str_n-strlen(al_str), p_ac, op_cnt, o_fnum);
	  o_fnum = fnum;
	  op_cnt=1;
	}
      }
      else {
	update_code(al_str,al_str_n-strlen(al_str),p_ac,op_cnt,o_fnum);
	op_cnt = 1; p_ac = 0; o_fnum = fnum = f_str->aa0ti[i0] + 1;
      }
      op = *rp++;
      lenc++;
      if (pst.pam2[0][f_str->aa0t[i0]][aa1p[i1]]>=0) {aln->nsim++;}
      sp0 = pst.sq[f_str->aa0t[i0++]];
      sp1 = pst.sq[aa1p[i1++]];
      if (toupper(sp0) == toupper(sp1)) aln->nident++;
    }
    else {
      if (op==0) op = *rp++;
      if (p_ac == 1) { op_cnt++;}
      else {
	update_code(al_str,al_str_n - strlen(al_str),p_ac,op_cnt,o_fnum);
	op_cnt = 1; p_ac = 1; fnum = f_str->aa0ti[i0] + 1;
      }
      op--; lenc++; i1++; len_gap++;
    }
  }
  update_code(al_str,al_str_n - strlen(al_str),p_ac,op_cnt,o_fnum);

  return lenc - len_gap;
}

/* update_code(): if "op" == 0, this is the end of a match of length
   "op_cnt" involving fragment "fnum"
   otherwise, this is an insertion (op==1) or deletion (op==2)
*/

update_code(char *al_str, int al_str_max, int op, int op_cnt, int fnum) {

  char op_char[4]={"=-+"};
  char tmp_cnt[20];

  if (op == 0)
    sprintf(tmp_cnt,"%c%d[%d]",op_char[op],op_cnt,fnum);
  else
    sprintf(tmp_cnt,"%c%d",op_char[op],op_cnt);

  strncat(al_str,tmp_cnt,al_str_max);
}

int calc_id(unsigned char *aa0, int n0,
	    unsigned char *aa1, int n1,
	    int *res, int nres,
	    struct a_struct *aln, struct pstruct pst,
	    struct f_struct *f_str)
{
  int i0, i1, nn1;
  int op, lenc, len_gap, nd, ns, itmp;
  unsigned char *aa1p;
  int sp0, sp1;
  int *rp;
  
#ifndef TFAST
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
  lenc=len_gap=aln->nident=aln->nsim=aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  i0 = aln->min0;
  i1 = aln->min1;
  
  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;

      if (pst.pam2[0][f_str->aa0t[i0]][aa1p[i1]]>=0) {aln->nsim++;}

      sp0 = pst.sq[f_str->aa0t[i0++]];
      sp1 = pst.sq[aa1p[i1++]];
      lenc++;
      if (toupper(sp0) == toupper(sp1)) aln->nident++;
    }
    else {
      if (op==0) { op = *rp++;}
      i1++;
      op--;
      len_gap++;
      lenc++;
    }
  }
  return lenc-len_gap;
}

#ifdef PCOMPLIB

#include "structs.h"
#include "p_mw.h"

void
update_params(struct qmng_str *qm_msg,
	      struct mngmsg *m_msg, struct pstruct *ppst)
{
  m_msg->n0 = ppst->n0 = qm_msg->n0;
  m_msg->nm0 = qm_msg->nm0;
  m_msg->escore_flg = qm_msg->escore_flg;
  m_msg->qshuffle = qm_msg->qshuffle;
}
#endif
