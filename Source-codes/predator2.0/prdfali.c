#include "prd.h"
#include "prdprot.h"


#define MIN(x,y) ((x)<=(y) ? (x) : (y))
#define EOSEQ 127
#define MAXSQ 32
#define MAXTST 10000
#define MAXLIB 50000l
#define MAXDIAG (MAXTST+MAXLIB)
#define NA 124
#define EL 125
#define ES 126
#define AAMASK 127
#define PAIRNULL (pairptr)NULL
static int  zero = 0;


int maxv, max0, max1;		/* best match values */
int min0, min1;			/* left end of match */
int minc, maxc;			/* consensus pointers */
int mins;
char *seqc0, *seqc1;	/* arrays for the consensus sequence */
int salloc=0;		/* flag for initseq */
int gdelval, ggapval;
char *sq;
char aa[MAXSQ] = {"ARNDCQEGHILKMFPSTWYVBZX"};
int naa = 23;
int nsq;
int haa[MAXSQ] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,6,0};
int *hsq;
char nt[MAXSQ]={"ACGTURYMWSKDHVBNX"};
int nnt = 17;
int hnt[MAXSQ] = {0,1,2,3,3,0,1,0,0,1,2,0,0,0,1,0,0};
int *pam;
int pam2[MAXSQ][MAXSQ];
int pamh1[MAXSQ];		/* used for kfact replacement */
int *sascii;
char *aa0=NULL, *aa1=NULL;	/* amino acid sequence data */
long sq0off=1, sq1off=1;
int bktup, bkfact, scfact, bestoff, bestscale, histint, bestmax;
int nseq=2;		/* same sequence twice */
int nn0, nn1;	/* length of aa0, length of aa1, nn0+nn1,	diagonal offset */
int match, mismh;
int smin0, smin1, smins;	/* set bounds for discons */
static int (*vv)[32];			/* substitution scores */
static int q, r;			/* gap penalties */
static int qr;				/* qr = q + r */
typedef struct ONE { int COL ;  struct ONE  *NEXT ; struct ONE *PREV; } pair, *pairptr;
pairptr *row, z, zz; /* for saving used aligned pairs */
static int FIRSTZ;

static int tt;
typedef struct NODE
	{ int  SCORE;
	  int  STARI;
	  int  STARJ;
	  int  ENDI;
	  int  ENDJ;
	  int  TOP;
	  int  BOT;
	  int  LEFT;
	  int  RIGHT; }  vertex,
*vertexptr;
vertexptr  *LIST;			/* an array for saving k best scores */
vertexptr  low = 0;			/* lowest score node in LIST */
vertexptr  most = 0;			/* latestly accessed node in LIST */
static int numnode;			/* the number of nodes in LIST */
static int *CC, *DD;			/* saving matrix scores */
static int *RR, *SS, *EE, *FF; 		/* saving start-points */
static int *HH, *WW;		 	/* saving matrix scores */
static int *II, *JJ, *XX, *YY; 		/* saving start-points */
static int  m1, mm, n1, nn;		/* boundaries of recomputed area */
static int  rl, cl;			/* left and top boundaries */
static int  lmin;			/* minimum score in LIST */
static int flag;			/* indicate if recomputation necessary*/



/* DIAG() assigns value to x if (ii,jj) is never used before */
#define DIAG(ii, jj, x, value)				\
{ for ( tt = 1, z = row[(ii)]; z != PAIRNULL; z = z->NEXT )	\
    if ( z->COL == (jj) )				\
      { tt = 0; break; }				\
  if ( tt )						\
    x = ( value );					\
}


/* replace (ss1, xx1, yy1) by (ss2, xx2, yy2) if the latter is large */
#define ORDER(ss1, xx1, yy1, ss2, xx2, yy2)		\
{ if ( ss1 < ss2 )					\
    { ss1 = ss2; xx1 = xx2; yy1 = yy2; }		\
  else							\
    if ( ss1 == ss2 )					\
      { if ( xx1 < xx2 )				\
	  { xx1 = xx2; yy1 = yy2; }			\
	else						\
	  if ( xx1 == xx2 && yy1 < yy2 )		\
	    yy1 = yy2;					\
      }							\
}


/* The following definitions are for function diff() */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

static int I, J;				/* current positions of A ,B */
static int no_mat; 				/* number of matches */ 
static int no_mis; 				/* number of mismatches */ 
static int al_len; 				/* length of alignment */
						/* Append "Delete k" op */
#define DEL(k)				\
{ I += k;				\
  al_len += k;				\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ J += k;				\
  al_len += k;				\
  if (last < 0)				\
    { sapp[-1] = (k); *sapp++ = last; }	\
  else					\
	 last = *sapp++ = (k);		\
}

						/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
  al_len += 1;				\
}

#define FCKALLOC ckalloc


static int abl50[] = {
  5,
 -2, 7,
 -1,-1, 7,
 -2,-2, 2, 8,
 -1,-4,-2,-4,13,
 -1, 1, 0, 0,-3, 7,
 -1, 0, 0, 2,-3, 2, 6,
  0,-3, 0,-1,-3,-2,-3, 8,
 -2, 0, 1,-1,-3, 1, 0,-2,10,
 -1,-4,-3,-4,-2,-3,-4,-4,-4, 5,
 -2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5,
 -1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6,
 -1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7,
 -3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8,
 -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,
  1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5,
  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5,
 -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15,
 -2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8,
  0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5,
 -2,-1, 4, 5,-3, 0, 1,-1, 0,-4,-4, 0,-3,-4,-2, 0, 0,-5,-3,-4, 5,
 -1, 0, 0, 1,-3, 4, 5,-2, 0,-3,-3, 1,-1,-4,-1, 0,-1,-2,-2,-3, 2, 5,
 -1,-1,-1,-1,-2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-1, 0,-3,-1,-1,-1,-1,-1};

/*      0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15	*/
int aascii[]={
	EL,NA,NA,NA,NA,NA,NA,NA,NA,NA,EL,NA,NA,EL,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ES,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	ES, 0,20, 4, 3, 6,13, 7, 8, 9,NA,11,10,12, 2,NA,
	14, 5, 1,15,16,NA,19,17,22,18,21,NA,NA,NA,NA,NA,
	NA, 0,20, 4, 3, 6,13, 7, 8, 9,NA,11,10,12, 2,NA,
	14, 5, 1,15,16,NA,19,17,22,18,21,NA,NA,NA,NA,NA};


void Align(char *Seq1, int Len1, char *Seq2, int Len2, ALN *Aln,
           int DelVal, int GapVal, int NAli, int MinAliLength, float Safety)
{
  static int FirstCall = TRUE;

  InitAll();

  initenv();

  aa0=(char *)ckalloc((int)(MAXTST+MAXLIB)*sizeof(char));

  strcpy(aa0,Seq1);
  nn0 = Len1;
  trick(aa0,nn0);

  aa1 = aa0 + nn0 + 2;
  
  strcpy(aa1,Seq2);
  nn1=Len2;
  trick(aa1,nn1); 
  
  initseq(MIN(nn0,nn1)*2);
  if( FirstCall ) {
    initpam2();			/* convert 1-d pam to 2-d pam2 */
    FirstCall = FALSE;
  }

  gdelval = DelVal;
  ggapval = GapVal;

  SIM(aa0-1,aa1-1,nn0,nn1,pam2,-(gdelval-ggapval),-ggapval,nseq,Aln,NAli,MinAliLength,Safety);

  if( Aln->NAli ) 
    Unlink(Aln->Beg1,Aln->Beg2,Aln->End1,Aln->End2,Aln->Len,Aln->NAli,Aln->Accepted);

  free(aa0);
  freesq();
}

void initenv(void)
{
  sascii = aascii;
  pam = abl50;
  sq = aa;
  hsq = haa;
  nsq = naa;
}


void initpam2(void)
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

int calcons(char *aa0, char *aa1, int *res, int *nc, int *nident)
{
  int i0, i1;
  int op, nid, lenc, nd;
  /* int ns */
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
      else if ( /*(dnaseq==1) && (*sp0=='T' && *sp1=='U') ||*/
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

void initseq(int seqsiz)		/* initialize arrays */
{
  seqc0=calloc((size_t)seqsiz,sizeof(char));
  seqc1=calloc((size_t)seqsiz,sizeof(char));
  if (seqc0==NULL || seqc1==NULL) {
    fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
    exit(1);
  }
}

void freesq(void)
{
  free(seqc0);
  free(seqc1);
}


void trick(char *Seq, int Len)
{
  register int i;

  for( i=0; i<Len ; i++ )
    Seq[i] = sascii[Seq[i]&AAMASK];
}

/* SIM(A,B,M,N,K,V,Q,R) reports K best non-intersecting alignments of
	the segments of A and B in order of similarity scores, where
	V[a][b] is the score of aligning a and b, and -(Q+R*i) is the score
	of an i-symbol indel.  						*/

void SIM(char *A, char *B, int M, int N, int V[][32], int Q, int R, int nseq, ALN *Aln,
         int NAli, int MinAliLength, float Safety)
{
  int endi, endj, stari, starj;	/* endpoint and startpoint */ 
  int count;				/* maximum size of list */	
  register  int  i, j, iii;			/* row and column indices */
  int  *S;				/* saving operations for diff */
  int ns, nident;		/* for display */
  vertexptr cur; 			/* temporary pointer */
  vertexptr findmax(void);	 		/* return the largest score node */
  float InfThreshold, INF;
  
  Aln->NAli = 0;

  /* allocate space for all vectors */
  j = (int)((N + 1) * sizeof(int));
  CC = ( int * ) ckalloc(j);
  DD = ( int * ) ckalloc(j);
  RR = ( int * ) ckalloc(j);
  SS = ( int * ) ckalloc(j);
  EE = ( int * ) ckalloc(j);
  FF = ( int * ) ckalloc(j);
  i = (int)((M + 1) * sizeof(int));
  HH = ( int * ) ckalloc(i);
  WW = ( int * ) ckalloc(i);
  II = ( int * ) ckalloc(i);
  JJ = ( int * ) ckalloc(i);
  XX = ( int * ) ckalloc(i);
  YY = ( int * ) ckalloc(i);
  S = ( int * ) ckalloc(MIN(i,j)*5/4);
  row = ( pairptr * ) ckalloc( (M + 1) * sizeof(pairptr));
  FIRSTZ = TRUE;
  
  /* set up list for each row */
  if (nseq == 2) for ( i = 1; i <= M; i++ ) row[i]= PAIRNULL;

  vv = V;
  q = Q;
  r = R;
  qr = q + r;

  LIST = ( vertexptr * ) ckalloc( NAli * sizeof(vertexptr));
  for ( i = 0; i < NAli ; i++ )
    LIST[i] = ( vertexptr ) FCKALLOC( (int) sizeof(vertex));
  
  numnode = lmin = 0;
  big_pass(A,B,M,N,NAli,nseq);
  
  /* Report the K best alignments one by one. After each alignment is
     output, recompute part of the matrix. First determine the size
     of the area to be recomputed, then do the recomputation         */
  
  for ( count = NAli - 1; count >= 0 ; count-- )
    { 
      if ( numnode == 0 ) {
	fprintf(stderr,"The number of alignments computed is too large\n");
        exit(0);
      }

      cur = findmax();	/* Return a pointer to a node with max score*/
      Aln->Score[Aln->NAli] = cur->SCORE;
      stari = ++cur->STARI;
      starj = ++cur->STARJ;
      endi = cur->ENDI;
      endj = cur->ENDJ;
      m1 = cur->TOP;
		mm = cur->BOT;
      n1 = cur->LEFT;
      nn = cur->RIGHT;
      rl = endi - stari + 1;
      cl = endj - starj + 1;
      I = stari - 1;
      J = starj - 1;
      sapp = S;
      last = 0;
      al_len = 0;
		no_mat = 0;
      no_mis = 0;
      diff(&A[stari]-1, &B[starj]-1,rl,cl,q,q);

		min0 = stari;
      min1 = starj;
      max0 = stari+rl-1;
      max1 = starj+cl-1;
      ns=calcons(A+1,B+1,S,&Aln->Nc[Aln->NAli],&nident);
      Aln->Percent[Aln->NAli] = nident*100.0/(float)(Aln->Nc[Aln->NAli]);
      INF = Inf(Aln->Percent[Aln->NAli]/100.0);
      if( Aln->Percent[Aln->NAli] > 50.0 )
        InfThreshold = 0.0;
      else
        InfThreshold = INFTHR;

      if( ns >= MinAliLength && INF > InfThreshold && 
/*      if( ns >= MinAliLength &&*/
          ( (Aln->HCri[Aln->NAli] = HsspCrit(seqc0,seqc1,ns)) + Safety) < Aln->Percent[Aln->NAli] ) {
	Aln->Beg1[Aln->NAli] = min0;
	Aln->Beg2[Aln->NAli] = min1;
	Aln->End1[Aln->NAli] = max0;
	Aln->End2[Aln->NAli] = max1;
	Aln->Len[Aln->NAli] = ns;
	for( iii=0; iii<Aln->Len[Aln->NAli]; iii++ ) {
	  Aln->A1[Aln->NAli][iii] = seqc0[iii];
	  Aln->A2[Aln->NAli][iii] = seqc1[iii];
	}
	Aln->NAli++;
      }
      if ( count )
	{ flag = 0;
	  locate(A,B,nseq);
	  if ( flag )
	    small_pass(A,B,count,nseq);
	}
    }

  free(CC);
  free(DD);
  free(RR);
  free(SS);
  free(EE);
  free(FF);
  free(HH);
  free(WW);
  free(II);
  free(JJ);
  free(XX);
  free(YY);
  free(S);
  free(row);
  for ( i = 0; i < NAli ; i++ )
	 free(LIST[i]);
  free(LIST);

  do {
    /* Update 2.1 */
    if( z->PREV != NULL )
    /* End of pdate 2.1 */
      zz = z->PREV;
    free(z);
    z=zz;
  } while( z->PREV != NULL );
  free(z);

}

/* A big pass to compute K best classes */

void big_pass(char *A, char *B, int M, int N, int K, int nseq)
/*char A[],B[]; int M,N,K,nseq;*/
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  *va;				/* pointer to vv(A[i], B[j]) */


	/* Compute the matrix and save the top K best scores in LIST
		CC : the scores of the current row
	   RR and EE : the starting point that leads to score CC
	   DD : the scores of the current row, ending with deletion
	   SS and FF : the starting point that leads to score DD        */
 	/* Initialize the 0 th row */
  for ( j = 1; j <= N ; j++ )
  {  CC[j] = 0;
  RR[j] = 0;
  EE[j] = j;
  DD[j] = - (q);
  SS[j] = 0;
  FF[j] = j;
  }
  for ( i = 1; i <= M; i++) 
  {  c = 0;				/* Initialize column 0 */
  f = - (q);
  ci = fi = i;
  va = vv[A[i]];
  if ( nseq == 2 )
  { p = 0;
  pi = i - 1;
  cj = fj = pj = 0;
  }
  else
  { p = CC[i];
  pi = RR[i];
  pj = EE[i];
  cj = fj = i;
  }
  for ( j = (nseq == 2 ? 1 : (i+1)) ; j <= N ; j++ )  
  {  f = f - r;
  c = c - qr;
  ORDER(f, fi, fj, c, ci, cj)
    c = CC[j] - qr; 
  ci = RR[j];
  cj = EE[j];
  d = DD[j] - r;
  di = SS[j];
  dj = FF[j];
  ORDER(d, di, dj, c, ci, cj)
    c = 0;
  DIAG(i, j, c, p+va[B[j]])		/* diagonal */
    if ( c <= 0 )
    { c = 0; ci = i; cj = j; }
    else
    { ci = pi; cj = pj; }
  ORDER(c, ci, cj, d, di, dj)
    ORDER(c, ci, cj, f, fi, fj)
    p = CC[j];
  CC[j] = c;

  pi = RR[j];
  pj = EE[j];
  RR[j] = ci;
  EE[j] = cj;
  DD[j] = d;
  SS[j] = di;
  FF[j] = dj;
  if ( c > lmin )	/* add the score into list */
    lmin = addnode(c, ci, cj, i, j, K, lmin);
  }
  }
}

/* Determine the left and top boundaries of the recomputed area */

void locate(char *A, char *B, int nseq)
/* char A[],B[]; int nseq; */
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  cflag, rflag;			/* for recomputation */
  int  *va;				/* pointer to vv(A[i], B[j]) */
  int  limit;				/* the bound on j */

	/* Reverse pass
		rows
		CC : the scores on the current row
		RR and EE : the endpoints that lead to CC
		DD : the deletion scores
		SS and FF : the endpoints that lead to DD

		columns
		HH : the scores on the current columns
		II and JJ : the endpoints that lead to HH
		WW : the deletion scores
		XX and YY : the endpoints that lead to WW
	*/
  for ( j = nn; j >= n1 ; j-- )
  {  CC[j] = 0;
  EE[j] = j;
  DD[j] = - (q);
  FF[j] = j;
  if ( nseq == 2 || j > mm )
    RR[j] = SS[j] = mm + 1;
  else
    RR[j] = SS[j] = j;
  }

  for ( i = mm; i >= m1; i-- )
  {  c = p = 0;
  f = - (q);
  ci = fi = i;
  pi = i + 1;
  cj = fj = pj = nn + 1;
  va = vv[A[i]];
  if ( nseq == 2 || n1 > i )
    limit = n1;
  else
    limit = i + 1;
  for ( j = nn; j >= limit ; j-- )
  {  f = f - r;
  c = c - qr;
  ORDER(f, fi, fj, c, ci, cj)
    c = CC[j] - qr;
  ci = RR[j];
  cj = EE[j];
  d = DD[j] - r;
  di = SS[j];
  dj = FF[j];
  ORDER(d, di, dj, c, ci, cj)
    c = 0;
  DIAG(i, j, c, p+va[B[j]])		/* diagonal */
    if ( c <= 0 )
    { c = 0; ci = i; cj = j; }
    else
    { ci = pi; cj = pj; }
  ORDER(c, ci, cj, d, di, dj)
    ORDER(c, ci, cj, f, fi, fj)
    p = CC[j];
  CC[j] = c;
  pi = RR[j];
  pj = EE[j];
  RR[j] = ci;
  EE[j] = cj;
  DD[j] = d;
  SS[j] = di;
  FF[j] = dj;
  if ( c > lmin )
    flag = 1;
  }
  if ( nseq == 2 || i < n1 )
  { HH[i] = CC[n1];
  II[i] = RR[n1];
  JJ[i] = EE[n1];
  WW[i] = DD[n1];
  XX[i] = SS[n1];
  YY[i] = FF[n1];
  }
  }

  for ( rl = m1, cl = n1; ; )
  { for ( rflag = cflag = 1; ( rflag && m1 > 1 ) || ( cflag && n1 > 1 ) ;  )
  { if ( rflag && m1 > 1 )	/* Compute one row */
  { rflag = 0;
  m1--;
  c = p = 0;
  f = - (q);
  ci = fi = m1;
  pi = m1 + 1;
  cj = fj = pj = nn + 1;
  va = vv[A[m1]];
  for ( j = nn; j >= n1 ; j-- )
  { f = f - r;
  c = c - qr;
  ORDER(f, fi, fj, c, ci, cj)
    c = CC[j] - qr;
  ci = RR[j];
  cj = EE[j];
  d = DD[j] - r;
  di = SS[j];
  dj = FF[j];
  ORDER(d, di, dj, c, ci, cj)
    c = 0;
  DIAG(m1, j, c, p+va[B[j]])		/* diagonal */
    if ( c <= 0 )
    { c = 0; ci = m1; cj = j; }
    else
    { ci = pi; cj = pj; }
  ORDER(c, ci, cj, d, di, dj)
    ORDER(c, ci, cj, f, fi, fj)
    p = CC[j];
  CC[j] = c;
  pi = RR[j];
  pj = EE[j];
  RR[j] = ci;
  EE[j] = cj;
  DD[j] = d;
  SS[j] = di;
  FF[j] = dj;
  if ( c > lmin )
    flag = 1;
  if ( ! rflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
                    || (fi > rl && fj > cl )) )
    rflag = 1;
  }
  HH[m1] = CC[n1];
  II[m1] = RR[n1];
  JJ[m1] = EE[n1];
  WW[m1] = DD[n1];
  XX[m1] = SS[n1];
  YY[m1] = FF[n1];
  if ( ! cflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
                    || (fi > rl && fj > cl) ) )
    cflag = 1;
  }

  if ( nseq == 1 && n1 == (m1 + 1) && ! rflag )
    cflag = 0;
  if ( cflag && n1 > 1 )	/* Compute one column */
  { cflag = 0;
  n1--;
  c = 0;
  f = - (q);
  cj = fj = n1;
  va = vv[B[n1]];
  if ( nseq == 2 || mm < n1 )
  { p = 0;
  ci = fi = pi = mm + 1;
  pj = n1 + 1;
  limit = mm;
  }
  else
  { p = HH[n1];
  pi = II[n1];
  pj = JJ[n1];
  ci = fi = n1;
  limit = n1 - 1;
  }
  for ( i = limit; i >= m1 ; i-- )
  { f = f - r;
  c = c - qr;
  ORDER(f, fi, fj, c, ci, cj)
    c = HH[i] - qr;
  ci = II[i];
  cj = JJ[i];
  d = WW[i] - r;
  di = XX[i];
  dj = YY[i];
  ORDER(d, di, dj, c, ci, cj)
    c = 0;
  DIAG(i, n1, c, p+va[A[i]])
    if ( c <= 0 )
    { c = 0; ci = i; cj = n1; }
    else
    { ci = pi; cj = pj; }
  ORDER(c, ci, cj, d, di, dj)
    ORDER(c, ci, cj, f, fi, fj)
    p = HH[i];
  HH[i] = c;
  pi = II[i];
  pj = JJ[i];
  II[i] = ci;
  JJ[i] = cj;
  WW[i] = d;
  XX[i] = di;
  YY[i] = dj;
  if ( c > lmin )
    flag = 1;
  if ( ! cflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
                    || (fi > rl && fj > cl) ) )
    cflag = 1;
  }
  CC[n1] = HH[m1];
  RR[n1] = II[m1];
  EE[n1] = JJ[m1];
  DD[n1] = WW[m1];
  SS[n1] = XX[m1];
  FF[n1] = YY[m1];
  if ( ! rflag && ( (ci > rl && cj > c) || (di > rl && dj > cl)
                    || (fi > rl && fj > cl) ) )
    rflag = 1;
  }
  }
  if ( (m1 == 1 && n1 == 1) || no_cross() )
    break;
  }
  m1--;
  n1--;
}

/* recompute the area on forward pass */
void small_pass(char *A, char *B, int count, int nseq)
/* char A[], B[]; int count, nseq;*/
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  *va;				/* pointer to vv(A[i], B[j]) */
  int  limit;				/* lower bound on j */

  for ( j = n1 + 1; j <= nn ; j++ )
  {  CC[j] = 0;
  RR[j] = m1;
  EE[j] = j;
  DD[j] = - (q);
  SS[j] = m1;
  FF[j] = j;
  }
  for ( i = m1 + 1; i <= mm; i++) 
  {  c = 0;				/* Initialize column 0 */
  f = - (q);
  ci = fi = i;
  va = vv[A[i]];
  if ( nseq == 2 || i <= n1 )
  { p = 0;
  pi = i - 1;
  cj = fj = pj = n1;
  limit = n1 + 1;
  }
  else
  { p = CC[i];
  pi = RR[i];
  pj = EE[i];
  cj = fj = i;
  limit = i + 1;
  }
  for ( j = limit ; j <= nn ; j++ )  
  {  f = f - r;
  c = c - qr;
  ORDER(f, fi, fj, c, ci, cj)
    c = CC[j] - qr;
  ci = RR[j];
  cj = EE[j];
  d = DD[j] - r;
  di = SS[j];
  dj = FF[j];
  ORDER(d, di, dj, c, ci, cj)
    c = 0;
  DIAG(i, j, c, p+va[B[j]])		/* diagonal */
    if ( c <= 0 )
    { c = 0; ci = i; cj = j; }
    else
    { ci = pi; cj = pj; }
  ORDER(c, ci, cj, d, di, dj)
    ORDER(c, ci, cj, f, fi, fj)
    p = CC[j];
  CC[j] = c;
  pi = RR[j];
  pj = EE[j];
  RR[j] = ci;
  EE[j] = cj;
  DD[j] = d;
  SS[j] = di;
  FF[j] = dj;
  if ( c > lmin )	/* add the score into list */
    lmin = addnode(c, ci, cj, i, j, count, lmin);
  }
  }
}

/* Add a new node into list.  */

int addnode(int c, int ci, int cj, int i, int j, int K, int cost)
/*int c, ci, cj, i, j, K, cost;*/
{ int found;				/* 1 if the node is in LIST */
  register int d;

  found = 0;
  if ( most != 0 && most->STARI == ci && most->STARJ == cj )
    found = 1;
  else
     for ( d = 0; d < numnode ; d++ )
	{ most = LIST[d];
	  if ( most->STARI == ci && most->STARJ == cj )
	    { found = 1;
	      break;
		 }
        }
  if ( found )
    { if ( most->SCORE < c )
        { most->SCORE = c;
          most->ENDI = i;
          most->ENDJ = j;
        }
      if ( most->TOP > i ) most->TOP = i;
      if ( most->BOT < i ) most->BOT = i;
		if ( most->LEFT > j ) most->LEFT = j;
      if ( most->RIGHT < j ) most->RIGHT = j;
    }
  else
	 { if ( numnode == K )	/* list full */
	 most = low;
      else
         most = LIST[numnode++];
      most->SCORE = c;
      most->STARI = ci;
      most->STARJ = cj;
      most->ENDI = i;
      most->ENDJ = j;
      most->TOP = most->BOT = i;
		most->LEFT = most->RIGHT = j;
    }
  if ( numnode == K )
    { if ( low == most || ! low ) 
		  { for ( low = LIST[0], d = 1; d < numnode ; d++ )
            if ( LIST[d]->SCORE < low->SCORE )
              low = LIST[d];
	}
      return ( low->SCORE ) ;
    }
  else
    return cost;
}

/* Find and remove the largest score in list */

vertexptr findmax(void)
{ vertexptr  cur;
  register int i, j;

  for ( j = 0, i = 1; i < numnode ; i++ )
    if ( LIST[i]->SCORE > LIST[j]->SCORE )
       j = i;
  cur = LIST[j];
  if ( j != --numnode )
    { LIST[j] = LIST[numnode];
		LIST[numnode] =  cur;
    }
  most = LIST[0];
  if ( low == cur ) low = LIST[0];
  return ( cur );
}

/* return 1 if no node in LIST share vertices with the area */

int no_cross(void)
{ vertexptr  cur;
  register int i;

      for ( i = 0; i < numnode; i++ )
	{ cur = LIST[i];
	  if ( cur->STARI <= mm && cur->STARJ <= nn && cur->BOT >= m1-1 && 
			 cur->RIGHT >= n1-1 && ( cur->STARI < rl || cur->STARJ < cl ))
	     { if ( cur->STARI < rl ) rl = cur->STARI;
	       if ( cur->STARJ < cl ) cl = cur->STARJ;
	       flag = 1;
			 break;
	     }
	}
		if ( i == numnode )
	return TRUE;
      else
	return FALSE;
}

/* diff(A,B,M,N,tb,te) returns the score of an optimum conversion between
	A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

int diff(char *A, char *B, int M, int N, int tb, int te)
/* char *A, *B; int M, N; int tb, te;*/

{ int   midi, midj, type;	/* Midpoint, type, and cost */
int midc;

{ register int   i, j;
register int c, e, d, s;
int t, *va;

/* Boundary cases: M <= 1 or N == 0 */

if (N <= 0)
{ if (M > 0) DEL(M)
               return - gap(M);
}
if (M <= 1)
{ if (M <= 0)
{ INS(N);
return - gap(N);
}
if (tb > te) tb = te;
midc = - (tb + r + gap(N) );
midj = 0;
va = vv[A[1]];
for (j = 1; j <= N; j++)
{  for ( tt = 1, z = row[I+1]; z != PAIRNULL; z = z->NEXT )	
  if ( z->COL == j+J )			
  { tt = 0; break; }		
if ( tt )			
{ c = va[B[j]] - ( gap(j-1) + gap(N-j) );
if (c > midc)
{ midc = c;
midj = j;
}
}
}
if (midj == 0)
{ INS(N) DEL(1) }
else
{ if (midj > 1) INS(midj-1)
                  REP
                  if ( A[1] == B[midj] )
                    no_mat += 1;
                  else
                    no_mis += 1;
/* mark (A[I],B[J]) as used: put J into list row[I] */	
I++; J++;
z = ( pairptr ) FCKALLOC( (int) sizeof(pair));
z->COL = J;			
z->NEXT = row[I];
if( FIRSTZ ) {
  zz = NULL;
  FIRSTZ = FALSE;
}
row[I] = z;
z->PREV = zz;
zz = z;
if (midj < N) INS(N-midj)
                }
return midc;
}

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

midi = M/2;			/* Forward phase:                          */
CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
t = -q;
for (j = 1; j <= N; j++)
{ CC[j] = t = t-r;
DD[j] = t-q;
}
t = -tb;
for (i = 1; i <= midi; i++)
{ s = CC[0];
CC[0] = c = t = t-r;
e = t-q;
va = vv[A[i]];
for (j = 1; j <= N; j++)
{ if ((c = c - qr) > (e = e - r)) e = c;
if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
DIAG(i+I, j+J, c, s+va[B[j]])
  if (c < d) c = d;
if (c < e) c = e;
s = CC[j];
CC[j] = c;
DD[j] = d;
}
}
DD[0] = CC[0];

RR[N] = 0;			/* Reverse phase:                          */
t = -q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
for (j = N-1; j >= 0; j--)
{ RR[j] = t = t-r;
SS[j] = t-q;
}
t = -te;
for (i = M-1; i >= midi; i--)
{ s = RR[N];
RR[N] = c = t = t-r;
e = t-q;
va = vv[A[i+1]];
for (j = N-1; j >= 0; j--)
{ if ((c = c - qr) > (e = e - r)) e = c;
if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
DIAG(i+1+I, j+1+J, c, s+va[B[j+1]])
  if (c < d) c = d;
if (c < e) c = e;
s = RR[j];
RR[j] = c;
SS[j] = d;
}
}
SS[N] = RR[N];

midc = CC[0]+RR[0];		/* Find optimal midpoint */
midj = 0;
type = 1;
for (j = 0; j <= N; j++)
if ((c = CC[j] + RR[j]) >= midc)
if (c > midc || (CC[j] != DD[j] && RR[j] == SS[j]))
{ midc = c;
midj = j;
}
for (j = N; j >= 0; j--)
if ((c = DD[j] + SS[j] + q) > midc)
{ midc = c;
midj = j;
type = 2;
}
}

/* Conquer: recursively around midpoint */

if (type == 1)
{ diff(A,B,midi,midj,tb,q);
diff(A+midi,B+midj,M-midi,N-midj,q,te);
}
else
{ diff(A,B,midi-1,midj,tb,zero);
DEL(2);
diff(A+midi+1,B+midj,M-midi-1,N-midj,zero,te);
}
return midc;
}

void InitAll(void)
{
  naa = 23;
  nnt = 17;
  sq0off=1;
  sq1off=1;
  nseq=2;
  low = 0;
  most = 0;
  maxv = max0 = max1 = min0 = min1 = minc = maxc = mins = salloc = gdelval = ggapval = 0;
  nsq = bktup = bkfact = scfact = bestoff = bestscale = histint = bestmax = nn0 = nn1 = 0;
  match = mismh = smin0 = smin1 = smins = 0;
  q = r = qr = tt = numnode = m1 = mm = n1 = nn = rl = cl = lmin = flag = I = J = last = 0;
  no_mat = no_mis = al_len = 0;
}

