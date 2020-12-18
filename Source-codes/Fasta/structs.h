
/* $Name: fa34t25d2 $ - $Id: structs.h,v 1.32 2004/11/10 21:06:04 wrp Exp $ */

#ifndef A_STRUCT
#define A_STRUCT
/* this structure keeps fixed alignment parameters */
/* must be identical to same structure in param.h. */
struct a_struct {
  int min0, max0;	/* boundaries of alignment in aa0 */
  int min1, max1;	/* boundaries of alignment in aa1 */
  int smin0, smin1, smins;
  int mins;		/* start of aligned region in seqc0,seqc1,seqca */
  int llen, llcntx, llcntx_flg, showall;
  int qlrev;		/* query is reversed */
  int qlfact;		/* query display increment/decrement (1/3) */
  int llrev;		/* library is reversed */
  int llfact;		/* library display increment/decrement (1/3) */
  int llmult;		/* library factor for coord. display (1/3) */
  int frame;		/* used for tfasta/s/f offsets */
  int a_len;		/* consensus alignment length */
  int nident, nsim, ngap_q, ngap_l, nfs;	/* number of identities, gaps in q, l */
  long d_start0,d_stop0;	/* actual coordinates for display */
  long d_start1,d_stop1;	/* which include various offsets */
};
#endif

struct hist_str {
  int histflg;
  int *hist_a;
  int histint, min_hist, max_hist, maxh;
  long entries;
  int z_calls;
  char stat_info[MAX_STR];
};

struct db_str {
  long entries;
  unsigned long length;
  int carry;
};

struct mngmsg 		/* Message from host to manager */
{
  int n0;		/* Integer returned by hgetseq */
  int nm0;		/* number of segments */
  int nmoff;		/* length of fastf segment */
  unsigned char *aa0a;	/* annotation array */
  char ann_arr[MAX_FN];	/* annotation characters */
  int ann_flg;		/* have annotation array, characters */
  char tname[MAX_FN];	/* Query sequence name */
  int tnamesize;	/* Query name size */
  int qsfnum[10];
  int nqsfnum;
  int qsfnum_n[10];
  int nqsfnum_n;
  char lname[MAX_FN];	/* Library  file  name */
  char *lbnames[MAX_LF]; /* list of library files */
  struct lmf_str *lb_mfd[MAX_LF];	/* list of opened file pointers */
  int lb_size[MAX_LF];	/* estimated size of database */
  /* int lb_types[MAX_LF]; */	/* library types */
  int maxn;		/* longest library sequence chunk */
  int dupn;		/* overlap to use when segmenting sequence (p_comp) */
  int qoff;		/* overlap when segmenting long query sequence */
  int loff;		/* overlap when segmenting long library sequences */
  int maxt3;		/* overlap for tranlated sequences */
  int qdnaseq;		/* query is protein (0)/dna (1) */
  int ldnaseq;		/* library is protein (0)/dna (1) */
  int qframe;		/* number of possible query frames */
  int nframe;		/* frame for TFASTA */
  int nitt1;		/* nframe-1 */
  int thr_fact;		/* fudge factor for threads */
  int s_int;		/* sampling interval for statistics */
  int ql_start;		/* starting query sequence */
  int ql_stop;		/* ending query sequence */
  int nln;		/* number of library names */
  int pbuf_siz;		/* buffer size for sequences send in p2_complib */
  char qtitle[MAX_FN];	/* query title */
  char ltitle[MAX_FN];	/* library title */
  char flstr[MAX_FN];	/* FASTLIBS string */
  char outfile[MAX_FN];
  char label [MAXLN];	/* Output label */
  char f_id0[4];	/* function id for markx==10 */
  char f_id1[4];	/* function id for markx==10 */
  char libstr[MAX_FN];	/* Title from query sequence */
  char sqnam[4];	/* "aa" or "nt" */ 
  char sqtype[10];	/* "DNA" or "protein" */
  int long_info;	/* long description flag*/
  long sq0off, sq1off;	/* offset into aa0, aa1 */
  int markx;		/* alignment display type */
  int seqnm;		/* query sequence number */
  int nbr_seq;		/* number of library sequences */
  int term_code;	/* add termination codes to proteins if absent */
  int n1_high;		/* upper limit on sequence length */
  int n1_low;		/* lower limit on sequence length */
  double e_cut;		/* e_value for display */
  double e_low;		/* e_value for display */
  int e_cut_set;	/* e_value deliberately set */
  int pamd1;		/* 1st dimension of pam matrix */
  int pamd2;		/* 2nd dimension of pam matrix */
  int revcomp;		/* flag to do reverse complement */
  int quiet;		/* quiet option */
  int nrelv;		/* number of interesting scores */
  int srelv;		/* number of scores to show in showbest */
  int arelv;		/* number of scores to show at alignment */
  int z_bits;		/* z_bits==1: show bit score, ==0 show z-score */
  char alab[3][24];	/* labels for alignment scores */
  int nohist;		/* no histogram option */
  int nshow;
  int mshow;		/* number of scores to show */
  int mshow_flg;
  int ashow;		/* number of alignments to show */
  int nmlen;		/* length of name label */
  int show_code;	/* show alignment code in -m 9;  ==1 => identity only, ==2 alignment code*/
  int self;		/* self comparison */
  int thold;		/* threshold */
  int last_calc_flg;	/* needs a last calculation stage */
  int qshuffle;	/* shuffle the query and do additional comparisons */
  int shuff_max;	/* number of shuffles to perform */
  int shuff_node;	/* number of shuffles/worker node */
  int shuff_wid;
  int stages;		/* number of stages */
  double Lambda, K, H;	/* Karlin-Altschul parameters */
  int escore_flg;	/* use escore calculated by do_work() */
  struct hist_str hist;
  struct db_str db;
  void *pstat_void;
  struct a_struct aln;	/* has llen, llnctx, llnctx_flg, showall */
  char dfile [MAX_FN];	/* file for dumping scores to */
};


