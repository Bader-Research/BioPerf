#include "structs.h"
#include <signal.h>
#include "labeltree.h"
#include "inittree.h"
#include "specialinit.h"
#include "gen_tree.h"
#include "circ_order.h"
#include "const_tree.h"
#include "condense.h"
#include "correction.h"
#include "circular_ordering.h"
#include "inversion_median_alberto.h"
#ifndef CYGWIN_NT_40
#include <unistd.h>
#include <time.h>
#endif
#include "randomBSD.h"
#include <errno.h>
#include "uf.h"
#include "binencode.h"
#include "invdist.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "neighborj.h"
#include "growTree.h"
#ifdef CYGWIN_NT_40
#include "getopt3.h"
#else
#include <sys/time.h>
#endif
#include "read_input.h"

/* this is the helper function to compute a tree score,
   which will init the tree, score it and iterate it*/
int init_score_tree(int COND,struct tNode *tree,int NUM_GENES,int NUM_GENOMES,
                      int tspsolver,int thresh, struct adj_struct *adj_list,
                      struct adj_struct *adj_pool, intpair_t *neighbors,
                      int *stack,int *incycle,int *outcycle,int **weights,
                      int *degree,int *otherEnd,edge_t *edges,
                      struct tNode *tpool,struct genome_struct *labels,
                      int *pred1,int *pred2,int *picked,int *decode,
                      int *genes, int CIRCULAR, int distmethod, 
                      distmem_t distmem, int CORRECTION, struct qNode *qpool, 
                      int *edgepool, struct genome_struct *genome_list,
                      smalledge_t *smalledges, int ** status, triple_t* triple,
                      int inittspsolver, int OPT, int initmethod);
/* remember, you must set tree->leaf = TRUE before you call this one */

#ifdef CYGWIN_NT_40
__int64 numTree(int num_genomes);
#else
long long numTree(int num_genomes);
#endif


#ifdef GMP
#include "gmp.h"
#endif

#ifdef MPBPA
#include "mpi.h"
int MYPROC, PROCS;
#endif

FILE *outfile;

void print_usage(char *progname) {
  fprintf(stderr,
	  "\nUsage: %s -f <datafile> [-r <constraintfile>] [-t tspsolver] [-T initsolver] [-n NJsolver] [-K threshold] [-i initmethod] [-o outfile] [-s step] [-N] [-L] [-d] [-S] [-b upperbound] [-B] [-C]",progname);
#ifdef TESTING
  fprintf(stderr," [-Y num]");
#endif
  fprintf(stderr,"\n\n");
  fprintf(stderr,
	  " -r: constraintfile: only generate trees that are\n");
  fprintf(stderr,
	  "     compatible with the parenthesized tree constraint\n");
#ifdef CONCORDE
  fprintf(stderr,
	  " -t: tspsolver: [%d==Greedy, %d==Exact (default), %d==SimpleLK, %d==ChainedLK, %d==InversionMedian, %d==Alberto's Inversion Median]\n",
	  TSP_COALESCED,TSP_BBTSP,TSP_GREEDYLK,TSP_CHLINKERN, INVERSION_MEDIAN, INVERSION_MEDIAN_FAST);
#else
  fprintf(stderr,
	  " -t: tspsolver: [%d==Greedy, %d==Exact (default), %d==InversionMedian, %d==Alberto's Inversion Median]\n",
	  TSP_COALESCED,TSP_BBTSP,INVERSION_MEDIAN, INVERSION_MEDIAN_FAST);
#endif
#ifdef CONCORDE
  fprintf(stderr,
	  " -T: tspsolver for initialization: [%d==Greedy, %d==Exact (default), %d==SimpleLK, %d==ChainedLK]\n",
	  TSP_COALESCED,TSP_BBTSP,TSP_GREEDYLK,TSP_CHLINKERN);
#else
  fprintf(stderr,
	  " -T: tspsolver for initialization: [%d==Greedy, %d==Exact (default)]\n",
	  TSP_COALESCED,TSP_BBTSP);
#endif
#ifdef CONCORDE
  fprintf(stderr,
	  " -n: tspsolver for NJ trees: [%d==Greedy, %d==Exact (default), %d==SimpleLK, %d==ChainedLK]\n",
	  TSP_COALESCED,TSP_BBTSP,TSP_GREEDYLK,TSP_CHLINKERN);
#else
  fprintf(stderr,
	  " -n: tspsolver for NJ trees: [%d==Greedy, %d==Exact (default)]\n",
	  TSP_COALESCED,TSP_BBTSP);
#endif
#ifdef CONCORDE
  fprintf(stderr,
	  " -K: threshold on size to use exact rather than LK solver (default 20)\n");
#endif
  fprintf(stderr,
	  " -i: initmethod: [%d==Random, %d==Nearest-Neighbors (small), %d==Nearest-Neighbors (S/B, default), %d==Nearest-Neighbors (large), %d==Propagate-Nearest, %d==Propagate-Median, %d==Uniform, %d==Adjacency-Parsimony]\n",
	  RAND,SMALLNN,SBNN,BIGNN,FASTPROP,MEDIANPROP,TRIV,ADJPARS);
  fprintf(stderr,
	  " -o: outfile: output filename (default to stdout)\n");
  fprintf(stderr,
	  " -s: step: next tree is <step> away from current (default: 1)\n");
  fprintf(stderr,
	  " -N: do not iterate labellings beyond initialization (default: off)\n");
  fprintf(stderr,
	  " -L: process linear (noncircular) genomes (default: circular)\n");
  fprintf(stderr,
	  " -d: Use Breakpoint distance (default: Inversion distance)\n");
  fprintf(stderr,
	  " -c: turn off tree counting (default: ON)\n");
  fprintf(stderr,
	  " -C: condense triples before calling TSP solver (default: OFF)\n");
  fprintf(stderr,
	  " -b: use the following integer as an upper bound in the search\n");
  fprintf(stderr,
	  " -B: do not use circular ordering lower bound; implied by -S; default: use the bound)\n");
  fprintf(stderr,
	  " -S: skip dist. matrix and NJ computations (for timings); also sets -B; default: compute matrix and NJ\n");
  fprintf(stderr,
	  " -m: tighten bound (slower)\n");
  fprintf(stderr,
	  " -e: EDE correction for the inversion distance\n");
  fprintf(stderr,
	  " -l: Layered solution\n");
  fprintf(stderr,
	  " -a: one level branch and bound method\n");

#ifdef TESTING
  fprintf(stderr,
	  " -Y: Developer's testing routine\n");
#endif
  return;
}

int contains(const char *s1, const char *s2) {
  int size1, size2, sdiff;
  if (!strcmp(s1, s2))
    return TRUE;

  size1 = strlen(s1);
  size2 = strlen(s2);
  sdiff = size1-size2;

  if (sdiff <= 0)
    return FALSE;

  if (!strcmp(s1+sdiff, s2))
    return TRUE;

  return FALSE;
}

void print_distmatrix(int **distmatrix,int num_genomes) {
  int i, j;

  fprintf(outfile,"\n");

  fprintf(outfile,"     ");
  for (j=0; j<num_genomes ; j++)
    fprintf(outfile,"%4d ",j+1); /* the +1 because arrays are 0-based */
  fprintf(outfile,"\n");

  fprintf(outfile,"\n");

  for (i=0 ; i<num_genomes ; i++) {
    fprintf(outfile,"%3d: ",i+1); /* the +1 because arrays are 0-based */
    for (j=0; j<=i ; j++)
      fprintf(outfile,"   . ");
    for (j=i+1 ; j<num_genomes ; j++) {
      fprintf(outfile,"%4d ",distmatrix[i][j]);
    }
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
    
  return;
}

/* Stacia added this so that the ouput can be directly fed into other
 * programs such as the tds suite. It expects # of taxa on the first line
 * by itself, and then a full symmetric matrix with the taxa name as the
 * first thing on a line.
 */
void print_full_distmatrix(int **distmatrix,int num_genomes, 
			   struct genome_struct *genome_list) {
  int i, j;

  fprintf(outfile,"\n   %d\n",num_genomes);

  for (i=0 ; i<num_genomes ; i++) {
    fprintf(outfile,"%s\t",genome_list[i].gnamePtr);
    for (j=0; j<i ; j++)
      fprintf(outfile,"%2d ",distmatrix[j][i]);
    fprintf(outfile," 0 ");
    for (j=i+1 ; j<num_genomes ; j++) {
      fprintf(outfile,"%2d ",distmatrix[i][j]);
    }
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  
  return;
}


/***********************************************************************/
int main(int argc, char *argv[])
{ 
  FILE *input; 
  struct genome_struct *genome_list;
  struct genome_struct *labels=(struct genome_struct*)NULL;
  /* genome_list stores leaf genomes;
     labels is a pool of genomes entries used to label internal nodes,
     plus 4 entries used for temporary condensing (in condense3) */
  struct tNode *tree;
  struct tNode *tpool=(struct tNode*)NULL; /* tpool is an array of tree nodes */
  struct qNode *qpool=(struct qNode*)NULL;
  struct adj_struct *adj_list=(struct adj_struct*)NULL;
  struct adj_struct *adj_pool=(struct adj_struct*)NULL;
  intpair_t *neighbors=(intpair_t*)NULL;
  triple_t *triple=(triple_t*)NULL;
  int i, j, k, step, currentLow, score, prev, newbest, upperLowerBound, lowerLowerBound;
  register int colb=0, firstg;
  int cutValue;
  int *condense_succ, *condense_decode;
  int orig_num_genes;
  int inittspsolver,tspsolver,njsolver,thresh,initmethod;
#ifdef GMP
  mpz_t stepping;
  mpz_t treenbr;
#else
  int stepping;
#endif
#ifdef MPBPA
#ifdef GMP
  mpz_t PROCS_stepping;
#else
  int PROCS_stepping;
#endif
#endif
  int *genes=(int*)NULL;
  int *incycle=(int*)NULL;
  int *outcycle=(int*)NULL;
  int *stack=(int*)NULL;
  int **distmatrix, **cpdist;
  int **weights=(int**)NULL;
  int **status=(int**)NULL;
  int *degree=(int*)NULL;
  int *picked=(int*)NULL;
  int *decode=(int*)NULL;
  int *pred1=(int*)NULL;
  int *pred2=(int*)NULL;
  int *edgepool=(int*)NULL;
  int *otherEnd=(int*)NULL;
  edge_t *edges=(edge_t*)NULL;
  smalledge_t *smalledges=(smalledge_t*)NULL;
  struct stat statbuf;
#ifdef TESTING
  struct timeval tm;
  double time0,time1;
  int not_done_count = 0;
#endif
  int distmethod;
  distmem_t distmem;
  int INVDIST, DISTMAT;
  int c, errflg;
  extern char *optarg;
  extern int optopt;
  char *inputfname;
  char *outputfname;
  char  *tmpoutname;
  char *constfname;
  int constrained;
  FILE *constFile;
  char constString[MAX_STR_LEN];
  int BOUND, NJ,  COND, CIRCULAR, OPT, YTEST, YTESTARG, CORRECTION,TIGHTBOUND, LAYERED;
  int GROWTREE, BRANCH;
  int NUM_GENES;   /* number of genes in input */
  int NUM_GENOMES; /* number of gene orders in input */

/* added for the branch and bound*/
  int flag = 0;
#ifdef MPBPA
#ifdef GMP
  mpz_t PROCS_treelevelcount;
  mpz_t PROCS_treelevelcount2;
  mpz_t PROCS_treelevel2count;
  mpz_t PROCS_treelevel2count2;
  
#else
  int PROCS_treelevelcount;
  int PROCS_treelevelcount2;
  int PROCS_treelevel2count;
  int PROCS_treelevel2count2;
  int PROCS_treelevelstepping;
#endif
#endif
  int treelevelcount;
  int treelevelcount2;
  int treelevel2count;
  int treelevel2count2;
  int minidist;
#ifdef CYGWIN_NT_40
__int64 treenbr, count, totalN, foundTree;
__int64** lobs;
__int64* lobn;
__int64* workLob;
__int64* tmpInt;
#else
long long treenbr, count, totalN, foundTree;
long long** lobs;
long long*  lobn;
long long* workLob;
long long*  tmpInt;
#endif
  
  int best_so_far = 0;
  int best_so_far_2; /* best_so_far_2 is the double of best_so_fa, to avoid the divide in
                        each lower bound calculation */

  bin_env_t  *bin_env=NULL;
  env_t    const_env;
  ConstraintTree_T const_tree;
  time_t start;
  time_t end;
  struct genome_struct* cpgenome=NULL;

  /* for layered method*/
  int*   lobt;
  FILE*  cacheFile;
  char   cacheName[128];

  struct tm *newtime;

#ifdef MPBPA
  double start_time,end_time;
#endif


#ifdef MPBPA
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MYPROC);
  MPI_Comm_size(MPI_COMM_WORLD, &PROCS);
  start_time=MPI_Wtime();
#endif
  outputfname = (char*) malloc(256*sizeof(char));
  INVDIST =  contains(argv[0], INVDISTEXE);
  DISTMAT =  contains(argv[0], DISTMATEXE);
  inputfname  = (char *)NULL;
  input       = (FILE *)NULL;
  tmpoutname = (char*)NULL;
  /*outputfname = (char *)NULL;*/
  constrained = FALSE;
  constfname  = (char *)NULL;
  outfile     = stdout;
  inittspsolver = TSP_BBTSP;
  tspsolver   = TSP_BBTSP;
  njsolver    = TSP_BBTSP;
  thresh      = 20;
  initmethod  = SBNN;
  cutValue = LARGENUM; /* v */

#ifdef GMP
  mpz_init_set_si(stepping, 1);
  mpz_init_set_si(treenbr,10);
#else
  stepping    = 1;
  treenbr     = 0;
#endif

#ifdef MPBPA
#ifdef GMP
  mpz_init(PROCS_stepping);
  mpz_mul_ui(PROCS_stepping, stepping, (unsigned long int)PROCS);
#else
  PROCS_stepping = PROCS * stepping;
#endif
#endif
  
  distmethod  = DIST_INV;
  COND        = FALSE;
  OPT         = TRUE; 
  CIRCULAR    = TRUE;
  NJ          = TRUE;
  BOUND	      = TRUE;
  YTEST       = FALSE;
  CORRECTION  = FALSE;
  TIGHTBOUND  = FALSE;
  YTESTARG    = 0;
  LAYERED     = FALSE;
  GROWTREE    = FALSE; /* g */
  BRANCH     = FALSE; /*u*/
  errflg = 0;
  while ((c = getopt(argc,argv,"v:i:f:r:t:T:n:K:o:s:b:aNLlCemdgPBSY:")) != -1) {
    switch (c) {
    case 'v':
      cutValue = atoi(optarg);
      break;
    case 'm':
      TIGHTBOUND = TRUE;
      break;
    case 'g':
      GROWTREE = TRUE;
      break;
    case 'e':
      CORRECTION = TRUE;
      break;
    case 'b':
      best_so_far = atoi(optarg);
      break;
    case 'i':
      initmethod = atoi(optarg);
      flag = 1;
      break;
    case 'f':
      inputfname = optarg;
      break;
    case 'r':
      constrained = TRUE;
      constfname  = optarg;
      break;
    case 't':
      tspsolver = atoi(optarg);
      break;
    case 'T':
      inittspsolver = atoi(optarg);
      break;
    case 'n':
      njsolver = atoi(optarg);
      break;
    case 'K':
      thresh = atoi(optarg);
      break;
    case 'o':
      tmpoutname = optarg;
#ifdef MPBPA
      sprintf(outputfname, "%s%d", optarg, MYPROC);
#else
      sprintf(outputfname, "%s", optarg);
#endif
      /*outputfname = optarg;*/
      break;
    case 'a':
      BRANCH = TRUE;
#ifdef MPBPA
#ifdef GMP
        mpz_set_si(PROCS_stepping, 1);
#else
      PROCS_stepping = 1;
#endif
#endif
      break;
    case 's':
#ifdef GMP
      if (!mpz_set_str(stepping, optarg, 10))
	fprintf(stderr,"ERROR: mpz_set_str of stepping failed\n");
#else
      stepping = atoi(optarg);
#endif
#ifdef MPBPA
#ifdef GMP
      mpz_mul_ui(PROCS_stepping, stepping, (unsigned long int)PROCS);
#else
      PROCS_stepping = PROCS * stepping;
#endif
#endif
      break;
    case 'L':
      CIRCULAR = FALSE;
      break;
    case 'l':
      LAYERED = TRUE;
      break;
    case 'N':
      OPT = FALSE;
      break;
    case 'C':
      COND = TRUE;
      break;
    case 'B':
      BOUND = FALSE;
      break;
    case 'S':
      NJ = FALSE;
      BOUND = FALSE;
      break;
#ifdef TESTING
    case 'Y':
      YTEST = TRUE;
      YTESTARG = atoi(optarg);
      break;
#endif
    case 'd':
      distmethod = DIST_BP;
      break;
    case ':':
      fprintf(stderr, "Option -%c requires an operand\n", optopt);
      errflg++;
      break;
    case '?' :
      fprintf(stderr, "Unrecognized option: -%c\n", optopt);
      errflg++;
    }
  }
  
  /* check we have usable parameters */
  if (errflg || (tspsolver < 1) || (tspsolver > SOLVERCHOICES) ||
      (inittspsolver < 1) || (inittspsolver > SOLVERCHOICES) ||
      (njsolver < 1) || (njsolver > SOLVERCHOICES) ||
      (initmethod < 1) || (initmethod > INITCHOICES) ||
#ifdef GMP
      (mpz_sgn(stepping) != 1)
#else
      (stepping <= 0)
#endif
      ) {
    print_usage(argv[0]);
    exit(-1);
  }

  /* set the default of init method for inversion median to 2 (NN small) */
  if (flag == 0 && (inittspsolver == 3 || inittspsolver == 4)) {
     initmethod = 2;
  }
  if (inputfname == (char *)NULL) {
    fprintf(stderr,"ERROR: input filename required\n");
    print_usage(argv[0]);
    exit(-1);
  }
  else {
    if (stat(inputfname, &statbuf)) {
      fprintf(stderr,"ERROR: Input file (%s): ",inputfname);
      perror("");
      print_usage(argv[0]);
      exit(-1);
    }
    if (statbuf.st_mode & S_IFDIR) {
      fprintf(stderr,"ERROR: Input file (%s) is a directory.\n",inputfname);
      print_usage(argv[0]);
      exit(-1);
    }
    input = fopen(inputfname,"r");
    if (input == (FILE*)NULL) {
      fprintf(stderr,"ERROR: Could not open input file (%s): ",inputfname);
      perror("");
      print_usage(argv[0]);
      exit(-1);
    }
  }

  if (/*outputfname*/ tmpoutname != (char *)NULL) {
    if (!stat(outputfname, &statbuf)) {
      fprintf(stderr,"ERROR: Output file (%s) already exists.\n",outputfname);
      print_usage(argv[0]);
      exit(-1);
    }
    outfile = fopen(outputfname,"a+");
    if (outfile == (FILE *)NULL) {
      fprintf(stderr,"ERROR: Could not open output file (%s): ",outputfname);
      perror("");
      print_usage(argv[0]);
      exit(-1);
    }
  }

  if (constrained) {
    if (stat(constfname, &statbuf)) {
      fprintf(stderr,"ERROR: Constraint file (%s): ",constfname);
      perror("");
      print_usage(argv[0]);
      exit(-1);
    }
    if (statbuf.st_mode & S_IFDIR) {
      fprintf(stderr,"ERROR: Constraint file (%s) is a directory.\n",
	      constfname);
      print_usage(argv[0]);
      exit(-1);
    }
    constFile = fopen(constfname,"r");
    if (constFile == (FILE*)NULL) {
      fprintf(stderr,"ERROR: Could not open constraint file (%s): ",
	      constfname);
      perror("");
      print_usage(argv[0]);
      exit(-1);
    }
    fscanf(constFile,"%[0123456789,()]",constString);
    if (strlen(constString)+1 > MAX_STR_LEN) {
      fprintf(stderr,"ERROR: Memory overrun. Increase MAX_STR_LEN (%d)\n",
	      MAX_STR_LEN);
      exit(-1);
    }
    fclose(constFile);
  }

  if (!(DISTMAT || INVDIST)) { /* print full info only if solving */
#ifdef MPBPA
    if (MYPROC == 0 || MYPROC != 0) {
#endif
      fprintf(outfile,"Running: \t\t%s\t",argv[0]);
#ifdef CONCORDE
      fprintf(outfile,"(with Concorde)\n");
#else
      fprintf(outfile,"(without Concorde)\n");
#endif
      fprintf(outfile,"Input: \t\t\t%s\n",inputfname);
      if (constrained) {
        fprintf(outfile,"Constraint File: \t%s\n",constfname);
        fprintf(outfile,"Constraint Tree: \t%s\n",constString);
      }
      if (CIRCULAR)
        fprintf(outfile,"Genome: \t\tCircular\n");
      else
        fprintf(outfile,"Genome: \t\tNoncircular\n"); 
      if (OPT) {
        fprintf(outfile,"TSP Solver: \t\t");
        switch (tspsolver) {
#ifdef CONCORDE
        case TSP_CHLINKERN : fprintf(outfile,"Chained Lin-Kernighan\n"); break;
#endif
        case TSP_BBTSP : fprintf(outfile,"Exact\n"); break;
#ifdef CONCORDE
        case TSP_GREEDYLK : fprintf(outfile,"Simple Lin-Kernighan\n"); break;
#endif
        case TSP_COALESCED : fprintf(outfile,"Coalesced simple paths\n"); break;
        case INVERSION_MEDIAN : fprintf(outfile,"Inversion median\n"); break;
        case INVERSION_MEDIAN_FAST : fprintf(outfile,"Inversion median fast\n"); break; 
        default: fprintf(stderr,"ERROR: Bad TSP Solver\n");
        }
      }
      else
        fprintf(outfile,"No iteration on labels beyond initialization\n");
      fprintf(outfile,"Initialization Method: \t");
      switch (initmethod) {
      case RAND : fprintf(outfile,"Random\n"); break;
      case SMALLNN   : fprintf(outfile,"Nearest-Neighbors (small)\n");   break;
      case SBNN   : fprintf(outfile,"Nearest-Neighbors (Sankoff/Blanchette)\n");   break;
      case BIGNN   : fprintf(outfile,"Nearest-Neighbors (big)\n");   break;
      case FASTPROP : fprintf(outfile,"Propagate-Nearest\n"); break;
      case MEDIANPROP : fprintf(outfile,"Propagate-Median\n"); break;
      case TRIV : fprintf(outfile,"Uniform (all identity)\n"); break;
      case ADJPARS: fprintf(outfile,"Adjacency-Parsimony\n"); break;
      default: fprintf(stderr,"ERROR: Bad Method\n");
      }
      fprintf(outfile,"Initial TSP Solver: \t");
      switch (inittspsolver) {
        case TSP_COALESCED : fprintf(outfile,"Coalesced simple paths\n"); break;
        case TSP_BBTSP : fprintf(outfile,"Exact\n"); break;
        case INVERSION_MEDIAN : fprintf(outfile,"Inversion median\n"); break;
        case INVERSION_MEDIAN_FAST : fprintf(outfile,"Inversion median fast\n"); break;
#ifdef CONCORDE
        case TSP_GREEDYLK : fprintf(outfile,"Simple Lin-Kernighan\n"); break;
        case TSP_CHLINKERN : fprintf(outfile,"Chained Lin-Kernighan\n"); break;
#endif
      default: fprintf(stderr,"ERROR: Bad Init TSP Solver\n");
      }
      fprintf(outfile,"NJ TSP Solver: \t\t");
      switch (njsolver) {
        case TSP_COALESCED : fprintf(outfile,"Coalesced simple paths\n"); break;
        case TSP_BBTSP : fprintf(outfile,"Exact\n"); break;
#ifdef CONCORDE
        case TSP_GREEDYLK : fprintf(outfile,"Simple Lin-Kernighan\n"); break;
        case TSP_CHLINKERN : fprintf(outfile,"Chained Lin-Kernighan\n"); break;
#endif
      default: fprintf(stderr,"ERROR: Bad Init TSP Solver\n");
      }
      fprintf(outfile,"Genome Distance: \t");
      switch (distmethod) {
      case DIST_BP  : fprintf(outfile,"Breakpoint\n"); break;
      case DIST_INV : fprintf(outfile,"Inversion\n");   break;
      default: fprintf(stderr,"ERROR: Bad Distance Method\n");
      }
      if (best_so_far > 0) {
	fprintf(outfile,"Lower Bound Provided: %d\n",best_so_far);
      }
#ifdef GMP
      fprintf(outfile,"Stepping: \t\t");
      mpz_out_str(outfile,10,stepping);
      fprintf(outfile,"\n");
#else
      fprintf(outfile,"Stepping: \t\t%d\n",stepping);
#endif
    
#ifdef MPBPA
      fprintf(outfile,"Processors: \t\t%d\n",PROCS);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  else { /* either INVDIST or DISTMAT -- simplify output */
#ifdef MPBPA
    if (MYPROC == 0) {
#endif
      fprintf(outfile,"Running: \t\t%s\t\n",argv[0]);
      fprintf(outfile,"Input: \t\t\t%s\n",inputfname);
      if (CIRCULAR)
        fprintf(outfile,"Genome: \t\tCircular\n");
      else
        fprintf(outfile,"Genome: \t\tNoncircular\n"); 
#ifdef MPBPA
      fprintf(outfile,"Processors: \t\t%d\n",PROCS);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  srandom(time(0));

  condense_succ = (int *)malloc((2*MAX_NUM_GENES+1)*sizeof(int));
  if (condense_succ == (int *)NULL)
    fprintf(stderr,"ERROR: condense_succ NULL\n");
  condense_decode = (int *)malloc((2*MAX_NUM_GENES+1)*sizeof(int));
  if (condense_decode == (int *)NULL)
    fprintf(stderr,"ERROR: condense_decode NULL\n");
  
  read_data(input,&genome_list,&NUM_GENES,&NUM_GENOMES, INVDIST, CIRCULAR,
	    condense_succ, condense_decode, &orig_num_genes);
  
  if (NUM_GENES < 2) {
    fprintf(stderr,"ERROR: Fewer than 2 gene fragments\n");
    exit(-1);
  }
  if (NUM_GENOMES < 3) {
    fprintf(stderr,"ERROR: Fewer than 3 genomes; use invdist\n");
    exit(-1);
  }

  /* ALLOCATE MEMORY */
  /* TABLE OF array ALLOCATIONS
    
        NAME             SIZE                     USED IN         INITIALIZED?
     -------------------------------------------------------------------------
     genes              genes           iterate_over_tree
     incycle 		2*genes+1	LK solvers, PROP
     outcycle		2*genes+1	all solvers, PROP
     degree		2*genes+1	bbtsp solver, PROP
     otherEnd		2*genes+1	bbtsp solver,
     stack		2*genes+1	bbtsp solver, PROP
     pred1		2*genes+1	condense3, MEDIANPROP
     pred2		2*genes+1	condense3, MEDIANPROP
     picked		2*genes+1	condense3 PROP
     decode		2*genes+1	condense3 only
     neighbors		2*genes+1	bbtsp solver,
     adj_pool		14*genes	bbtsp solver only
     adj_list		2*genes+1	bbtsp solver only
     edges		7*genes		bbtsp solver,
     smalledges		7*genes		adjpars and avtsp init + solvers
     weights		(2*genes)^2	LK solvers & AP init		NO
     status		(2*genes)^2	AP init only     		NO
     tpool		2*GENOMES+6	everywhere
     labels		2*GENOMES+6	everywhere
     edgepool		2*GENOMES+1	all solvers
     qpool		2*GENOMES+1	FASTPROP only			YES
     triple		2*GENOMES	bigNN only
     distmatrix		GENOMES^2	circular bounding, print matrix
     distmem.hammingArr genes+1		everywhere
     distmem.perm1	2*genes+2	invdist only
     distmem.perm2	2*genes+2	invdist only
     distmem.perm	2*genes+2	invdist only
     distmem.done	2*genes+2	invdist only
     distmem.greyEdges	2*genes+2	invdist only
     distmem.stack	2*genes+2	invdist only
     distmem.oriented	2*genes+2	invdist only
     distmem.cc		2*genes+2	invdist only
     distmem.labeled	2*genes+2	invdist only
     distmem.components	2*genes+2	invdist only
     distmem.uf		2*genes+2	invdist only
  */

if (!(INVDIST || DISTMAT)) {
  outcycle = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
  if (outcycle == (int*)NULL)
    fprintf(stderr,"ERROR: outcycle NULL\n");

  genes = (int *)malloc(NUM_GENES*sizeof(int));
  if (genes == (int*)NULL)
    fprintf(stderr,"ERROR: genes NULL\n");
  
  incycle = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
  if (incycle == (int*)NULL)
    fprintf(stderr,"ERROR: incycle NULL\n");
	
  if (
#ifdef CONCORDE
      (tspsolver==TSP_CHLINKERN) || (tspsolver==TSP_GREEDYLK) ||
      (njsolver==TSP_CHLINKERN) || (njsolver==TSP_GREEDYLK) ||
      (inittspsolver==TSP_CHLINKERN) || (inittspsolver==TSP_GREEDYLK) ||
#endif
      initmethod==ADJPARS)
    {
      weights = (int **)malloc((2*NUM_GENES)*sizeof(int *));
      if (weights == (int**)NULL)
	fprintf(stderr,"ERROR: weights NULL\n");
      for (i=0; i < 2*NUM_GENES; i++) {
	weights[i] = (int *)malloc((2*NUM_GENES)*sizeof(int));
	if (weights[i] == (int*)NULL)
	  fprintf(stderr,"ERROR: weights[i] NULL\n");
      }
    }
  else {
    weights = NULL;
  }
  distmatrix = NULL;

  if (initmethod==ADJPARS) {
    status = (int **)malloc((2*NUM_GENES)*sizeof(int *));
    if (status == (int**)NULL)
      fprintf(stderr,"ERROR: status NULL\n");
    for (i=0; i < 2*NUM_GENES; i++) {
      status[i] = (int *)malloc((2*NUM_GENES)*sizeof(int));
      if (status[i] == (int*)NULL)
	fprintf(stderr,"ERROR: status[i] NULL\n");
    }
  }
  else {
    status = NULL;
  }
  
  degree = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
  if (degree == (int*)NULL)
    fprintf(stderr,"ERROR: degree NULL\n");

  otherEnd = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
  if (otherEnd == (int *)NULL)
    fprintf(stderr,"ERROR: otherEnd NULL\n");

  stack = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
  if (stack == (int*)NULL)
    fprintf(stderr,"ERROR: stack NULL\n");

  pred1 = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
  if (pred1 == (int*)NULL)
    fprintf(stderr,"ERROR: pred1 NULL\n");

  pred2 = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
  if (pred2 == (int*)NULL)
    fprintf(stderr,"ERROR: pred2 NULL\n");

  edges = (edge_t *)malloc((7*NUM_GENES)*sizeof(edge_t));
  if (edges == (edge_t *)NULL)
    fprintf(stderr,"ERROR: edges NULL\n");

  smalledges = (smalledge_t *)NULL;
  if (initmethod==ADJPARS) {
    smalledges =
      (smalledge_t *)malloc((NUM_GENES*(2*NUM_GENES-1))*sizeof(smalledge_t));
    if (smalledges == (smalledge_t *)NULL)
      fprintf(stderr,"ERROR: edges NULL\n");
  }

  /* tree pool, edgepool, labels, triple, and qpool are only used for tree
   * -- so size is determined by NUM_GENOMES */
  /* the extra 4 locations are used in condense3.c */
  tpool = (struct tNode *)malloc((2*NUM_GENOMES+6)*sizeof(struct tNode));
  if (tpool == (struct tNode *)NULL)
    fprintf(stderr,"ERROR: tpool NULL\n");
  
  labels = (struct genome_struct *)
    malloc((2*NUM_GENOMES+6)*sizeof(struct genome_struct));
  if (labels == (struct genome_struct *)NULL)
    fprintf(stderr,"ERROR: labels NULL\n");

  for (i=0 ; i< (2*NUM_GENOMES+6) ; i++) {
    labels[i].genes = (int *)malloc(NUM_GENES * sizeof(int));
    if (labels[i].genes == NULL) {
      fprintf(stderr,"ERROR: labels genes  NULL\n");
    }
  }
  

  edgepool = (int *)malloc((2*NUM_GENOMES+1)*sizeof(int));
  if (edgepool == (int *)NULL)
    fprintf(stderr,"ERROR: edgepool NULL\n");

  if ((initmethod==FASTPROP) || (initmethod==MEDIANPROP) ||
      (initmethod==SMALLNN) || (initmethod==SBNN) || (initmethod==BIGNN)) {
     qpool = (struct qNode *)malloc((2*NUM_GENOMES+1)*sizeof(struct qNode));
     if (qpool == (struct qNode *)NULL)
       fprintf(stderr,"ERROR: qpool NULL\n");
     /* now thread qpool circularly */
     for (i=0; i<(2*NUM_GENOMES); i++)
       qpool[i].next = &qpool[i+1];
     qpool[2*NUM_GENOMES].next = &qpool[0];
  }
  else {
     qpool = NULL;
  }
  
  if (initmethod==BIGNN) {
     triple = (triple_t *)malloc((2*NUM_GENOMES)*sizeof(triple_t));
     if (triple == (triple_t *)NULL)
       fprintf(stderr,"ERROR: triple NULL\n");
  }
  else {
     triple = NULL;
  }
  
  adj_pool =
    (struct adj_struct *)malloc((14*NUM_GENES)*sizeof(struct adj_struct));
  if (adj_pool == (struct adj_struct *)NULL)
    fprintf(stderr,"ERROR: adj_pool NULL\n");
  
  adj_list =
    (struct adj_struct *)malloc((2*NUM_GENES+1)*sizeof(struct adj_struct));
  if (adj_list == (struct adj_struct *)NULL)
    fprintf(stderr,"ERROR: adj_list NULL\n");

  picked = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
  if (picked == (int*)NULL)
    fprintf(stderr,"ERROR: picked NULL\n");

  if (COND) {
    decode = (int *)malloc((2*NUM_GENES+1)*sizeof(int));
    if (decode == (int*)NULL)
      fprintf(stderr,"ERROR: decode NULL\n");
    decode += NUM_GENES;
  }
  else {
     decode = NULL;
  }

  neighbors = (intpair_t *)malloc((2*NUM_GENES+1)*sizeof(intpair_t));
  if (neighbors == (intpair_t *)NULL)
    fprintf(stderr,"ERROR: neighbors NULL\n");
}  /* end allocation for solvers only */

  distmatrix = (int **)malloc((NUM_GENOMES)*sizeof(int *));
  if (distmatrix == (int**)NULL)
    fprintf(stderr,"ERROR: distmatrix NULL\n");
  for (i=0; i < NUM_GENOMES; i++) {
    distmatrix[i] = (int *)malloc((NUM_GENOMES)*sizeof(int));
    if (distmatrix[i] == (int*)NULL)
      fprintf(stderr,"ERROR: distmatrix[i] NULL\n");
  }
  cpdist = (int **)malloc((NUM_GENOMES)*sizeof(int *));
  if (cpdist == (int**)NULL)
    fprintf(stderr,"ERROR: cpdist NULL\n");
  for (i=0; i < NUM_GENOMES; i++) {
    cpdist[i] = (int *)malloc((NUM_GENOMES)*sizeof(int));
    if (cpdist[i] == (int*)NULL)
      fprintf(stderr,"ERROR: cpdist[i] NULL\n");
  }

  distmem.hammingArr = (int *)malloc((NUM_GENES+1)*2*sizeof(int));
  if (distmem.hammingArr == (int *)NULL)
    fprintf(stderr,"ERROR: hammingArr NULL\n");
  distmem.perm1 = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.perm1 == (int *)NULL)
    fprintf(stderr,"ERROR: perm1 NULL\n");
  distmem.perm2 = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.perm2 == (int *)NULL)
    fprintf(stderr,"ERROR: perm2 NULL\n");
  distmem.perm = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.perm == (int *)NULL)
    fprintf(stderr,"ERROR: perm NULL\n");
  distmem.done = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.done == (int *)NULL)
    fprintf(stderr,"ERROR: done NULL\n");
  distmem.greyEdges = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.greyEdges == (int *)NULL)
    fprintf(stderr,"ERROR: greyEdges NULL\n");
  distmem.stack = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.stack == (int *)NULL)
    fprintf(stderr,"ERROR: stack NULL\n");
  distmem.oriented = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.oriented == (int *)NULL)
    fprintf(stderr,"ERROR: oriented NULL\n");
  distmem.cc = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.cc == (int *)NULL)
    fprintf(stderr,"ERROR: cc NULL\n");
  distmem.labeled = (int *)malloc((2*NUM_GENES+2)*sizeof(int));
  if (distmem.labeled == (int *)NULL)
    fprintf(stderr,"ERROR: labeled NULL\n");
  distmem.components = (component_t *)
    malloc((2*NUM_GENES+2)*sizeof(component_t));
  if (distmem.components == (component_t *)NULL)
    fprintf(stderr,"ERROR: components NULL\n");
  distmem.uf = UFalloc(2*NUM_GENES+2);
  /* end of memory allocation */

  /*allocate memory for tree generator*/
  if (!constrained) {
    bin_env=alloc4_bin_tree_env(NUM_GENOMES);
    /*const_tree is used by NJ*/
    alloc4_const_tree_env(NUM_GENOMES,&const_env); 
  }
  else {
    init_const(&const_env,constString,const_tree,tpool);
    alloc4_const_tree_env(const_env.num_leaves,&const_env);
  }

  if (DISTMAT) {
#ifdef MPBPA
    if (MYPROC == 0) {
#endif
      fprintf(outfile,"\n");
      setBPmatrix(distmatrix,genome_list,NUM_GENES,NUM_GENOMES,&distmem,CIRCULAR);
      fprintf(outfile,"Breakpoint Distance Matrix:\n");
      print_full_distmatrix(distmatrix,NUM_GENOMES,genome_list);
      setinvmatrix(distmatrix,genome_list,NUM_GENES,NUM_GENOMES,&distmem,CIRCULAR);
      fprintf(outfile,"Inversion Distance Matrix:\n");
      print_full_distmatrix(distmatrix,NUM_GENOMES,genome_list);
#ifdef MPBPA
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    fflush(outfile);
    fclose(outfile);

    UFfree(distmem.uf);
    free(distmem.components);
    free(distmem.labeled);
    free(distmem.cc);
    free(distmem.oriented);
    free(distmem.stack);
    free(distmem.greyEdges);
    free(distmem.done);
    free(distmem.perm1);
    free(distmem.perm2);
    free(distmem.perm);
    free(distmem.hammingArr);
    for (i=0; i < NUM_GENOMES; i++)
      free(cpdist[i]);
    free(cpdist);
    for (i=0; i < NUM_GENOMES; i++)
      free(distmatrix[i]);
    free(distmatrix);
    free(condense_decode);
    free(condense_succ);
    for (i=0 ; i<NUM_GENOMES ; i++) {
      free(genome_list[i].genes);
      free(genome_list[i].gnamePtr);
    }
    free(genome_list);
#ifdef MPBPA
    MPI_Finalize();
#endif
    exit(1);
  }

  if (NJ) {
#ifdef MPBPA
    if (MYPROC == 0) {
#endif
      neighborj_SN(genome_list, NUM_GENES, NUM_GENOMES, &distmem, CIRCULAR, 1);
      neighborj_SK(genome_list, NUM_GENES, NUM_GENOMES, &distmem, CIRCULAR, 1);
#ifdef MPBPA
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    newbest = neighborj_score(genome_list, NUM_GENES, NUM_GENOMES,
				  &distmem,CIRCULAR,
				  tpool,edgepool,labels,initmethod,COND,qpool,
				  adj_list, adj_pool, stack,degree,otherEnd,
				  neighbors,smalledges,edges,outcycle,
				  pred1,pred2,picked,decode,njsolver,
				  triple,incycle,njsolver,distmethod,thresh,
				  weights,status,&const_env,genes,
				  condense_succ, condense_decode,
				  orig_num_genes, CORRECTION);
#ifdef MPBPA
    MPI_Allreduce(&newbest, &i, 1, MPI_INT, MPI_MIN,
		  MPI_COMM_WORLD);
    newbest = i;
    if (MYPROC == 0) {
      fprintf(outfile,"Best Neighbor-joining score: %d\n",newbest);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif    
  }
  else {
    newbest = 3*NUM_GENES*NUM_GENOMES;
  }
  
  if (YTEST) {
#ifdef TESTING
    time_linear=0.0;
    time_BH=0.0;
    gettimeofday(&tm,NULL);
    time0=(double)tm.tv_sec+(double)tm.tv_usec/1000000.0;
    for (i=0 ; i<YTESTARG ; i++)
      calc_invmatrix(genome_list,NUM_GENES,NUM_GENOMES,&distmem,CIRCULAR);
    gettimeofday(&tm,NULL);
    time1=(double)tm.tv_sec+(double)tm.tv_usec/1000000.0;
    fprintf(outfile,"running time of calc_invmatrix:%9.6f\n",time1-time0);
    fprintf(outfile,"running time of connected_component:%9.6f\n",
	    time_linear);
    fflush(outfile);

    gettimeofday(&tm,NULL);
    time0=(double)tm.tv_sec+(double)tm.tv_usec/1000000.0;
    for (i=0 ; i<YTESTARG ; i++)
      calc_invmatrix_BH(genome_list,NUM_GENES,NUM_GENOMES,&distmem,CIRCULAR);
    gettimeofday(&tm,NULL);
    time1=(double)tm.tv_sec+(double)tm.tv_usec/1000000.0;
    fprintf(outfile,"running time of calc_invmatrix_BH:%9.6f\n",
	    time1-time0);
    fprintf(outfile,"running time of connected_component_BH:%9.6f\n",
	    time_BH);
    fflush(outfile);
#endif
    fclose(outfile);
    exit(0);
  }
  

    if (distmethod == DIST_BP) { /* fill matrix with breakpoint distances */
      setBPmatrix(distmatrix,genome_list,NUM_GENES,NUM_GENOMES,&distmem,CIRCULAR);
    }
    else { /* fill matrix with inversion distances */
      setinvmatrix(distmatrix,genome_list,NUM_GENES,NUM_GENOMES,&distmem,CIRCULAR);

      if (CORRECTION) {
          for (i = 0; i < NUM_GENOMES; i++) {
                  for (j = 0; j < NUM_GENOMES; j++) {
                      if (i != j) distmatrix[i][j]=ede(distmatrix[i][j], NUM_GENES);
                  }
              }
      }
    }

/* are used for the branch and bound */
  minidist = 1000000;
  for (j = 0; j < NUM_GENOMES-1; j++) {
      for (k = j+1; k < NUM_GENOMES-1; k++) {
          if ((flag=distmatrix[j][NUM_GENOMES-1]+distmatrix[NUM_GENOMES-1][k]-distmatrix[j][k]) < minidist) 
              minidist = flag;
      }
  }
#ifdef MPBPA
#ifdef GMP
      mpz_mul_ui(PROCS_treelevelcount, 2, (unsigned long int)NUM_GENOMES);
      mpz_add(PROCS_treelevelcount, PROCS_treelevelcount, -5);
      mpz_mul_ui(PROCS_treelevelcount, PROCS_treelevelcount, (unsigned long int)PROCS);
      mpz_add(PROCS_treelevelcount2, PROCS_treelevelcount, -1);
#else
      PROCS_treelevelcount = (2*NUM_GENOMES -5)*PROCS;
      PROCS_treelevelcount2 = PROCS_treelevelcount - 1;
      PROCS_treelevelstepping = (2*NUM_GENOMES -5)*(MYPROC+1);
#endif
#endif
  treelevelcount = 2*NUM_GENOMES -5;
  treelevelcount2 = 2*NUM_GENOMES -6;
/* these two are reserved for the second level */
  treelevel2count = (2*NUM_GENOMES -5)*(2*NUM_GENOMES-7);
  treelevel2count2 = (2*NUM_GENOMES -5)*(2*NUM_GENOMES-7) -1;

  if ((newbest < best_so_far) || (best_so_far == 0)) best_so_far = newbest;

  count = 0;

  time(&start);
  
  newtime = gmtime( &start );
  fprintf(outfile, "Start time is %s\n", asctime( newtime ) );


  if (GROWTREE) {
      cutValue = MST(NUM_GENOMES, distmatrix);

      for (i = 0; i < NUM_GENOMES; i++) {
          for (j = 0; j < NUM_GENOMES; j++) {
              cpdist[i][j] = distmatrix[i][j];
              if (distmatrix[i][j] > cutValue) cpdist[i][j] = -1;
          }
      }
      triangulate(cpdist, distmatrix, 0.0, NUM_GENOMES);

      for (i = 0; i < NUM_GENOMES; i++) {
          orderGenome(genome_list, cpgenome, cpdist, distmatrix, NUM_GENES, NUM_GENOMES,cutValue, i);
          tree = NULL;
          score = growTree(COND,&tree,NUM_GENES,NUM_GENOMES,tspsolver,thresh,
				                      adj_list,adj_pool,neighbors,stack,incycle,
				                      outcycle,weights,degree,otherEnd,edges,
				                      tpool,labels,pred1,pred2,picked,decode,
				                      genes,CIRCULAR, distmethod,distmem, CORRECTION, qpool,
				                      edgepool, cpgenome, smalledges, status, triple,
				                      inittspsolver, OPT, initmethod, bin_env);
          print_tree_nexus(tree);
          if (best_so_far > score) {
              best_so_far = score;
          }
          fprintf(outfile, "Grow Tree %d\n", score);

      }
      fflush(outfile);

  }
  
  if (tspsolver == INVERSION_MEDIAN_FAST) init_global_variables(NUM_GENES);
  upperLowerBound = circular_ordering(NUM_GENOMES,NUM_GENES,
          0,neighbors, distmatrix,
          stack,outcycle,degree,otherEnd,edges,
          CIRCULAR)/2;
  lowerLowerBound = circular_ordering(NUM_GENOMES,NUM_GENES,
          1,neighbors, distmatrix,
          stack,outcycle,degree,otherEnd,edges,
          CIRCULAR)/2;


  /* Major division: layered method or not? */
  if (LAYERED) {

      /* start layered approach */
      best_so_far_2 = best_so_far*2;
      totalN = 1024*1024;
      count = best_so_far-lowerLowerBound+1;
      if (count < 1) {
         count = 1;
         lowerLowerBound = best_so_far;
      }
#ifndef CYGWIN_NT_40
      lobs = (long long**)malloc(sizeof(long long*)*count);
      lobn = (long long*)malloc(sizeof(long long)*count);
#else
      lobs = (__int64**)malloc(sizeof(__int64*)*count);
      lobn = (__int64*)malloc(sizeof(__int64)*count);
#endif
      lobt = (int*)malloc(sizeof(int)*count);
      for (i = 0; i < count; i++) {
#ifndef CYGWIN_NT_40
          lobs[i] = (long long*)malloc(sizeof(long long)*(totalN+1024));
#else
          lobs[i] = (__int64*)malloc(sizeof(__int64)*(totalN+1024));
#endif
          lobn[i] = 0;
          lobt[i] = 0;
      }
      foundTree = 0;
      currentLow = LARGENUM;
      step = 0;
      treenbr = 1;

      init_gen_bin(NUM_GENOMES,bin_env);
      tree = first(NUM_GENOMES,tpool,0, bin_env);

      while (tree != NULL) {
          tree->leaf = TRUE;
          prev = tree->tag;
          /* tree traversal */
          firstg = prev;
          colb = compute_co_score(distmatrix,tree,&prev);
          colb += distmatrix[firstg-1][prev-1];
          if (colb <= best_so_far_2 && TIGHTBOUND) {
              colb = test_all_co_score_fast(best_so_far_2, colb, distmatrix,tree,
                          tree->lChild, firstg, &prev);
          }
          if ((colb <= best_so_far_2) && TIGHTBOUND) {
			  colb = test_all_co_score(best_so_far_2, colb, distmatrix,tree,
                          tree->lChild, firstg, &prev, 0);
          }
          if (colb > best_so_far_2) goto CONTINUE;
		  colb = colb / 2;
          if (currentLow > colb) currentLow = colb;
          colb = colb - lowerLowerBound;
          if (colb < 0) colb = 0;

          tmpInt = &(lobn[colb]);
          lobs[colb][*tmpInt] = treenbr;
          if (*tmpInt >= totalN) {
             fprintf(outfile,"cache one\n");
             *tmpInt = 0;
             sprintf(cacheName, "%s_%d_%d", inputfname, colb+lowerLowerBound, lobt[colb]);
             /* possible memory overflow*/
             cacheFile = fopen(cacheName, "w"); /* need to compress*/
             for (i = 0; i < totalN; i++) {
#ifndef CYGWIN_NT_40
                 fprintf(cacheFile, "%lld\n", lobs[colb][i]);
#else
                 fprintf(cacheFile, "%I64d\n", lobs[colb][i]);
#endif
             }
             fprintf(cacheFile, "%d\n", -1);
             fclose(cacheFile);
             lobt[colb]++;
          }
          else {
             (*tmpInt)++;
          }

          CONTINUE:
          tree = next(NUM_GENOMES,tpool,stepping, bin_env);
          treenbr += stepping;
      }

      for (i = 0; i < count; i++) {
          fprintf(outfile,"cache one %d\n", colb);
          sprintf(cacheName, "%s_%d_%d", inputfname, i+lowerLowerBound, lobt[i]);
          /* possible memory overflow*/
          cacheFile = fopen(cacheName, "w"); /* need to compress*/
          foundTree = lobn[i];
          for (j = 0; j < foundTree; j++) {
#ifndef CYGWIN_NT_40
              fprintf(cacheFile, "%lld\n", lobs[i][j]);
#else
              fprintf(cacheFile, "%I64d\n", lobs[i][j]);
#endif
          }
          fprintf(cacheFile, "%d\n", -1);
          fclose(cacheFile);
      }

      free(lobn);
      for (i = 0; i < count; i++) {
         free(lobs[i]);
      }
      free(lobs);
#ifndef CYGWIN_NT_40
      workLob = (long long*)malloc(sizeof(long long)*(totalN+1024));
#else
      workLob = (__int64*)malloc(sizeof(__int64)*(totalN+1024));
#endif
      fprintf(outfile,"Finished caching\n");

      count = 0;
      while (best_so_far >= currentLow) {
          fprintf(outfile,"Layer %d, best %d, current %d\n", step, best_so_far,
                 currentLow);
          foundTree = 1;
          treenbr = 1;
          /* now iterate over all trees */
          init_gen_bin(NUM_GENOMES,bin_env);
          tree = first(NUM_GENOMES,tpool,0, bin_env);
          i = 0;
          j = lobt[currentLow-lowerLowerBound];
          while(i <= j) {
             sprintf(cacheName, "%s_%d_%d",inputfname, currentLow, i);
             cacheFile = fopen(cacheName, "r");
             if (cacheFile == NULL) {
                 fprintf(outfile, "No file %s\n", cacheName);
                 break;
             }
             for (k = 0; k < totalN; k++) {
#ifndef CYGWIN_NT_40
                fscanf(cacheFile, "%lld", &workLob[k]);
#else
                fscanf(cacheFile, "%I64d", &workLob[k]);
#endif 
                if (workLob[k] < 0) break;
             }
             fclose(cacheFile);
             for (k = 0; k < totalN; k++) {
                count++;
                treenbr = workLob[k];
                if (treenbr < 0) break;
                foundTree = treenbr-foundTree;
                if (foundTree>0) {
                    tree = next(NUM_GENOMES,tpool, foundTree, bin_env);
                }
                foundTree = treenbr;
                tree->leaf = TRUE;
                score = init_score_tree(COND,tree,NUM_GENES,NUM_GENOMES,tspsolver,thresh,
 				              adj_list,adj_pool,neighbors,stack,incycle,
 				              outcycle,weights,degree,otherEnd,edges,
 				              tpool,labels,pred1,pred2,picked,decode,
 				              genes,CIRCULAR, distmethod,distmem, CORRECTION, qpool,
 				              edgepool, genome_list, smalledges, status, triple,
 				              inittspsolver, OPT, initmethod);
                if (score <= best_so_far) {
#ifndef CYGWIN_NT_40
                       fprintf(outfile,"With Tree# %lld\n", foundTree);
#else
                       fprintf(outfile,"With Tree# %I64d\n", foundTree);
#endif
                        print_tree_uncondensed(tree,NUM_GENES,condense_succ,condense_decode,
                                               orig_num_genes);
                        best_so_far = score;
 	                switch (distmethod) {
 	                case DIST_BP:
 	                  fprintf(outfile,"breakpoint score = %12d\n",score);
 	                  fprintf(outfile,"inversion score  = %12d\n",
                           score_tree(tree,NUM_GENES,CIRCULAR,DIST_INV,&distmem, CORRECTION));
 	                  fflush(outfile);
 	                  break;
 	                case DIST_INV:
 	                  fprintf(outfile,"inversion score  = %12d\n",score);
 	                  fprintf(outfile,"breakpoint score = %12d\n",
 		                  score_tree(tree,NUM_GENES,CIRCULAR,DIST_BP,&distmem, CORRECTION));
 	                  fflush(outfile);
 	                  break;
 	                }
                }
             }
             i++;
          }
          step++;
          currentLow++;
      }
      free(lobt);
      free(workLob);
      treenbr = numTree(NUM_GENOMES);
      count = treenbr-count;
      /* end layered approach */
  }

  else {
    /* no layering */
    totalN = numTree(NUM_GENOMES);
    #ifdef MPBPA
    #ifdef GMP
      if (BRANCH) mpz_set_si(treenbr, 1);
      else  mpz_set_si(treenbr, MYPROC+1);
    #else
      if (BRANCH) treenbr = 1;  
      else treenbr = MYPROC+1;

    #endif
    #else
    #ifdef GMP
        mpz_set_si(treenbr, 1);
    #else
        treenbr = 1;
    #endif
    #endif
      /* now iterate over all trees */
      if (!constrained) {
        init_gen_bin(NUM_GENOMES,bin_env);
    #ifdef MPBPA
       if (!BRANCH)
           tree = first(NUM_GENOMES,tpool,MYPROC, bin_env);
       else 
           tree = first(NUM_GENOMES,tpool,0, bin_env);
    #else
        tree = first(NUM_GENOMES,tpool,0, bin_env);
    #endif
      }    
      else {
    #ifdef MPBPA
        init_const(&const_env,constString,const_tree,tpool);
        tree = first_const(const_tree,tpool,MYPROC,&const_env);
    #else
        init_const(&const_env,constString,const_tree,tpool);
        tree = first_const(const_tree,tpool,0,&const_env);
    #endif
      }

#ifdef MPBPA
      if (BRANCH && MYPROC!= 0) {
         tree = next(NUM_GENOMES,tpool,PROCS_treelevelstepping-1-treelevelcount, 
                     bin_env);
     #ifdef GMP
         mpz_add(treenbr, treenbr, PROCS_treelevelstepping);
         mpz_add(treenbr, treenbr, 1);
         mpz_sub(treenbr, treenbr, treelevelcount);
     #else
         treenbr += (PROCS_treelevelstepping-1-treelevelcount);
     #endif
      }
      if (BRANCH && MYPROC == 0) {
         tree = next(NUM_GENOMES,tpool,PROCS_treelevelcount2,bin_env);
     #ifdef GMP
         mpz_add(treenbr, treenbr, PROCS_treelevelcount2);
     #else
         treenbr += (PROCS_treelevelcount2); 
     #endif
      }
#endif
      tree->leaf = TRUE;

      while (tree != NULL && treenbr <= totalN) {

        best_so_far_2 = best_so_far*2; /* to avoid the divide in lower bound */
        tree->leaf = TRUE;
        prev = tree->tag;
        if (BOUND && upperLowerBound > best_so_far) {
          /* tree traversal */
          firstg = prev;
          colb = compute_co_score(distmatrix,tree,&prev);
          colb += distmatrix[firstg-1][prev-1];
          if (colb <= best_so_far_2 && TIGHTBOUND) {
              colb = test_all_co_score_fast(best_so_far_2, colb,distmatrix,tree,
                          tree->lChild, firstg, &prev);
          }
          if ((colb <= best_so_far_2) && TIGHTBOUND) {
              colb = test_all_co_score(best_so_far_2, colb, distmatrix,tree,
                          tree->lChild, firstg, &prev, 0);
          }

        }
        else {
          colb = 0;
        }

        
        if (best_so_far_2 < colb) {
            count++;

          /* if lower bound exceeds upper bound, no need to score this tree */
    #ifdef TESTING
          not_done_count++;
          /* fprintf(outfile,"Entire tree pruned: %d out of %d\n",not_done_count,treenbr); */
    #endif
        }
        else  { /* score it */
          score = init_score_tree(COND,tree,NUM_GENES,NUM_GENOMES,tspsolver,thresh,
				              adj_list,adj_pool,neighbors,stack,incycle,
				              outcycle,weights,degree,otherEnd,edges,
				              tpool,labels,pred1,pred2,picked,decode,
				              genes,CIRCULAR, distmethod,distmem, CORRECTION, qpool,
				              edgepool, genome_list, smalledges, status, triple,
				              inittspsolver, OPT, initmethod);
              /* print next best found -- will the first best tree that way */
        #ifdef GMP
	        fprintf(outfile,"Tree #");
	        mpz_out_str(outfile,10,treenbr);
	        fprintf(outfile,", score %5d\n",score);
        #else
        #endif
	        fflush(outfile);
              if (score <= best_so_far) {
#ifndef CYGWIN_NT_40
                    fprintf(outfile,"With Tree# %lld\n", treenbr);     
#else
                    fprintf(outfile,"With Tree# %I64d\n", treenbr);
#endif
                    print_tree_uncondensed(tree,NUM_GENES,condense_succ,condense_decode,
                                               orig_num_genes);
                    best_so_far = score;

	            switch (distmethod) {
	            case DIST_BP:
	              fprintf(outfile,"breakpoint score = %12d\n",score);
	              fprintf(outfile,"inversion score  = %12d\n",
                      score_tree(tree,NUM_GENES,CIRCULAR,DIST_INV,&distmem, CORRECTION));
	              fflush(outfile);
	              break;
	            case DIST_INV:
	              fprintf(outfile,"inversion score  = %12d\n",score);
	              fprintf(outfile,"breakpoint score = %12d\n",
		              score_tree(tree,NUM_GENES,CIRCULAR,DIST_BP,&distmem, CORRECTION));
	              fflush(outfile);
	              break;
	            }
          }
        } /* if ub < lb */
        if (!constrained) {
#ifdef MPBPA
          if (BRANCH && treenbr%(PROCS_treelevelcount) == 
             (PROCS_treelevelstepping-treelevelcount) && 
             treenbr > 1)
#else
          if (BRANCH && treenbr%treelevelcount == 0 && treenbr > 1)
#endif
          {

              flag = 0;
              tree = testonelevel(NUM_GENOMES,tpool,0, bin_env,
                  &flag, distmatrix,
                  best_so_far_2, minidist);
               while (flag == 1) {
                  flag = 0;

            #ifdef MPBPA
                  tree = next(NUM_GENOMES,tpool,PROCS_treelevelcount2, bin_env);
            #else
                  tree = next(NUM_GENOMES,tpool,treelevelcount2, bin_env);
            #endif
                  tree = testonelevel(NUM_GENOMES,tpool,stepping, bin_env, &flag, distmatrix,
                                   best_so_far_2, minidist);
            #ifdef MPBPA
            #ifdef GMP
                  mpz_add(treenbr, treenbr, PROCS_treelevelcount);
                  mpz_add(count, count, PROCS_treelevelcount);
            #else
                  treenbr += PROCS_treelevelcount;
                  count += PROCS_treelevelcount;
            #endif
            #else
                  treenbr += treelevelcount;
                  count += treelevelcount;
            #endif
               }
           }
#ifdef MPBPA
           else if (treenbr%treelevelcount == 0 && treenbr > 1 
             && BRANCH)
           {
               tree = next(NUM_GENOMES,tpool,PROCS_treelevelcount2-treelevelcount-1, bin_env);
               treenbr += (PROCS_treelevelcount-treelevelcount-1);
           }
#endif
           else {
	
    #ifdef MPBPA
          tree = next(NUM_GENOMES,tpool,PROCS_stepping, bin_env);
    #else
          tree = next(NUM_GENOMES,tpool,stepping, bin_env);
    #endif
           }
        }
        else {
    #ifdef MPBPA
          tree = next_const(const_tree,tpool,PROCS_stepping,&const_env);
    #else
          tree = next_const(const_tree,tpool,stepping,&const_env);
    #endif
        }
    #ifdef MPBPA
       #ifdef GMP
          mpz_add(treenbr, treenbr, PROCS_stepping);
       #else
          treenbr += PROCS_stepping;
       #endif
    #else
       #ifdef GMP
          mpz_add(treenbr, treenbr, stepping);
       #else
          treenbr += stepping;
       #endif

    #endif
          
    }
  }
  /* all done; gather stats */
  time(&end);
  newtime = gmtime( &end );
  fprintf(outfile, "End time is %s\n", asctime( newtime ) );
  fflush(outfile);

#ifdef MPBPA
  MPI_Barrier(MPI_COMM_WORLD);
  end_time=MPI_Wtime();
  MPI_Reduce(&best_so_far, &i, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  if (MYPROC == 0) {
    fprintf(outfile,"Running time: %9.6f\n",end_time-start_time);
    fprintf(outfile,"The best tree has score: %d\n",i);
    fflush(outfile);
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (tspsolver == INVERSION_MEDIAN_FAST) free_variables();

#ifndef CYGWIN_NT_40
  fprintf(outfile, "Pruned %lld of %lld (%f), used %ds, best score = %d\n",
          count, treenbr, count/(treenbr-0.0),
          (int)(end-start), best_so_far);
#else
  fprintf(outfile, "Pruned %I64d of %I64d (%f), used %ds ,best score = %d\n",
          count, treenbr, count/(treenbr-0.0),
          (int)(end-start), best_so_far);
#endif
  
  if (outfile != stdout)
    fclose(outfile);
  
  /* free memory before stopping */
#ifdef GMP
  mpz_clear(treenbr);
  mpz_clear(stepping);
#ifdef MPBPA
  mpz_clear(PROCS_stepping);
  mpz_clear(PROCS_treelevelcount);
  mpz_clear(PROCS_treelevelcount2);
  mpz_clear(PROCS_treelevel2count);
  mpz_clear(PROCS_treelevel2count2);
#endif
#endif
  
  /*free memory used by tree generator*/
  if (!constrained) {
    free_bin_tree_env(NUM_GENOMES,bin_env);
    free_const_tree_env(NUM_GENOMES,&const_env);
  }
  else {
     free_const_tree_env(const_env.num_leaves,&const_env);
  } 
    
  /*free memory used by inversion distance calculation*/

  UFfree(distmem.uf);
  free(distmem.components);
  free(distmem.labeled);
  free(distmem.cc);
  free(distmem.oriented);
  free(distmem.stack);
  free(distmem.greyEdges);
  free(distmem.done);
  free(distmem.perm1);
  free(distmem.perm2);
  free(distmem.perm);
  free(distmem.hammingArr);

  if (COND) {
     decode += -NUM_GENES;
     free(decode);
  }
  free(adj_list);
  free(adj_pool);
  free(picked);
  for (i=0 ; i< (2*NUM_GENOMES+6) ; i++)
    free(labels[i].genes);
  free(labels);
  free(neighbors);
  free(tpool);
  if (initmethod==BIGNN) {
     free(triple);
  }
  if ((initmethod==FASTPROP) || (initmethod==MEDIANPROP) ||
      (initmethod==SMALLNN) || (initmethod==SBNN) || (initmethod==BIGNN)) {
     free(qpool);
  }
  free(edgepool);
  free(edges);
  free(otherEnd);
  free(degree);

  
  if (initmethod==ADJPARS) {
    for (i=0; i < 2*NUM_GENES; i++)
      free(status[i]);
    free(status);
  }

  if (
#ifdef CONCORDE
      (tspsolver==TSP_CHLINKERN) || (tspsolver==TSP_GREEDYLK) ||
      (inittspsolver==TSP_CHLINKERN) || (inittspsolver==TSP_GREEDYLK) ||
#endif
      initmethod==ADJPARS) {
    for (i=0; i < 2*NUM_GENES; i++)
      free(weights[i]);
    free(weights);
  }
  free(stack);
  free(pred1);
  free(pred2);
  free(outcycle);
  free(incycle);
  free(genes);

  if (GROWTREE) {
    for (i = 0; i < NUM_GENOMES+1; i++){ 
      free(cpgenome[i].genes);
    }
    free(cpgenome);
  }
  for (i=0; i < NUM_GENOMES; i++) 
    free(cpdist[i]);
  free(cpdist);

  for (i=0; i < NUM_GENOMES; i++) 
    free(distmatrix[i]);
  free(distmatrix);

  free(condense_decode);
  free(condense_succ);
  for (i=0 ; i<NUM_GENOMES ; i++) {
    free(genome_list[i].genes);
    free(genome_list[i].gnamePtr);
  }
  free(genome_list);
#ifdef MPBPA
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  return(0);
}

#ifndef CYGWIN_NT_40
long long numTree(int num_genomes) {
#else 
__int64 numTree(int num_genomes) {
#endif
    if (num_genomes <= 3) return 1;
    return (2*num_genomes-5)*numTree(num_genomes-1);
}


/* this is the helper function to compute the score of a tree*/				
int init_score_tree(int COND,struct tNode *tree,int NUM_GENES,int NUM_GENOMES,
		      int tspsolver,int thresh,
		      struct adj_struct *adj_list,struct adj_struct *adj_pool,
		      intpair_t *neighbors,
		      int *stack,int *incycle,int *outcycle,int **weights,
		      int *degree,int *otherEnd,edge_t *edges,
		      struct tNode *tpool,struct genome_struct *labels,
		      int *pred1,int *pred2,int *picked,int *decode,
		      int *genes, int CIRCULAR,
		      int distmethod, distmem_t distmem, int CORRECTION,
            struct qNode *qpool, int *edgepool, struct genome_struct *genome_list,
            smalledge_t *smalledges, int ** status, triple_t* triple,
            int inittspsolver, int OPT, int initmethod)
            /* remember, you must set tree->leaf = TRUE before you call this one*/
{
    int index, i, score;
    index = 0;
    SetTreeEdgePtrs(tree,edgepool,&index);
    for (i=0; i<2*NUM_GENOMES+6; i++) {
        tpool[i].genome = NULL;
    }
    add_genomes_to_tree(tree,labels,genome_list,NUM_GENES);

    switch (initmethod) {
        case RAND:
	        initialize_tree_random(tree,labels,NUM_GENES,NUM_GENOMES);
	        break;
        case SMALLNN:
	        initialize_tree_SNN(COND,tree,tpool,labels,qpool,adj_list,adj_pool,
			            weights,stack,degree,otherEnd,neighbors,edges,
			            incycle,outcycle,
			            pred1,pred2,picked,decode,NUM_GENES,NUM_GENOMES,
			            inittspsolver,thresh,CIRCULAR);
	        /* reuses pred1, pred2, picked */
	        break;
        case SBNN:
	        initialize_tree_SBNN(COND,tree,tpool,labels,qpool,adj_list,adj_pool,
			         weights,stack,degree,otherEnd,neighbors,edges,
	         incycle,outcycle,
	         pred1,pred2,picked,decode,NUM_GENES,NUM_GENOMES,
	         inittspsolver,thresh,CIRCULAR);
            /* reuses pred1, pred2, picked */
           break;
        case BIGNN:
           initialize_tree_BNN(COND,tree,tpool,labels,triple,adj_list,adj_pool,
	            weights,stack,degree,otherEnd,neighbors,edges,
	            incycle,outcycle,
	            pred1,pred2,picked,decode,NUM_GENES,NUM_GENOMES,
	            inittspsolver,thresh,CIRCULAR);
           /* reuses pred1, pred2, picked */
           break;
        case FASTPROP:
           initialize_tree_propagate(COND,tree,tpool,labels,qpool,adj_list,adj_pool,
		          weights,stack,degree,otherEnd,neighbors,edges,
		          incycle,outcycle,NUM_GENES,NUM_GENOMES,
		          pred1,pred2,picked,decode,FAST,inittspsolver,
		          thresh,CIRCULAR);
           break;
        case MEDIANPROP:
           initialize_tree_propagate(COND,tree,tpool,labels,qpool,adj_list,adj_pool,
		          weights,stack,degree,otherEnd,neighbors,edges,
		          incycle,outcycle,NUM_GENES,NUM_GENOMES,
		          pred1,pred2,picked,decode,MEDIAN,inittspsolver,
		          thresh,CIRCULAR);
           break;
        case TRIV:
           initialize_tree_trivial(tree,labels,NUM_GENES,NUM_GENOMES);
           break;
        case ADJPARS:
           initialize_tree_adjpars(tree,tpool,labels,stack,degree,otherEnd,
		        neighbors,smalledges,incycle,outcycle,
		        NUM_GENES,NUM_GENOMES,
		        weights,status,inittspsolver,thresh,CIRCULAR);
           /* reuses pred1, pred2, picked */
           break;
        default:
           fprintf(stderr,"ERROR: no initialization given\n");
    }

    score = score_tree(tree,NUM_GENES,CIRCULAR,distmethod,&distmem, CORRECTION);
    if (OPT) {
          score = iterate_over_tree(COND,tree,NUM_GENES,NUM_GENOMES,tspsolver,thresh,
		          score,adj_list,adj_pool,neighbors,stack,incycle,
		          outcycle,weights,degree,otherEnd,edges,
		          tpool,labels,pred1,pred2,picked,decode,
		          genes,CIRCULAR, distmethod,&distmem, CORRECTION);
    }

    return score;
}
