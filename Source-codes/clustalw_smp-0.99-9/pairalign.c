/*
 * pairalign.c
 *
 * SMP (0.99-9) version of clustalw 1.82
 * 
 * author: Ognen Duzlevski, http://www.aquaos.org/index.php
 * Plant Biotechnology Institute, National Research Council of Canada
 *
 * Changes from the original (non-SMP version):
 *  - added parallel() which implements a single thread
 *  - added envget_maxthreads() which reads an environment variable and sets
 *    the number of supported threads and associated queues accordingly
 *  - added setup_threads() which sets up the whole mechanism for threads
 *  - localized the global memory shared variables into heap allocated
 *      structs accessible only to th thread who allocated the structure
 *  - added some header files and reorganized things to suit my model
 */
/* Change int h to int gh everywhere  DES June 1994 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "clustalw.h"

/* SMP */
#include <pthread.h>
#include <errno.h>
/*
 *   Global variables
 */
#ifdef MAC
#define pwint   short
#else
#define pwint   int
#endif
sint     int_scale;

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define gap(k)  ((k) <= 0 ? 0 : g + gh * (k))
#define tbgap(tb,k,gh)  ((k) <= 0 ? 0 : tb + gh * (k))
#define tegap(te,k,gh)  ((k) <= 0 ? 0 : te + gh * (k))

extern double   **tmat;
extern float    pw_go_penalty;
extern float    pw_ge_penalty;
extern float    transition_weight;
extern sint     nseqs;
extern sint     max_aa;
extern sint     gap_pos1,gap_pos2;
extern sint     max_aln_length;
extern sint     *seqlen_array;
extern sint     mat_avscore;
extern short    blosum30mt[],pam350mt[],idmat[],pw_usermat[],pw_userdnamat[];
extern short    clustalvdnamt[],swgapdnamt[];
extern short    gon250mt[];
extern short    def_dna_xref[],def_aa_xref[],pw_dna_xref[],pw_aa_xref[];
extern Boolean  dnaflag;
extern char     **seq_array;
extern char     *amino_acid_codes;
extern char     pw_mtrxname[];
extern char     pw_dnamtrxname[];

sint     matrix[NUMRES][NUMRES];
float   gscale,ghscale;

#include "parallel.h"
    
/*
*   Prototypes
*/
void add(sint v, global_t *gmem);
sint calc_score(sint iat, sint jat, sint v1, sint v2, global_t *gmem);
float tracepath(sint tsb1,sint tsb2, global_t *gmem);
void forward_pass(char *ia, char *ib, sint n, sint m, global_t *gmem);
void reverse_pass(char *ia, char *ib, global_t *gmem);
sint diff(sint A, sint B, sint M, sint N, sint tb, sint te, global_t *gmem);
void del(sint k, global_t *gmem);

/*
 * this is the thread function that does the work
 * each thread function has its own queue
 * the thread will take work requests off queue's start
 * and the requests are added to the queue's end
 * a NULL entry means the end of the queue
 */
void parallel(qentry *start) {
    int len2, m, i, j;
    /*
     * There is a rather lengthy story behind this volatile declaration:
     * on x86 platforms using the gcc compiler, the code in the threaded
     * version is left to operate in higher precision (this is the default)
     * however, either because C standard does not address this issue or due
     * to gcc non-conformance, the mm_score is kept on a higher precision
     * than is needed. The volatile keyword will force this rounding
     * behavior and thus equal the results of the threaded app to the
     * sequential one. To see what I am talking about, on a Linux on an x86 
     * platform, remove the volatile keyword and recompile and run the
     * program. On certain sequences the results will be off by 0.00001 in
     * the .dnd file but this should not affect the .aln file. Since we are
     * using volatile, we are killing some optimizations especially on
     * platforms that do not exhibit this behavior. On these platforms you
     * are advised to remove volatile, recompile and go about your business.
     */
    volatile float mm_score;
    char c;
    
    /* we are isolating thread variables in a heap allocated struct */
    global_t *gmem;
    
    /* work request from the qentry in the queue */
    parallel_t *pt;
    /* a qentry in the queue */
    qentry *st;
    
    /* wait until the thread has been told that it is OK to start working */
    pthread_mutex_lock(&thread_start_mutex);
        while (job_start == 0)
            pthread_cond_wait(&thread_start, &thread_start_mutex);
    pthread_mutex_unlock(&thread_start_mutex);

loop:
    if (start->next == NULL)
        goto loop_exit;
    /* get a request */
    pt = start->wr;
        
    /* start the work */
    gmem = malloc(sizeof(global_t));
    if (!gmem)
        fatal("parallel(...): Memory allocation error.");
    
    /* the allocations are thread specific */
    gmem->displ=(sint *)ckalloc(((max_aln_length+1) << 1) * sizeof(sint));
    gmem->HH = (pwint *)ckalloc((max_aln_length) * sizeof(pwint));
    gmem->DD = (pwint *)ckalloc((max_aln_length) * sizeof(pwint));
    gmem->RR = (pwint *)ckalloc((max_aln_length) * sizeof(pwint));
    gmem->SS = (pwint *)ckalloc((max_aln_length) * sizeof(pwint));
    
    if ((!gmem->HH) || (!gmem->DD) || (!gmem->RR) || (!gmem->SS))
        fatal("parallel(...): Memory allocation error.");
    m = seqlen_array[pt->sj+1];
    len2 = 0;
    for (i=1;i<=m;i++) {
        c = seq_array[pt->sj+1][i];
        if ((c !=gap_pos1) && (c != gap_pos2)) 
            len2++;
    }
    if (dnaflag) {
        gmem->g = 2 * (float)pw_go_penalty * int_scale*gscale;
        gmem->gh = pw_ge_penalty * int_scale*ghscale;
    } else {
        if (mat_avscore <= 0)
            gmem->g = 2 * (float)(pw_go_penalty + log((double)(MIN(pt->n,m))))*int_scale;
        else 
            gmem->g = 2 * mat_avscore * (float)(pw_go_penalty + log((double)(MIN(pt->n,m))))*gscale;
        gmem->gh = pw_ge_penalty * int_scale;
    }

    /*
        align the sequences
    */
    gmem->seq1 = pt->si+1;
    gmem->seq2 = pt->sj+1;

    forward_pass(&seq_array[gmem->seq1][0], &seq_array[gmem->seq2][0], pt->n, m, gmem);
    reverse_pass(&seq_array[gmem->seq1][0], &seq_array[gmem->seq2][0], gmem);

    gmem->last_print = 0;
    gmem->print_ptr = 1;

    /* use Myers and Miller to align two sequences -> includes SGI fix */
    
    gmem->maxscore = diff(gmem->sb1-1, gmem->sb2-1, gmem->se1-gmem->sb1+1,
                          gmem->se2-gmem->sb2+1, (sint)0, (sint)0, gmem);
 
    /* calculate percentage residue identity */
    mm_score = tracepath(gmem->sb1, gmem->sb2, gmem);

    if(pt->len1==0 || len2==0)
        mm_score=0;
    else
        mm_score /= (float)MIN(pt->len1,len2);

    tmat[pt->si+1][pt->sj+1] = ((float)100.0 - mm_score)/(float)100.0;
    tmat[pt->sj+1][pt->si+1] = ((float)100.0 - mm_score)/(float)100.0;
    
    /* do not use info() which uses vfprintf() family of functions */
    fprintf(stdout, "\nSequences (%d:%d) Aligned. Score:  %d", (pint)pt->si+1,(pint)pt->sj+1,(pint)mm_score);
    /* free memory and exit gracefully */
    gmem->displ=ckfree((void *)gmem->displ);
    gmem->HH=ckfree((void *)gmem->HH);
    gmem->DD=ckfree((void *)gmem->DD);
    gmem->RR=ckfree((void *)gmem->RR);
    gmem->SS=ckfree((void *)gmem->SS);

    /* free memory */
    free(gmem);
        
    /* destroy the start of the queue */
    st = start;
    start = start->next;
    free(st);
    /*  free up some more memory (dmalloc reported leak) */
    free(pt);

    goto loop;
loop_exit:
    /* memory-leak, albeit a small one (0.99-3) */
    free(start);
    
    /* signal that we are done */
    pthread_mutex_lock(&barrier_mutex);
        threads_num--;
        pthread_cond_signal(&barrier_cond);
    pthread_mutex_unlock(&barrier_mutex);
    
    pthread_exit(NULL);
}

/* 
 * following function will get the number of threads from the THREADS
 * environment variable
 */
int getenv_var() {
    char *s;
    int res;
    
    s = getenv("THREADS");
    if (!s) {
        fprintf(stdout, "\n\nSet the THREADS variable to the number of processors.\n");
        fprintf(stdout, "On a Linux system try export THREADS=<number>\n");
        fprintf(stdout, "On Solaris/Tru64/IRIX try setenv THREADS <number>\n\n");
        exit(-1);
    }
    res = atoi(s);
    if (res == 0) {
        fprintf(stdout, "\n\nTHREADS variable must be a positive number.\n");
        exit(-1);
    }
    return res;
}

/*
 * this function starts up the thread pool, initializes the individual
 * thread work request queues
 */
void bootstrap_threads(int thread_max) {
    int i;
    pthread_t pt;
    int res;
    pthread_attr_t pattr;
    
    pthread_attr_init(&pattr);
    pthread_attr_setdetachstate(&pattr, PTHREAD_CREATE_DETACHED);
    pthread_attr_setscope(&pattr, PTHREAD_SCOPE_SYSTEM);
    
    /* make the queues - some of the variables are global! */
    queue = malloc(sizeof(qinfo) * thread_max);
    if (!queue)
        fatal("%s.%s: %s\n", "pairalign.c", "setup_threads()", strerror(errno));
    
    /* prepare the threads and their queues */
    for (i=0; i < thread_max; i++) {
        queue[i].start = malloc(sizeof(qentry));
        if (!queue[i].start)
            fatal("%s.%s: %s\n", "pairalign.c", "setup_threads()", strerror(errno));
        queue[i].end = queue[i].start;
        queue[i].start->next = NULL;
        
        /* create the threads and feed them the queue starts */
        res = pthread_create(&pt, &pattr, (void *)parallel, queue[i].start);
        if (res != 0)
            fatal("%s.%s: %s\n", "pairalign.c", "setup_threads()", strerror(errno));
    }
}
        
sint pairalign(sint istart, sint iend, sint jstart, sint jend) {
    short *mat_xref;
    sint si, sj, i;
    sint n,m,len1;
    sint maxres;
    short *matptr;
    char c;
    int res;
    
    int thread_max, next;
    int thrcnt = 0;
    
    qentry *q;
    parallel_t *pt;

    /* obtain maximum num of threads and setup everything */
    thread_max = getenv_var();
    bootstrap_threads(thread_max);
    threads_num = thread_max;
    
#ifdef MAC
    int_scale = 10;
#else
    int_scale = 100;
#endif
    
    gscale=ghscale=1.0;
    if (dnaflag) {
        if (strcmp(pw_dnamtrxname, "iub") == 0) {
            matptr = swgapdnamt;
            mat_xref = def_dna_xref;
      } else if (strcmp(pw_dnamtrxname, "clustalw") == 0) { 
            matptr = clustalvdnamt;
            mat_xref = def_dna_xref;
            gscale=0.6667;
            ghscale=0.751;
      } else {
            matptr = pw_userdnamat;
            mat_xref = pw_dna_xref;
        }
        maxres = get_matrix(matptr, mat_xref, matrix, TRUE, int_scale);
        if (maxres == 0) 
            return((sint)-1);

        matrix[0][4]=transition_weight*matrix[0][0];
        matrix[4][0]=transition_weight*matrix[0][0];
        matrix[2][11]=transition_weight*matrix[0][0];
        matrix[11][2]=transition_weight*matrix[0][0];
        matrix[2][12]=transition_weight*matrix[0][0];
        matrix[12][2]=transition_weight*matrix[0][0];
    } else {
        if (strcmp(pw_mtrxname, "blosum") == 0) {
            matptr = blosum30mt;
            mat_xref = def_aa_xref;
        } else if (strcmp(pw_mtrxname, "pam") == 0) {
            matptr = pam350mt;
            mat_xref = def_aa_xref;
        } else if (strcmp(pw_mtrxname, "gonnet") == 0) {
            matptr = gon250mt;
            int_scale /= 10;
            mat_xref = def_aa_xref;
        } else if (strcmp(pw_mtrxname, "id") == 0) {
            matptr = idmat;
            mat_xref = def_aa_xref;
        } else {
            matptr = pw_usermat;
            mat_xref = pw_aa_xref;
        }

        maxres = get_matrix(matptr, mat_xref, matrix, TRUE, int_scale);
        if (maxres == 0) 
            return((sint)-1);
    }
    
    for (si=MAX(0,istart);si<nseqs && si<iend;si++) {
        n = seqlen_array[si+1];
        len1 = 0;
        for (i=1;i<=n;i++) {
            c = seq_array[si+1][i];
            if ((c!=gap_pos1) && (c != gap_pos2)) 
                len1++;
        }

        for (sj=MAX(si+1,jstart+1);sj<nseqs && sj<jend;sj++) {
            m = seqlen_array[sj+1];
            if(n==0 || m==0) {
                tmat[si+1][sj+1]=1.0;
                tmat[sj+1][si+1]=1.0;
                continue;
            }
            
            /* add a request to the next thread queue */
            q = malloc(sizeof(qentry));
            if (!q)
                fatal("%s.%s: %s\n", "pairalign.c", "pairalign()", "q malloc failed");
            queue[thrcnt].end->next = q;
            /* fill out the structure */
            queue[thrcnt].end->wr = malloc(sizeof(parallel_t));
            pt = queue[thrcnt].end->wr;
            if (!pt)
                fatal("%s.%s: %s\n", "pairalign.c", "pairalign()", "thrcnt malloc failed");
            pt->n = n;
            pt->len1 = len1;
            pt->si = si;
            pt->sj = sj;
            /* link up to the end of queue */
            queue[thrcnt].end = queue[thrcnt].end->next;
            queue[thrcnt].end->next = NULL;
            /* move on to the next thread */
            thrcnt++;
            /* and check if we need to roll-over to thread zero */
            if (thrcnt >= thread_max)
                thrcnt = 0;
        }           
    }
    /* this is the moment where we start our work! */
    pthread_mutex_lock(&thread_start_mutex);
        job_start++;
        pthread_cond_broadcast(&thread_start);
    pthread_mutex_unlock(&thread_start_mutex);

    /* we need a barrier here to wait for all threads to finish */
    pthread_mutex_lock(&barrier_mutex);
        while (threads_num > 0)
            pthread_cond_wait(&barrier_cond, &barrier_mutex);
    pthread_mutex_unlock(&barrier_mutex);
        
    /*  release some memory (dmalloc reported leak) */
    free(queue);
    return((sint)1);
}

void add(sint v, global_t *gmem)
{
    if(gmem->last_print<0) {
        gmem->displ[gmem->print_ptr-1] = v;
        gmem->displ[gmem->print_ptr++] = gmem->last_print;
    } else
        gmem->last_print = gmem->displ[gmem->print_ptr++] = v;
}

sint calc_score(sint iat,sint jat,sint v1,sint v2, global_t *gmem)
{
    sint ipos,jpos;
    sint ret;

    ipos = v1 + iat;
    jpos = v2 + jat;

    ret=matrix[(int)seq_array[gmem->seq1][ipos]][(int)seq_array[gmem->seq2][jpos]];

    return(ret);
}

float tracepath(sint tsb1,sint tsb2, global_t *gmem)
{
    char c1,c2;
    sint  i1,i2,r;
    sint i,k,pos,to_do;
    sint count;
    float score;
    char s1[100], s2[100];

    to_do=gmem->print_ptr-1;
    i1 = tsb1;
    i2 = tsb2;

    pos = 0;
    count = 0;
    for(i=1;i<=to_do;++i) {
        if(gmem->displ[i]==0) {
            c1 = seq_array[gmem->seq1][i1];
            c2 = seq_array[gmem->seq2][i2];

            if ((c1!=gap_pos1) && (c1 != gap_pos2) && (c1 == c2))
                count++;
            ++i1;
            ++i2;
            ++pos;
        } else {
            if((k=gmem->displ[i])>0) {
                i2 += k;
                pos += k;
            } else {
                i1 -= k;
                pos -= k;
            }
        }
    }

    score = 100.0 * (float)count;
    
    return(score);
}

void forward_pass(char *ia, char *ib, sint n, sint m, global_t *gmem)
{

    sint i,j;
    pwint f,hh,p,t;

    gmem->maxscore = 0;
    gmem->se1 = 0;
    gmem->se2 = 0;
    for (i=0;i<=m;i++) {
        gmem->HH[i] = 0;
        gmem->DD[i] = -(gmem->g);
    }

    for (i=1;i<=n;i++) {
        hh = p = 0;
        f = -(gmem->g);

        for (j=1;j<=m;j++) {
            f -= gmem->gh; 
            t = hh - gmem->g - gmem->gh;
            if (f<t) 
                f = t;

            gmem->DD[j] -= gmem->gh;
            t = gmem->HH[j] - gmem->g - gmem->gh;
            if (gmem->DD[j]<t) 
                gmem->DD[j] = t;

            hh = p + matrix[(int)ia[i]][(int)ib[j]];
            if (hh<f) 
                hh = f;
            if (hh<gmem->DD[j]) 
                hh = gmem->DD[j];
            if (hh<0) 
                hh = 0;

            p = gmem->HH[j];
            gmem->HH[j] = hh;

            if (hh > gmem->maxscore) {
                gmem->maxscore = hh;
                gmem->se1 = i;
                gmem->se2 = j;
            }
        }
    }
}

void reverse_pass(char *ia, char *ib, global_t *gmem)
{

    sint i,j;
    pwint f,hh,p,t;
    pwint cost;

    cost = 0;
    gmem->sb1 = 1;
    gmem->sb2 = 1;
    for (i=gmem->se2;i>0;i--) {
        gmem->HH[i] = -1;
        gmem->DD[i] = -1;
    }

    for (i=gmem->se1;i>0;i--) {
        hh = f = -1;
        if (i == gmem->se1) 
            p = 0;
        else 
            p = -1;

        for (j=gmem->se2;j>0;j--) {
            f -= gmem->gh; 
            t = hh - gmem->g - gmem->gh;
            if (f<t) 
                f = t;
            gmem->DD[j] -= gmem->gh;
            t = gmem->HH[j] - gmem->g - gmem->gh;
            if (gmem->DD[j]<t) 
                gmem->DD[j] = t;

            hh = p + matrix[(int)ia[i]][(int)ib[j]];
            if (hh<f) 
                hh = f;
            if (hh<gmem->DD[j]) 
                hh = gmem->DD[j];

            p = gmem->HH[j];
            gmem->HH[j] = hh;

            if (hh > cost) {
                cost = hh;
                gmem->sb1 = i;
                gmem->sb2 = j;
                if (cost >= gmem->maxscore) 
                    break;
            }
        }
        if (cost >= gmem->maxscore) 
            break;
    }
}

int diff(sint A, sint B, sint M, sint N, sint tb, sint te, global_t *gmem)
{
    sint type;
    sint midi,midj,i,j;
    int midh;
    pwint f, hh, e, s, t;
    
    if(N<=0)  {
        if(M>0) {
            del(M, gmem);
        }
        return(-(int)tbgap(tb, M,gmem->gh));
    }

    if (M<=1) {
        if(M<=0) {
            add(N, gmem);
            return(-(int)tbgap(tb, N,gmem->gh));
        }

        midh = -(tb+gmem->gh) - tegap(te, N,gmem->gh);
        hh = -(te+gmem->gh) - tbgap(tb, N, gmem->gh);
        if (hh>midh) 
            midh = hh;
            midj = 0;
            for(j=1;j<=N;j++) {
                hh = calc_score(1,j,A,B,gmem) - tegap(te, N-j,gmem->gh) - tbgap(tb, j-1,gmem->gh);
                if(hh>midh) {
                    midh = hh;
                    midj = j;
                }
            }

        if(midj==0) {
            del(1, gmem);
            add(N, gmem);
        } else {
            if(midj>1)
                add(midj-1, gmem);
            gmem->displ[gmem->print_ptr++] = gmem->last_print = 0;
            if(midj<N)
                add(N-midj, gmem);
        }
        return midh;
    }

    /* Divide: Find optimum midpoint (midi,midj) of cost midh */

    midi = M / 2;
    gmem->HH[0] = 0.0;
    t = -tb;
    for(j=1;j<=N;j++) {
        gmem->HH[j] = t = t - gmem->gh;
        gmem->DD[j] = t - gmem->g;
    }

    t = -tb;
    for(i=1;i<=midi;i++) {
        s=gmem->HH[0];
        gmem->HH[0] = hh = t = t - gmem->gh;
        f = t - gmem->g;
        for(j=1;j<=N;j++) {
            if ((hh = hh - gmem->g - gmem->gh) > (f = f - gmem->gh)) 
                f=hh;
            if ((hh = gmem->HH[j] - gmem->g - gmem->gh) > (e = gmem->DD[j] - gmem->gh)) 
                e=hh;
            hh = s + calc_score(i,j,A,B,gmem);
            if (f>hh) 
                hh = f;
            if (e>hh) 
                hh = e;

            s = gmem->HH[j];
            gmem->HH[j] = hh;
            gmem->DD[j] = e;
        }
    }   

    gmem->DD[0]=gmem->HH[0];
    gmem->RR[N]=0;
    t = -te;
    for(j=N-1;j>=0;j--) {
        gmem->RR[j] = t = t-gmem->gh;
        gmem->SS[j] = t - gmem->g;
    }

    t = -te;
    for(i=M-1;i >= midi;i--) {
        s = gmem->RR[N];
        gmem->RR[N] = hh = t = t - gmem->gh;
        f = t - gmem->g;

        for(j=N-1;j>=0;j--) {
            if ((hh = hh - gmem->g - gmem->gh) > (f = f - gmem->gh)) 
                f=hh;
            if ((hh = gmem->RR[j] - gmem->g - gmem->gh) > (e = gmem->SS[j] - gmem->gh)) 
                e=hh;
            hh = s + calc_score(i+1,j+1,A,B,gmem);
            if (f>hh) 
                hh = f;
            if (e>hh) 
                hh = e;

            s = gmem->RR[j];
            gmem->RR[j] = hh;
            gmem->SS[j] = e;
        }
    }

    gmem->SS[N]=gmem->RR[N];
    midh=gmem->HH[0] + gmem->RR[0];
    midj=0;
    type=1;
    for(j=0;j<=N;j++) {
        hh = gmem->HH[j] + gmem->RR[j];
        if(hh>=midh) {
            if(hh>midh || (gmem->HH[j]!=gmem->DD[j] && gmem->RR[j]==gmem->SS[j])) {
                midh=hh;
                midj=j;
            }
        }
    }

    for(j=N;j>=0;j--) {
        hh = gmem->DD[j] + gmem->SS[j] + gmem->g;
        if(hh>midh) {
            midh=hh;
            midj=j;
            type=2;
        }
    }

    /* Conquer recursively around midpoint  */

    if(type==1) {             /* Type 1 gaps  */
        diff(A,B,midi,midj,tb,gmem->g, gmem);
        diff(A+midi,B+midj,M-midi,N-midj,gmem->g,te, gmem);
    } else {
        diff(A,B,midi-1,midj,tb,0.0, gmem);
        del(2, gmem);
        diff(A+midi+1,B+midj,M-midi-1,N-midj,0.0,te, gmem);
    }

    return midh;       /* Return the score of the best alignment */
}

void del(sint k, global_t *gmem)
{
    if(gmem->last_print<0)
        gmem->last_print = gmem->displ[gmem->print_ptr-1] -= k;
    else
        gmem->last_print = gmem->displ[gmem->print_ptr++] = -k;
}
