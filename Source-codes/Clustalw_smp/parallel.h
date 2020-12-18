/*
 * parallel.h
 *
 * SMP (0.99-9) version of clustalw 1.82
 * Plant Biotechnology Institute, National Research Council of Canada (Saskatoon, SK) 
 * author: Ognen Duzlevski, http://www.aquaos.org/index.php
 *
 * This file contains all the structures and variables needed for the SMP
 * version
 */

/* attribute used to setup the threads */
pthread_attr_t pa_attr;

/*
 * The following structure is used within parallel to simulate the
 * global memory approach used in clustalw
 */
typedef struct {
	sint *displ;
	pwint *HH, *DD, *RR, *SS;
	pwint maxscore;
	sint last_print, print_ptr, se1, se2, sb1, sb2, g, gh, seq1, seq2;
} global_t;

/*
 *	The following structure is sufficient to bootstrap parallel(...)
 */ 
typedef struct {
	int n, len1, si, sj;
} parallel_t;

/*
 * Entry within a thread queue follows. Notice the mutex which is used
 * to control access to the linked list of requests.
 */
struct q_entry {
	parallel_t *wr;
	struct q_entry *next;
};

typedef struct q_entry qentry;

typedef struct {
	qentry *start;
	qentry *end;
} qinfo;

/* global variable for the thread queues */
qinfo *queue;

/* needed to signal threads to start working */
pthread_mutex_t thread_start_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t thread_start = PTHREAD_COND_INITIALIZER;

/* needed to control the barrier for all threads to finish */
pthread_mutex_t barrier_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t barrier_cond = PTHREAD_COND_INITIALIZER;

static threads_num = 0;		/* the number of currently active threads */
static job_start = 0;
