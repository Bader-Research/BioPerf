/*
 * clustalw-smp 0.99-9
 *
 * This file measures thread creation performance on a given machine
 *
 * Plant Biotechnology Institute, National Research Council of Canada
 *
 * Author: Ognen Duzlevski, http://www.aquaos.org/index.php
 * ognen@gene.pbi.nrc.ca
 */

#include <pthread.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#define THREADS 64	/* every system should be able to create 64 threads */
#define LOOPS 15
#define DEFAULT_USEC 150

/* a dummy thread */
void mythread() {
	pthread_exit(NULL);
}

float thr_create_usec() {
	pthread_t pt;
	pthread_attr_t attr;
	int i,j,res;
	long elapsed;
	float q;
	struct timeval tvstart, tvend;
	long overall = 0;
	
	if ((pthread_attr_init(&attr)) != 0)
		return DEFAULT_USEC;
	if ((pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM)) != 0)
		return DEFAULT_USEC;
	if ((pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE)) != 0)
		return DEFAULT_USEC;
	/* pthread_attr_setstacksize(&attr, 256); */
	for (j=0; j<LOOPS; j++) {
		res = gettimeofday(&tvstart, NULL);
		for (i=0; i<THREADS; i++) {
			res = pthread_create(&pt, &attr, (void *)mythread, NULL);
			if (res != 0)
				return DEFAULT_USEC;
		}
		res = gettimeofday(&tvend, NULL);
		if (tvend.tv_sec > tvstart.tv_sec)
			elapsed = tvend.tv_usec+(1000000*(tvend.tv_sec-tvstart.tv_sec)) - tvstart.tv_usec;
		else
			elapsed = tvend.tv_usec - tvstart.tv_usec;
		overall += elapsed;
	}
	q = (float)overall / LOOPS;
	return q/THREADS;
}