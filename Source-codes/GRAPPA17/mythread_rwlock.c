/* $Id: mythread_rwlock.c,v 1.3 2001/12/20 03:17:10 acs Exp $
   Written by Adam Siepel, Spring 2001 
   Copyright 2001, Adam Siepel */

/* Simple implementation of reader-writer locks that uses pthread
   mutexes and POSIX semaphores.  Interface is modeled on that of the
   "pthread_rwlock" functions, which seem not to be widely supported
   (otherwise I would have used them directly).  Note that the
   implementation here gives priority to readers (i.e., there is a
   possibility of starvation of writers).  */


#include "mythread_rwlock.h"
#include <pthread.h>

int mythread_rwlock_init(mythread_rwlock_t* rwlock) {
  rwlock->readcount = 0;
  pthread_mutex_init(&rwlock->countermutex, NULL);
  sem_init(&rwlock->writesem, 0, 1);

  return 0;
}

int mythread_rwlock_wrlock(mythread_rwlock_t* rwlock) {
  sem_wait(&rwlock->writesem);
  return 0;
}

int mythread_rwlock_rdlock(mythread_rwlock_t* rwlock) {
  pthread_mutex_lock(&rwlock->countermutex);
  rwlock->readcount++;
  if (rwlock->readcount == 1) sem_wait(&rwlock->writesem);
  pthread_mutex_unlock(&rwlock->countermutex);

  return 0;
}

int mythread_rwlock_rdunlock(mythread_rwlock_t* rwlock) {
  pthread_mutex_lock(&rwlock->countermutex);
  rwlock->readcount--;
  if (rwlock->readcount == 0) sem_post(&rwlock->writesem);
  pthread_mutex_unlock(&rwlock->countermutex);        

  return 0;
}

int mythread_rwlock_wrunlock(mythread_rwlock_t* rwlock) {
  sem_post(&rwlock->writesem);
  return 0;
}
