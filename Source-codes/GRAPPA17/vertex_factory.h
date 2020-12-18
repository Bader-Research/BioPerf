/* $Id: vertex_factory.h,v 1.8 2001/12/19 22:32:05 acs Exp $
   Written by Adam Siepel, Spring 2001 
   Copyright 2001, Adam Siepel */

/* "Factory" for vertices, allowing more efficient memory allocation
   and deallocation.  Compile with -DTHREADSAFE for concurrent access
   from multiple threads. */

#ifndef VFACT_H
#define VFACT_H

#include <pthread.h>

typedef struct vertex Vertex;
struct vertex {
  int *perm;
  int distance;                 /* distance from home */
  int best_possible_score;      /* score of best possible median */
  int worst_possible_score;     /* score of worst possible median */
  int memidx;                   /* used by vertex factory */
  int d1, d2;
};

typedef struct vertex_factory VertexFactory;
struct vertex_factory {
  Vertex* vertices;
  int* perms;
  int nalloc;
  int ngenes;
  int capacity;
  short* allocarray;
  int cursoridx;
  void (*clear_mem)(void*);
  void *clear_mem_arg;
#ifdef THREADSAFE
  pthread_mutex_t mutex;
#endif
};

VertexFactory* new_vf(int startsize, int ngenes, void (*clear_mem)(void*), 
                      void *arg);
Vertex* get_vertex(VertexFactory *vf);
void return_vertex(VertexFactory *vf, Vertex* v);
void vf_free(VertexFactory *vf);

#endif
