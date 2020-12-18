/* $Id: vertex_factory.c,v 1.5 2001/12/19 22:32:05 acs Exp $
   Written by Adam Siepel, Spring 2001 
   Copyright 2001, Adam Siepel */

/* "Factory" for vertices, allowing more efficient memory allocation
   and deallocation.  Compile with -DTHREADSAFE for concurrent access
   from multiple threads. */

#include <stdlib.h>
#include "vertex_factory.h"
#include <stdio.h>

/* Last two arguments pertain to a function to be called when capacity
   has been reached, before reallocating more memory.  Pass NULL for
   both if you do not wish to specify such a function */
VertexFactory* new_vf(int startsize, int ngenes, void (*clear_mem)(void*), 
                      void* arg) {
  int i;
  VertexFactory *vf;
  vf = (VertexFactory*)malloc(sizeof(VertexFactory));
  vf->vertices = (Vertex*)malloc(startsize * sizeof(Vertex));
  vf->perms = (int*)malloc(startsize * ngenes * sizeof(int));
  vf->nalloc = 0;
  vf->ngenes = ngenes;
  vf->capacity = startsize;
  vf->allocarray = (short*)malloc(startsize * sizeof(short));
  vf->cursoridx = 0;
  vf->clear_mem = clear_mem;
  vf->clear_mem_arg = arg;
  for (i = 0; i < startsize; i++) {
    vf->vertices[i].perm = &(vf->perms[i*ngenes]);
    vf->allocarray[i] = 0;
  }

#ifdef THREADSAFE
  pthread_mutex_init(&vf->mutex, NULL);
#endif

  return vf;
}

Vertex* get_vertex(VertexFactory *vf) {
  int oldcap, i;
  Vertex* retval;

#ifdef THREADSAFE
  pthread_mutex_lock(&vf->mutex);
#endif

  if (vf->nalloc == vf->capacity) {
    if (vf->clear_mem != NULL) {
#ifdef THREADSAFE
      pthread_mutex_unlock(&vf->mutex);
#endif
      fprintf(stderr, "Invoking clear_mem\n");
      vf->clear_mem(vf->clear_mem_arg);
#ifdef THREADSAFE
      pthread_mutex_lock(&vf->mutex);
#endif
    }

    if (vf->nalloc == vf->capacity) {
      oldcap = vf->capacity;
      vf->capacity *= 2;
      vf->vertices = (Vertex*)realloc(vf->vertices, 
                                      vf->capacity * 
                                      sizeof(Vertex));
      vf->perms = (int*)realloc(vf->perms, 
                                vf->capacity*vf->ngenes*sizeof(int));
      vf->allocarray = (short*)realloc(vf->allocarray, 
                                       vf->capacity * sizeof(short));
      for (i = oldcap; i < vf->capacity; i++) {
        vf->vertices[i].perm = &(vf->perms[i*vf->ngenes]);
        vf->allocarray[i] = 0;
      }
    }
  }

  if (vf->cursoridx == vf->capacity) vf->cursoridx = 0;
  while (vf->allocarray[vf->cursoridx] == 1) {
    vf->cursoridx++;
    if (vf->cursoridx == vf->capacity) vf->cursoridx = 0;
  }

  vf->nalloc++;
  vf->vertices[vf->cursoridx].memidx = vf->cursoridx;
  vf->allocarray[vf->cursoridx] = 1;

  retval = &(vf->vertices[vf->cursoridx++]);

#ifdef THREADSAFE
  pthread_mutex_unlock(&vf->mutex);
#endif

  return retval;
}

void return_vertex(VertexFactory *vf, Vertex* v) {
#ifdef THREADSAFE
  pthread_mutex_lock(&vf->mutex);
#endif

/*   printf("%d\n", v->memidx); */
  vf->nalloc--;
  vf->allocarray[v->memidx] = 0;

#ifdef THREADSAFE
  pthread_mutex_unlock(&vf->mutex);
#endif
}

void vf_free(VertexFactory *vf) {
  free(vf->vertices);
  free(vf->perms);
  free(vf->allocarray);
  free(vf);
}

