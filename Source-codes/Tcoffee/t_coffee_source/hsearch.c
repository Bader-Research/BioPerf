
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

float HASCH_K;

void hstat(HaschT*T)
  {
    fprintf ( stderr, "\nSIZE: %d ALLOCATED: %d FREED: %d TOTAL: %d", T->ne, T->allocated, T->freed, T->total);
  }

HaschT * hcreate ( int n_elements)
       {
	 HaschT *T;
	 int a;
	 
	 n_elements=n_elements*2+1;

	 if ( HASCH_K==0)HASCH_K=(sqrt((double)5)-1)/2;
	 
	 T=vcalloc ( 1, sizeof (HaschT));
	 T->ne=n_elements;	 
	 T->p=vcalloc (n_elements,sizeof ( Hasch_entry*));
	 for ( a=0; a<n_elements; a++)
	   {
	     T->p[a]=allocate_hasch_entry(NULL,ADD);
	     (T->p[a])->k=-1;
	   }
	 return T;
       }
HaschT *hdestroy (HaschT *T)
       {
	 int a;
	 Hasch_entry *p, *pp;
	 
	 if ( T==NULL)return NULL;

	 for (a=0; a< T->ne; a++)
	   {
	     p=T->p[a];
	     while (p)
	       {
		 pp=p;
		 p=p->next;
		 allocate_hasch_entry(pp,REMOVE);
	       }
	   }
	 free ( T);
	 return NULL;
       }


double hsearch (HaschT *T,double data, int k, int action)
       {
	 int h;
	 /*int low;*/
	 double a;
	 Hasch_entry *p, *pp;
	

	 

	 /* find the key: k->h*/
	 
	 h=k%T->ne;
	 

	 if ( action==ENTER)
	   {	
	    
	     if (data==UNDEFINED)
	       {
		 hsearch (T, data, k, REMOVE);
	       }
	     else
	       {
		 p=pp=(T->p[h]);
		 while (p && pp->k!=k)
		   {
		     pp=p;
		     p=p->next;		 
		   } 
		 
		 
		 if (pp->k==k)p=pp;
		 else
		   {
		     
		     T->allocated++;
		     T->total++;
		     pp->next=allocate_hasch_entry(NULL,ADD);
		     p=pp->next;
		     p->previous=pp;
		     
		   }
		 
		 p->k=k;
		 p->data=data;
	       }
	     
	     return data;
	   }
	 else if (action==ADD)
	   {
	     
	     a=hsearch ( T,0, k, FIND);
	     
	     if ( a==UNDEFINED)return hsearch ( T,data, k,ENTER);
	     else return hsearch (T,a+data, k,ENTER);
	   }
	 else if (action==MULTIPLY)
	   {
	     a=hsearch ( T,0, k, FIND);
	     if ( a==UNDEFINED)return hsearch ( T,data, k,ENTER);
	     else return hsearch (T,a*data, k,ENTER);
	   }
	 else if ( action==FIND)
	   {
	     p=T->p[h];
	     while (p && p->k!=k)p=p->next;
	     if (!p)return UNDEFINED;
	     else return p->data;
	   } 
	 else if ( action==REMOVE)
	   {
	    
	     if ( hsearch ( T, 0, k, FIND)==UNDEFINED)return UNDEFINED;
	     p=T->p[h];
	     while (p && p->k!=k)p=p->next;
	     if ( !p)return UNDEFINED;
	     else 
	       {
		 
		
		 pp=p->previous;
		 pp->next=p->next;
		 if (pp->next)(pp->next)->previous=pp;		 		 
		 allocate_hasch_entry(p, REMOVE);
		
		 T->freed++;
		 T->total--;

	       }
	     return UNDEFINED;
	     
	   }
	 else
	   {
	     fprintf ( stderr, "\nERROR: Unknown action in hsearch [FATAL:%s]", PROGRAM);	     
	     exit(0);
	     return UNDEFINED;
	     
	   }
       }


	 
Hasch_entry * allocate_hasch_entry ( Hasch_entry * p, int action)
   {
     static Hasch_entry *first;
     static Hasch_entry *last;
     static int tot;
     int a;

     if ( action==DESTROY_STACK)allocate_hasch_entry(last->previous, DESTROY);
     else if ( action==DESTROY)
       {
	 free (p);
	 tot--;
       }
     else if (action==REMOVE)
       {
	 p->previous=last;
	 p->next=NULL;
	 last->next=p;
	 	 
	 p->k=0;
	 p->data=UNDEFINED;
	 
	 last=p;
	 return NULL;
       }
     else if (action==ADD )
       {
	 
	 if (!first)
	   {
	     first=vcalloc ( 1, sizeof (Hasch_entry));
	     first->previous=first;
	     last=first;
	   }

	 if (last->previous==first)
	   {
	     tot++;
	     for ( a=0; a< 1; a++)
	       {
		 p=vcalloc ( 1, sizeof (Hasch_entry));
		 p->previous=last;
		 last->next=p;
		 last=p;
	       }
	   }

	
	 p=last;

	 
	 last=last->previous;
	 p->previous=NULL;
	 p->next=NULL;
	 return p;
       }
     else 
       {
	 fprintf ( stderr, "ERROR: allocate_hasch_entry, unknown mode [FATA:%s]", PROGRAM);
	 crash ("");
	 
       }
     return NULL;
   }
