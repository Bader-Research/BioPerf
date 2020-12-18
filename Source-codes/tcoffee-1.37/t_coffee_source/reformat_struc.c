#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h" 
#include "define_header.h"



#define FATAL "fatal:reformat_struc"


char *     normalize_pdb_file  (char *name, char *out_file)
{
  char command[1000];
  if ( !is_pdb(name))
    {
      fprintf(stderr, "\nERROR %s is not a pdb file[FATAL:%s]\n", name, PROGRAM);
    }
 
  sprintf ( command, "extract_from_pdb -infile %s -atom ALL -chain FIRST -nodiagnostic > %s", name, out_file);
  system ( command);
  
  return out_file;
  }

Ca_trace * trim_ca_trace (Ca_trace *T, char *seq )
{
  /*This function removes from Ca trace all the residues that are not in st*/
  Alignment *ALN;
  Atom *A;  
  int a,l, s, r, is_r, is_s;
  int *seq_cache, *struc_cache;
  
  if ( strm ( T->seq,seq))return T;
  else
    {
     ALN=align_two_sequences (T->seq,seq, "est_idmat",-1, 0,"fasta_pair_wise");  
          
     struc_cache=vcalloc (ALN->len_aln+1, sizeof (int));
     seq_cache  =vcalloc (ALN->len_aln+1, sizeof (int));
     
      for ( r=0, s=0,a=0; a< ALN->len_aln; a++)
       {
	 is_r=!is_gap(ALN->seq_al[0][a]);
	 is_s=!is_gap(ALN->seq_al[1][a]);
	 
	 r+=!is_gap(ALN->seq_al[0][a]);
	 s+=!is_gap(ALN->seq_al[1][a]);
	 
	 if ( is_s && is_r)
	   {
	     struc_cache[r-1]=s-1;
	     seq_cache[s-1]=r-1;
	   }
	 else if ( is_s && !is_r)
	   {
	     seq_cache[s-1]=-1;
	   }
	 else if ( !is_s && is_r)
	   {
	     struc_cache[r-1]=-1;
	   }
       }
     

      T->ca=vrealloc ( T->ca, sizeof (Atom*)*(ALN->len_aln));
      T->peptide_chain=vrealloc ( T->peptide_chain, sizeof (Amino_acid*)*(ALN->len_aln));

      for ( a=0; a< T->n_atom; a++)
	{
	  
	  A=(T->structure[a]);
	  if ( struc_cache[A->res_num-1]==-1)continue;
	  A->res_num=struc_cache[A->res_num-1]+1;
	  if (strm (A->type, "CA")) T->ca[A->res_num-1]=A;
	  if ( strm (A->type, "CA"))(T->peptide_chain[A->res_num-1])->CA=A;
	  if ( strm (A->type, "C"))(T->peptide_chain[A->res_num-1] )->C=A;
	  if ( strm (A->type, "CB"))(T->peptide_chain[A->res_num-1])->CB=A;
	  if ( strm (A->type, "N"))(T->peptide_chain[A->res_num-1] )->N=A;
	}
     l=strlen(seq);
     for ( a=0;a< l; a++)
       {
	 if ( seq_cache[a]==-1)
	   {
	     T->ca[a]=NULL;
	     
	     T->peptide_chain[a]->CA=NULL;
	     T->peptide_chain[a]->C =NULL;
	     T->peptide_chain[a]->CB=NULL;
	     T->peptide_chain[a]->N=NULL;
	   }
       }
     free_aln (ALN);
     free(seq_cache);
     free(struc_cache);
     
	  
    }
  return T;
}

Ca_trace * read_ca_trace (char *name )
    {
	/*This function reads a pdb file into a Ca_trace structure*/
	 
	 

	int a, c, n;
	FILE *fp;
	char *tp_name;
	char *command;
	Atom *A;
	char res;
	char *buf;
	Ca_trace *T=NULL;
	
	
	
	tp_name=vtmpnam (NULL);
	command=vcalloc (strlen (tp_name)+strlen (name)+200, sizeof (char));
	
	sprintf ( command, "extract_from_pdb -infile %s -atom C CA CB N -chain FIRST -mode simple> %s", name, tp_name);
	system ( command);
	
	buf=vcalloc ( VERY_LONG_STRING, sizeof (char));
	n=count_n_line_in_file (tp_name );
	
	if ( !T)
	    {
	    T=vcalloc ( 1, sizeof ( Ca_trace));
	    declare_name (T->name);
	    }
	

	fp=vfopen (tp_name, "r");
	while ( (c=fgetc(fp))!='>');
	fscanf ( fp, "%s",T->name );
	while ( (c=fgetc(fp))!='\n');
	fscanf ( fp, "%s",buf );
	while ( (c=fgetc(fp))!='\n');
	
	T->len=strlen (buf);
	T->seq=vcalloc ( T->len+1, sizeof (char));
	buf=lower_string (buf);
	sprintf ( T->seq, "%s", buf);
	
	
	n+=T->len;
	T->structure=vcalloc ( n, sizeof (Atom*));
	for ( a=0; a< n; a++)T->structure[a]=vcalloc ( 1, sizeof (Atom));
	
	
	T->ca=vcalloc ( T->len+1, sizeof ( Atom*));
	
	a=0;

		
	while ((c=fgetc(fp))!=EOF)
	   {
	       ungetc(c, fp);
	       A=T->structure[a];
	       A->num=a;
	       fscanf (fp, "ATOM %s %c %*c %d %f %f %f\n",A->type, &res,&A->res_num, &A->x, &A->y, &A->z);
	       res=tolower (res);
	       if ( res!=T->seq[A->res_num-1])
	           {
		       fprintf ( stderr, "\nPROBLEM IN YOUR PDB FILE: STRUC != SEQ on RES %d [%c.. %c]\n", A->res_num, res, T->seq[A->res_num-1]);
		       fprintf ( stderr, "\n%s\n", T->seq);
		       exit (1);
		   }
	       
	       if ( strm ( A->type, "CA")) T->ca[A->res_num-1]=A;
	       a++;
	   }
	T->n_atom=a;



	T->peptide_chain=vcalloc (T->len, sizeof (Amino_acid*));
	for ( a=0; a< T->len; a++) T->peptide_chain[a]=vcalloc (T->len, sizeof (Amino_acid));				   
	for ( a=0; a< T->n_atom; a++)
	  {
	    A=T->structure[a];

	    if ( strm (A->type, "CA"))(T->peptide_chain[A->res_num-1])->CA=A;
	    if ( strm (A->type, "C"))(T->peptide_chain[A->res_num-1] )->C=A;
	    if ( strm (A->type, "CB"))(T->peptide_chain[A->res_num-1])->CB=A;
	    if ( strm (A->type, "N"))(T->peptide_chain[A->res_num-1] )->N=A;
	  }
	    

	vfclose (fp);	
	free (buf);
	free ( command);
	remove ( tp_name);
	free ( tp_name);
	
	return T;
    }
Ca_trace * hasch_ca_trace    ( Ca_trace *T)
    {
      
      T=hasch_ca_trace_nb (T);
      T=hasch_ca_trace_bubble (T);
      T=hasch_ca_trace_transversal (T);
      return T;
    }
Ca_trace * hasch_ca_trace_transversal ( Ca_trace *TRACE)
    {
	/*This function gets the Coordinates of a protein and computes the distance of each Ca to its 
	  
	  given a Ca,
	  Compute the distance between, CA-x and CA+x with x=[1-N_ca]
	  T->nb[a][0]-->Number of distances.
	  T->nb[a][1... T->nb[a][0]]-->ngb index with respect to the Ca chain
	  T->d_nb[a][1... T->d_nb[a][0]]-->ngb index with respect to the Ca chain
	*/
	
	int a, b, d;
	float dist;
	Atom *A, *B;
		
	Struct_nb *T;
	Pdb_param *PP;


	TRACE->Transversal=vcalloc ( 1, sizeof (Struct_nb));
	
	T=TRACE->Transversal;
	PP=TRACE->pdb_param;

	if ( !T->nb)T->nb=declare_int (TRACE->len+1, 1);
	if ( !T->d_nb)T->d_nb=declare_float (TRACE->len+1, 1);

	for (d=0,a=0; a< TRACE->len; a++)
	    {

	        for ( b=1; b<=PP->N_ca; b++)
		    {
		    if ( (a-b)<0 || (a+b)>=TRACE->len)continue;
		    A=TRACE->ca[a-b];
		    B=TRACE->ca[a+b];
		    dist=get_atomic_distance ( A, B);
		    
		    T->nb[a]=vrealloc ( T->nb[a], (++T->nb[a][0]+1)*sizeof (int));
		    T->nb[a][T->nb[a][0]]=b;

		    T->d_nb[a]=vrealloc ( T->d_nb[a], (T->nb[a][0]+1)*sizeof (float));
		    T->d_nb[a][T->nb[a][0]]=dist;
		    
		    d++;
		    }
		T->max_nb=MAX (T->max_nb, T->nb[a][0]);
	    }
	return TRACE;

    }			       

Ca_trace * hasch_ca_trace_nb ( Ca_trace *TRACE)
    {
	/*This function gets the Coordinates of a protein and computes the distance of each Ca to its 
	  T->N_ca Ngb.
	  The first Ngb to the left and to the right are excluded
	  Ngd to the left get negative distances
	  Ngb to the right receive positive distances
	  T->nb[a][0]-->Number of ngb.
	  T->nb[a][1... T->nb[a][0]]-->ngb index with respect to the Ca chain
	  T->d_nb[a][1... T->d_nb[a][0]]-->ngb index with respect to the Ca chain
	*/
	
	int a, b, d;
	float dist;
	Atom *A, *B;
		
	Struct_nb *T;
	Pdb_param *PP;


	TRACE->Chain=vcalloc ( 1, sizeof (Struct_nb));
	
	T=TRACE->Chain;
	PP=TRACE->pdb_param;

	if ( !T->nb)T->nb=declare_int (TRACE->len+1, 1);
	if ( !T->d_nb)T->d_nb=declare_float (TRACE->len+1, 1);

	for (d=0,a=0; a< TRACE->len; a++)
	    {
		for ( b=MAX(0,a-PP->N_ca); b< MIN( a+PP->N_ca, TRACE->len); b++)
		        {
			if (FABS(a-b)<2)continue;
			A=TRACE->ca[a];
			B=TRACE->ca[b];
			if ( !A || !B)continue;
			dist=get_atomic_distance ( A, B);
			if (b<a) dist=-dist;
			
			T->nb[a]=vrealloc ( T->nb[a], (++T->nb[a][0]+1)*sizeof (int));
			T->nb[a][T->nb[a][0]]=b;
			
			T->d_nb[a]=vrealloc ( T->d_nb[a], (T->nb[a][0]+1)*sizeof (float));
			T->d_nb[a][T->nb[a][0]]=dist;
			d++;
			}
		    
		T->max_nb=MAX (T->max_nb, T->nb[a][0]);
	    }
	return TRACE;

    }	
Ca_trace * hasch_ca_trace_bubble ( Ca_trace *TRACE)
    {

	int a, b;
	float dist;
	Atom *A, *B;
	float **list;
	Pdb_param *PP;
	Struct_nb *T;

	PP=TRACE->pdb_param;
	TRACE->Bubble=vcalloc ( 1, sizeof (Struct_nb));
	T=TRACE->Bubble;
	
		
	
	if ( !T->nb)T->nb=declare_int (TRACE->len+1, 1);
	if ( !T->d_nb)T->d_nb=declare_float (TRACE->len+1, 1);
	list=declare_float ( TRACE->n_atom, 3);
		
	
	for (a=0; a< TRACE->len; a++)
	    {
		for ( b=0; b< TRACE->len; b++)
		    {
			A=TRACE->ca[a];
			B=TRACE->ca[b];
			if ( !A || !B)continue;
			dist=get_atomic_distance ( A, B);
			
			if ( dist<PP->maximum_distance && FABS((A->res_num-B->res_num))>2)
			   {
			     T->nb[a][0]++;			       			       
			     T->nb[a]=vrealloc ( T->nb[a], (T->nb[a][0]+1)*sizeof (int));
			     T->nb[a][T->nb[a][0]]=(TRACE->ca[b])->num;
			     
			     T->d_nb[a]=vrealloc ( T->d_nb[a], (T->nb[a][0]+1)*sizeof (float));
			     T->d_nb[a][T->nb[a][0]]= ((a<b)?dist:-dist);
			   }						
		    }
		T->max_nb=MAX (T->max_nb, T->nb[a][0]);
	
	    }
	
	for ( a=0; a< TRACE->len; a++)
	    {
	    for ( b=0; b< T->nb[a][0]; b++)
	        {
		    list[b][0]=T->nb[a][b+1];
		    list[b][1]=T->d_nb[a][b+1];
		    list[b][2]=(TRACE->structure[T->nb[a][b+1]])->res_num;
		}
	    
	    sort_float ( list, 3,2, 0, T->nb[a][0]-1);
	    for ( b=0; b< T->nb[a][0]; b++)
	        {
		    T->nb[a][b+1]=list[b][0];
		    T->d_nb[a][b+1]=list[b][1];
		}
	    }

	free_float ( list, -1);
	return TRACE;

    }

       
float ** measure_ca_distances(Ca_trace *T)
   {
       int a, b;
       Atom *A, *B;
       float **dist;
       
       dist=declare_float ( T->len, T->len);
       
       for (a=0; a< T->len-1; a++)
	    {
		for ( b=a+1; b< T->len; b++)
		    {
			A=T->ca[a];
			B=T->ca[b];
			dist[a][b]=dist[b][a]=get_atomic_distance ( A,      B);
		    }
	    }
       return dist;
   }
float    get_atomic_distance ( Atom *A, Atom*B)
   {
       float dx, dy, dz;
       
       if ( !A || !B)return UNDEFINED;

       
       dx=A->x - B->x;
       dy=A->y - B->y;
       dz=A->z - B->z;
       
       return (float) sqrt ( (double) ( dx*dx +dy*dy +dz*dz));
   }
