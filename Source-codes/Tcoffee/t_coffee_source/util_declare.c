#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

/************************************************************************/
/*                                                                      */
/*            CONSTRAINT_LIST                                           */
/*                                                                      */
/*                                                                      */
/************************************************************************/ 

Constraint_list * declare_constraint_list ( Sequence *S, char *name, int **L, int ne,FILE *fp, int **M)
    {
    Constraint_list *CL;

    CL=vcalloc (1, sizeof ( Constraint_list));
  
    
    CL->S=S;
    CL->M=M;
    CL->list_name=vcalloc ( STRING, sizeof (char));
    if ( name!=NULL)
	{
	sprintf ( CL->list_name, "%s", name);
	
	}
    CL->cpu=1;
    CL->fp=fp;
    CL->L=L;
    CL->ne=ne;
    CL->entry_len=LIST_N_FIELDS;
    CL->el_size=sizeof (CLIST_TYPE);
    CL->matrices_list=declare_char(20,20);
    
    CL->chunk=1000;
    CL->weight_field=WE;
    CL->seq_for_quadruplet=vcalloc ( S->nseq, sizeof (int));
    
    return CL;
    }
Constraint_list *cache_dp_value4constraint_list ( char mode[],Constraint_list *CL)
    {
	static char dp_mode[100];
	static int gop;
	static int f_gop;
	static int f_gep;
	static int gep;
	static int TG_MODE;
	static int F_TG_MODE;
	static int maximise;
	static char matrix_for_aa_group[100];
	static int diagonal_threshold;
	static int ktup;
	static int fasta_step;
	static int use_fragments;
	static int extend_jit;
	static int (*get_dp_cost)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *);
        static int (*evaluate_residue_pair)(struct Constraint_list *, int, int, int, int);


	if ( strm (mode, "cache"))
	   {
	   sprintf (dp_mode, "%s", CL->dp_mode);
	   gop=CL->gop;
	   f_gop=CL->f_gop;
	   gep=CL->gep;
	   f_gep=CL->f_gep;
	   TG_MODE=CL->TG_MODE;
	   F_TG_MODE=CL->F_TG_MODE;
	   maximise=CL->maximise;
	   sprintf (matrix_for_aa_group, "%s", CL->matrix_for_aa_group);
	   diagonal_threshold=CL->diagonal_threshold;
	   ktup=CL->ktup;
	   fasta_step=CL->fasta_step;
	   use_fragments=CL->use_fragments;
	   extend_jit=CL->extend_jit;
	   get_dp_cost=CL->get_dp_cost;
	   evaluate_residue_pair=CL->evaluate_residue_pair;
	   return CL;
	   }
	else if ( strm (mode, "restaure"))
	   {
	   sprintf ( CL->dp_mode, "%s",dp_mode);
	   CL->gop=gop;
	   CL->f_gop=f_gop;
	   CL->gep=gep;
	   CL->f_gep=f_gep;
	   
	   CL->TG_MODE=TG_MODE;
	   CL->F_TG_MODE=F_TG_MODE;
	   CL->maximise=maximise;
	   sprintf ( CL->matrix_for_aa_group, "%s", matrix_for_aa_group);
	   CL->diagonal_threshold=diagonal_threshold;
	   CL->fasta_step=fasta_step;
	   CL->ktup=ktup;
	   CL->use_fragments=use_fragments;
	   CL->extend_jit=extend_jit;
	   CL->get_dp_cost=get_dp_cost;
	   CL->evaluate_residue_pair=evaluate_residue_pair;
	   return CL;
	   
	   }
	else
	   {
	       fprintf (stderr, "\nUNKNOWN MODE %s for cache_dp_value4comstraint_list[FATAL]", mode);
	       crash ( "");
	       return 0;
	   }
    }	       
Constraint_list *duplicate_constraint_list (Constraint_list *CL)
    {   
    Constraint_list *NCL;
    Sequence *S;
    int a, b;
    
    NCL=vcalloc ( 1, sizeof ( Constraint_list*));
    /*Sequences*/
      S=duplicate_sequence (CL->S);
      NCL=declare_constraint_list (S, NULL, NULL,0, NULL, NULL);
      NCL->DO_S=duplicate_sequence (CL->DO_S);
      NCL->out_aln_format=CL->out_aln_format;
      NCL->n_out_aln_format=CL->n_out_aln_format;
      NCL->W=duplicate_weights (CL->W);
    /*DATA*/
      if (CL->fp)NCL->fp=vtmpfile();
    
    /*copy of L*/
      for ( a=0; a< CL->ne; a++)
	for ( b=0; b< CL->entry_len; b++)
	    {
	    vwrite_clist(NCL, a, b, vread_clist(CL, a, b));
	    }
      if (CL->M)NCL->M=copy_int ( CL->M,NCL->M,-1, -1);
    
    /*List Information*/   
      NCL->ne=CL->ne;
      if ( CL->list_name)sprintf ( NCL->list_name, "%s", CL->list_name);
      NCL->entry_len=CL->entry_len;
      NCL->el_size=CL->el_size;

    /*Normalisation information*/
      NCL->normalise=CL->normalise;
      NCL->max_ext_value=CL->max_ext_value;
      NCL->max_value=CL->max_value;
    
    /*Pair wise alignment method*/
      NCL->pw_parameters_set=CL->pw_parameters_set;
      NCL->gop=CL->gop;
      NCL->f_gop=CL->f_gop;
      NCL->gep=CL->gep;
      NCL->f_gep=CL->f_gep;
      
      NCL->nomatch=CL->nomatch;
      NCL->TG_MODE=CL->TG_MODE;
      NCL->F_TG_MODE=CL->F_TG_MODE;
      sprintf ( NCL->dp_mode, "%s", CL->dp_mode);
      NCL->maximise=CL->maximise;
      sprintf ( NCL->matrix_for_aa_group, "%s", CL->matrix_for_aa_group);
      NCL->diagonal_threshold=CL->diagonal_threshold;
      NCL->ktup=CL->ktup;
      NCL->use_fragments=CL->use_fragments;
      NCL->fasta_step=CL->fasta_step;
      NCL->lalign_n_top=CL->lalign_n_top;
      NCL->sw_min_dist=CL->sw_min_dist;
      NCL->matrices_list=duplicate_char (CL->matrices_list, -1, -1);
      NCL->n_matrices=CL->n_matrices;
    /*Functions used for dynamic programming and Evaluation*/ 
      NCL->get_dp_cost=CL->get_dp_cost;
      NCL->evaluate_residue_pair=CL->evaluate_residue_pair;
      NCL->pair_wise=CL->pair_wise;

      NCL->weight_field=CL->weight_field;

    /*Parameters for domain extraction*/
      if ( CL->moca)NCL->moca=duplicate_moca ( CL->moca);
   
      

    /*Functions for hiding forbiden pairs of residues*/
      /* Not To be copied yet */
    
   /*extention properties:*/
      NCL->nseq_for_quadruplet=CL->nseq_for_quadruplet;
      NCL->seq_for_quadruplet=vcalloc ( S->nseq, sizeof(int));
      for ( a=0; a< S->nseq; a++)NCL->seq_for_quadruplet[a]=CL->seq_for_quadruplet[a];
      
   /*extention properties: Do Not copy*/
      /* Not To be copied yet */

    /*Lookup table parameteres*/
     /* Not To be copied yet */

    /*PDB STRUCTURE ALIGNMENTS*/
      /* Not To be copied yet */
    
    /*MISC*/  
       NCL->cpu=CL->cpu;
       NCL->local_stderr=CL->local_stderr;
       
    
    return NCL;
    }
Sequence *free_constraint_list (Constraint_list *CL)
    {

    
    Sequence *S;
    int a, b;

    if  (CL->list_name)vfree(CL->list_name);
    S=CL->S;
    if ( CL->L)free_int(CL->L,-1);
    
    if ( CL->fp)vfclose (CL->fp);	 
    
    
    if ( CL->M)free_int(CL->M, -1);
    if ( CL->matrices_list)free_char (CL->matrices_list, -1);
    if ( CL->start_index)free_int ( CL->start_index,-1);
    if ( CL->end_index)free_int ( CL->end_index,-1);
    if ( CL->residue_index)
      {
	for ( a=0; a< (CL->S)->nseq; a++)
	  {
	    for ( b=0; b<=(CL->S)->len[a]; b++)
	      free(CL->residue_index[a][b]);
	    free (CL->residue_index[a]);
	  }
	free(CL->residue_index);
      }
    
    if ( CL->DO_S && CL->DO_S!=CL->S)free_sequence ( CL->DO_S, (CL->DO_S)->nseq);
    if ( CL->W)free_weights (CL->W);
    
    if ( CL->translation)vfree(CL->translation);
    if ( CL->moca)free_moca (CL->moca);
    vfree (CL->seq_for_quadruplet);
    
    vfree(CL);
    return S;
    }

/************************************************************************/
/*                                                                      */
/*            MOCA Functions                                            */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Moca * duplicate_moca ( Moca *m)
      {
	Moca *nm;

	nm=vcalloc ( 1, sizeof (Moca));

	nm->moca_scale=m->moca_scale;	
	nm->evaluate_domain=m->evaluate_domain;
	nm->moca_threshold=m->moca_threshold;
	nm->cache_cl_with_domain=m->cache_cl_with_domain;
	if ( m->forbiden_residues)nm->forbiden_residues=copy_int  (m->forbiden_residues,nm->forbiden_residues,    -1, -1);
	nm->make_nol_aln=m->make_nol_aln;

	
	return nm;
      }
Moca * free_moca ( Moca *m)
      {
	if ( m->forbiden_residues)free_int ( m->forbiden_residues, -1);
	vfree ( m);
	return NULL;
      }

/************************************************************************/
/*                                                                      */
/*            PDB Functions                                             */
/*                                                                      */
/*                                                                      */
/************************************************************************/	     
Structure* declare_structure ( int n, char **array)
    {
    Structure *S;
    int a;

    S=vcalloc (1, sizeof (Structure));
    S->n_fields=1;
    S->nseq=n;
    
    S->struc=vcalloc ( n, sizeof (int**));
    S->len=vcalloc ( n, sizeof (int));
    for ( a=0; a< n; a++)
        {
	S->len[a]=strlen(array[a]);
	S->struc[a]=declare_int ( strlen ( array[a])+2, 1);	
	}
    return S;
    }

Structure *extend_structure ( Structure *S)
    {
    int a, b;


    for ( a=0; a< S->nseq; a++)
        {
	for ( b=0; b< S->len[a]; b++)
		S->struc[a][b]=realloc ( S->struc[a][b],( S->n_fields+1)*sizeof (int));		   
	}
    S->n_fields++;
    return S;
    }

Sequence * declare_sequence ( int min, int max, int nseq)
    {
    Sequence *LS;
    
    

    LS=vcalloc (1, sizeof ( Sequence));
    LS->comment=declare_char ( nseq,LONG_STRING+1);
    
    LS->file=declare_char( nseq,STRING+1);
    LS->seq=declare_char ( nseq, max+1);
    LS->name=declare_char( nseq,STRING+1);
    LS->len=vcalloc ( nseq, sizeof (int));
    LS->max_len=max;
    LS->min_len=min;
    LS->nseq=nseq;
    LS->max_nseq=nseq;
    LS->type=vcalloc(30, sizeof (char));
    LS->Profile=vcalloc (nseq, sizeof (Alignment*));
    return LS;
    }
Sequence * realloc_sequence   (Sequence *OUT, int new_nseq, int max_len)
    {
      int a;

      if ( new_nseq<OUT->max_nseq)return OUT;
      else 
      OUT->min_len=MIN(OUT->min_len,max_len);
      OUT->max_len=MAX(OUT->max_len,max_len);
      OUT->comment=new_realloc_char ( OUT->comment, new_nseq,LONG_STRING+1);
      OUT->seq    =new_realloc_char ( OUT->seq,     new_nseq,OUT->max_len+1);
      OUT->name   =new_realloc_char ( OUT->name,    new_nseq,STRING+1);
      OUT->file   =new_realloc_char ( OUT->file,    new_nseq,STRING+1);
      OUT->len    =vrealloc     ( OUT->len,     (new_nseq+1)*sizeof (int));
      OUT->Profile=vrealloc     ( OUT->Profile,     new_nseq*sizeof (Alignment *));
      for (a=OUT->max_nseq; a<new_nseq; a++)OUT->Profile[a]=NULL;
      OUT->max_nseq=new_nseq;
      return OUT;
    }
	
Sequence * duplicate_sequence (Sequence *S )
    {
    Sequence *LS;
    int a;
    
    LS=declare_sequence (S->min_len, S->max_len, S->nseq);
    
    for ( a=0; a<S->nseq; a++)
        {
	sprintf ( LS->file[a], "%s", S->file[a]);
	sprintf ( LS->comment[a], "%s", S->comment[a]);
	sprintf ( LS->seq[a], "%s", S->seq[a]);
	sprintf ( LS->name[a], "%s", S->name[a]);
	LS->len[a]=S->len[a];
	}
    LS->max_len=S->max_len;
    LS->min_len=S->min_len;
    LS->nseq=S->nseq;
    LS->Profile=vcalloc ( S->nseq, sizeof (Alignment*));
    for ( a=0; a< S->nseq; a++)
      {
	if ( S->Profile[a])LS->Profile[a]=copy_aln (S->Profile[a], NULL);
      }
    LS->max_nseq=S->nseq;
    return LS;
    }
void free_sequence ( Sequence *LS, int nseq)
	{
	int a;
	free_char ( LS->file, -1);
	free_char ( LS->comment, -1);
	free_char ( LS->seq, -1);
	free_char ( LS->name,-1);
	vfree (LS->type);
	vfree (LS->len);
	for ( a=0; a< LS->nseq; a++)free_aln(LS->Profile[a]);
	free (LS->Profile);
	free (LS);
	}
/************************************************************************/
/*                                                                      */
/*            Weights Functions                                         */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Weights* declare_weights ( int nseq)
	{
	Weights *W;
	
	W=calloc ( 1, sizeof ( Weights));
	W->comments=calloc ( 1000, sizeof (char));
	W->nseq=nseq;
	W->mode=vcalloc (FILENAMELEN, sizeof (char));
	W->seq_name= declare_char ( W->nseq*2, 200);
	W->PW_SD=declare_float ( W->nseq, W->nseq);
	W->PW_ID=declare_float ( W->nseq, W->nseq);
	W->SEQ_W=calloc ( W->nseq, sizeof ( float));
	return W;
	}  
Weights* duplicate_weights (Weights *W)
     {
       Weights *NW;
       int a, b;

       NW=declare_weights (W->nseq);
       sprintf ( NW->comments, "%s", W->comments);
       sprintf ( NW->mode, "%s", W->mode);
       for (a=0; a< W->nseq; a++)
	 {
	   sprintf ( NW->seq_name[a], "%s", W->seq_name[a]);
	   NW->SEQ_W[a]=W->SEQ_W[a];
	   for(b=0; b< W->nseq; b++)
	     {
	       	NW->PW_SD[a][b]=W->PW_SD[a][b];
		NW->PW_ID[a][b]=W->PW_ID[a][b];
	     }
	 }
       return NW;
     }
Weights* free_weights ( Weights* W)
	{

	
	
	free(W->comments);
	
	
	free(W->mode);
	
	free_char(W->seq_name, -1);	
	free_float(W->PW_SD,-1);	
	free_float(W->PW_ID, -1);	
	vfree(W->SEQ_W);
	return NULL;
	} 
Alignment* copy_aln ( Alignment *A, Alignment *B)
        {
	  int a, b;
	  if ( A==NULL){free_aln(B); return NULL;}
	    
	    free_aln (B);
	    if ( A->S && A->nseq>(A->S)->nseq)
	      {
		B=declare_aln2(A->nseq, MAX((A->S)->max_len+1, A->len_aln+1));
		B->S=A->S;
	      }
	    else 
	      B=declare_aln ((A->S));
	    
	   
	    
/*SIZES*/    
	    B->max_len=A->max_len;
	    B->min_len=A->min_len;
	    B->declared_len=A->declared_len;
	    B->max_n_seq=A->max_n_seq;

	    B->nseq=A->nseq;	
	    B->len_aln=A->len_aln;
	   
	   
/*sequence Information*/	

	    if ( (A->S)==NULL){vfree (B->len); B->len=vcalloc ( A->max_n_seq, sizeof (int));}
	    ga_memcpy_int ( A->len, B->len, B->nseq);
	    
	    B->comment=copy_char ( A->comment,  B->comment,  -1,-1);
	    B->name=copy_char ( A->name,     B->name,     -1,-1);
	    B->file=copy_char ( A->file,     B->file,     -1,-1);
	    B->tree_order=copy_char ( A->tree_order,     B->tree_order,     -1,-1);
	    free_char ( B->seq_al, -1);
	    B->seq_al=declare_char(B->max_n_seq, B->declared_len);
	    for ( a=0; a< A->max_n_seq; a++)
	      {
		for ( b=0; b< A->declared_len; b++)
		  B->seq_al[a][b]=A->seq_al[a][b];
	      }
	    
	    
	    
	    B->order=copy_int  ( A->order,    B->order,    -1, -1);
	    B->S=A->S;
	    if (A->seq_cache)
	        {
		B->seq_cache=copy_int  ( A->seq_cache,    B->seq_cache,-1,-1);
		}
	    
	    if (A->cdna_cache)
	        {
		B->cdna_cache=copy_int  ( A->cdna_cache,    B->cdna_cache,-1,-1);
		}
	    B->S=A->S;
	    B->P=A->P;
	    B->Dp_result=A->Dp_result;
	    
/*Score*/

	    if ( (A->S)==NULL){vfree (B->score_seq); B->score_seq=vcalloc ( A->max_n_seq, sizeof (int));}
	    ga_memcpy_int(  A->score_seq,B->score_seq,B->nseq);
	    
	    B->pw_score_seq=copy_int  ( A->pw_score_seq, B->pw_score_seq,-1,-1);
	    B->score_aln=A->score_aln;
	    B->score=A->score;
	    B->cpu=A->cpu;
	    B->finished=A->finished;

/*Output Options*/    
	    B->output_res_num=A->output_res_num;
	    B->residue_case=A->residue_case;
	    
	    B->CL=A->CL;
	    
	return B;
	}
Alignment* shrink_aln ( Alignment *A, int nseq, int *list)
        {
	Alignment *B=NULL;
	int a,seq;
	
	B=copy_aln (A, B);
	for ( a=0; a< nseq; a++)
	    {
	    seq=list[a];
	    sprintf ( A->comment[a], "%s",B->comment[seq]);
	    sprintf ( A->seq_al [a], "%s",B->seq_al [seq]);
	    A->order[a][0]=B->order[seq][0];
	    A->order[a][1]=B->order[seq][1];
	    A->score_seq[a]=B->score_seq[seq];
	    A->len[a]=B->len[seq];
	    }
	A->nseq=nseq;
	A->len_aln=strlen (A->seq_al[0]);
	free_aln (B);
	return A;
	}
Alignment* extract_sub_aln ( Alignment *B, int nseq, int *list)
        {
	Alignment *A=NULL;
	int a,b,seq;
	
	A=declare_aln2(nseq, B->len_aln+1);
	for ( a=0; a< nseq; a++)
	    {
	    seq=list[a];
	    sprintf ( A->comment[a], "%s",B->comment[seq]);
	    sprintf ( A->name[a], "%s",B->name[seq]);
	    for (b=0; b<=B->len_aln; b++)A->seq_al [a][b]=B->seq_al [seq][b];
	    A->order[a][0]=B->order[seq][0];
	    A->order[a][1]=B->order[seq][1];
	    A->score_seq[a]=B->score_seq[seq];
	    A->len[a]=B->len[seq];
	    }
	A->nseq=nseq;
	A->len_aln=B->len_aln;
	return A;
	}
	    
Alignment *declare_aln2 ( int nseq, int len)
        {
	  Sequence *S;
	  Alignment *A;
	  
	  S=vcalloc ( 1, sizeof ( Sequence));
	  S->nseq=nseq;
	  S->max_len=len;
	  
	  A=declare_aln (S);
	  A->S=NULL;
	  free(S);
	  return A;
	}
  


Alignment *declare_aln ( Sequence *S){return declare_Alignment(S);}

Alignment *declare_Alignment ( Sequence *S)
	{
	Alignment *LA;
	int a;
	
	/*ordre:
	  [x][0]= which is the xth seq of aln
	  [x][1]= how many deleted residues before the first one
	*/

	LA=calloc (1, sizeof ( Alignment));
	if ( S==NULL)
	    {
	    LA->declared_len=MAX_LEN_ALN;
	    LA->max_n_seq=MAX_N_SEQ;
	    }
	else
	    {
	    LA->declared_len=2*S->max_len+1;
	    LA->max_n_seq=S->nseq+1;
	    }
	

	LA->comment=declare_char (LA->max_n_seq, LONG_STRING);
	LA->seq_al=declare_char ( LA->max_n_seq,LA->declared_len );
	LA->name=declare_char (LA->max_n_seq, STRING);
	LA->file=declare_char (LA->max_n_seq, STRING);
	LA->tree_order=declare_char (LA->max_n_seq, STRING);
	LA->order= declare_int (LA->max_n_seq , 2);
	LA->score_seq= vcalloc (LA->max_n_seq, sizeof (int));
	LA->pw_score_seq= declare_int (LA->max_n_seq,LA->max_n_seq);

	for ( a=0; a< LA->max_n_seq; a++)LA->order[a][0]=a;
	
	LA->len_aln=0;
	LA->score_aln=0;
	LA->len=vcalloc (LA->max_n_seq, sizeof (int));		
	
	LA->S=S;	
	if (S && S->name)for ( a=0; a<S->nseq; a++)sprintf ( LA->name[a], "%s", S->name[a]);
	return LA;
		
	}
Alignment * realloc_aln ( Alignment *A, int new_len){return realloc_alignment(A, new_len);}
Alignment * realloc_alignment ( Alignment *A, int new_len)
	{
	if (A==NULL)A=declare_Alignment (NULL);

	return realloc_alignment2( A, A->max_n_seq,new_len);

	}

Alignment * realloc_aln2 ( Alignment *A, int n_nseq, int n_len){return realloc_alignment2(A, n_nseq, n_len);}



Alignment * realloc_alignment2 ( Alignment *A, int n_nseq, int n_len)
	{
	int a;
	int len, nseq;
	int delta_len, delta_nseq;
	
        if ( A==NULL) A=declare_Alignment(NULL);
	  
	n_len++;
	n_nseq++;

	len=A->declared_len;
	nseq=A->max_n_seq;
	
	n_len=MAX(len, n_len);
	n_nseq=MAX(nseq,n_nseq);
	delta_nseq=MAX(0,n_nseq-nseq);
	delta_len =MAX(0,n_len-len);

	if ( delta_nseq<=0 && delta_len<=0)return A;
	
	
	else
	    {	    
	        A->len          =vrealloc( A->len      , sizeof (int)*n_nseq);
		for (a=nseq; a< n_nseq; a++)A->len[a]=0;
		A->declared_len =n_len;
		A->max_n_seq    =n_nseq;

		
		A->comment=new_realloc_char ( A->comment, n_nseq, -1);
		A->name   =new_realloc_char ( A->name, n_nseq, -1);
		A->file   =new_realloc_char ( A->file, n_nseq, -1);
		A->tree_order   =new_realloc_char ( A->tree_order, n_nseq, -1);
		A->seq_al =new_realloc_char ( A->seq_al, n_nseq, n_len);
		A->order  =new_realloc_int  ( A->order, n_nseq, -1);
		if ( A->seq_cache) A->seq_cache=new_realloc_int  ( A->seq_cache, n_nseq,n_len);
		if ( A->cdna_cache)A->cdna_cache=new_realloc_int  ( A->cdna_cache, n_nseq,n_len);
	

		A->score_seq    =vrealloc( A->score_seq, sizeof (int)*(n_nseq));
		for ( a=nseq; a< n_nseq; a++)A->score_seq[a]=0;
		A->pw_score_seq =new_realloc_int  ( A->pw_score_seq,  n_nseq,n_nseq );
			   
	    }	  
	return A;
	}

Alignment* free_aln ( Alignment *LA){ return free_Alignment(LA);}
Alignment* free_Alignment ( Alignment *LA)
	{

	if ( LA==NULL){return NULL;}


	free_char ( LA->file, -1);
	free_char ( LA->seq_al, -1);
	free_int  ( LA->seq_cache, -1);
	free_int  ( LA->cdna_cache, -1);
	free_char ( LA->name,-1);	
	free_char ( LA->tree_order,-1);	
	free_char ( LA->comment, -1);  
	free_int  ( LA->order, -1);
	free_int  ( LA->pw_score_seq, -1);
	vfree ( LA->score_seq);
	vfree ( LA->len);
        vfree ( LA);
	return NULL;
	}

Profile   *declare_profile(char *alphabet, int len)
{
  Profile *P;
  P=vcalloc ( 1, sizeof ( Profile));
  P->alp_size=strlen(alphabet);
  P->max_len=len;
  P->alphabet=vcalloc ( strlen (alphabet)+1, sizeof (char));
  sprintf ( P->alphabet, "%s", alphabet);
  
  P->count=declare_int( P->alp_size, len);
  return P;
}
Profile * free_profile ( Profile *P)
{
  vfree (P->alphabet);
  free_int ( P->count, -1);
  vfree (P);
  return NULL;
}


/************************************************************************/
/*                                                                      */
/*             ALLOCATION                                               */
/*                                                                      */
/*                                                                      */
/************************************************************************/   
void * vmalloc ( size_t size)
	{
	void * x;
	
	
	if ( size==0)
	    return NULL; /*crash ("\n0 bytes in vmalloc\n");*/
	else
	    {
	    x= malloc (size);
	    if ( x==NULL)
		{
		crash ( "\nFAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)\n");
		return NULL;
		}
	    else
	        {
		
		return x;
		}
	    }
	}
void *vcalloc ( size_t nobj, size_t size)
	{
	void *x;

	if ( nobj<=0 || size<=0)return NULL;/*crash ("\n0 bytes in vmalloc\n");*/
	else x= calloc ( nobj, size);
	if ( x==NULL)
		{
		crash ( "\nFAILED TO ALLOCATE REQUIRED MEMORY (vcalloc)\n");
		return NULL;
		}
	else
	        {	
		return x;
		}
	}

void *vrealloc ( void *p, size_t size)
	{
	void *x;
	if ( size<=0)return p;
	else if ( p==NULL)
	    {
	    x=vmalloc (size);
	    memset(x, 0, size);
	    }
	else
		x=realloc ( p, size);
	if ( x==NULL)
		{
		crash ( "\nFAILED TO ALLOCATE REQUIRED MEMORY (realloc)\n");		
		return 0;
		}
	else
		return x;
	}
	
void vfree ( void *p)
     {
     if (p)free(p);
     }
/************************************************************************/
/*                                                                      */
/*             DECLARE 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/   


#define DECLARE_ARRAY(type,wf,rf,function)\
type**  function (int first, int second)\
  {type **array;\
   int a;\
\
   if ( first<0)return NULL;\
   array=vcalloc (first+1, sizeof (type*));\
   array++;\
   array[-1]=vcalloc ( 2*SIZE_OF_INT, sizeof(type));\
   wf(first,array[-1],0);\
   wf(MAX(second,0),array[-1],1);\
   for ( a=0; a< first; a++)\
       array[a]=vcalloc ( MAX(0,(second)), sizeof (type));\
   return array;\
   }

DECLARE_ARRAY(short,write_size_short,read_size_short,declare_short)
DECLARE_ARRAY(char,write_size_char,read_size_char,declare_char)
DECLARE_ARRAY(int,write_size_int,read_size_int,declare_int)
DECLARE_ARRAY(float,write_size_float,read_size_float,declare_float)
DECLARE_ARRAY(double,write_size_double,read_size_double,declare_double)


Alignment ** declare_aln_array ( int first)
    {
    Alignment ** array;
    int a;
    
    array=vcalloc (first+1, sizeof (Alignment*));
    array++;
    array[-1]=vcalloc ( 1, sizeof(Alignment));
    (array[-1])->nseq=first;
    
    for ( a=0; a< first; a++)
	array[a]=declare_Alignment (NULL);
    
    return array;
    }
/************************************************************************/
/*                                                                      */
/*             Realloc 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/   

#define REALLOC_ARRAY(type,wf,rf,function1,function2,function3) \
type ** function1 ( type **array, int first, int second, int ext1, int ext2)\
    {int a, b;\
     int delta1;\
     int delta2;\
\
    if (first==-1 && second==-1)\
       {\
       first=(array==NULL)?0:rf(array[-1],0);\
       second=(array==NULL)?0:rf(array[-1],1);\
       if ( ext1==-1)ext1=first;\
       if ( ext2==-1)ext2=second;\
       ext1=ext1-first;\
       ext2=ext2-second;\
       }\
     if (GIVE_MEMORY_BACK==1)\
	{\
	ext1=MAX(0,ext1);\
	ext2=MAX(0,ext2);\
	}\
\
    if ( array!=NULL){array--;first++;}\
    delta1=first +ext1;\
    delta2=second+ext2;\
\
    if ( delta1<=0)\
       {\
       function3(array, -1);\
       return NULL;\
       }\
    else if ( ext1==0 && ext2==0);\
    else if ( array==NULL)return function2(delta1, delta2);\
    else if ( ext1>=0 && ext2>=0)\
	{\
	if ( ext1>0)array=vrealloc ( array, (sizeof (type*))*delta1);\
	for ( a=1; a<first; a++)\
	        {\
		if ( ext2>0)array[a]=vrealloc(array[a],(sizeof(type))*(delta2));\
                for ( b=second; b< delta2; b++)\
                         array[a][b]=0;\
                }\
	for ( a=first; a< (delta1); a++)\
	    array[a]=vcalloc( (delta2), sizeof (type));\
	}\
    else if (ext1<=0 && ext2<=0)\
        {\
	if ( ext1<0)for (a=first-1; a>=(delta1); a--)vfree ( array[a]);\
	if (ext1<0)array=vrealloc ( array,(sizeof (type*))*(delta1));\
	for (a=1; a<(delta1);a++)\
		if (ext2<0)array[a]=vrealloc ( array[a],(sizeof(type))*(delta2));\
	}\
    else if ( ext1<=0 && ext2>=0)\
        {\
	if ( ext1<0)for (a=first-1; a>=(delta1); a++)vfree ( array[a]);\
	if ( ext1<0)array=vrealloc ( array,(sizeof (type*))*(delta1));\
	for (a=1; a<(delta1);a++)\
	    {\
	    if (ext2>0)array[a]=vrealloc ( array[a],(sizeof(type))*(delta2));\
	    for ( b=second; b<delta2; b++)array[a][b]=0;\
            }\
	}\
    else if ( ext1>=0 && ext2<=0)\
        {\
	if (ext1>0)array=vrealloc ( array, (sizeof (type*))*(delta1));\
        for ( a=1; a<first; a++)\
	    if(ext2<0)array[a]=vrealloc(array[a],(sizeof(type))*(delta2));\
	for ( a=first; a< (delta1); a++)\
	    array[a]=vcalloc( (delta2), sizeof (type));\
	}\
    array++;\
    wf(delta1-1,array[-1],0);\
    wf(delta2, array[-1],1);\
    return array;\
    }
REALLOC_ARRAY(short,write_size_short,read_size_short,realloc_short,declare_short,free_short)
REALLOC_ARRAY(char,write_size_char,read_size_char,realloc_char,declare_char,free_char)
REALLOC_ARRAY(int,write_size_int,read_size_int,realloc_int,declare_int,free_int)
REALLOC_ARRAY(float,write_size_float,read_size_float,realloc_float,declare_float,free_float)
REALLOC_ARRAY(double,write_size_double,read_size_double,realloc_double,declare_double,free_double)

#define NEW_REALLOC_ARRAY(type,wf,rf,function1,function2,function3)\
type ** function1 ( type **array, int ext1, int ext2)\
    {int a, b;\
     int first, l1;\
     int second, l2;\
     type **new_array;\
\
     first=rf(array[-1],0);\
     second=rf(array[-1],1);\
     \
     if ( ext1==-1)ext1=first;\
     if ( ext2==-1)ext2=second;\
     l1=MIN(ext1, first);\
     l2=MIN(ext2, second);\
     new_array=function2(ext1, ext2);\
     for ( a=0; a<l1; a++)\
       for ( b=0; b<l2; b++)new_array[a][b]=array[a][b];\
     function3(array, -1);\
     return new_array;\
    }
NEW_REALLOC_ARRAY(short,write_size_short,read_size_short,new_realloc_short,declare_short,free_short)
NEW_REALLOC_ARRAY(char,write_size_char,read_size_char,new_realloc_char,declare_char,free_char)
NEW_REALLOC_ARRAY(int,write_size_int,read_size_int,new_realloc_int,declare_int,free_int)
NEW_REALLOC_ARRAY(float,write_size_float,read_size_float,new_realloc_float,declare_float,free_float)
NEW_REALLOC_ARRAY(double,write_size_double,read_size_double,new_realloc_double,declare_double,free_double)
Alignment ** realloc_aln_array ( Alignment **array, int ext1)
    {
    int a;
    int first;
    
    
     if ( array==NULL)
	 {	 
	 array=declare_aln_array(ext1);	
	 return array;
	 }
    first=(array[-1])->nseq+1;
    array--;
    
    if ( ext1>0)
	{
	array=vrealloc ( array, (sizeof (Alignment*))*(first+ext1));
        for ( a=first; a<first+ext1; a++)array[a]=declare_Alignment (NULL);
	}
    else if ( ext1==0);
    else if ( ext1<0)
         {
	 for ( a=first-1; a>=(first+ext1);a--)free_Alignment (array[a]);
	 array=vrealloc ( array, (sizeof (Alignment*))*(first+ext1));
	 }
    
   
    array++;
    (array[-1])->nseq+=ext1;
    
   
    if ( ext1<0)exit (1);
    return array;    
    }

/************************************************************************/
/*                                                                      */
/*            free 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/   
#define FREE_ARRAY(type,wf,rf,function) \
type ** function (type **array, int first)\
    {int a, len;\
    if ( array==NULL)return(type **) NULL; \
    else \
	{array--; \
	len=rf(array[0],0);\
	for ( a=0; a<=len;a++)vfree ( array[a]);\
	vfree (array);\
	}\
    return (type **)NULL;\
    }
FREE_ARRAY(short,write_size_short,read_size_short,free_short)
FREE_ARRAY(char,write_size_char,read_size_char,free_char)
FREE_ARRAY(int,write_size_int,read_size_int,free_int)
FREE_ARRAY(float,write_size_float,read_size_float,free_float)
FREE_ARRAY(double,write_size_double,read_size_double,free_double)

Alignment ** free_aln_array (Alignment **array)
   {
   int a;
   int len;


   if ( array==NULL)return NULL;
   array--;
   len=array[0]->nseq;
   free_Alignment (array[0]);

   for ( a=1; a< len; a++)free_Alignment(array[a]);
   vfree ( array);
   return  NULL;
   }


Fname *declare_fname ()
   {
   Fname *F;

   F=vcalloc ( 1, sizeof (Fname));
   F->name  =vcalloc ( FILENAMELEN, sizeof (char));
   F->path  =vcalloc ( FILENAMELEN, sizeof (char));
   F->suffix=vcalloc ( FILENAMELEN, sizeof (char));
   return F;
   }

Fname *free_fname ( Fname *F)
   {
   vfree (F->name);
   vfree (F->path);
   vfree (F->suffix);
   return NULL;
   }
