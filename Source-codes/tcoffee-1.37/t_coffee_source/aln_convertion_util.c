#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
 
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"
int ** minimise_repeat_coor (int **coor, int nseq, Sequence *S)
    {
    int **new_coor;
    int a, min;
    new_coor=declare_int ( nseq, 3);
    min=return_min_int (coor, nseq, 2);
    for ( a=0; a< nseq; a++)
        {
	new_coor[a][0]=coor[a][0];
	new_coor[a][1]=coor[a][1];
	new_coor[a][2]=min;
	}
    return new_coor;
    }
int ** get_nol_seq ( Constraint_list *CL, int **coor, int nseq, Sequence *S)
    {
    int a, s, p, l, nl;
    int **buf;
    int **new_coor;
    
    new_coor=declare_int ( nseq+1, 3);
    
    
    buf=get_undefined_list ( CL);
    
    

    for ( a=0; a< nseq; a++)buf[coor[a][0]][coor[a][1]]=1;

    
    for ( a=0; a< nseq; a++)
        {
	s=coor[a][0];
	p=coor[a][1]+1;
	l=strlen(S->seq[s]);
	nl=0;
	while ( p<=l && !buf[s][p++])nl++;
	new_coor[a][0]=s;
	new_coor[a][1]=coor[a][1];
	new_coor[a][2]=nl;
	}
    free_int ( buf, -1);
    return new_coor;
    }
int compare_pos_column( int **pos1,int p1, int **pos2,int p2, int nseq)
    {
    int a,v1, v2;
    int identical=0;


    
    for ( a=0; a< nseq; a++)
	{
	
	v1=pos1[a][p1];
	v2=pos2[a][p2];
	
	if (v1>0 || v2>0) 
	    {
	    if ( v1!=v2)return 0;
	    else identical=1;
	    }
	}
    
    return identical;
    }

  

int ** seq2aln_pos      (Alignment *A, int *ns, int **l_s)
    {
    int **code;
    int a, b,c, d,l, p , g;
    

    l=MAX(strlen (A->seq_al[l_s[0][0]]), strlen (A->seq_al[l_s[1][0]]));
    code=declare_int ((A->S)->nseq,l+1);
    
    for (c=0; c<2; c++)
        {
	l=strlen (A->seq_al[l_s[c][0]]);
	for (d=0; d<ns[c]; d++)
	    {
	    a=A->order[l_s[c][d]][0];
	    for (p=0, b=0; b<l; b++)
	        {
		    g=is_gap (A->seq_al[l_s[c][d]][b]);
		    if (!g){p++; code[a][p]=b+1;}
		}
	    }
	}
    return code;
    }

int ** aln2pos_simple_2 (Alignment *A)
    {
    int **pos1;
    int **pos2;
    pos1=aln2pos_simple (A, A->nseq);
    pos2=duplicate_int  (pos1, A->nseq,read_size_int (pos1[-1], 1));
    pos1=aln2pos_simple (NULL, 0);
    return pos2;
    }
int ** aln2pos_simple (Alignment *A, int n_nseq, ...)
    {
    /*
    function documentation: start
    int ** aln2pos_simple (Alignment *A, int n_nseq, ...)

####with two parameter only: Alignment *A, int n_nseq

    this function turns A into pos, a matrix where each residue is replace by its inice acvcording to the complete sequence.
    the indices in pos are computed using A->order[x][1] that contains the indice of the first residue of seq x of A
    
    n_nseq MUST not be null

####with more than two param:
     int ** aln2pos_simple (Alignment *A, int n_nseq, int *ns, int **ls)
     n_nseq must be set to 0 for the param 3 and four to be read
     
     ns[x]=number seq in group 
     ls[x]=list of the sequences in group x ( size=ns[x])
     
    The computation of the indices is only carried out on the scpecified residues

####IMPORTANT
      in pos, the numbering of the residues goes from 1 to L:
        pos[0][0]=3, means that the first position of the first sequence 
	in the alignmnet contains residue #3 from sequence A->order[0][0];
	
    function documentation: end
    */

    int a, b,c, p, g,l;
    static int **T;
    int nseq;
    int max_nseq;
    int len;
    int n_len=0;

    int *list=NULL;
    int *ns=NULL;
    int **ls=NULL;



    va_list ap;
    
    
    if ( A==NULL)
       {
       nseq=len=0;
       free_int (T, -1);
       T=NULL;   
       return T;
       }
    else
       {
       if ( n_nseq>0)
          {
	  list=vcalloc(n_nseq, sizeof (int));
	  for ( a=0; a< n_nseq; a++)list[a]=a;
	  }
       else
          { 
	  va_start (ap, n_nseq);
	  ns=va_arg(ap, int * );
	  ls=va_arg(ap, int **);
	  va_end(ap);
	  list=vcalloc ( ns[0]+ns[1], sizeof (int));
	  n_nseq=0;
	  for ( a=0; a< ns[0]; a++)list[n_nseq++]=ls[0][a];
	  for ( a=0; a< ns[1]; a++)list[n_nseq++]=ls[1][a];

	  }
       
    
       nseq=(T==NULL)?0:(read_size_int ( T[-1],0));
       len =(T==NULL)?0:(read_size_int ( T[-1],1));

       max_nseq=MAX(read_size_int(A->order[-1],0),return_max_int (A->order, read_size_int(A->order[-1],0),0))+1;
       n_len=get_longest_string ( A->seq_al,A->max_n_seq, NULL, NULL)+1;  

       T=realloc_int ( T,nseq, len,MAX(0,max_nseq-nseq), MAX(0,n_len-len));
       
       for ( c=0; c< n_nseq; c++)
           {
	   a=list[c];	       
	   l=strlen ( A->seq_al[a]);
	   
	   for ( p=A->order[a][1],b=0; b<l; b++)
	       {
	       g=1-is_gap(A->seq_al[a][b]);
	       p+=g;
	       T[a][b]=(g==1)?p:-(1+p);
	       if ( A->seq_al[a][b]==UNDEFINED_RESIDUE)T[a][b]=0;
	       if ( A->seq_cache && T[a][b]>0)T[a][b]=A->seq_cache[A->order[a][0]][T[a][b]];
	       } 
	   }
       free (list);
       }
   
   return T;
   }
Alignment ** split_seq_in_aln_list ( Alignment **aln, Sequence *S, int n_seq, char **seq_list)
        {
	int a, b, c;
	char * long_seq=NULL;
	int    len,l;
	int  **translation;
	int  **table;


	

	if ( aln==NULL)return NULL;
	translation=declare_int ( S->nseq,2);
	
	for (len=0,a=0; a< S->nseq; a++)
	    {
	    if((b=name_is_in_list (S->name[a],seq_list, n_seq, 100))!=-1)
	       {
	       l=strlen(S->seq[a])+1;
	       long_seq=vrealloc(long_seq,(len+l+1)*sizeof(char));
	       long_seq=strcat(long_seq, S->seq[a]);
	       long_seq=strcat(long_seq, "*");   
	       
	       translation[a][0]=b;
	       translation[a][1]=len;
	       len+=l;
	       }
	    else translation[a][0]=-1;
	    }

	long_seq[len-1]='\0';
	len--;

	table=declare_int ( len+1, 2);

	for ( b=0,a=0; a< S->nseq; a++)
	    {
	    if ( translation[a][0]!=-1)
	       {
	       c=1;
	       while (long_seq[b]!='\0' && long_seq[b]!='*')
		   {
		   table[b+1][1]=c++;
		   table[b+1][0]=translation[a][0];
		   b++;
		   }
	       table[b][1]=c++;
	       table[b][0]=translation[a][0];
	       b++;
	       }
	    }

	for ( a=0; a< (aln[-1])->nseq; a++)
	    {
	    for ( b=0; b< (aln[a])->nseq; b++)
	        {
		
		(aln[a])->order[b][0]=table[(aln[a])->order[b][1]][0];
		(aln[a])->order[b][1]=table[(aln[a])->order[b][1]][1];
		sprintf ( (aln[a])->name[b],"%s_%d_%d", S->name[(aln[a])->order[b][0]],a+1,b+1); 
		}
	    }
	free_int (translation, -1);
	free_int (table,       -1);
	return aln;
	}



Sequence  *  fill_sequence_struc ( int nseq, char **sequences, char **seq_name)
	{
	int a;
	Sequence *S;
	
	S=declare_sequence (get_shortest_string( sequences, nseq, NULL, NULL), get_longest_string (sequences, nseq, NULL, NULL),nseq);
	S->nseq=nseq;

	S->seq=copy_char ( sequences, S->seq, nseq, -1);
	S->name=copy_char ( seq_name, S->name,nseq, -1);
	
	ungap_array (S->seq,nseq);
	for ( a=0; a< S->nseq; a++)S->len[a]=strlen(S->seq[a]);
	return S;
	}

Alignment * expand_aln (Alignment *A)
  {
  /*This function expands the profiles within an alignment*/
  
  
  int a, b, d, e;
  Alignment *MAIN=NULL, *SUB=NULL;
  int n_sub_seq=0;
  int new_nseq=0;
  int *list;
  
  
  if ( !A)return A;


  
  list=vcalloc (A->nseq, sizeof (int)); 
  for ( a=0; a< A->nseq; a++)
    {
      if ( (A->S)->Profile && (A->S)->Profile[A->order[a][0]])new_nseq+=((A->S)->Profile[A->order[a][0]])->nseq;
      else 
	{
	  new_nseq++;
	  list[n_sub_seq++]=a;
	}      
    }
  
  if ( n_sub_seq==A->nseq){vfree(list);return A;}
  else if (n_sub_seq==0){MAIN=copy_aln (A, MAIN);MAIN->nseq=0;}
  else
    {
      MAIN=extract_sub_aln (A, n_sub_seq, list);
    }
  vfree(list);
  
  
  for ( a=0; a< A->nseq; a++)
    {
      if ( (A->S)->Profile && (A->S)->Profile[A->order[a][0]])
	{
	  SUB=copy_aln ( (A->S)->Profile[A->order[a][0]],SUB);
	  SUB=realloc_aln2(SUB, SUB->nseq, A->len_aln+1);
	  
	  for ( e=0,b=0; b<= A->len_aln; b++)
	    {
	      if ( is_gap(A->seq_al[a][b]))
		{for (d=0; d< SUB->nseq; d++)SUB->seq_al[d][b]='-';}
	      else 
	         {
		   for(d=0; d<SUB->nseq; d++)SUB->seq_al[d][b]=((A->S)->Profile[A->order[a][0]])->seq_al[d][e];
		   e++;
		 }
		   
	    }
	  MAIN=stack_aln(MAIN, SUB);
	}
    }
  free_aln (SUB);
  free_aln (A);
  return MAIN;
  }
Alignment * expand_number_aln (Alignment *A)
  {
  /*This function expands the profiles within an alignment*/
  

  int a, b, d, e;
  Alignment *MAIN=NULL, *SUB=NULL, *C=NULL;
  int n_sub_seq=0;
  int new_nseq=0;
  int *list;
  

  if ( !A)return A;
  
  list=vcalloc (A->nseq, sizeof (int)); 
  for ( a=0; a< A->nseq; a++)
    {
      if ( (A->S)->Profile[A->order[a][0]])new_nseq+=((A->S)->Profile[A->order[a][0]])->nseq;
      else 
	{
	  new_nseq++;
	  list[n_sub_seq++]=a;
	}      
    }
  
  if ( n_sub_seq==A->nseq){vfree(list);return A;}
  else if (n_sub_seq==0){MAIN=copy_aln (A, MAIN);MAIN->nseq=0;}
  else
    {
      MAIN=extract_sub_aln (A, n_sub_seq, list);      
    }
  
  
  list[0]=A->nseq;
  C=extract_sub_aln (A,1, list);  
  vfree(list);
  
  
  
  for ( a=0; a< A->nseq; a++)
    {
      if ( (A->S)->Profile[A->order[a][0]])
	{
	  SUB=copy_aln ( (A->S)->Profile[A->order[a][0]],SUB);
	  SUB=realloc_aln2(SUB, SUB->nseq, A->len_aln+1);
	  
	  for ( e=0,b=0; b<= A->len_aln; b++)
	    {

	      
		  
	      if (is_gap(A->seq_al[a][b]))for ( d=0; d<SUB->nseq; d++)SUB->seq_al[d][b]='-';
	      else
		{
		  for ( d=0; d<SUB->nseq; d++)
		    {
		      if ( is_gap (((A->S)->Profile[A->order[a][0]])->seq_al[d][e]))SUB->seq_al[d][b]='-';
		      else SUB->seq_al[d][b]=A->seq_al[a][b];
		    }
		  e++;
		}
	    }
	  for (d=0; d< SUB->nseq; d++)SUB->score_seq[d]=A->score_seq[a];
	  
	  MAIN=stack_aln(MAIN, SUB);
	}
    }
  
  MAIN=stack_aln(MAIN, C);
  MAIN->nseq--;
  MAIN->score=MAIN->score_aln=A->score_aln;
  free_aln (SUB);
  free_aln (A);
  
  free_aln (C);
  
  return MAIN;
  }
Alignment * ungap_sub_aln (Alignment *A, int ns, int *ls)
        {

	int a, b, c,t;
	int len;

	len=strlen ( A->seq_al[ls[0]]);

	for ( c=0,a=0; a<len; a++)
		{
		for ( t=0,b=0; b<ns; b++)
			t+=is_gap(A->seq_al[ls[b]][a]);
		if (t==ns);
		else
		    {
		    for ( b=0; b<ns; b++)
		    	A->seq_al[ls[b]][c]=A->seq_al[ls[b]][a];
		    c++;
		    }
	 	}
	 for ( b=0; b<ns; b++)A->seq_al[ls[b]][c]='\0';  
	 return A;
	}
Sequence * ungap_seq ( Sequence *S)
        {
	  int a;
	  
	  if ( !S)return NULL;
	  ungap(S->seq[0]);
	  S->max_len=S->min_len=strlen (S->seq[0]);
	  for ( a=0; a< S->nseq; a++)
	    {
	      ungap(S->seq[a]);
	      S->len[a]=strlen (S->seq[a]);
	      S->max_len=MAX(S->max_len,S->len[a]);
	      S->min_len=MAX(S->min_len,S->len[a]);
	    }
	  return S;
	  
	}

Alignment *ungap_aln_n ( Alignment *A, int p)
	{
	int a, b, c;
	int t;
	

	for ( c=0,a=0; a< A->len_aln; a++)
		{
		for ( t=0,b=0; b<A->nseq; b++)
			t+=is_gap(A->seq_al[b][a]);
		if ((t*100)/A->nseq>=p || (t==A->nseq && p==100) || (t && p==1));
		else
		    {
		    for ( b=0; b<A->nseq; b++)
		    	A->seq_al[b][c]=A->seq_al[b][a];
		    c++;
		    }
	 	}
	 for ( b=0; b<A->nseq; b++)A->seq_al[b][c]='\0';
	 A->len_aln=c; 
	 return A;

	 }

Alignment *ungap_aln ( Alignment *A)
{
  return ungap_aln_n (A, 100);
}
/*
Alignment *ungap_aln ( Alignment *A)
	{
	int a, b, c,t;
	
	for ( c=0,a=0; a< A->len_aln; a++)
		{
		for ( t=0,b=0; b<A->nseq; b++)
			t+=is_gap(A->seq_al[b][a]);
		if (t==A->nseq);
		else
		    {
		    for ( b=0; b<A->nseq; b++)
		    	A->seq_al[b][c]=A->seq_al[b][a];
		    c++;
		    }
	 	}
	 for ( b=0; b<A->nseq; b++)A->seq_al[b][c]='\0';
	 A->len_aln=c; 
	 return A;

	 }
*/


Alignment *remove_end (Alignment *A)
        {
	int a, b, d;
	int left, right;

	for (a=0; a< A->len_aln; a++)
	    {
	    for ( b=0, d=0; b< A->nseq; b++)
		if ( !is_gap(A->seq_al[b][a]))d++;
	    if ( d>1)break;
	    }
	left=a;
	for (a=A->len_aln-1; a>0; a--)
	    {
	    for ( b=0, d=0; b< A->nseq; b++)
		if ( !is_gap(A->seq_al[b][a]))d++;
	    if ( d>1)break;
	    }
	right=a;

	return extract_aln(A, left, right+1);
	}
void compress_aln ( Alignment *A)
        {
	int a, b, c, d;
        
        

	for (c=0, a=0; a< A->len_aln; a++)
	  {
	    for ( b=0, d=0; b< A->nseq; b++)
		if ( A->seq_al[b][a]!='-'){d=1; break;}
	    if ( d==0);
	    else
	      {
		for (b=0; b< A->nseq; b++)
		  A->seq_al[b][c]=A->seq_al[b][a];
	      c++;
	      }
	  }
	 A->len_aln=c;
	
	for ( a=0; a< A->nseq; a++)
	  A->seq_al[a][c]='\0';
	}

Alignment *seq_coor2aln ( Sequence *S, Alignment *A, int **coor, int nseq)
        {
	int a;
	char *buf;

	A=realloc_alignment2(A, nseq, return_maxlen ( S->seq, S->nseq)+1);
	for ( a=0; a< S->nseq; a++)sprintf ( A->file[a], "%s", S->file[a]);
	for ( a=0; a< nseq; a++)
	    {
	    sprintf (A->name[a], "Repeat_%d_%d", a, coor[a][0]);
	    buf=extract_char ( S->seq[coor[a][0]], coor[a][1]-1, coor[a][2]);
	    sprintf ( A->seq_al[a],"%s", buf);
	    free(buf);
	    A->order[a][0]=0;
	    A->order[a][1]=coor[a][1]-1;
	    }
	A->nseq=nseq;
	return A;
	}

Alignment *strings2aln (int nseq,...)
        {
	  /*strings2aln(nseq, <name1>, <seq1>, <name2>, <seq2>....)*/
	  va_list ap;
	  char **list, **list2;
	  char **name, **name2;
	  Sequence *S;
	  Alignment *A;
	  int a, max;

	  va_start(ap, nseq);
	  list=vcalloc (nseq, sizeof (char*));
	  name=vcalloc (nseq, sizeof (char*));
	  for ( a=0; a< nseq; a++)
	    {
	      name[a]=va_arg(ap,char*);
	      list[a]=va_arg(ap,char*);
	      
	    }
	  va_end(ap);
	  
	  for ( max=0,a=0; a< nseq; a++)
	    {
	      max=(strlen (list[a])>max)?strlen(list[a]):max;
	    }
	  list2=declare_char (nseq, max+1);
	  name2=declare_char (nseq, MAXNAMES+1);
	  
	  for ( a=0; a< nseq; a++)
	    {
	      sprintf ( list2[a], "%s", list[a]);
	      sprintf ( name2[a], "%s", name[a]);
	    }

	  
	  S=fill_sequence_struc(nseq,list2,name2);
	  
	  free_char (list2, -1);
	  free_char (name2, -1);
	  free (list);
	  free(name);
	  A=seq2aln(S,NULL, 1); 
	  return A;
	}
	  
Alignment *seq2aln ( Sequence *S, Alignment *A,int rm_gap)
	{
	int a, max_len;
	char *buf2;

	

	A=realloc_alignment2(A, S->nseq, S->max_len+1); 	
	for ( a=0; a< S->nseq; a++)sprintf ( A->file[a], "%s", S->file[a]);
	A->nseq=S->nseq;
	A->max_len=S->max_len;
	A->min_len=S->min_len;

	for ( a=0; a< S->nseq; a++)
		{
		A->order[a][0]=a;
		A->order[a][1]=0;
		sprintf ( A->comment[a], "%s", S->comment[a]);
		
		sprintf ( A->name[a], "%s", S->name[a]);
		sprintf ( A->seq_al[a], "%s", S->seq[a]);
	
		ungap ( A->seq_al[a]);
		A->len[a]=strlen ( A->seq_al[a]); 
		
		if ( rm_gap==0)sprintf ( A->seq_al[a], "%s", S->seq[a]);
		
		}
	
	max_len=get_longest_string  (A->seq_al,A->nseq, NULL, NULL);
	
	for (a=0; a<A->nseq; a++)
	    {
	    buf2=generate_null (max_len-strlen ( A->seq_al[a]));
	    strcat ( A->seq_al[a], buf2);	
	    vfree (buf2);	    
	    }	

	A->len_aln=max_len;
	A->S=S;
	return A;
	}

Alignment * add_align_seq2aln ( Alignment *A, char *seq, char *seq_name)
        {
	  if ( strlen (seq)!=A->len_aln)
	    {
	      fprintf ( stderr, "\nError: Attempt to stack incompatible aln and aligned sequence[FATAL]\n");
	      exit (1);
	      A=NULL;
	    }
	  else
	    {

	      A=realloc_aln2 ( A, A->nseq+1, A->len_aln+1);
	      sprintf ( A->name[A->nseq], "%s", seq_name);
	      sprintf ( A->seq_al[A->nseq], "%s", seq);
	      A->nseq++;
	    }
	  return A;
	}
  
Alignment *aln2number (Alignment *A)
        {
	A->seq_al=char_array2number(A->seq_al, A->nseq);
	return A;
	}
Sequence *seq2number (Sequence *A)
        {
	A->seq=char_array2number(A->seq, A->nseq);
	return A;
	}


Sequence * aln2seq (Alignment *A)
	{
	Sequence *LS;
	int a;
	
	LS=declare_sequence ( A->len_aln+1, A->len_aln+1, A->nseq);
	LS->nseq=A->nseq;
	for ( a=0; a< LS->nseq; a++)
		{
		sprintf (LS->file[a], A->file[a]); 
		LS->len[a]=A->len[a];
		sprintf ( LS->seq[a], A->seq_al[a]);
		if ( LS->contains_gap==0)ungap ( LS->seq[a]);
		
		LS->len[a]=strlen ( LS->seq[a]);
		sprintf ( LS->comment[a], A->comment[a]);
		sprintf ( LS->name[a], "%s", A->name[a]);
		}
	return LS;
	}
Alignment * filter_aln ( Alignment *A, Alignment *ST, int value)
        {
	    int a, b, c;
	    
	    for ( a=0; a< A->nseq; a++)
	        {		
		for ( c=0; c< ST->nseq; c++)if ( strm(ST->name[c], A->name[a]))break;
		for ( b=0; b< A->len_aln; b++)
		    {
			if ( ST->seq_al[c][b]<=value || ST->seq_al[c][b]==NO_COLOR_RESIDUE)A->seq_al[a][b]='-';
		    }
		}
	    return A;

	}
Alignment * filter_aln_upper_lower ( Alignment *A, Alignment *ST, int value)
        {
	    int a, b, c, st;
	    

	    /*Turns to lowercase everything below or equal to value*/
	    A->residue_case=2;
	    for ( a=0; a< A->nseq; a++)
	        {
		if(value!=10)for ( c=0; c< ST->nseq; c++)if ( strm(ST->name[c], A->name[a]))break;
		for ( b=0; b< A->len_aln; b++)
		    {
		        if ( value==10)st=0;
		        else if ( isdigit( ST->seq_al[c][b]))st=ST->seq_al[c][b]-'0';
			else st=ST->seq_al[c][b];
			if ( st<=value || st==NO_COLOR_RESIDUE)A->seq_al[a][b]=tolower (A->seq_al[a][b]);
		
		    }
		}
	    return A;

	}
Alignment * filter_aln_lower_upper ( Alignment *A, Alignment *ST, int value)
        {
	  int a, b, c, st;
	  /*Turns to uppercase everything above or equal to value*/
	    A->residue_case=2;
	    for ( a=0; a< A->nseq; a++)
	        {
		if(value!=10)for ( c=0; c< ST->nseq; c++)if ( strm(ST->name[c], A->name[a]))break;
		for ( b=0; b< A->len_aln; b++)
		    {
		      if ( value==10)st=10;
		      else if ( isdigit( ST->seq_al[c][b]))st=ST->seq_al[c][b]-'0';
		      else st=ST->seq_al[c][b];
		      if ( st>=value || st==NO_COLOR_RESIDUE)A->seq_al[a][b]=toupper (A->seq_al[a][b]);
			
		    }
		}
	    return A;  
	}
Alignment * filter_aln_convert ( Alignment *A, Alignment *ST, int value, int n_symbol,char **symbol_list)
        {
	  int a, b, c;
	  int st;

	  A->residue_case=2;
	  for ( a=0; a< A->nseq; a++)
	        {
		if(value!=10)for ( c=0; c< ST->nseq; c++)if ( strm(ST->name[c], A->name[a]))break;
		for ( b=0; b< A->len_aln; b++)
		    {
		      if ( value==10)st=11;
		      else st=(isdigit(ST->seq_al[c][b]))?ST->seq_al[c][b]-'0':ST->seq_al[c][b];
		      if ( st>=value) 
			{
			  A->seq_al[a][b]=convert(A->seq_al[a][b],n_symbol,symbol_list);
			}
		    }
		}
	  return A;
	} 

char *dna_aln2cons_seq ( Alignment *A)
        {
	int a, b, best;
	static int **column_count;
	static int **old_tot_count;
	static int **new_tot_count;
	static char *string1, *string2;
	int **count_buf;
	char r1, r2,*seq;
	int NA=0, NG=1, NC=2, NT=3, IGAP=4;
	static int   MAX_EST_SIZE=10000;
	static int   size_increment=1000;
	static int first;
	int overlap=0, best_overlap=0;
	

	seq=vcalloc ( A->len_aln+1, sizeof (char));

	if (!column_count )
	  {
	    column_count=vcalloc(MAX_EST_SIZE, sizeof (int*));
	    for ( a=0; a< MAX_EST_SIZE; a++)
	      column_count[a]=vcalloc (5, sizeof (int));
	  
	    old_tot_count=vcalloc(MAX_EST_SIZE, sizeof (int*));
	    new_tot_count=vcalloc(MAX_EST_SIZE, sizeof (int*));
	    A->P=declare_profile( "agct-",MAX_EST_SIZE);
	    string1=vcalloc (MAX_EST_SIZE, sizeof (char));
	    string2=vcalloc (MAX_EST_SIZE, sizeof (char));
	  }
	else if (A->len_aln>MAX_EST_SIZE)
	  {
	    if ( column_count)
	      {
		for ( a=0; a< MAX_EST_SIZE; a++)
		  free(column_count[a]);
		free(column_count);
		free(old_tot_count);
		free(new_tot_count);
		free(string1);
		free(string2);
	      }
	    
	  column_count=vcalloc(MAX_EST_SIZE+ size_increment, sizeof (int*));
	  for ( a=0; a< MAX_EST_SIZE+ size_increment; a++)
	      column_count[a]=vcalloc (5, sizeof (int));
	  
	  old_tot_count=vcalloc(MAX_EST_SIZE+ size_increment, sizeof (int*));
	  new_tot_count=vcalloc(MAX_EST_SIZE+ size_increment, sizeof (int*));
	  
	  for (a=0; a< MAX_EST_SIZE; a++)
	    {
	      old_tot_count[a]=*(column_count++);
	      for ( b=0; b<5; b++)old_tot_count[a][b]=(A->P)->count[b][a];
	    }
	  free_int ( (A->P)->count, -1);
	  
	  (A->P)->count=declare_int (5, MAX_EST_SIZE+ size_increment);
	  (A->P)->max_len=MAX_EST_SIZE+ size_increment;
	  MAX_EST_SIZE+= size_increment;
	  string1=vcalloc (MAX_EST_SIZE, sizeof (char));
	  string2=vcalloc (MAX_EST_SIZE, sizeof (char));
	  }
	
	
	sprintf ( string1, "%s",A->seq_al[0]);
	sprintf ( string2, "%s",A->seq_al[1]);
	

	string1=mark_internal_gaps(string1,'.');
	string2=mark_internal_gaps(string2,'.');

	
	
	for (b=0,a=0; a< A->len_aln; a++)
	  {
	    r1=string1[a];
	    r2=string2[a];
	    
	    if ( r1==r2)
	      {
		overlap++;
	      }
	    else
	      {
		best_overlap=MAX(overlap, best_overlap);
		overlap=0;
	      }


	    if (!is_gap(r1) && first==1)new_tot_count[a]=old_tot_count[b++]; 
	    else if (is_gap(r1) || first==0){new_tot_count[a]=*column_count;column_count++;};
	    
	    if ( first==0)
	      {
		if(r1=='a')       new_tot_count[a][NA]++;
		else if ( r1=='g')new_tot_count[a][NG]++;
		else if ( r1=='c')new_tot_count[a][NC]++;
		else if ( r1=='t')new_tot_count[a][NT]++;	
		else if (is_gap(r1));
		else
		  {
		   new_tot_count[a][NA]++;
		   new_tot_count[a][NG]++;
		   new_tot_count[a][NC]++;
		   new_tot_count[a][NT]++;
		  }
	      }
	    if ( a> 0 && a<A->len_aln-1 && r1=='.')
	      {
		new_tot_count[a][IGAP]+=((new_tot_count[a-1][NA]+new_tot_count[a-1][NG]+new_tot_count[a-1][NC]+new_tot_count[a-1][NT]));
	      }
	    

	    if(r2=='a')       new_tot_count[a][NA]++;
	    else if ( r2=='g')new_tot_count[a][NG]++;
	    else if ( r2=='c')new_tot_count[a][NC]++;
	    else if ( r2=='t')new_tot_count[a][NT]++;
	    else if ( r2=='.')new_tot_count[a][IGAP]++;
	    else if ( r2=='-');
	    else 
	      {
		new_tot_count[a][NA]++;
		new_tot_count[a][NG]++;
		new_tot_count[a][NC]++;
		new_tot_count[a][NT]++;	
	      }
	    (A->P)->count[0][a]=new_tot_count[a][NA];
	    (A->P)->count[1][a]=new_tot_count[a][NG];
	    (A->P)->count[2][a]=new_tot_count[a][NC];
	    (A->P)->count[3][a]=new_tot_count[a][NT];
	    (A->P)->count[4][a]=new_tot_count[a][IGAP];

	    best_int(4,1, &best,new_tot_count[a][NA], new_tot_count[a][NG],new_tot_count[a][NC],new_tot_count[a][NT]); 
	    if( best==0)      seq[a]='a';
	    else if ( best==1)seq[a]='g';
	    else if ( best==2)seq[a]='c';
	    else if ( best==3)seq[a]='t';
	  }

	first=1;

	seq[a]='\0';
	fprintf ( stderr, "[Best Overlap: %d Residues]", best_overlap);
	count_buf=old_tot_count;
	old_tot_count=new_tot_count;
	new_tot_count=count_buf;

	return seq;
	
	}
char *aln2cons_seq ( Alignment *A, int ns, int *ls, int n_groups, char **group_list)
        {
	char *seq;
	int a, b, c;
	int best_group=0;
	int aa_group=0;
	int *group;
	int len;

	len=strlen  (A->seq_al[ls[0]]);
	seq=vcalloc (len+1, sizeof (char));



	if ( !group_list)
	   {
	       group_list=declare_char ( 26, 2);
	       for ( a=0; a<26; a++)group_list[a][0]=a+'a';
	       n_groups=26;
	       aa_group=1;
	   }
	
	
	for ( a=0; a<len; a++)
	    {
		group=vcalloc (n_groups+1, sizeof (int));
		for (best_group=0,b=0; b< ns; b++)
		    {
		    if ( !is_gap(A->seq_al[ls[b]][a]))
			 {
			 for (c=0; c< n_groups; c++)
			     if ( is_in_set (tolower(A->seq_al[ls[b]][a]), group_list[c]))
			                 {group[c]++;
					  best_group=(group[c]>group[best_group])?c:best_group;
					 }
			 }
		    seq[a]=group_list[best_group][0];
		    }
		vfree (group);
	    }
	seq[a]='\0';
	if ( aa_group) free_char (group_list, -1);


	return seq;
	}
		
char *aln2cons_seq_mat ( Alignment *A, char *mat_name)
        {
	  int a, b, c;
	  char *seq, r1, r2;
	  int **mat;
	  int score=0, best_score, best_r;
	  
	  seq=vcalloc ( A->len_aln+1, sizeof (char));

	  mat=read_matrice (mat_name);
	  for ( a=0; a< A->len_aln; a++)	    
	    {
	      for (b=0; b<20; b++)
		{
		  r1=AA_ALPHABET[b];
		  for ( score=0,c=0; c< A->nseq; c++)
		    {
		      if (is_gap(A->seq_al[c][a]))continue;
		      else 
			{
			  r2=A->seq_al[c][a];
			  score+=mat[r1-'a'][r2-'a'];
			}
		    }
		  if ( b==0 || score>best_score){best_score=score; best_r=r1;}
		}
	      seq[a]=best_r;
	    }
	  free_int (mat, -1);
	  return seq;
	}
			  
	  
	  
  
int  seq_list2fasta_file( Sequence *S,  char *list, char *file)
        {
	FILE *fp;
	int n, a, s;
	
	fp=vfopen ( file, "w");
	n= atoi(strtok (list,SEPARATORS));
	for ( a=0; a< n; a++)
	    {
	    s=atoi(strtok (NULL, SEPARATORS));
	    fprintf ( fp, ">%s\n%s\n", S->name[s], S->seq[s]);
	    }
	vfclose (fp);
	return 1;
	}
Structure * seq2struc ( Sequence *S, Structure *ST)
        {
	int a, b;
	
	for ( a=0; a< S->nseq; a++)
	    for ( b=0; b< S->len[a]; b++)
		ST->struc[a][b+1][ST->n_fields-1]=S->seq[a][b];
	return ST;
	}

void aln2struc (Alignment *A, Structure *ST) 
        {
	int a, b, c;

	for ( a=0; a< A->nseq; a++)
	    for (c=0, b=0; b< A->len_aln; b++)
	        {
		if ( !is_gap (A->seq_al[a][b]))
		     {
		     ST->struc[a][c][ST->n_fields-1]=A->seq_al[a][b];
		     c++;
		     }
		}
	}
Alignment *stack_aln (Alignment *A, Alignment *B)
        {
	int a,b;

	if ( B==NULL)return A;


	A=realloc_aln2 ( A, ((A==NULL)?0:A->nseq)+B->nseq, MAX((A==NULL)?0:A->len_aln,B->len_aln)+1);
	
	for (a=A->nseq,b=0; b< B->nseq; b++, a++)
	    {
	    sprintf ( A->comment[a] , "%s", B->comment[b]);
	    sprintf ( A->seq_al [a] , "%s", B->seq_al [b]);
	    sprintf ( A->name   [a] , "%s", B->name[b]);
	    sprintf ( A->file   [a], "%s" , B->file[b]);
	    A->order[a][0]=B->order[b][0];
	    A->order[a][1]=B->order[b][1];
	    A->score_seq[a]=B->score_seq[b];
	    A->len[a]=B->len[b];
	    }
	
	A->len_aln=MAX(A->len_aln, B->len_aln);
	A->nseq=A->nseq+B->nseq;
	A->score_aln=A->score_aln+B->score_aln;
	
	A->finished=A->finished+B->finished;
	return A;
	}
	    
Alignment *chseqIaln(char *name, int seq_n, int start,int len,Sequence *S, int seqIaln, Alignment *A)
        {
	char *seq;

	seq=extract_char ( S->seq[seq_n], start, len);
	A=realloc_aln2 (A, (A==NULL)?(seqIaln+1):MAX(A->nseq,seqIaln+1), ((A==NULL)?(strlen (seq)):MAX(strlen (seq),A->len_aln))+1);
	
	
	sprintf ( A->seq_al[seqIaln], "%s",seq);

	
	A->order[seqIaln][0]=seq_n;
	A->order[seqIaln][1]=start;
	sprintf ( A->name[seqIaln], "%s", name);
	A->nseq=MAX(A->nseq, seqIaln+1);
	A->len_aln=return_maxlen(A->seq_al, A->nseq);
	A->S=S;
	free (seq);
	return A;
	}

Alignment * aln_gap2random_aa(Alignment *A)
        {
	 int a, b,l;
	 char alp[200];
	 
	 if (strm ( (A->S)->type, "PROTEIN"))
	   sprintf ( alp, "acefghiklmnpqrstuvwy");
	 else if ( strm ( (A->S)->type, "DNA"))
	   sprintf ( alp, "agct");
	 l=strlen (alp);
	 
	 
	 for (a=0; a<A->nseq; a++)
	    for ( b=0; b<A->len_aln; b++)
	      if ( is_gap (A->seq_al[a][b]))A->seq_al[a][b]=alp[(int)rand()%(l)];
	  return A;
	}

Alignment * make_random_aln(Alignment *A,int nseq, int len, char *alphabet)
        {
	int a;
	    
	A=realloc_aln2(A, nseq, len+1);
	A->nseq=0;
	A->len_aln=len;
	for ( a=0; a< A->nseq; a++)sprintf ( A->file[a], "random alignment");
	for ( a=0; a< nseq; a++)
	    A=add_random_sequence2aln(A,alphabet);
	return A;
	}
Alignment * add_random_sequence2aln( Alignment *A, char *alphabet)
        {
	static int init;
	int a, n;
	if (init==0)addrandinit((init=100));
	n=strlen(alphabet);
	A=realloc_alignment2 (A, A->nseq+1, A->len_aln+1);
	    
	for ( a=0; a< A->len_aln; a++)A->seq_al[A->nseq][a]=alphabet[addrand((unsigned long)n)];
	A->nseq++;
	return A;
	}

Sequence *get_defined_residues( Alignment *A)
        {
	    char *buf;
	    Sequence *S;
	    int a, b, s, l, r;
	    if ( !A || !A->S) return NULL;

	    S=duplicate_sequence (A->S);
	    for ( a=0; a< S->nseq; a++)
		for ( b=0; b< S->len[a]; b++)S->seq[a][b]=UNDEFINED_RESIDUE;
	    buf=vcalloc(A->len_aln+1,sizeof (char));
	    for ( a=0; a< A->nseq; a++)
	        {
		    sprintf ( buf, "%s",A->seq_al[a]);
		    ungap(buf);
		    l=strlen (buf);
		    s=A->order[a][0];
		    
		    for ( b=1; b<= l; b++)
			{
			r=A->seq_cache[s][b];
			S->seq[s][r-1]=(A->S)->seq[s][r-1];
			}
		}
	    free(buf);
	    return S;
	}
Alignment *thread_defined_residues_on_aln ( Alignment *A, Sequence *S1)
	{
	int a, b;
	int gap, r,s, r2;
	for ( a=0; a< A->nseq; a++)
	    {
		s=A->order[a][0];
		r=A->order[a][1];
		for (b=0;b< A->len_aln; b++)
		    {
		    gap=is_gap(A->seq_al[a][b]);
		    
		    if (!gap)
			{
			r+=!gap;
			r2=A->seq_cache[s][r]-1;
			if ( S1->seq[s][r2]==UNDEFINED_RESIDUE)
			    A->seq_al[a][b]=UNDEFINED_RESIDUE;
			}
		    }
	    }
	return A;
	}

int ** trim_aln_borders (char **seq1, char **seq2, int nseq)
        {
	int a, b, c,l1,l2;
	char *buf1;
	char *buf2;
	int max;

	

	
	max=MAX(get_longest_string (seq1,-1, NULL, NULL),get_longest_string (seq2,-1, NULL, NULL))+1;
	buf1=vcalloc ( max, sizeof(char));
	buf2=vcalloc ( max, sizeof(char));
	
	for ( a=0; a< nseq; a++)
	    {
	    sprintf ( buf1, "%s", seq1[a]);
	    sprintf ( buf2, "%s", seq2[a]);


	   
	    ungap (buf1);
	    ungap (buf2);

	    if (str_overlap ( buf1, buf2,'*')!=0)
		{		      
		l1=strlen ( seq1[a]);
	        l2=strlen ( seq2[a]);
		for ( b=0,c=0; c< l1; c++)
		    if ( !is_gap(seq1[a][c]))seq1[a][c]=buf1[b++];
		seq1[a][c]='\0';
		for ( b=0,c=0; c< l2; c++)
		    if ( !is_gap(seq2[a][c]))seq2[a][c]=buf2[b++]; 
		seq2[a][c]='\0';
      		}
	    }
	free (buf1);
        free (buf2);
	return NULL;

	}
Sequence * merge_seq    ( Sequence *IN, Sequence *OUT)
        {
	int a;
	
	
	
	if ( OUT==NULL)return duplicate_sequence (IN);
	else
	    {
	     if ( IN && check_list_for_dup( IN->name, IN->nseq))
	          {
		      fprintf ( stderr, "\nERROR: %s is duplicated in file %s[FATAL]\n", check_list_for_dup( IN->name, IN->nseq), IN->file[0]);
		      exit (1);
		  }
	    for ( a=0; a< IN->nseq; a++)
		if ((OUT=add_sequence ( IN, OUT, a))==NULL)return NULL;
	    return OUT;
	    }
	}
char * seq_name2coor ( char *s, int *start, int *end, char sep)
{
  /*name|start|end */
  char n1[100], n2[100];
  int a=0, b=0, c=0;
  
  n1[0]=n2[0]='\0';
  start[0]=end[0]=0;
  
  while ( s[a]!=sep && s[a]!='\0')a++;
  if ( s[a]=='\0')return s;
  else 
    s[a++]='\0';

 
  
  while ( s[a]!=sep && s[a]!='\0')n1[b++]=s[a++];
  
  if ( s[a]=='\0'){n1[b]='\0';if ( n1[0])start[0]=atoi(n1);return s;}
  else s[a++]=n1[b]='\0';
 
  
  while ( s[a]!=sep && s[a]!='\0')n2[c++]=s[a++];
  n2[c]='\0';

  
  if ( n1[0])start[0]=atoi(n1);
  if ( n2[0])end[0]=atoi(n2);


  return s;
}
  
Sequence *extract_one_seq(char *n,int start, int end, Alignment *S, int keep_name)
       {
	
	 int seq, a;
	 FILE*fp;
	 char *name;
	 Sequence *OUT_S;
	 
	 if ( n[0]=='#')seq=S->nseq;
	 else if (is_number (n) && (seq=atoi(n))!=0) seq=atoi(n);
	 else if ( (seq=name_is_in_list (n, S->name, S->nseq, 100)+1)!=0);
	 else
	   {
	     fprintf ( stderr, "\nCould not find Sequence %s [FATAL]", n);
	     exit (1);
	   }
	 seq--;
	 
	 name=vtmpnam ( NULL);
	 fp=vfopen ( name, "w");
	 if ( start && end &&!keep_name)fprintf (fp, ">%s_%d_%d\n",S->name[seq],start, end);
	 else if ( start && end==0)fprintf (fp, ">%s_%d_%d\n",S->name[seq],start,strlen ( S->seq_al[seq]));
	 else fprintf (fp, ">%s\n", S->name[seq]);
	 
	 if ( start==0 && end==0)fprintf (fp, "%s\n", S->seq_al[seq]);
	 else if (end==0)fprintf (fp, "%s\n", S->seq_al[seq]+start-1);
	 else
	   {
	     for ( a=start-1; a<end; a++)fprintf ( fp, "%c", S->seq_al[seq][a]);
	     fprintf ( fp, "\n");
	   }
	 
	 
	 vfclose (fp);
	 OUT_S=get_fasta_sequence (name, NULL);
	 return OUT_S;
       }
	 

	 
Sequence * extract_sub_seq( Sequence  *COOR, Sequence *S)
        {
	int a, b, c,s;
	int start, end;
	
	for ( a=0; a< S->nseq; a++)
	    {
	    if ( (s=name_is_in_list ( S->name[a], COOR->name, COOR->nseq, 100))!=-1)
		 {
		 
		 sscanf ( COOR->comment[s], "%d %d", &start, &end);
		 for (c=0,b=start-1; b< end; b++, c++)S->seq[a][c]=S->seq[a][b];
		 S->seq[a][c]='\0';
		 sprintf ( S->comment[a], "%s",COOR->comment[s]);
		 
		 }
	    }
	S=reorder_seq ( S, COOR->name, COOR->nseq);
	return S;
	}
		 
	     
Alignment * fix_aln_seq  ( Alignment *A, Sequence *S)
        {
	int a, b, c, d;
	char *buf;
	static char *mat;
	static char *dp_mode;


	
	Alignment *B;
	

	if ( !mat)mat=vcalloc (STRING, sizeof (char));
	if ( !dp_mode) dp_mode=vcalloc (STRING, sizeof (char));

	if ( S==NULL)return A;
	buf=vcalloc ( A->len_aln+1, sizeof (char));
	
	A->seq_cache=declare_int ( S->nseq, S->max_len+2);
	
	for ( a=0; a< S->nseq; a++)
	    {
	    for (b=0; b< A->nseq; b++) 
	        {
		if (strm ( S->name[a], A->name[b]))
		   {
		   A->order[b][0]=a;
		   sprintf (buf, "%s", A->seq_al[b]);
		   ungap (buf);
		   
		   if ( strm (buf, S->seq[a]))
		       {
			   
			   for ( c=0; c<S->len[a]; c++)A->seq_cache[a][c+1]=c+1;
		       }
		   else
		       {
			   sprintf ( mat, "idmat");
			   sprintf ( dp_mode, "fasta_pair_wise");
			   B=align_two_sequences (S->seq[a],buf,mat, 0, 0, dp_mode);
			   if ( S->len[a]!=B->len_aln)
			      {
				  fprintf (stderr, "\nERROR: IMPOSSIBLE TO RECONCILIATE SEQ %s in file %s\nCHECK THAT TWO != SEQUENCES DO NOT HAVE THE SAME NAME\n",S->name[a], A->file[0]);
				  print_aln (B);
				  exit (1);
			      }
			   else 
			      {
				  /*
				    if ( !flag){fprintf ( stderr, "\nWARNING: RECONCILIATION(S) IN FILE %s\n",A->file);flag=1;}
				    fprintf (stderr, "\tSEQ %s [%d-%d]\n", S->name[a], B->len_aln, strlen (buf));
				  */
				  
				   for (d=0,c=0; c<B->len_aln; c++)
				       {
					   if (!is_gap(B->seq_al[1][c]))
					      {
						  d++; 
						  A->seq_cache[a][d]=c+1;
						  
					      }
				       }
			       }
			   free_aln (B);
		       }
		   
		   }
		}
	    }
	return A;
	}

Sequence * add_prf2seq  ( Alignment *A, Sequence *S)
    {
      char **new_seq;
      Sequence *NS;
      new_seq=declare_char (1, A->len_aln+1);
      sprintf ( new_seq[0], "%s",aln2cons_seq_mat(A, "blosum62mt"));
      
      NS=fill_sequence_struc(1, new_seq,A->file);
      
      S=add_sequence (NS, S, 0);
      free_sequence (NS, NS->nseq);
      S->Profile[S->nseq-1]=copy_aln (A, NULL);
      free_char( new_seq, -1);
      return S;
    }
      
      
Sequence * add_sequence ( Sequence *IN, Sequence *OUT, int i)
    {
    int s;
   
    char *buf;
    if (OUT==NULL)
      {
	OUT=duplicate_sequence (IN);
	return OUT;
      }

    /*Adds sequence i of IN at the end of OUT*/

    if ((s=name_is_in_list ( IN->name[i], OUT->name, OUT->nseq,STRING))==-1 )
       {
	 OUT=realloc_sequence (OUT, OUT->nseq+1, IN->len[i]);	 
	 sprintf ( OUT->name[OUT->nseq],"%s",IN->name[i]);
	 sprintf ( OUT->file[OUT->nseq],"%s",IN->file[i]);
	 sprintf ( OUT->comment[OUT->nseq],"%s",IN->comment[i]);
	 sprintf ( OUT->seq[OUT->nseq],"%s",IN->seq[i]);
	 OUT->len[OUT->nseq]=IN->len[i];
	 OUT->nseq++;
	 return OUT;
       }
    else if ( s!=-1 && strcmp ( IN->seq[i], OUT->seq[s])!=0)
       {
       	 
       fprintf ( stderr, "\nWARNING: DISCREPANCY:%s in [%s] and  [%s]", IN->name[i], IN->file[i], OUT->file[s]);
       
       if (((buf=build_consensus(IN->seq[i], OUT->seq[s],"fasta_pair_wise" ))!=NULL)||((buf=build_consensus(IN->seq[i], OUT->seq[s],"myers_miller_pair_wise" ))!=NULL))
	   {
	   OUT->max_len=MAX(OUT->max_len, strlen(buf));
	   OUT->min_len=MIN(OUT->min_len, strlen(buf));
	   OUT->seq    =realloc_char ( OUT->seq, -1, -1,OUT->nseq,OUT->max_len+1);
	 
	   sprintf ( OUT->seq[s],"%s",buf);
	   OUT->len[s]=strlen (buf);
	   free (buf);
	   return OUT;
	   }
       else
           {
	       fprintf ( stderr, "IMPOSSIBLE TO RECONCILIATE SOME SEQUENCES[FATAL:%s]\n", PROGRAM);
	       print_aln ( align_two_sequences (IN->seq[i], OUT->seq[s], "idmat", 0, 0, "fasta_pair_wise"));
	       return NULL;
	   }
       
       }
    else
       {
       return OUT;
       }
    }
			   

Sequence  * trim_seq       ( Sequence  *A, Sequence  *B)
        {
	int a;
	static char **name_list;
	int n;
	
	name_list=realloc_char (name_list, -1, -1, MAX(A->nseq, B->nseq),STRING+1);

	for (n=0, a=0; a< A->nseq; a++)
	    {  
	    if ( name_is_in_list ( A->name[a], B->name, B->nseq,STRING+1)!=-1)
	        {
		sprintf ( name_list[n++], "%s", A->name[a]);
		}
	    }
	reorder_seq ( A, name_list, n);
	reorder_seq ( B, name_list, n);
	return B;
	}
Sequence * trim_aln_seq ( Alignment *A, Alignment *B)
        {
	int a;
	static char **name_list;
	int n=0;
	Sequence *SA, *SB;
	
	/*This function inputs two alignments A and B
	  It removes sequences that are not common to both of them
	  It rearange the sequences so that they are in the same order
	  A decides on the order
	  The Sequences (A->S) and (B->S) are treated the same way
	  Sequences are also merged in order to detects discrepencies.
	  A pointer to S is returned
	*/

        name_list=realloc_char (name_list, -1, -1, MAX(A->nseq, B->nseq), STRING+1);
	for ( a=0; a< A->nseq; a++)
	    {  
	    if ( name_is_in_list ( A->name[a], B->name, B->nseq,STRING)!=-1)
	        {
		sprintf ( name_list[n++], "%s", A->name[a]);
		}
	    }
	
	
	
	reorder_aln ( A, name_list, n);
	reorder_aln ( B, name_list, n);
	for ( a=0; a< n; a++)A->order[a][0]=B->order[a][0]=a;
	

	
	SA=aln2seq(A);
	SB=aln2seq(B);
	
	A->S=B->S=merge_seq (SA, SB);
	

	return A->S;
	}






Alignment * reorder_aln ( Alignment *A, char **name, int nseq)
	{
	int a, b,sn, sn2;
	static Alignment *BUF;
	int  n=0;
	int *tpp_int;

	
	BUF=copy_aln ( A, BUF);


	for ( a=0; a<nseq; a++)
		{
		sn =name_is_in_list ( name[a],BUF->name, A->nseq,STRING);
		sn2=name_is_in_list ( name[a],(A->S)->name, (A->S)->nseq,STRING);
		if ( sn==-1)
			{
			    ;
			}
		else
		    {
		    SWAPP(A->order[n], BUF->order[sn], tpp_int);
		    
		    sprintf ( A->name[n], "%s", BUF->name[sn]);		    
		    sprintf ( A->seq_al[n], "%s",BUF->seq_al[sn]);
		    if ( A->seq_cache)for ( b=0; b< (A->S)->max_len+2; b++)A->seq_cache[n][b]=BUF->seq_cache[sn2][b];
		    sprintf ( A->comment[n], "%s",  BUF->comment[sn]);
		    n++;
		    }
		}
	for ( a=n; a< A->nseq; a++)A->name[a][0]=A->seq_al[a][0]='\0';
	A->nseq=n;
	
	return A;
	} 
Sequence * reorder_seq ( Sequence *A, char **name, int nseq)
	{
	int a,sn,n;
	static char **name_buf;
	static char **seq_buf;
	static char **comment_buf;

	name_buf=realloc_char (name_buf, -1, -1, A->nseq, STRING+1);
	seq_buf =realloc_char (seq_buf , -1, -1, A->nseq,get_longest_string(A->seq,-1, NULL, NULL)+1);
	comment_buf=realloc_char (comment_buf , -1, -1, A->nseq,STRING+1);

	
	for ( a=0; a< A->nseq; a++)
		{
		sprintf ( seq_buf[a], "%s", A->seq[a]);
		sprintf ( name_buf[a], "%s",A->name[a]);
		sprintf ( comment_buf[a], "%s",A->comment[a]);
		}
	
	
	for ( n=0,a=0; a< A->nseq; a++)
		{
		sn=name_is_in_list ( name_buf[a],name, nseq, 100);
		if ( sn==-1)
			{
			
			;
			}
		else
		    {
		    n++;
		    sprintf ( A->name[sn], "%s", name_buf[a]);
		    sprintf ( A->seq[sn], "%s",seq_buf[a]);
		    sprintf ( A->comment[sn], "%s",comment_buf[a]);
		    }
		}
	A->nseq=n;
	return A;
	} 
char * concatenate_seq ( Sequence *S, char *conc, int *order)
        {
	    int a;
	    
	    vfree (conc);
	    conc=vcalloc ( S->nseq*S->max_len, sizeof (char));

	    for ( a=0; a< S->nseq; a++)
	        {
		    conc=strcat ( conc, S->seq[order[a]]);
		}
	    return conc;

	}



Alignment * extract_nol_local_aln(Alignment *A, int start, int max_end)
     {
     A=extract_aln ( A, start, max_end);
     A=trunkate_local_aln (A);
     return A;
     }
Alignment * extract_aln ( Alignment *A, int start, int end)
     {
     int a, b;
     
     if ( start>=end)	
	 {
	 for ( b=0; b<start; b++)A->seq_al[b][0]=0;
	 }
     else
         {
	 for ( a=0; a< A->nseq; a++)
	     {
	     for ( b=0; b<start; b++)A->order[a][1]+=1-is_gap ( A->seq_al[a][b]);
	     crop_string ( A->seq_al[a], start, end); 	 	 
	     }
	 }
     A->len_aln=strlen ( A->seq_al[0]);
     return A;
     }

Alignment * trunkate_local_aln ( Alignment *A)
     {
     int a, b;
     int **pos;
     int **cache;
     int seq;
   
     
     cache=declare_int (return_max_int (A->order,read_size_int ( A->order[-1],0),0)+1,return_max_int (A->order,read_size_int ( A->order[-1],0),1)+A->len_aln+1);    
     pos=aln2pos_simple(A,A->nseq);
     
     for ( b=0; b<A->len_aln; b++)
	 for ( a=0; a< A->nseq; a++)	 
	     {
	     seq=A->order[a][0];
	     if ( pos[a][b]<=0);
	     else if ( pos[a][b]>0)
		 {
		 
		 if (cache[seq][pos[a][b]]==0)cache[seq][pos[a][b]]++;
		 else if ( cache[seq][pos[a][b]]>=1)
		      {	     
		      cache[seq][pos[a][b]]++;
		      A->seq_al[a][b]='\0';
		      }
		 }
	     }
     
     A->len_aln=get_shortest_string ( A->seq_al, A->nseq, NULL, NULL);
     pad_string_array ( A->seq_al, A->nseq, A->len_aln, '-');
     
     free_int ( cache,-1);
     
     
     return A;
     }

int get_nol_aln_border ( Alignment *A, int start, int direction)
     {
     int a, b;
     int **pos;
     int **cache;
     int seq,end;
     
     /*This Function Returns the limit position for a non overlaping alignment*/
     
     cache=declare_int (return_max_int (A->order,read_size_int ( A->order[-1],0),0)+1,return_max_int (A->order,read_size_int ( A->order[-1],0),1)+A->len_aln+1);
     pos=aln2pos_simple(A,A->nseq);
     end=(direction==GO_RIGHT)?A->len_aln:-1;
     

     for ( b=start; b!=end;b+=direction)
	 for ( a=0; a< A->nseq; a++)	 
	     {
	     seq=A->order[a][0];
	     if ( pos[a][b]<=0);
	     else if ( pos[a][b]>0)
		 {
		 
		 if (cache[seq][pos[a][b]]==0)cache[seq][pos[a][b]]++;
		 else if ( cache[seq][pos[a][b]]>=1)
		      {	     
		      cache[seq][pos[a][b]]++;
		      free_int(cache, -1);
		      return b-direction;
		      }
		 }
	     }
     
     free_int ( cache,-1);
     return end-direction;
     }





char * extract_defined_seq ( char *in, int in_of, int in_start, int *aa_def, int dir, int *out_start, char *out)
     {
     int start, end,l;
     int b, c, d;

     

     if ( dir==GO_LEFT){start=in_start-1;}
     else if ( dir==GO_RIGHT){start=in_start+1;}        
	
     end=start;
     while (aa_def[end]!=UNDEFINED)
	 {
	 end+=dir;
	 }
     end-=dir;
     
     if (end<start)SWAP(end,start);
     
     l=strlen ( in);
     out_start[0]=-1;
     for (b=0,d=0,c=in_of;b<l; b++)
         {
	 c+=1-is_gap(in[b]);
	 if ( c>=start && c<=end)
	     {
	     if ( out_start[0]==-1)out_start[0]=c-!is_gap(in[b]);
	     out[d++]=in[b];
	     }
	 }
     out[d]='\0';
     
    
     return out;
     }
Alignment * aln_cat ( Alignment *A, Alignment *B)
     { 
     int a;
     
     if ( A->nseq!=B->nseq) 
	 {
	 fprintf ( stderr, "\nERROR IN ALN CAT: DIFFERENT NSEQ\n");
	 exit(1);
	 }

     A=realloc_alignment2(A, A->nseq,A->len_aln+B->len_aln+1);
    
     for ( a=0;a< A->nseq; a++)
         {	
	 strcat ( A->seq_al[a], B->seq_al[a]);
	 }
     A->len_aln+=B->len_aln;
     return A;
     }
int verify_aln ( Alignment *A, Sequence *S, char *message)
     {
     int a, b, c,s,r;


     for ( a=0;a< A->nseq; a++)
         {
	 s=A->order[a][0];
	 r=A->order[a][1];
	 for ( b=0, c=0; b< A->len_aln; b++)
	     {
	     if ( !is_gap(A->seq_al[a][b]))
		  {
		  if (tolower(A->seq_al[a][b])!=tolower(S->seq[s][c+r]))
		      {
		      fprintf ( stderr, "\n%s\nResidue [%c %d, %c %d] line %d seq %d",message,A->seq_al[a][b], b,S->seq[s][c+r], c+r,a,s);  
		      output_Alignment_with_res_number(A, stderr);
		      exit(1);
		      return 0;
		      }
		  c++;
		  }
	     }
	 }
     return 1;
     }

Alignment *adjust_est_aln ( Alignment *PW, Alignment *M, int s)
{
  /*This function reajusts M, threading M onto PW
   two seqences in PW
   s+1 seq in M
   
   seq 0 PW ----> 0->s-1 in M
   seq 1 PW ----> 1->s   in M;
   
   */
  int a, b;
  static char **array;

  
  int top_M=0;
  int bottom_M=0;
  
  
  if ( array==NULL)
    {
      array=declare_char (500, 100000);
    }

  for ( a=0; a< PW->len_aln; a++)
    {
      if ( is_gap(PW->seq_al[0][a]))
	{
	  for ( b=0; b< s; b++)
	    array[b][a]='-';
	}
      else
	{
	  for ( b=0; b< s; b++)
	    array[b][a]=M->seq_al[b][top_M];
	top_M++;
	}
      
      if ( is_gap(PW->seq_al[1][a]))
	{
	  array[s][a]='-';
	}
      else
	{
	  
	  array[s][a]=M->seq_al[s][bottom_M];
	bottom_M++;
	} 
    }
  
  M->len_aln=PW->len_aln;
  for (a=0; a<s; a++)
    {
      for (b=0; b<PW->len_aln; b++)
	M->seq_al[a][b]=array[a][b];
      M->seq_al[a][b]='\0';
    }


  M->nseq=s+1;
  
  return M;
}
/********************************************************************/
/*                                                                  */
/*                   ALIGNMENT ANALYSES                             */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
int aln_is_aligned ( Alignment *A)
{
  int a, b;
  
  if ( !A)return 0;
  for (a=0; a< A->nseq; a++)
    for ( b=A->len_aln-1; b>0; b--)
      {
	if (!is_gap(A->seq_al[a][b]) && is_gap(A->seq_al[a][b-1]))return 1;
      }
  return 0;
}
	
int* get_cdna_seq_winsim ( int *cache, char *string1, char *string2, char *ignore, char *mode,int *w )
	{
	int len1, len2;
	int a, x;

	
	len1=strlen (string1);
	len2=strlen (string2);
		
	if ( len1!=len2)
	  {
	    fprintf ( stderr, "\nTHE TWO cDNAs DO NOT HAVE THE SAME LENGTH [FATAL:get_cdna_seq_sim:%s", PROGRAM);
	    crash("");
	  }
	
	x=get_cdna_seq_sim(cache, string1, string2, ignore, "");
	for ( a=0; a< len1; a++)
	  w[a]=x;

	fprintf (stderr, "\nWARNING: winsim not implemented for cDNA");
	return w;
	}

int get_cdna_seq_sim ( int *cache, char *string1, char *string2, char *ignore, char *mode)
	{
	int len1;
	int len2;
	int a;
	int pos=0;
	int sim=0;
	char r1, r2;
	
	len1=strlen (string1);
	len2=strlen (string2);


	
	if ( len1!=len2)
	  {
	    fprintf ( stderr, "\nTHE TWO cDNAs DO NOT HAVE THE SAME LENGTH [FATAL:get_cdna_seq_sim:%s", PROGRAM);
	    crash("");
	  }
	    
	for ( a=0; a< len1;)
		{
		 
		  if ( cache[a]==0){a++;continue;}
		  else if ( cache[a]==1)
		    {
		      
		      r1=translate_dna_codon (string1+a, 'x');
		      r2=translate_dna_codon (string2+a, 'x');
		      a+=3;
		    }
		  
		  if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			pos++;
			if (is_in_same_group_aa(r1,r2,0, NULL,mode+4))
				{
				sim++;
				}
			}
		}



	if (pos==0)
		 return 0;
	else	
		return (int) (sim*100)/pos;
	
	}	

int* get_seq_winsim ( char *string1, char *string2, char *ignore, char *mode, int*w)
	{
	int len1, len2, len;
	int left, right;
	int a,b;
	int sim=0;
	int window;
	char r1, r2;

	len1=strlen (string1);
	len2=strlen (string2);
	window=atoi(mode);
	len=2*window+1;
	
	if ( len1!=len2)return 0;
	if (window==0 || (window*2+1)>=len1)
	  {
	    sim=get_seq_sim (string1, string2, ignore, "");
	    for (a=0; a<len1; a++)w[a]=sim;
	    return w;
	  }
	

	for ( a=0; a< len1; a++)
		{
		  
		  left =MAX(0, a-window);
		  right=MIN(len1, left+len);
		  for (sim=0,b=left; b<right; b++)
		    {
		      r1=string1[b];
		      r2=string2[b];
		      if (  !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			  if (r1==r2)sim++;
			}
		    }
		  w[a]=(sim*100)/len;
		}
	return w;
	}


int get_seq_sim ( char *string1, char *string2, char *ignore, char *mode)
	{
	int len1;
	int len2;
	int a;
	int pos=0;
	int sim=0;
	char r1, r2;

	len1=strlen (string1);
	len2=strlen (string2);
	
	if ( len1!=len2)return 0;
	
	for ( a=0; a< len1; a++)
		{
		  
		  r1=string1[a];
		  r2=string2[a];
		
		  if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			pos++;
			if (is_in_same_group_aa(r1,r2,0, NULL, mode))
				{
				sim++;
				}
			}
		}

	
	if (pos==0)
		 return 0;
	else	
		return (int) (sim*100)/pos;
	
	}	
int get_seq_sim_2 ( char *string1, char *string2, char *ignore, char **gr, int ng)
	{
	int len1;
	int len2;
	int a;
	int pos=0;
	int sim=0;
	char r1, r2;
	
	
	len1=strlen (string1);
	len2=strlen (string2);
	
	if ( len1!=len2)return 0;
	
	for ( a=0; a< len1; a++)
		{
		r1=string1[a];
		r2=string2[a];
		if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			pos++;
			if (is_in_same_group_aa(r1,r2,ng, gr, NULL))
				{
				sim++;
				}
			}
		}
	
	if (pos==0)
		 return 0;
	else	
		return (int) (sim*100)/pos;
	
	}

int get_seq_sim_3 ( char *string1, char *string2, char *ignore, int **mat)
	{
	int len1;
	int len2;
	int a;
	
	int sim=0;
	char r1, r2;
	
	
	len1=strlen (string1);
	len2=strlen (string2);
	
	if ( len1!=len2)return 0;
	
	for ( a=0; a< len1; a++)
		{
		r1=string1[a];
		r2=string2[a];
		if ( !is_in_set (r1, ignore) && !is_in_set (r2, ignore))
			{
			sim+=mat[r1-'a'][r2-'a'];
			}
		}
	return sim;
	
	}	
int * get_aln_col_weight ( Alignment *A, char *mode)
	{
	int a, b;
	char *col;
	int *weight;
	
	col=vcalloc ( A->nseq, sizeof (int));
	weight=vcalloc (A->len_aln, sizeof (int));
	
	for (a=0; a< A->len_aln; a++)
		{
		for ( b=0; b< A->nseq; b++)
			col[b]=A->seq_al[b][a];
		weight[a]=(find_group_aa_distribution (col, A->nseq,0,NULL,NULL, mode )*100)/A->nseq;		 
		}
	free (col);
	return weight;
	
	}	
			
int analyse_aln_column ( Alignment *B, int col)
    {

    char r=' ';
    int a, b, c;
    static char *mat;
    static int ng_cw_star;
    static char **cw_star;
    int *cw_star_count;
    
    static int ng_cw_col;
    static char **cw_col;
    int *cw_col_count;
    
    static int ng_cw_dot;
    static char **cw_dot;
    int *cw_dot_count;
    
    if ( !B->S)B= get_aln_type (B);
    if ( !B->S || !(B->S)->type)B= get_aln_type (B);
  
    if ( !mat)mat=vcalloc ( STRING, sizeof (char));

    if ( !ng_cw_star)
       {
	   cw_star=make_group_aa ( &ng_cw_star, strcpy ( mat,"idmat"));
	   cw_col=make_group_aa ( &ng_cw_col, strcpy (mat,"clustalw_col"));
	   cw_dot=make_group_aa ( &ng_cw_dot, strcpy (mat, "clustalw_dot"));
       }

    cw_star_count=vcalloc (ng_cw_star, sizeof (int));
    cw_col_count=vcalloc ( ng_cw_col, sizeof (int));
    cw_dot_count=vcalloc (ng_cw_dot, sizeof (int));
    
    for ( a=0; a< B->nseq; a++)
        {
	c=tolower (B->seq_al[a][col]);
	if (is_gap(c)){r=' ';break;}
	  
	for ( b=0; b< ng_cw_star; b++)
	    cw_star_count[b]+=is_in_set (c, cw_star[b]);	
	for ( b=0; b< ng_cw_col; b++)
	    cw_col_count[b]+=is_in_set  (c, cw_col[b]);
	for ( b=0; b< ng_cw_dot; b++)
	    cw_dot_count[b]+=is_in_set  (c, cw_dot[b]);
	}
    
   
   
   
    
    if ( !is_gap(c) && r==' ')
	for ( b=0; b< ng_cw_star; b++)if ( cw_star_count[b]==B->nseq){r='*'; break;}
    if ( !is_gap(c) && r==' ' && !(strm((B->S)->type, "DNA")))
	for ( b=0; b< ng_cw_col ; b++)if ( cw_col_count [b]==B->nseq){r=':'; break;}
    if ( !is_gap(c) && r==' ' && !(strm((B->S)->type, "DNA")))
	for ( b=0; b< ng_cw_dot ; b++)if ( cw_dot_count [b]==B->nseq){r='.'; break;}
    
    
    
    vfree(cw_star_count);
    vfree(cw_col_count);
    vfree(cw_dot_count);   
	
    return r;
    }
int ** get_sim_aln_array ( Alignment *A, char *mode)
	{
	int **w;
	int a, b;
	
	w=declare_int ( A->nseq, A->nseq);
	
	
	for ( a=0; a< A->nseq-1; a++)
		for ( b=a+1; b< A->nseq; b++)
			{
			  if ( strm (mode, "cdna"))
			    w[a][b]=w[b][a]=get_cdna_seq_sim ( A->cdna_cache[0], A->seq_al[a], A->seq_al[b], "-", mode);			
			  else
			    w[a][b]=w[b][a]=get_seq_sim ( A->seq_al[a], A->seq_al[b], "-", mode);
			}
	return w;
	}

int *** get_winsim_aln_array ( Alignment *A,char *mode, int ***w)
	{
	int a, b;
	for ( a=0; a< A->nseq; a++)
		for ( b=0; b< A->nseq; b++)
			{
			  if ( strm (mode, "cdna"))
			    w[a][b]=get_cdna_seq_winsim ( A->cdna_cache[0], A->seq_al[a], A->seq_al[b], "-", mode, w[a][b]);			
			  else
			    w[a][b]=get_seq_winsim ( A->seq_al[a], A->seq_al[b], "-", mode, w[a][b]);
			}
	return w;
	}



int** aln2count_mat ( Alignment *A)
    { /*
	function documentation: start
	
	int output_freq_mat ( char *outfile, Aligmnent *A)

	This function counts the number of residues in each column of an alignment (Prot/NA)
	It outputs these values in the following format

	This format can be piped into:
	The routine used for computing the p-value  gmat-inf-gc-v2c
	
	function documentation: end
      */
      
    int a, b,x;
    int **freq_mat;
    
    freq_mat=declare_int (27, A->len_aln);
      
    for ( a=0; a<A->len_aln; a++)
      {
	for ( b=0; b< A->nseq; b++)
	  {
	    if ( is_gap ( A->seq_al[b][a]))freq_mat[26][a]++;
	    else
	      {
		x=tolower(A->seq_al[b][a]);
		freq_mat[x-'a'][a]++;
	      }
	  }
      }
    return freq_mat;
    }
char *aln2random_seq (Alignment *A, int pn1, int pn2, int pn3, int gn)
    {

      /* 

	 
	 Given the frequencies in A ( read as total counts of each Residue in
	 freq[A->nseq][A->len_aln], and pn1, pn2 and pn3:
	  
	                  1-Generate a new amino-acid at each position
			  2-Insert Gaps, using a HMM.

	 
	 pn3=Weight of the noise induced with sub mat.

	 pn1=% noise type 1 ( Varies with enthropia)
	  n1=Ration noise type 1
	  
	 T =Nseq 
	 t1=Noise 1 expressed in Nseq
	 al=alphabet size;
	 ncat=number of non 0 cat for a given position
	 ICi initial count for residue i

	 Ci=freq[seq][AA]
	 t1=T*n1*(1-1/ncat);
	 t2=T*n2;
	 
	 Ci= ICi*(T-(t1+t2))/T +(t1)/al+(t2)/al
	 
      */
      
      int **freq;
      int **count;
      float T, tot_t1, tot_t2,tot_t3, n1, n2, n3;
      
      
      float ncat;
      
      double gf;
      double *init_freq;
      double *blur_freq;
      double *t1, *t2,*t3;
      int a, b, c, x;
      char *seq;
      int tot;
      /*Viterbi  Parameters */
      
      int p;
      int AL=0;        /*Allowed Transition*/
      int F=-100000; /*Forbiden Transition*/
      int PENALTY=100;
      int GAP_TRANSITION;
      int IGAP=0, IAA=1;
      
      int state,best_state, score, best_score;
      int p_state;
      int e;
      int **score_tab;
      int **state_tab;
      int nstate=2;
      int **transitions;
      
      int max;

      seq=vcalloc ( A->len_aln+1, sizeof (char));     
      count=aln2count_mat(A);
      freq=aln2count_mat(A);

      T=100;

      n1=(float)pn1/100;
      n2=(float)pn2/100;
      n3=(float)pn3/100;
      
      for ( a=0; a< A->len_aln; a++)
	{
	  for ( b=0; b<26; b++)
	    freq[b][a]=freq[b][a]*((T)/(A->nseq-freq[26][a]));
	  freq[26][a]= (freq[26][a]*T)/A->nseq;
	}

      
      init_freq=vcalloc ( 26, sizeof (double));
      blur_freq=vcalloc ( 26, sizeof (double));
      
      tot_t1=tot_t2=tot_t3=0;
      
      t1=vcalloc ( 27, sizeof (double));
      t2=vcalloc ( 27, sizeof (double));
      t3=vcalloc ( 27, sizeof (double));
      for (a=0; a< A->len_aln; a++)
	{

        /*Compute Frequencies*/
	  for (tot=0, b=0; b<26; b++)
		{
		  if ( is_aa(b+'a'))
		    {
		      init_freq[b]=freq[b][a];
		      tot+=freq[b][a];
		    }
		}
        /*Count the number of  different amino acids*/
	for ( ncat=0, b=0; b<=26; b++)
	  {
	    ncat+=(freq[b][a]!=0)?1:0;	    
	  }
	/*Blurr the distribution using */
       	blur_freq=compute_matrix_p (init_freq,tot);	
		
	
	/*compute noise 1: biased with blurred content * enthropy--> keeps prosite motifs*/
	tot_t1=T*n1*(1-1/ncat);
	for (  b=0; b< 26; b++)if ( is_aa(b+'a')){t1[b]=blur_freq[b]*(1-1/ncat)*n1;}
	
	/*Compute noise 2: completely random*/
	tot_t2=T*n2;
	for (  b=0; b< 26; b++)if ( is_aa(b+'a')){t2[b]=tot_t2/21;}
	
	/*compute noise 3: biased with the sole content(pam250mt)*/
	tot_t3=T*n3;
	for (  b=0; b<26; b++)if ( is_aa(b+'a')){t3[b]=blur_freq[b]*n3;}
	
	for ( b=0; b<26; b++)
	  {
	    if ( is_aa('a'+b))
	      freq[b][a]=freq[b][a]*(T-(tot_t1+tot_t2+(tot_t3)))/T+t1[b]+t2[b]+t3[b];
	  }
	
	/*end of the loop that mutates position a*/
	}
     
        vfree (blur_freq);
	vfree (init_freq);
	vfree ( t3);
	
      /*1-Generate the amino acids of the new sequence new*/
      srand(time(NULL));
      
      for ( a=0; a< A->len_aln; a++)
	{

	  for (T=0,b=0; b<26; b++)T+=freq[b][a];
	  x=rand ()%((int)T);
	  for (c=0,b=0; b<26; b++)
	    {
	     c+=freq[b][a];
	     if ( c>=x)
	       {
		 seq[a]='a'+b;
		 c=-1;
		 break;
	       }
	    }
	  if ( c!=-1)seq[a]='-';
	}
      seq[a]='\0';
      

      /*2 Generate the gaps in the new sequence*/
      
      if ( gn<0);
      else
	{

	  transitions=declare_int ( nstate, nstate);
	  score_tab=declare_int ( A->len_aln+2, nstate       );
	  state_tab=declare_int ( A->len_aln+2, nstate       );
	  
	  
	  
	  for (a=0; a<nstate;a++)
	    for (b=0; b<nstate;b++)
	      {transitions[a][b]=F;}
	  
	  GAP_TRANSITION=AL-PENALTY*(rand()%((int)gn+1));
	  
	 



	  transitions[IGAP ][IGAP ]=AL;
	  transitions[IAA][IAA]=AL;
	  transitions[IAA ][IGAP]=GAP_TRANSITION;
	  transitions[IGAP][IAA ]=GAP_TRANSITION; 
	  
	  
	  for ( p=1; p<=A->len_aln; p++){for (state=0; state< nstate; state++){score_tab[p][state]=F;state_tab[p][state]=-1;} }
	  
	  for (p=1; p<= A->len_aln; p++)
	    {
	      for (max=0,a=0; a<26; a++)max=MAX(max, freq[a][p-1]);
	      max=(max*(A->nseq-count[26][p-1]))/A->nseq;
	      
	      for (state=0; state< nstate; state++)
		{
		  
		  
		  gf=freq[26][p-1];
		  if      ( state==IGAP)  e=gf-50;
		  else if ( state==IAA )  e=max-50;
		  for (p_state=0; p_state<nstate; p_state++)
		    {
		      score=(score_tab[p-1][p_state]==F)?F:(e+transitions[p_state][state]+score_tab[p-1][p_state]);
		      if(p_state==0 || score>best_score){ best_score=score;best_state=p_state;}
		    }
		  score_tab[p][state]=best_score;
		  state_tab[p][state]=best_state;
		}
	    }
	  
	  for (state=0; state<nstate; state++)
	    {
	      if (state==0 || score_tab[p-1][state]>best_score){best_score=score_tab[p-1][state]; best_state=state;}
	    }
	  
	  for (p=A->len_aln; p>0;)
	    {
	      if ( best_state==IGAP)
		{
		  seq[p-1]='-';
		}
	      else if ( best_state==IAA)
		{
		  seq[p-1]=seq[p-1];
		}
	      best_state=state_tab[p][best_state];
	      p--;
	    }
        }
      
      free_int (freq, -1);
      return seq;
    }

/********************************************************************/
/*                                                                  */
/*			Weighting functions                         */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/

Alignment * trimseq( Alignment *A, Sequence *S,char *mode)
   {
     Alignment *NA;
     char *p;
     int use_aln=0, percent=0, max_nseq=0, min_sim=0;
     char weight_mode[1000];
     char method[1000];
     int statistics=0;
     
     /*
       mode: 
           (trim)_<seq or aln>_%<percentage of tot weight to keep>_n<number of seq to keep>_w<weight mode>
     */
     
    
     use_aln=aln_is_aligned(A);
     
     
     if ( mode[0]=='\0')
       {
	
	 percent=50;
	 min_sim=0;
	 max_nseq=0;
	 sprintf (weight_mode, "pwsim");
	 sprintf ( method, "clustering2");
       }
     else 
       {
	
	 percent=min_sim=max_nseq;
	 sprintf (weight_mode, "pwsim");
	 sprintf ( method, "clustering2");
       }
    

     while ( (p=strtok(mode, "_")))
	   {	     
	     mode=NULL;
	     if (strm (p, "seq"))use_aln=0;
	     else if ( strm(p,"aln"))use_aln=1;
	     else if  (p[0]=='s')statistics=1;
	     else if  (p[0]=='%')percent=atoi(p+1);
	     else if  (p[0]=='n')max_nseq=atoi(p+1);
	     else if  (p[0]=='N')max_nseq=atoi(p+1)*-1;
	     else if  (p[0]=='w')sprintf (weight_mode, "%s", p+1);
	     else if  (p[0]=='M')sprintf (method, "%s", p+1);
	     else if  (p[0]=='m')min_sim=atoi(p+1);
	   }
     
     if ( !percent && !max_nseq && !min_sim)percent=50;
     
     

     if  (!S)
       {
	 fprintf ( stderr, "\ntrimseq requires a set of sequences[FATAL:%s]\n", PROGRAM);
	 crash("");
       }
     
     else if ( percent >100)
       {
	 percent=100;
	 fprintf ( stderr, "\nWARNING: trimseq: percent(%%) max_val=100%% [Automatic reset]\n");
	 
       }
      else if ( min_sim >100)
       {
	 min_sim=100;
	 fprintf ( stderr, "\nWARNING: trimseq: min_id(m)  max_val=100%% [Automatic reset]\n");
	 
       }
     else if ( max_nseq> S->nseq)
       {
	 max_nseq=S->nseq;
       }
     else if ( max_nseq<0)
       {
	 if ( max_nseq<-100)
	   {
	     fprintf ( stderr, "\nWARNING: trimseq: Nseq(N)  max_val=100%% [Automatic reset]\n");
	     max_nseq=-100;
	   }
	 
	 max_nseq=(int)((float)S->nseq*((float)max_nseq/100)*-1);
       }

    
     if ( strm (method, "clustering1"))
       NA=seq2subseq1 (A, S,use_aln,percent,max_nseq,min_sim,weight_mode);
     else if ( strm (method, "clustering2"))
       NA=seq2subseq2 (A, S,use_aln,percent,max_nseq,min_sim, weight_mode);
     
       
     if ( statistics)
       {
	 fprintf ( stderr, "\nTRIM Informations:\n");
	 fprintf ( stderr, "\tUse...........: %s\n",(use_aln)?"multiple_aln":"pairwise_aln");
	 fprintf ( stderr, "\tcluster_mode..: %s\n"  ,method);
	 fprintf ( stderr, "\tsim_mode......: %s\n"  ,weight_mode);
	 fprintf ( stderr, "\tmin_id........: %d%%\n"  ,(min_sim==0)?-1:min_sim);
	 fprintf ( stderr, "\tpercent.......: %d%%\n",(percent==0)?-1:percent);
	 fprintf ( stderr, "\tn_seq.........: %d (out of %d)\n"  ,NA->nseq, S->nseq);
	 fprintf ( stderr, "\treduction.....: %d%% of original set\n"  ,(NA->nseq*100)/S->nseq);
       }
     
     return NA;
   }
     

Alignment* seq2subseq2( Alignment *A, Sequence *S,int use_aln, int percent,int max_nseq, int ms, char *weight_mode)
{
  float **sim_weight, **max_sim;
  char **seq, **name;
  int *unused_seq_list;
  float sim, small;
  int a, b, best_a, best_b, nchosen, condition1, condition2;
  Sequence *NS;
  Alignment *NA;
  
 
  sim_weight=get_weight   ((use_aln)?A:NULL, S, weight_mode);
  unused_seq_list=vcalloc ( S->nseq, sizeof (int));
  name=declare_char (S->nseq, (MAXNAMES+1));
  seq= declare_char (S->nseq, S->max_len+1);
  max_sim=declare_float (S->nseq, 2);
  
  nchosen=S->nseq;
  /*1 Remove outlayers*/
  if (ms)
    {
      condition1=1;
      condition2=1;
      

      while (condition1 && condition2)
	{
	  for (a=0; a< A->nseq; a++)
	    {
	      max_sim[a][0]=a;
	      max_sim[a][1]=0;
	      if (unused_seq_list[a])
		{
		  max_sim[a][1]=100;
		  continue;
		}

	      for ( b=0; b<A->nseq; b++)
		{
		  if (unused_seq_list[b] || a==b)continue;
	      
		  sim=100-sim_weight[a][b];
		  max_sim[a][1]=MAX(max_sim[a][1],sim);
		}
	    }
	  sort_float (max_sim,2,1,0, S->nseq-1);

	  if (max_sim[0][1]<ms)
	    {
	      best_a=(int)max_sim[0][0];
	    }
	  else
	    {
	      best_a=-1;
	    }
	  if (best_a==-1)condition1=0;
	  else condition1=1;
	  if ((max_nseq && (nchosen-1)<max_nseq) ||(nchosen-1==0)||(best_a==-1))condition2=0;
	  else condition2=1;
	   
	  if ( condition1 && condition2)
	     {	
	       nchosen--;
	       unused_seq_list[best_a]=1;
	     }
	}
    }
	  
   
  if ( !percent && !max_nseq){condition1=condition2=0;}
  else {condition1=condition2=1;}

  
  while (condition1 && condition2)
    {
      best_a=best_b=-1;
      for (small=0, a=0; a< A->nseq-1; a++)
	{
	  if (unused_seq_list[a])continue;
	  
	  for ( b=a+1; b<A->nseq; b++)
	    {
	      if (unused_seq_list[b])continue;
	      
	      sim=100-sim_weight[a][b];
	      if ( sim>small)
		{
		  small=sim;
		  best_a=a;
		  best_b=b;
		}
	    }
	}
      
     
      if ( (percent && small<percent)||(best_a==-1 || best_b==-1))condition1=0;
      else condition1=1;
      
      if ((max_nseq && (nchosen-1)<max_nseq)||(best_a==-1 || best_b==-1))condition2=0;
      else condition2=1;
           
      if ( condition1 && condition2)
	{
	nchosen--;
	unused_seq_list[best_a]=1;
	}
    }

  for (nchosen=0, a=0; a<S->nseq; a++) 
	  {
	    if ( !unused_seq_list[a])
	      {
		sprintf ( name[nchosen], "%s", S->name[a]);
		sprintf ( seq[nchosen] , "%s",(use_aln)?A->seq_al[a]: S->seq[a] );
		nchosen++;
	      }
	  }
  
  
  NS=fill_sequence_struc (nchosen,seq,name);
  NA=seq2aln(NS,NULL,1);	 
  
  if ( use_aln && A)
    {
      NA=realloc_aln2  ( NA,A->max_n_seq,A->len_aln+1);
      for (nchosen=0, a=0; a<S->nseq; a++) 
	{
	  if ( !unused_seq_list[a])
	      {
		sprintf ( NA->seq_al[nchosen] , "%s",A->seq_al[a]);
		nchosen++;
	      }
	}
      NA->len_aln=A->len_aln;
      ungap_aln(NA);
    }
  
  return NA;
}
  
  
      
      
Alignment* seq2subseq1( Alignment *A, Sequence *S,int use_aln, int percent,int max_nseq, int ms,char *weight_mode)
   {
     float **pw_weight,**sim_weight, **seq_weight;
     int a,b,c,d;
     float sum, chosen,last_chosen, last_nchosen,nchosen;
     int condition1, condition2;
     Sequence *NS;
     Alignment *NA;
     char **name, **seq;
     float score, best_score;
     int best_seq;
     int *seq_list, *used_seq_list;
     
     /*
       mode: 
           (trim)_<seq or aln>_%<percentage of tot weight to keep>_n<number of seq to keep>_w<weight mode>
     */
     
     sim_weight=get_weight   ((use_aln)?A:NULL, S, weight_mode);
     pw_weight=declare_float (S->nseq, S->nseq);
     seq_weight=declare_float ( S->nseq, 2);

     
     for (best_score=0,a=0; a<S->nseq; a++)
       {
	 for ( b=0; b<S->nseq; b++)
	   {
	     if ( a==b)continue;
	     seq_weight[a][0]+=sim_weight[a][b];
	   }
	 seq_weight[a][0]=seq_weight[a][0]/(S->nseq-1);
	 score=seq_weight[a][0]=100-seq_weight[a][0];
	 
	 if ( score>best_score)
	   {
	     best_seq=a;
	     best_score=score;
	   }

       }
      for (a=0; a<S->nseq; a++)
       {
	 for ( b=0; b<S->nseq; b++)
	   {
	     if ( a==b)continue;
	     pw_weight[a][b]=sim_weight[a][b]*seq_weight[a][0]*seq_weight[b][0]/(100*100);
	     
	   }
       }
      
 
     seq_list=vcalloc ( S->nseq, sizeof (int));
     used_seq_list=vcalloc ( S->nseq, sizeof (int));

    

     name=declare_char (S->nseq, (MAXNAMES+1));
     seq= declare_char (S->nseq, S->max_len+1);
     
     /*compute the normalization factor*/
     for (sum=0,d=0; d< S->nseq; d++)
	   {
	     for (score=0,c=0; c<S->nseq; c++)
	       {
		 if ( c!=d)
		   score=MAX(score, 100-sim_weight[c][d]);
	       }
	     sum+=score;
	   }
     sum=sum/S->nseq;
     /*chose the first sequence */
     for ( best_score=0,a=0; a< S->nseq; a++)
       {
	 for (score=0, b=0; b< S->nseq; b++)
	   {
	     score+=100-sim_weight[a][b];
	   }
	 if ( score>best_score)
	   {
	     best_seq=a;
	     best_score=score;
	   }
	 
       }


     last_chosen=chosen=((best_score/S->nseq)*100)/sum;
     nchosen=last_nchosen=1;
     seq_list[0]=best_seq;
     used_seq_list[best_seq]=1;

     sprintf ( name[0],"%s", S->name[seq_list[0]]);
     sprintf ( seq[0],"%s", S->seq[seq_list[0]]);
     nchosen=last_nchosen=1;
     

     fprintf ( stderr, "\nTRIM:\n");
     fprintf ( stderr, "\n1-Chosen Sequences\n");
     /*Assemble the list of sequences*/
     for (a=1; a< S->nseq; a++)
       {
	 for (best_score=0,b=0; b< S->nseq; b++)
	   {
	     if (used_seq_list[b]);
	     else
	       {
		 score=pw_weight[seq_list[0]][b]+1;
		 for (c=0; c<a; c++)
		   score=MIN(score,pw_weight[seq_list[c]][b]);
		   
		 if ( score>=best_score)
		   {
		   best_seq=b;
		   best_score=score;
		   }
	       
	       }
	   }
	 seq_list[a]=best_seq;
	 used_seq_list[best_seq]=1;
	 
	 

	 for ( chosen=0,d=0; d< S->nseq; d++)
	   {
	     for (score=0, c=0; c<=a; c++)
	       {
		 if ( seq_list[c]!=d)
		   score=MAX(score, 100-sim_weight[seq_list[c]][d]);
	       }
	     chosen+=score;
	     
	   }
	 
	 chosen=((chosen/S->nseq)*100)/sum;
	 nchosen=a+1;
	 	 
	 condition1= (int)chosen<=(int)percent || !percent;
	 condition2=(nchosen)<=max_nseq     || !max_nseq;
	 
	 if (condition1 && condition2)
	   {
	     fprintf ( stderr, "\tADD %s (set score: %.2f %%)\n", S->name[seq_list[a]], chosen);
	     sprintf ( name[a],"%s", S->name[seq_list[a]]);
	     sprintf ( seq[a],"%s", S->seq[seq_list[a]]);
	     
	   }
	 else
	   {
	     break;	 
	   }
	 last_chosen=chosen;
	 last_nchosen=nchosen;
       }
       
     NS=fill_sequence_struc (last_nchosen,seq,name);
     NA=seq2aln(NS,NULL,1);	 
     fprintf ( stderr, "\n2-Informations:\n");
     fprintf ( stderr, "\tUse...........: %s\n",(use_aln)?"multiple_aln":"pairwise_aln");
     fprintf ( stderr, "\tweight_mode...: %s\n"  ,weight_mode);
     fprintf ( stderr, "\tpercent_weight: %.2f%% (max=%d%%)\n",last_chosen,percent);
     fprintf ( stderr, "\tn_seq.........: %d\n"  ,NS->nseq);
     fprintf ( stderr, "\treduction.....: %d%% of original set\n"  ,(NS->nseq*100)/S->nseq);
     
     return NA;
   }   
float ** get_weight ( Alignment *A, Sequence *S, char *mode)
{
 char *aln_name;
 char *weight_name;
 char *seq_name;
 char command[LONG_STRING];
 char program[LONG_STRING];
 float **weight;
 FILE *fp;
 int c;
  
 if ( !mode || !mode[0] || strm (mode, "msa"))
      {
	if ( getenv ( "SEQ2MSA_WEIGHT")==NULL)sprintf (program, "%s",SEQ2MSA_WEIGHT);
	else sprintf ( program, "%s", (getenv ( "SEQ2MSA_WEIGHT")));   
      }
 else if ( strm(mode, "pwsim"))
      {
	return seq2pwsim (A, S, mode);
      }
 else
     {
       if (getenv (mode))sprintf ( program, "%s", (getenv (mode)));
       else fprintf ( stderr, "\nERROR: %s is not a valid mode for weight computation [FATAL:%s]", mode, PROGRAM);
     }

 /*MSA weights*/
 seq_name=vtmpnam(NULL);
 aln_name=vtmpnam(NULL);
 weight_name=vtmpnam(NULL);
 weight=declare_float (S->nseq+1, 2);
 

  
  if (A)
    {
      output_clustal_aln (seq_name,A);
      output_fasta_seq   (aln_name,A);
      sprintf ( command, "%s %s -i %s -w %s", program, seq_name, aln_name, weight_name);
    }
  else
    {
      A=seq2aln(S,A,1);
      output_fasta_seq   (seq_name,A);
      sprintf ( command, "%s %s -w %s", program, seq_name, weight_name);
    }
  
 
  system ( command);
  
  fp=vfopen( weight_name, "r");
  while ( (c=fgetc(fp))!='$');
  c=fgetc(fp);
  c=0;
  while ( (fscanf (fp, "%*s %f\n",&(weight[c][1])))==1)
    {weight[c][0]=c;c++;}
  vfclose (fp);
  
  
  return weight;
}
float **seq2pwsim (	   Alignment *A, Sequence *S, char *mode)
{
  int a, b, c;
  float d,t;
  float  **W;
  Alignment *B;
  W=declare_float (S->nseq, S->nseq);

  for (a=0; a< S->nseq; a++)
	for ( b=a; b<S->nseq; b++)
	  {
	    if ( a==b){d=1;}
	    else if (!A)
	      {
		B=align_two_sequences ((S)->seq[a], (S)->seq[b],"pam250mt", -10, -1, "fasta_pair_wise");
		
		for (t=0,d=0,c=0; c<B->len_aln; c++)
		  {
		    d+=(B->seq_al[0][c]==B->seq_al[1][c] && !is_gap(B->seq_al[0][c]));
		    t+=(!is_gap(B->seq_al[0][c]) && !is_gap(B->seq_al[1][c]));
		  }
		
		d=d/(t==0)?1:t;
		free_aln(B);
	      }
	    else
	      {
		for (t=0,d=0,c=0; c<A->len_aln; c++)
		  {
		    d+=(A->seq_al[a][c]==A->seq_al[b][c] && !is_gap(A->seq_al[a][c]));
		    t+=(!is_gap(A->seq_al[a][c]) && !is_gap(A->seq_al[b][c]));
		  }
		d/=(t==0)?1:t;
	      }
	  
	    W[a][b]=W[b][a]=(1-d)*100;
	  }
 
  
  return W;

}
  

/********************************************************************/
/*                                                                  */
/*			AMINO ACID FUNCTIONS                        */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/

char** make_group_aa (int *ngroup, char *mode)
	{
/*mode:         indicates which matrix will be used for the grouping*/
/*n_group:      pointer to the number of groups                     */
/*return value: an array of strings containing the AA of each group */


	int **matrix;
	int a, b,c,is_in;
	char buf[28];
	char **group_list;
	char *matrix_name;
	matrix_name=vcalloc ( 100, sizeof (char));

	
	ngroup[0]=0;
	group_list=declare_char ( 100, 27);
		
	if ( mode[0]=='_'){mode++;sprintf ( matrix_name, "%s", mode);}

	if (mode==NULL || mode[0]=='\0')sprintf ( matrix_name, "idmat");
	else if ( strm (mode, "simple"))
	     {
		 ngroup[0]=6;
		 sprintf ( group_list[0], "avilm");
		 sprintf ( group_list[1], "dekr");
		 sprintf ( group_list[2], "stcnqh");
		 sprintf ( group_list[3], "wfy");
		 sprintf ( group_list[4], "g");
		 sprintf ( group_list[5], "p");
		 vfree (matrix_name);
		 return group_list;
	     }
	else if ( strm (mode, "clustalw"))
	     {
		 ngroup[0]=11;
		 sprintf ( group_list[0] ,"asta");
		 sprintf ( group_list[1] ,"bneqk");
		 sprintf ( group_list[2] ,"cnhqk");
		 sprintf ( group_list[3] ,"dndeq");
		 sprintf ( group_list[4] ,"eqhrk");
		 sprintf ( group_list[5] ,"fmilv");
		 sprintf ( group_list[6] ,"gmilf");
		 sprintf ( group_list[7] ,"hhy");
		 sprintf ( group_list[8] ,"ifyw");
		 sprintf ( group_list[9] ,"jc");
		 sprintf ( group_list[10],"kp");
		 vfree (matrix_name);
		 return group_list;
	     }
	else if ( strm (mode, "polarity"))
	     {
		 ngroup[0]=6;
		 sprintf ( group_list[0], "eqrsdnkht");
		 sprintf ( group_list[1], "p");
		 sprintf ( group_list[2], "g");
		 sprintf ( group_list[3], "c");
		 sprintf ( group_list[4], "fyw");
		 sprintf ( group_list[5], "iavlm");
		 vfree (matrix_name);
		 return group_list;
	     }
	else if ( strm (mode, "vasiliky"))
	     {
		 ngroup[0]=0;
		 sprintf ( group_list[ngroup[0]++], "rk");
		 sprintf ( group_list[ngroup[0]++], "de");
		 sprintf ( group_list[ngroup[0]++], "qh");
		 sprintf ( group_list[ngroup[0]++], "vilm");
		 sprintf ( group_list[ngroup[0]++], "fy");
		 sprintf ( group_list[ngroup[0]++], "s");
		 sprintf ( group_list[ngroup[0]++], "w");
		 sprintf ( group_list[ngroup[0]++], "a");
		 sprintf ( group_list[ngroup[0]++], "c");
		 sprintf ( group_list[ngroup[0]++], "g");
		 sprintf ( group_list[ngroup[0]++], "n");
		 sprintf ( group_list[ngroup[0]++], "p");
		 sprintf ( group_list[ngroup[0]++], "t");
		 vfree (matrix_name);
		 return group_list;
	     }
	else if ( strm (mode, "clustalw_col"))
	     {
		 sprintf ( group_list[ngroup[0]++], "sta");
		 sprintf ( group_list[ngroup[0]++], "neqk");
		 sprintf ( group_list[ngroup[0]++], "nhqk");
		 sprintf ( group_list[ngroup[0]++], "ndeq");
		 sprintf ( group_list[ngroup[0]++], "qhrk");
		 sprintf ( group_list[ngroup[0]++], "milv");
		 sprintf ( group_list[ngroup[0]++], "milf");
		 sprintf ( group_list[ngroup[0]++], "hy");
		 sprintf ( group_list[ngroup[0]++], "fyw");
		 sprintf ( group_list[ngroup[0]++], "g");
		 sprintf ( group_list[ngroup[0]++], "p");
		 sprintf ( group_list[ngroup[0]++], "c");
		 vfree (matrix_name);
		 
		 return group_list;
	     }
	else if ( strm (mode, "clustalw_dot"))
	     {
		 sprintf ( group_list[ngroup[0]++], "csa");
		 sprintf ( group_list[ngroup[0]++], "atv");
		 sprintf ( group_list[ngroup[0]++], "sag");
		 sprintf ( group_list[ngroup[0]++], "stnk");
		 sprintf ( group_list[ngroup[0]++], "stpa");
		 sprintf ( group_list[ngroup[0]++], "sgnd");
		 sprintf ( group_list[ngroup[0]++], "sndeqk");
		 sprintf ( group_list[ngroup[0]++], "ndeqhk");
		 sprintf ( group_list[ngroup[0]++], "neqhrk");
		 sprintf ( group_list[ngroup[0]++], "fvlim");
		 sprintf ( group_list[ngroup[0]++], "hfy");
		 vfree (matrix_name);
		 return group_list;
	     }
	else if ( strm (mode, "make_all"))
	     {
		 ngroup[0]=1;
		 sprintf ( group_list[0], "abcdefghijklmnopqrstuvwxyz");
		 vfree (matrix_name);
		 return group_list;
	     }
	else sprintf ( matrix_name, "%s", mode);


	matrix=read_matrice ( matrix_name); 
	
	for ( a=0;a< 26; a++)
		{
		if ( matrix[a][a]>0)
			{
			for ( c=0,b=0;b< 26; b++)
				{
				
				if ( matrix[a][b]>0 && matrix[b][b]>0)
					{
					buf[c++]=b+'a';
					}
				}
			buf[c]='\0';
			for ( is_in=0,b=0; b< ngroup[0]; b++)if ( strcmp (buf, group_list[b])==0)is_in=1;
			if (is_in==0)sprintf ( group_list[ngroup[0]++], "%s", buf);
					
			}
		}
	free_int (matrix, -1); 
	vfree (matrix_name);
		 
	return group_list;
	}

int find_group_aa_distribution (char *col, int nseq,int n_group, char **gl,  int *distrib, char *mode )
	{
	static int *distribution;
	static char **lgl;
	static int ln_group;
	int a, b, c;
	int *d;
	char **gl2;
	int n_group2;
	
	
	
	if ( lgl==NULL)
		lgl=make_group_aa ( &ln_group, mode);
		
	if ( gl==NULL)
		{
		gl2=lgl;
		n_group2=ln_group;
		}
	else
		{
		gl2=gl;
		n_group2=n_group;
		}
	
	if ( distribution==NULL || ln_group<n_group)distribution=vcalloc ( n_group2, sizeof (int));
	if ( distrib==NULL)d=distribution;
	else d=distrib;
	
	
	for ( a=0; a< n_group2; a++)d[a]=0;
	
	for ( a=0; a< nseq; a++)
		{
		for ( b=0; b< n_group2; b++)
			d[b]+=is_in_set (col[a], gl2[b]);
		}
	c=d[0];
	for ( a=0; a< n_group2; a++)
		c=(d[a]>c)?d[a]:c;
	return c;
	}



int is_in_same_group_aa ( char r1, char r2, int n_group, char **gl, char *mode)
	{
	int a;
	static char **lgl;
	static int ln_group;
	
	char **gl2;
	int n_group2;
	
	/*use mode=idmat for similarity based on id*/

	if ( strm (mode, "clean"))
	     {
	     free_char (lgl, -1);
	     lgl=NULL;
	     ln_group=0;
	     return 0;
	     }
	
	if ( lgl==NULL)
		lgl=make_group_aa ( &ln_group, mode);
		
	if ( gl==NULL)
		{
		gl2=lgl;
		n_group2=ln_group;
		}
	else
		{
		gl2=gl;
		n_group2=n_group;
		}
	
	for ( a=0; a< n_group2; a++)
		if ( is_in_set ( r1, gl2[a]) && is_in_set ( r2, gl2[a]))return 1;
		return 0;
	}
			

Alignment * gene2prot (Alignment *A){return A; }
char * test_gene2prot (Constraint_list *CL, int s1)
       {
	   int a, b,q, nal;
	   int F=-10000000; /*FORBIDEN STATE*/
	   int AL=0;       /*ALLOWED STATE*/
	   int SPLICE_PENALTY=1000;
	   int FRAME_PENALTY=1000;
	   

	   int START,  ORF1,  ORF2, ORF3, s5NC; 
	   int s3NC,ORF3_G1, ORF3_T2, ORF3_NC, ORF3_A3, ORF3_T4;
	   int U1_G1,   U1_T2,   U1_NC,   U1_A3,   U1_T4;
	   int U2_G1,   U2_T2,   U2_NC,   U2_A3,   U2_T4;
	   int U1,     U2, U3,  U4, U5, END;
	   
	   int nstate=0;
	   int **transitions;
	   int **v_tab;
	   int **v_tab_p;
	   int **last_coding;
	   int **last_t4;
	   int *potential;
	   int v;

	   int orf1, orf2, orf3, ncp, p, state, pstate, e, best_state_p, best_state_v, best_pstate_p, best_pstate_v;
	   char *seq, *seq2, *seq3;
	   int l;
	   int  *is_coding;
	   int *is_t4;
	   char *codon;

	   static int *entry;
	   int tot=0;
	   
	   seq=vcalloc ( strlen ((CL->S)->seq[s1])+1, sizeof (char));
	   seq2=vcalloc ( strlen ((CL->S)->seq[s1])+1, sizeof (char));
	   seq3=vcalloc ( strlen ((CL->S)->seq[s1])+1, sizeof (char));
	   sprintf ( seq, "%s", (CL->S)->seq[s1]);
	   ungap (seq);

	   l=strlen (seq);
	   for ( a=0; a< l; a++) seq[a]=tolower ( seq[a]);
	   for ( a=0; a< l; a++) seq[a]=(seq[a]=='t')?'u': seq[a];
	   
	   fprintf ( stderr, "\n%s", seq);
	   potential=vcalloc (l+1, sizeof (int));
	   CL=index_constraint_list ( CL);
	   for (nal=0, a=0; a<(CL->S)->nseq; a++)
	       for ( b=CL->start_index[s1][a]; b< CL->end_index[s1][a];b++)
	           {
		       entry=extract_entry(entry, b, CL);
		       if ( entry[SEQ1]==s1)potential[entry[R1]-1]+=entry[WE];
		       else if (  entry[SEQ2]==s1)potential[entry[R2]-1]+=entry[WE];
		       tot+=entry[WE];
		       nal++;
		   }


	   SPLICE_PENALTY=10000;
	   FRAME_PENALTY=1000;

	   
	   nstate=0;
	   START=nstate++;  ORF1=nstate++;  ORF2=nstate++; ORF3=nstate++; s5NC=nstate++; 
	   s3NC=nstate++;
	   ORF3_G1=nstate++;U1_G1=nstate++;U2_G1=nstate++; 
	   ORF3_T2=nstate++;U1_T2=nstate++;U2_T2=nstate++;
	   ORF3_NC=nstate++;U1_NC=nstate++;U2_NC=nstate++; 
	   ORF3_A3=nstate++;U1_A3=nstate++;U2_A3=nstate++; 
	   ORF3_T4=nstate++;U1_T4=nstate++;U2_T4=nstate++;
	   
	   
	   U1=nstate++;     U2=nstate++; U3=nstate++;  U4=nstate++; U5=nstate++;  
	   END=nstate++;
	   
	   is_coding=vcalloc ( nstate, sizeof (int));
	   is_coding[ORF1]=is_coding[ORF2]=is_coding[ORF3]=is_coding[U1]=is_coding[U2]=1;
	   is_coding[U3]=is_coding[U4]=is_coding[U5]=1;
	   
	   is_t4=vcalloc ( nstate, sizeof (int));
	   is_t4[ORF3_T4]=is_t4[U1_T4]=is_t4[U2_T4]=1;
	   transitions=declare_int ( nstate, nstate);
	   for (a=0; a< nstate; a++)
		   for ( b=0; b< nstate; b++)transitions[a][b]=F;
	   
	   transitions[START][ORF1]=AL;
	   transitions[START][s5NC]=AL-FRAME_PENALTY;
	   transitions[s5NC][s5NC]=AL;

	   transitions[s5NC][ORF1]=AL-FRAME_PENALTY;

	   transitions[ORF1][ORF2]=AL;
	   transitions[ORF2][ORF3]=AL;
	   transitions[ORF3][U1]=AL;
	   transitions[ORF3][ORF1]=AL;
	   transitions[ORF3][ORF3_G1]=AL-SPLICE_PENALTY;
	   
	   
	   transitions[ORF3_G1][ORF3_T2]=AL;
	   transitions[ORF3_T2][ORF3_NC]=AL;
	   transitions[ORF3_NC][ORF3_NC]=AL;
	   transitions[ORF3_NC][ORF3_A3]=AL;
	   transitions[ORF3_A3][ORF3_T4]=AL;
	   transitions[ORF3_T4][ORF1]=AL-SPLICE_PENALTY;

	   transitions[U1][U2]=AL;
	   transitions[U1][U1_G1]=AL-SPLICE_PENALTY;
	   transitions[U1_G1][U1_T2]=AL;
	   transitions[U1_T2][U1_NC]=AL;
	   transitions[U1_NC][U1_NC]=AL;
	   transitions[U1_NC][U1_A3]=AL;
	   transitions[U1_A3][U1_T4]=AL;
	   transitions[U1_T4][U3]=AL-SPLICE_PENALTY;
	   transitions[U3][U4]=AL;
	   transitions[U4][ORF1]=AL;
	   
	   transitions[U2][U2_G1]=AL-SPLICE_PENALTY;
	   transitions[U2_G1][U2_T2]=AL;
	   transitions[U2_T2][U2_NC]=AL;
	   transitions[U2_NC][U2_NC]=AL;
	   transitions[U2_NC][U2_A3]=AL;
	   transitions[U2_A3][U2_T4]=AL;
	   transitions[U2_T4][U5]=AL-SPLICE_PENALTY;
	   transitions[U5][ORF1]=AL;
	   
	   transitions[ORF3][s3NC]=AL-FRAME_PENALTY;
	   transitions[ORF3][END]=AL;
	   transitions[s3NC][END]=AL;

	   	   
	   v_tab=declare_int ( l+1,nstate);
	   v_tab_p=declare_int ( l+1,nstate);
	   last_coding=declare_int ( l+1,nstate);
	   last_t4=declare_int ( l+1,nstate);

	   for (a=0; a< l; a++) potential[a]-=200;

	   codon=vcalloc ( 4, sizeof (char));
	   best_pstate_p=START;
	   best_pstate_v=0;
	   nal=0;
	   for ( p=1; p<=l; p++)
	       {
	       if  (translate_dna_codon (seq+(p-1), 'x')=='x' || p>(l-2))orf1=F;
	       else orf1=potential[p-1];

	       if  (p<2 || translate_dna_codon (seq+(p-2), 'x')=='x' || p>(l-1))orf2=F;
	       else orf2=potential[p-1];

	       
	       if  (p<3 || translate_dna_codon (seq+(p-3), 'x')=='x' || p>l)orf3=F;
	       else orf3=potential[p-1];
	       
	       if ( best_int (3, 1, &a, orf1, orf2, orf3)!=F)ncp=-best_int (3, 1, &a, orf1, orf2, orf3);
	       else ncp=1000;
	       
	       for ( state=0; state< nstate; state++)
	           {
		       
		       if      ( state==ORF1)e=orf1;
		       else if ( state==ORF2)e=orf2;
		       else if ( state==ORF3)e=orf3;
		       else if ( state>=U1 && state<=U3)
			   {
			   e=0;
			   }
		       else if ( state==U4)
		          {
			      codon[2]=seq[p-1];
			      codon[1]=seq[last_coding[p-1][U3]-1];
			      codon[0]=seq[last_coding[p-2][U1_T4]-1];
			      if ( translate_dna_codon (codon, 'x')=='x')e=F;
			      else e=0;
			  }
		       else if ( state==U5)
		          {
			      codon[2]=seq[p-1];
			      codon[1]=seq[last_coding[p-1][U2_T4]-1];
			      q=seq[last_coding[p-1][U2_T4]];
			      codon[0]=seq[last_coding[q-1][U1]-1];
			      if ( translate_dna_codon (codon, 'x')=='x')e=F;
			      else e=0;
			  }

		       else if (state>=ORF3_G1 && state<=U2_G1)e=(p<l-1 && seq[p-1]=='g' && seq[p]=='u')?ncp:F;
		       else if ( state>=ORF3_T2 && state<=U2_T2)
			   {
			   e=(p>1 && seq[p-2]=='g' && seq[p-1]=='u')?ncp:F;
			   }
		       else if ( state>=ORF3_A3 && state<=U2_A3)e=(seq[p-1]=='a')?ncp:F;
		       else if ( state>=ORF3_T4 && state<=U2_T4)e=(seq[p-1]=='u')?ncp:F;
		       else e=ncp;
		       
		       for ( pstate=0; pstate<nstate; pstate++)
		           {
			       if (e==F ||  transitions[pstate][state]==F || v_tab[p-1][pstate]==F)v=F;
			       else v=e+transitions[pstate][state]+v_tab[p-1][pstate];
				    
			       if ( pstate==0 || v>best_pstate_v)
			          {best_pstate_v=v;best_pstate_p=pstate;}
			   }
		      v_tab[p][state]=best_pstate_v;
		      v_tab_p[p][state]=best_pstate_p; 
		      
		      if (!is_coding[state])last_coding[p][state]=last_coding[p-1][best_pstate_p];
		      else if (is_coding[state])last_coding[p][state]=p;
		      
		      if (!is_t4[state])
		         {
			     if (is_coding[state] && last_t4[p-1][best_pstate_p]==0)last_t4[p][state]=p;
			     else last_t4[p][state]=last_t4[p-1][best_pstate_p];
			 }
		      else if (is_t4[state])last_t4[p][state]=p;
		      
		      if (state==0 ||best_pstate_v>best_state_v ){best_state_p=state; best_state_v=best_pstate_v;}
		   }
	       }
	   tot=0;
	   for ( p=l; p>0; p--)
	           {
		       if ( best_state_p>=ORF1 &&  best_state_p<=ORF3){seq2[tot++]=tolower (seq[p-1]);}
		       else if ( best_state_p>=U1 && best_state_p<=U5){seq2[tot++]=tolower (seq[p-1]);}
		       if (best_state_p==ORF1)seq[p-1]=toupper (seq[p-1]);
		       else if (best_state_p==ORF2 || best_state_p==ORF3)seq[p-1]=tolower (seq[p-1]);
		       else if ( best_state_p==ORF3_NC || best_state_p==U1_NC ||  best_state_p==U2_NC) seq[p-1]='.';
		       else if ( best_state_p==U1 || best_state_p==U2 || best_state_p==U3 || best_state_p==U4 || best_state_p==U5) seq[p-1]=best_state_p-U1+'1';
		       else seq[p-1]=toupper (seq[p-1]);
		       best_state_p=v_tab_p[p][best_state_p];
		   }

	   for ( a=0, b=tot-1; b>=0; b--, a++)
	       seq3[a]=seq2[b];
	   
	   fprintf ( stderr, "\n%s\n", seq);
	   fprintf ( stderr, "\nN coding=%d\n", tot);
	   for ( a=0; a< tot; a+=3)
	        {
		b=translate_dna_codon (seq3+a, 'x');
		fprintf ( stderr, "%c",b);
		if ( b=='x'){fprintf ( stderr, "\n");exit (1);}
		}
	    
	    fprintf ( stderr, "\n");	  
	    exit (1);
	    return 0;
	    
	       
	    
       }
Alignment * dna_aln2_3frame_cdna_aln(Alignment *A,int *ns,int **l_s)
{
  Alignment *B;
  int a;
  B=realloc_aln2 (NULL,6,strlen(A->seq_al[l_s[0][0]])+strlen(A->seq_al[l_s[1][0]]));
  for ( a=0; a< 3; a++)  
    {	 
      B->seq_al[a]=translate_dna_seq (A->seq_al[l_s[0][0]]+a, 0, 'o',B->seq_al[a]);
      B->seq_al[a+3]=translate_dna_seq (A->seq_al[l_s[1][0]]+a, 0, 'o',B->seq_al[a+3]);
    }
  for ( a=1; a<3; a++)
    {
      if ( strlen(B->seq_al[a])<strlen(B->seq_al[0])) B->seq_al[a]=strcat ( B->seq_al[a], "x");
      if ( strlen(B->seq_al[a+3])<strlen(B->seq_al[3])) B->seq_al[a+3]=strcat ( B->seq_al[a+3], "x");
    }
  
  B->nseq=6;
  B->len_aln=strlen (B->seq_al[0]);
  return B;
}
