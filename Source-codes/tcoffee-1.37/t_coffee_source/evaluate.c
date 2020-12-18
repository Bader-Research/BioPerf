#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

#include "dp_lib_header.h"
float compute_lambda (int **matrix,char *alphabet);
/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR EVALUATING THE CONSISTENCY BETWEEN ALN AND CL                       */
/*                                                                                           */
/*********************************************************************************************/

/*Fast:         score= extended_res/max_extended_residue_for the whole aln
  slow:         score= extended_res/sum all extended score for that residue
  non_extended  score= non_ext     /sum all non extended  score for that residue
  heuristic     score= extended    /sum of extended score of all pairs in the library
                                    (i.e. Not ALL the possible pairs)
*/				
Alignment * main_coffee_evaluate_output ( Alignment *IN,Constraint_list *CL, const char *mode )
{


   if ( CL->evaluate_residue_pair==evaluate_matrix_score)
     {
       return copy_aln (IN, NULL);
     }
   else if ( strm4 ( mode, "t_coffee_fast","tcoffee_fast","fast_tcoffee", "fast_t_coffee"))
     {
       return fast_coffee_evaluate_output ( IN,CL);
     }
   else if ( strm4 ( mode, "t_coffee_slow","tcoffee_slow","slow_tcoffee","slow_t_coffee" ))
     {
       return slow_coffee_evaluate_output ( IN,CL);
     }
   
   else if ( strm4 ( mode, "tcoffee_non_extended","t_coffee_non_extended","non_extended_tcoffee","non_extended_t_coffee"))
     {
       return non_extended_t_coffee_evaluate_output ( IN,CL);
     }
   else if ( strm5 ( mode, "tcoffee_heuristic","t_coffee_heuristic","heuristic_tcoffee","heuristic_t_coffee", "dali"))
     {
       return heuristic_coffee_evaluate_output ( IN,CL);
     }
   else 
     {
       fprintf ( stderr, "\nUNKNOWN MODE FOR ALIGNMENT EVALUATION: *%s* [FATAL:%s]",mode, PROGRAM);
       crash ("");
       return NULL;
     }
   
}
Alignment * coffee_evaluate_output ( Alignment *IN,Constraint_list *CL)
    {
    fprintf ( stderr, "\n[WARNING:%s]THE FUNCTION coffee_evaluate_output IS NOT ANYMORE SUPPORTED\n", PROGRAM);
    fprintf ( stderr, "\n[WARNING]fast_coffee_evaluate_output WILL BE USED INSTEAD\n");
    
    return fast_coffee_evaluate_output (IN,CL);
    }

Alignment * fast_coffee_evaluate_output ( Alignment *IN,Constraint_list *CL)
    {
    int a,b, c, m,res, s, s1, s2, r1, r2;	
    Alignment *OUT=NULL;
    int **pos;

    double score_col=0, score_aln=0, score_res=0;
    int max_col, max_aln;
    int **score_pw_seq, *max_seq, *score_seq;
    int local_m;
    int local_nseq;
    int **seq_count;
    int *entry=NULL;

    

    seq_count=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
    

    if ( !CL->evaluate_residue_pair){fprintf ( stderr, "\nWARNING: CL->evaluate_residue_pair Not set\nSet to: extend_residue_pair\n");CL->evaluate_residue_pair= extend_residue_pair; }
	
    OUT=copy_aln (IN, OUT);
    pos=aln2pos_simple(IN, IN->nseq);


    score_pw_seq=declare_int ( IN->nseq, IN->nseq);
    max_seq=vcalloc ( IN->nseq, sizeof (int));
    score_seq=vcalloc ( IN->nseq, sizeof (int));
    
    for (a=0; a< CL->ne; a++)
      {
	entry=extract_entry(entry, a, CL);
	seq_count[entry[SEQ1]][entry[R1]]++;
	seq_count[entry[SEQ2]][entry[R2]]++;
      }
    
    /*1: Identify the highest scoring pair within the alignment*/
 
    for ( m=0, a=0; a< IN->len_aln; a++)
        {
	    for ( b=0; b< IN->nseq; b++)
		{
		s1=IN->order[b][0];
		r1=pos[b][a];
	

		for ( c=0; c< IN->nseq; c++)
		    {
		    s2=IN->order[c][0];
		    r2=pos[c][a];
		    if ( s1==s2 && !CL->do_self)continue;
	
		    if ( s1< s2)s=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
		    else        s=(CL->evaluate_residue_pair)( CL, s2, r2, s1, r1);
		    s=(s!=UNDEFINED)?s:0;
		    m=MAX(m, s);
		    }
		}
	}
    
    local_m=m;
    

    sprintf ( OUT->name[IN->nseq], "Cons");
    for ( max_aln=0,score_aln=0,a=0; a< IN->len_aln; a++)
	{
	OUT->seq_al[IN->nseq][a]=NO_COLOR_RESIDUE;
	for ( local_nseq=0,b=0; b<IN->nseq; b++){local_nseq+=(pos[b][a]>0)?1:0;}
	local_m=m*(local_nseq-1);
	
	for ( max_col=0, score_col=0,b=0; b< IN->nseq; b++)
	    {
	    OUT->seq_al[b][a]=NO_COLOR_RESIDUE;
	    s1=IN->order[b][0];
	    r1=pos[b][a];
	    
	    
	    if (r1<=0 || seq_count[s1][r1]==0)continue;
	    
	    for ( score_res=0,c=0; c< IN->nseq; c++)
	        {
		    s2=IN->order[c][0];
		    r2=pos[c][a];		    
		    
		    if ((s1==s2 && !CL->do_self) || r2<=0)continue;	
		    max_col   +=m;
		    max_seq[b]+=m;
		    max_aln   +=m;
		   
		    if ( s1< s2)s=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
		    else        s=(CL->evaluate_residue_pair)( CL, s2, r2, s1, r1);
		    s=(s!=UNDEFINED)?s:0;
		    
		    score_res+=s;	
		    score_col+=s;		    
		    score_seq[b]+=s;
		    score_pw_seq[b][c]+=s;
		    score_aln+=s;		    
		}
	    
	    res=(local_m==0)?NO_COLOR_RESIDUE:((score_res*10)/local_m);
	    
	    (OUT)->seq_al[b][a]=(res==NO_COLOR_RESIDUE)?res:(MIN(res, 9));	    
	    
	    }
       
	res=(max_col==0)?NO_COLOR_RESIDUE:((score_col*10)/max_col);	
	OUT->seq_al[IN->nseq][a]=(res==NO_COLOR_RESIDUE)?res:(MIN(res,9));
	}    
   
    IN->score_aln=OUT->score_aln=(max_aln==0)?0:((score_aln*100)/max_aln);
    
    for ( a=0; a< OUT->nseq; a++)
	{
	OUT->score_seq[a]=(max_seq[a]==0)?0:((score_seq[a]*100)/max_seq[a]);
	for ( b=0; b< OUT->nseq; b++)
	    {
		if ( a==b)continue;
		OUT->pw_score_seq[a][b]=((max_seq[a]+max_seq[b])==0)?0:((score_pw_seq[a][b]*200)/(max_seq[a]+max_seq[b]));
	    }
	}
	
    free_int (score_pw_seq, -1);
    free ( score_seq);
    free ( max_seq);
    free(entry);
    free_int(seq_count, -1);
   

    return OUT;
    }
Alignment * slow_coffee_evaluate_output ( Alignment *IN,Constraint_list *CL)
    {
    int a,b, c,res, s, s1, s2, r1, r2;	
    Alignment *OUT=NULL;
    int **pos;
    int max_score_r, score_r;
    double score_col=0, score_aln=0;
    int max_score_col, max_score_aln;
    int *max_score_seq, *score_seq;
    int **tot_extended_weight;
    int **res_extended_weight;
    int n_res_in_col;

    /*
      Residue x: sum of observed extended X.. /sum of possible X..
    */

    

    if ( !CL->evaluate_residue_pair){fprintf ( stderr, "\nWARNING: CL->evaluate_residue_pair Not set\nSet to: extend_residue_pair\n");CL->evaluate_residue_pair= extend_residue_pair; }
	
    OUT=copy_aln (IN, OUT);
    pos=aln2pos_simple(IN, IN->nseq);


    max_score_seq=vcalloc ( IN->nseq, sizeof (int));
    score_seq=vcalloc ( IN->nseq, sizeof (int));
    
    tot_extended_weight=list2residue_total_extended_weight(CL);
    res_extended_weight=declare_int ((CL->S)->nseq, (CL->S)->max_len+1);
         
    for (a=0; a< IN->len_aln; a++)
        {
	    for ( b=0; b< IN->nseq-1; b++)
		{
		s1=IN->order[b][0];
		r1=pos[b][a];
		for ( c=b+1; c< IN->nseq; c++)
		    {
		    s2=IN->order[c][0];
		    r2=pos[c][a];	
		    if ( s1==s2 && !CL->do_self)continue;
		    else if ( r1<=0 || r2<=0)   continue;		    
		    else 
		      {
			if ( s1< s2)s=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
			else        s=(CL->evaluate_residue_pair)( CL, s2, r2, s1, r1);
			res_extended_weight[s1][r1]+=s;
			res_extended_weight[s2][r2]+=s;
		      }
		    }
		}
	}
        
  
    sprintf ( OUT->name[IN->nseq], "Cons");
    for ( max_score_aln=0,score_aln=0,a=0; a< IN->len_aln; a++)
	{
	OUT->seq_al[IN->nseq][a]=NO_COLOR_RESIDUE;
	for ( n_res_in_col=0,b=0; b<IN->nseq; b++){n_res_in_col+=(pos[b][a]>0)?1:0;}
	for ( max_score_col=0, score_col=0,b=0; b< IN->nseq; b++)
	    {
	    OUT->seq_al[b][a]=NO_COLOR_RESIDUE;
	    s1=IN->order[b][0];
	    r1=pos[b][a];
	    if (r1<=0)continue;
	    else
	      {
		max_score_r  =tot_extended_weight[s1][r1];
		score_r=res_extended_weight[s1][r1];
		res=(max_score_r==0 ||n_res_in_col<2 )?NO_COLOR_RESIDUE:((score_r*10)/max_score_r);
		(OUT)->seq_al[b][a]=(res==NO_COLOR_RESIDUE)?res:(MIN(res, 9));
		max_score_col+=max_score_r;
		    score_col+=score_r;
		max_score_seq[b]+=max_score_r;
		    score_seq[b]+=score_r;
		max_score_aln+=max_score_r;
		    score_aln+=score_r;
	      }
	    res=(max_score_col==0 || n_res_in_col<2)?NO_COLOR_RESIDUE:((score_col*10)/max_score_col);	
	    OUT->seq_al[IN->nseq][a]=(res==NO_COLOR_RESIDUE)?res:(MIN(res,9));
	    }
	}
    IN->score_aln=OUT->score_aln=(max_score_aln==0)?0:((score_aln*100)/max_score_aln);
    for ( a=0; a< OUT->nseq; a++)
      {
	OUT->score_seq[a]=(max_score_seq[a]==0)?0:((score_seq[a]*100)/max_score_seq[a]);
	for ( b=0; b< OUT->nseq; b++)
	  {
	    OUT->pw_score_seq[a][b]+=OUT->score_seq[b];
	  }
	OUT->pw_score_seq[a][b]=OUT->pw_score_seq[a][b]/OUT->nseq;
	
      }

    
    free ( score_seq);
    free ( max_score_seq);
    
    free_int (tot_extended_weight, -1);
    free_int (res_extended_weight, -1);
    
    
    return OUT;
    }

Alignment * heuristic_coffee_evaluate_output ( Alignment *IN,Constraint_list *CL)
    {
    int a,b, c,res, s, s1, s2, r1, r2;	
    Alignment *OUT=NULL;
    int **pos;
    int max_score_r, score_r;
    double score_col=0, score_aln=0;
    int max_score_col, max_score_aln;
    int *max_score_seq, *score_seq;
    int **tot_extended_weight;
    int **res_extended_weight;
    int n_res_in_col;

    /*
      Residue x: sum of observed extended X.. /sum of possible X..
    */

    if ( !CL->evaluate_residue_pair){fprintf ( stderr, "\nWARNING: CL->evaluate_residue_pair Not set\nSet to: extend_residue_pair\n");CL->evaluate_residue_pair= extend_residue_pair; }
	
    OUT=copy_aln (IN, OUT);
    pos=aln2pos_simple(IN, IN->nseq);


    max_score_seq=vcalloc ( IN->nseq, sizeof (int));
    score_seq=vcalloc ( IN->nseq, sizeof (int));
    
    tot_extended_weight=list2residue_partial_extended_weight(CL);
    res_extended_weight=declare_int ((CL->S)->nseq, (CL->S)->max_len+1);
         
    for (a=0; a< IN->len_aln; a++)
        {
	    for ( b=0; b< IN->nseq-1; b++)
		{
		s1=IN->order[b][0];
		r1=pos[b][a];
		for ( c=b+1; c< IN->nseq; c++)
		    {
		    s2=IN->order[c][0];
		    r2=pos[c][a];	
		    if ( s1==s2 && !CL->do_self)continue;
		    else if ( r1<=0 || r2<=0)   continue;		    
		    else 
		      {
			if ( s1< s2)s=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
			else        s=(CL->evaluate_residue_pair)( CL, s2, r2, s1, r1);
			res_extended_weight[s1][r1]+=s;
			res_extended_weight[s2][r2]+=s;
		      }
		    }
		}
	}
        
  
    sprintf ( OUT->name[IN->nseq], "Cons");
    for ( max_score_aln=0,score_aln=0,a=0; a< IN->len_aln; a++)
	{
	OUT->seq_al[IN->nseq][a]=NO_COLOR_RESIDUE;
	for ( n_res_in_col=0,b=0; b<IN->nseq; b++){n_res_in_col+=(pos[b][a]>0)?1:0;}
	for ( max_score_col=0, score_col=0,b=0; b< IN->nseq; b++)
	    {
	    OUT->seq_al[b][a]=NO_COLOR_RESIDUE;
	    s1=IN->order[b][0];
	    r1=pos[b][a];
	    if (r1<=0)continue;
	    else
	      {
		max_score_r  =tot_extended_weight[s1][r1];
		score_r=res_extended_weight[s1][r1];
		res=(max_score_r==0 ||n_res_in_col<2 )?NO_COLOR_RESIDUE:((score_r*10)/max_score_r);
		(OUT)->seq_al[b][a]=(res==NO_COLOR_RESIDUE)?res:(MIN(res, 9));
		max_score_col+=max_score_r;
		    score_col+=score_r;
		max_score_seq[b]+=max_score_r;
		    score_seq[b]+=score_r;
		max_score_aln+=max_score_r;
		    score_aln+=score_r;
	      }
	    res=(max_score_col==0 || n_res_in_col<2)?NO_COLOR_RESIDUE:((score_col*10)/max_score_col);	
	    OUT->seq_al[IN->nseq][a]=(res==NO_COLOR_RESIDUE)?res:(MIN(res,9));
	    }
	}
    IN->score_aln=OUT->score_aln=MIN(100,((max_score_aln==0)?0:((score_aln*100)/max_score_aln)));
    for ( a=0; a< OUT->nseq; a++)
      {
	OUT->score_seq[a]=MIN(100,((max_score_seq[a]==0)?0:((score_seq[a]*100)/max_score_seq[a])));
	for ( b=0; b< OUT->nseq; b++)
	  {
	    OUT->pw_score_seq[a][b]+=OUT->score_seq[b];
	  }
	OUT->pw_score_seq[a][b]=MIN(100,(OUT->pw_score_seq[a][b]/OUT->nseq));
	
      }

    
    free ( score_seq);
    free ( max_score_seq);
    
    free_int (tot_extended_weight, -1);
    free_int (res_extended_weight, -1);
    
    
    return OUT;
    }
Alignment * non_extended_t_coffee_evaluate_output ( Alignment *IN,Constraint_list *CL)
    {
    int a,b, c,res, s1, s2, r1, r2;	
    Alignment *OUT=NULL;
    int **pos;
    int max_score_r, score_r;
    double score_col=0, score_aln=0;
    int max_score_col, max_score_aln;
    int *max_score_seq, *score_seq;
    int local_nseq;
    int **tot_non_extended_weight;
    int **res_non_extended_weight;
    int **l;
    CLIST_TYPE  *entry=NULL;
    int p;
    
    entry=vcalloc (CL->entry_len, CL->el_size);
    if ( !CL->evaluate_residue_pair){fprintf ( stderr, "\nWARNING: CL->evaluate_residue_pair Not set\nSet to: extend_residue_pair\n");CL->evaluate_residue_pair= extend_residue_pair; }
	
    OUT=copy_aln (IN, OUT);
    pos=aln2pos_simple(IN, IN->nseq);


    max_score_seq=vcalloc ( IN->nseq, sizeof (int));
    score_seq=vcalloc ( IN->nseq, sizeof (int));
    
    tot_non_extended_weight=list2residue_total_weight(CL);
    res_non_extended_weight=declare_int ((CL->S)->nseq, (CL->S)->max_len+1);
         
    for (a=0; a< IN->len_aln; a++)
        {
	    for ( b=0; b< IN->nseq-1; b++)
		{
		s1=IN->order[b][0];
		r1=pos[b][a];
		for ( c=b+1; c< IN->nseq; c++)
		    {
		    s2=IN->order[c][0];
		    r2=pos[c][a];	
		    if ( s1==s2 && !CL->do_self)continue;
		    else if ( r1<=0 || r2<=0)   continue;		    
		    else 
		      {
			entry[SEQ1]=(s1< s2)?s1:s2;
			entry[SEQ2]=(s1< s2)?s2:s1;
			entry[R1]=(s1< s2)?r1:r2;
			entry[R2]=(s1< s2)?r2:r1;
			if ((l=main_search_in_list_constraint (entry,&p,4,CL))!=NULL)
			  {
			    	res_non_extended_weight[s1][r1]+=l[0][WE];
				res_non_extended_weight[s2][r2]+=l[0][WE];
			  }
		      }
		    }
		}
	}
  
    sprintf ( OUT->name[IN->nseq], "Cons");
    for ( max_score_aln=0,score_aln=0,a=0; a< IN->len_aln; a++)
	{
	OUT->seq_al[IN->nseq][a]=NO_COLOR_RESIDUE;
	for ( local_nseq=0,b=0; b<IN->nseq; b++){local_nseq+=(pos[b][a]>0)?1:0;}
	
       	for ( max_score_col=0, score_col=0,b=0; b< IN->nseq; b++)
	    {
	    OUT->seq_al[b][a]=NO_COLOR_RESIDUE;
	    s1=IN->order[b][0];
	    r1=pos[b][a];
	    if (r1<=0)continue;
	    else
	      {
		max_score_r  =tot_non_extended_weight[s1][r1];
		score_r=res_non_extended_weight[s1][r1];
		res=(max_score_r==0 || local_nseq<2 )?NO_COLOR_RESIDUE:((score_r*10)/max_score_r);
	
		(OUT)->seq_al[b][a]=(res==NO_COLOR_RESIDUE)?res:(MIN(res, 9));
		max_score_col+=max_score_r;
		    score_col+=score_r;
		max_score_seq[b]+=max_score_r;
		    score_seq[b]+=score_r;
		max_score_aln+=max_score_r;
		    score_aln+=score_r;
	      }
	    res=(max_score_col==0 || local_nseq<2)?NO_COLOR_RESIDUE:((score_col*10)/max_score_col);	
	    OUT->seq_al[IN->nseq][a]=(res==NO_COLOR_RESIDUE)?res:(MIN(res,9));
	    }
	}
    IN->score_aln=OUT->score_aln=(max_score_aln==0)?0:((score_aln*100)/max_score_aln);
    for ( a=0; a< OUT->nseq; a++)
      {
	OUT->score_seq[a]=(max_score_seq[a]==0)?0:((score_seq[a]*100)/max_score_seq[a]);
	OUT->score_seq[a]=(OUT->score_seq[a]>100)?100:OUT->score_seq[a];
	
	for ( b=0; b< OUT->nseq; b++)
	  {
	    OUT->pw_score_seq[a][b]+=OUT->score_seq[b];
	  }
	OUT->pw_score_seq[a][b]=OUT->pw_score_seq[a][b]/OUT->nseq;
	OUT->score_seq[a]=(OUT->score_seq[a]>100)?100:OUT->score_seq[a];
      }
    OUT->score_aln=(OUT->score_aln>100)?100:OUT->score_aln;
    
    free ( score_seq);
    free ( max_score_seq);
    
    free_int (tot_non_extended_weight, -1);
    free_int (res_non_extended_weight, -1);
    vfree(entry);
    
    return OUT;
    }




/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR GETING THE COST : (Sequences) ->evaluate_residue_pair               */
/*                                                                                           */
/*********************************************************************************************/
int evaluate_cdna_matrix_score (Constraint_list *CL, int s1, int r1, int s2, int r2)
    {
      char a1, a2;
      if (r1>0 && r2>0) 
       {
	 r1--;
	 r2--;
	 
	 a1=translate_dna_codon((CL->S)->seq[s1]+r1,'x');
	 a2=translate_dna_codon((CL->S)->seq[s2]+r2,'x');
	 
	 
	 
	 if (a1=='x' || a2=='x')return 0;
	 else return CL->M[a1-'a'][a2-'a'];
       }
      else
	{
	  return 0;
	}
    }
int evaluate_matrix_score ( Constraint_list *CL, int s1, int r1, int s2, int r2)
    {
      /*
	function documentation: start
	
	int evaluate_matrix_score ( Constraint_list *CL, int s1, int s2, int r1, int r2)
	
	this function evaluates the score for matching residue S1(r1) wit residue S2(r2)
	using Matrix CL->M;
	
	function documentation: end
      */
	
      if (r1>0 && r2>0)
	    {
	    r1--;
	    r2--;
	    
	  
	    return CL->M[(CL->S)->seq[s1][r1]-'a'][(CL->S)->seq[s2][r2]-'a'];
	    }
	else
	    return 0;
    }
	    

int residue_pair_non_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2 )
        {
		
	/*
	  This is the generic Function->works with everything
	  	  
	  int residue_pair_non_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2, int field );
	  
	  Computes the non extended score for aligning residue seq1(r1) Vs seq2(r2)
	  
	  This function can compare a sequence with itself.
	  
	  Associated functions: See util constraint list, list extention functions.
	  function documentation: end
	*/

	int p;


	static int *entry;
	int **r;
	int field;

	field = CL->weight_field;


	if ( r1<=0 || r2<=0)return 0;
	else if ( !CL->extend_jit)
	   {
	    if ( !entry) entry=vcalloc (LIST_N_FIELDS , sizeof (int));
	    entry[SEQ1]=s1;
	    entry[SEQ2]=s2;
	    entry[R1]=r1;
	    entry[R2]=r2;
	    if ( r1==r2 && s1==s2) return UNDEFINED;
	    r=main_search_in_list_constraint( entry,&p,4,CL);
	    if (r==NULL)return 0;
	    else return r[0][field];
	   }
	else
	  return UNDEFINED;/*ERROR*/

	
	}



int residue_pair_extended_list_mixt (Constraint_list *CL, int s1, int r1, int s2, int r2 )
        {
	  int score=0;
	  
	  score+= residue_pair_extended_list_quadruplet(CL, s1, r1, s2, r2);
	  score+= residue_pair_extended_list (CL, s1, r1, s2, r2);
	  
	  return score;
	}

int residue_pair_extended_list_quadruplet (Constraint_list *CL, int s1, int r1, int s2, int r2 )
        {	
	  int t_s, t_r, t_w, q_s, q_r, q_w;
	  int a, b;
	  static int **hasch;
	  int score=0;
	  
	  int field;
	  /* This measure the quadruplets cost on a pair of residues*/
	  
	  field=CL->weight_field;
	  
	  if ( r1<=0 || r2<=0)return 0;
	  if ( !hasch)
	    {
	      hasch=vcalloc ( (CL->S)->nseq, sizeof (int*));
	      for ( a=0; a< (CL->S)->nseq; a++)hasch[a]=vcalloc ( (CL->S)->len[a]+1, sizeof (int));
	    }
	  
	  CL=index_res_constraint_list ( CL, field);	  
	  hasch[s1][r1]=100000;
	  for (a=1; a< CL->residue_index[s1][r1][0]; a+=3)
	    {
	      t_s=CL->residue_index[s1][r1][a];
	      t_r=CL->residue_index[s1][r1][a+1];
	      t_w=CL->residue_index[s1][r1][a+2];
	      if ( CL->seq_for_quadruplet[t_s])
		{
		  for ( b=1; b<CL->residue_index[t_s][t_r][0]; b+=3)
		    {
		      q_s=CL->residue_index[t_s][t_r][b];		  
		      q_r=CL->residue_index[t_s][t_r][b+1];		 
		      if (CL-> seq_for_quadruplet[q_s])
			  hasch[q_s][q_r]=MIN(CL->residue_index[t_s][t_r][b+2],t_w);
		      
		    }
		}
	    }
	  
	  
	  for (a=1; a< CL->residue_index[s2][r2][0]; a+=3) 
	    {
	      t_s=CL->residue_index[s2][r2][a];
	      t_r=CL->residue_index[s2][r2][a+1];
	      t_w=CL->residue_index[s2][r2][a+2];
	      if ( CL->seq_for_quadruplet[t_s])
		{
		  for ( b=1; b<CL->residue_index[t_s][t_r][0]; b+=3)
		    {
		      q_s=CL->residue_index[t_s][t_r][b];		  
		      q_r=CL->residue_index[t_s][t_r][b+1];	
		      q_w=CL->residue_index[t_s][t_r][b+2];	
		      if (hasch[q_s][q_r] && CL->seq_for_quadruplet[q_s])
			score+=MIN(hasch[q_s][q_r],MIN(CL->residue_index[t_s][t_r][b+2],q_w));
		    }
		}
	    }
	  
	  score=(CL->normalise)?((score*CL->normalise)/CL->max_ext_value):score;
	  
	  for (a=1; a< CL->residue_index[s1][r1][0]; a+=3)
	    {
	      t_s=CL->residue_index[s1][r1][a];
	      t_r=CL->residue_index[s1][r1][a+1];
	      t_w=CL->residue_index[s1][r1][a+2];
	      if ( CL->seq_for_quadruplet[t_s])
		{
		  for ( b=1; b<CL->residue_index[t_s][t_r][0]; b+=3)
		    {
		      q_s=CL->residue_index[t_s][t_r][b];		  
		      q_r=CL->residue_index[t_s][t_r][b+1];		 
		      hasch[q_s][q_r]=0;
		    }
		}
	    }
	  
	  return score;
	}  


  
int residue_pair_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2 )
        {
	int a, t_s, t_r;
	static int **hasch;
	int score=0;

	int field;
	/*
	  function documentation: start

	  int residue_pair_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2);
	  
	  Computes the extended score for aligning residue seq1(r1) Vs seq2(r2)
	  Computes: matrix_score
	            non extended score
		    extended score

	  The extended score depends on the function index_res_constraint_list.
	  This function can compare a sequence with itself.
	  
	  Associated functions: See util constraint list, list extention functions.
	  
	  function documentation: end
	*/



	field=CL->weight_field;

	if ( r1<=0 || r2<=0)return 0;
	if ( !hasch)
	       {
	       hasch=vcalloc ( (CL->S)->nseq, sizeof (int*));
	       for ( a=0; a< (CL->S)->nseq; a++)hasch[a]=vcalloc ( (CL->S)->len[a]+1, sizeof (int));
	       }
	
	CL=index_res_constraint_list ( CL, field);
	
	hasch[s1][r1]=100000;
	for (a=1; a< CL->residue_index[s1][r1][0]; a+=3)
	  {
	    t_s=CL->residue_index[s1][r1][a];
	    t_r=CL->residue_index[s1][r1][a+1];
	    hasch[t_s][t_r]=CL->residue_index[s1][r1][a+2];		
	  }
	
	
	for (a=1; a< CL->residue_index[s2][r2][0]; a+=3) 
	  {
	    t_s=CL->residue_index[s2][r2][a];
	    t_r=CL->residue_index[s2][r2][a+1];
	    if (hasch[t_s][t_r])
	      {
		if (field==WE){score+=MIN(hasch[t_s][t_r],CL->residue_index[s2][r2][a+2]);}
	      }
	  }

	score=(CL->normalise)?((score*CL->normalise)/CL->max_ext_value):score;
	
	for (a=1; a< CL->residue_index[s1][r1][0]; a+=3)
	  {
	    t_s=CL->residue_index[s1][r1][a];
	    t_r=CL->residue_index[s1][r1][a+1];
	    hasch[t_s][t_r]=0;	      
	  }
	hasch[s1][r1]=hasch[s2][r2]=0;

	return score;
	}

int extend_residue_pair ( Constraint_list *CL, int s1, int r1, int s2, int r2)
        {
	int a, t_s, t_r, p;
	static int **hasch;
	int score=0;
	static int *entry;
	int **r;
	int field;


	
	/*
	  This is the generic Function->works with everything
	  should be gradually phased out

	  
	  int extend_residue_pair ( Constraint_list *CL, int s1, int r1, int s2, int r2, int field )
	  
	  Computes the extended score for aligning residue seq1(r1) Vs seq2(r2)
	  Computes: matrix_score
	            non extended score
		    extended score

	  The extended score depends on the function index_res_constraint_list.
	  This function can compare a sequence with itself.
	  
	  Associated functions: See util constraint list, list extention functions.
	  function documentation: end
	*/


	field=CL->weight_field;

	if ( r1<=0 || r2<=0)return 0;
	else if ( !CL->L && CL->M)
	   {
	    return evaluate_matrix_score (CL, s1,r1, s2, r2);	
	   }
	
	else if ( !CL->extend_jit)
	   {
	    if ( !entry) entry=vcalloc (LIST_N_FIELDS , sizeof (int));
	    entry[SEQ1]=s1;
	    entry[SEQ2]=s2;
	    entry[R1]=r1;
	    entry[R2]=r2;
	    r=main_search_in_list_constraint( entry,&p,4,CL);
	    if (r==NULL)return 0;
	    else return r[0][field];
	   }
	else
	   {
	   if ( !hasch)
	       {
	       hasch=vcalloc ( (CL->S)->nseq, sizeof (int*));
	       for ( a=0; a< (CL->S)->nseq; a++)hasch[a]=vcalloc ( (CL->S)->len[a]+1, sizeof (int));
	       }
	
	   CL=index_res_constraint_list ( CL, field);
	
	   hasch[s1][r1]=100000;
	   for (a=1; a< CL->residue_index[s1][r1][0]; a+=3)
	        {
		t_s=CL->residue_index[s1][r1][a];
		t_r=CL->residue_index[s1][r1][a+1];
		hasch[t_s][t_r]=CL->residue_index[s1][r1][a+2];		
		}
	
	
	   for (a=1; a< CL->residue_index[s2][r2][0]; a+=3) 
	       {
	       t_s=CL->residue_index[s2][r2][a];
	       t_r=CL->residue_index[s2][r2][a+1];
	       if (hasch[t_s][t_r])
		   {
		   if (field==WE)score+=MIN(hasch[t_s][t_r],CL->residue_index[s2][r2][a+2] );

		   }
	       }
	   score=(CL->normalise)?((score*CL->normalise)/CL->max_ext_value):score;
	   for (a=1; a< CL->residue_index[s1][r1][0]; a+=3)
	       {
	       t_s=CL->residue_index[s1][r1][a];
	       t_r=CL->residue_index[s1][r1][a+1];
	       hasch[t_s][t_r]=0;	      
	       }
	   hasch[s1][r1]=hasch[s2][r2]=0;

	   return score;
	   }
	}
/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR GETTING THE COST :  CL->get_dp_cost                                     */
/*                                                                                           */
/*********************************************************************************************/



int get_cdna_best_frame_dp_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
    {
	int a, b;
	int n=4;
	int s;
	char a1, a2;
	static int l1, l2;
	static Alignment *B;
	static int **score;
	
	if ( !score)score=declare_int(3, 2);
	
	if (!A)
	   {
	       free_aln(B);
	       B=NULL;
	       return UNDEFINED;
	   }
	if (!B)
	   {
	       if (ns1+ns2>2){fprintf ( stderr, "\nERROR: get_cdna_dp_cost mode is only for pair-wise ALN [FATAL]\n");crash("");} 
	       free_aln (B);
	       B=copy_aln (A, NULL);
	       
	       l1=strlen ( A->seq_al[list1[0]]);
	       for ( b=0; b<l1; b++)
		       B->seq_al[list1[0]][b]=translate_dna_codon (A->seq_al[list1[0]]+b, 'x');
	       l2=strlen ( A->seq_al[list2[0]]);
	       for ( b=0; b<l2; b++)
		       B->seq_al[list2[0]][b]=translate_dna_codon (A->seq_al[list2[0]]+b, 'x');
	   }
	
/*Set the frame*/

	for ( a=0; a< 3; a++)score[a][0]=score[a][1]=0;
	for ( a=col1-(n*3),b=col2-(n*3); a<col1+(n*3) ; a++, b++)
	        {
		    if ( a<0 || b<0 || a>=l1 || b>=l2)continue;
		    
		    a1=B->seq_al[list1[0]][a];
		    a2=B->seq_al[list2[0]][b];
		    
		    score[a%3][0]+=(a1=='x' || a2=='x')?0:CL->M[a1-'a'][a2-'a'];
		    score[a%3][1]++;
		 }
	
	for ( a=0; a< 3; a++)score[a][0]=(score[a][1]>0)?(score[a][0]/score[a][1]):0;
	if ( score[0][0]>score[1][0] &&  score[0][0]>score[2][0])
	    s=score[0][0];
	else if ( score[1][0]>score[0][0] &&  score[1][0]>score[2][0])
	    s=score[1][0];
	else s=score[2][0];
	
	return s*SCORE_K;
	
	} 

int get_dp_cost_quadruplet ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
        {
	  int score;
	  
	 
	  if ( ns1==1 || ns2==1)
	     score=slow_get_dp_cost ( A, pos1, ns1, list1,col1, pos2, ns2, list2, col2, CL);
	  else
	    score=fast_get_dp_cost_quadruplet ( A, pos1, ns1, list1,col1, pos2, ns2, list2, col2, CL);
	 
	  return (score==UNDEFINED)?UNDEFINED:(score-SCORE_K*CL->nomatch);
	}  
int get_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
        {
	int MODE=1;
	int score;
	


	if (A==NULL)return 0;
	
	if (MODE==0 || (!CL->L && CL->M) || (!CL->L && CL->T)|| ns1==1 || ns2==1)
	  score=slow_get_dp_cost ( A, pos1, ns1, list1,col1, pos2, ns2, list2, col2, CL);
	else if (MODE==1)
	  score=fast_get_dp_cost ( A, pos1, ns1, list1,col1, pos2, ns2, list2, col2, CL);
	else
	  score=0;


	return (score==UNDEFINED)?UNDEFINED:(score-SCORE_K*CL->nomatch);
	}


int very_fast_get_dp_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
    {
    /*WARNING: WORKS ONLY WITH A SUBSTITUTION MATRIX!!!!!!!!*/	
      int res, a, b, score=0;
      int  *used_res1,  *used_res2;
      static int **list_res1, **list_res2;
      int n_used_res1, n_used_res2;

      
      if ( !A)return UNDEFINED;
      if ( !CL->M)crash("\nERROR: Matrix Needed [FATAL:very_fast_get_dp_cost]\n");      
      if ( CL->evaluate_residue_pair!=evaluate_matrix_score)crash("\nERROR: pair function not correct [FATAL]\n");
      n_used_res1=n_used_res2=0;
      used_res1=vcalloc ( 1000, sizeof (int));
      used_res2=vcalloc ( 1000, sizeof (int));
      
      if ( !list_res1)
	{
	  list_res1=declare_int (30, 2);
	  list_res2=declare_int (30, 2);
	}
      


     

      for (a=0; a<ns1; a++)
	  {
	    res=A->seq_al[list1[a]][col1];
	    /*if (is_gap(res))continue;*/
	    if (!used_res1[res])
	      {
		used_res1[res]=n_used_res1;
		list_res1[n_used_res1][0]=res;
		list_res1[n_used_res1][1]=1;
		n_used_res1++;
	      }
	    else 
	      {
		list_res1[used_res1[res]][1]++;
	      }
	  }

	for (a=0; a<ns2; a++)
	  {
	    res=A->seq_al[list2[a]][col2];
	    /*if (is_gap(res))continue;*/
	    if (!used_res2[res])
	      {
		used_res2[res]=n_used_res2;
		list_res2[n_used_res2][0]=res;
		list_res2[n_used_res2][1]=1;
		n_used_res2++;
	      }
	    else 
	      {
		list_res2[used_res2[res]][1]++;
	      }
	  }
	
	if ( n_used_res1==0 || n_used_res2==0)
	  {
	    vfree (used_res1);
	    vfree (used_res2);
	    return 0;
	  }

	for ( a=0; a< n_used_res1; a++)
	  for ( b=0; b< n_used_res2; b++)
	    {
	      if ( is_gap(list_res1[a][0]) && is_gap (list_res2[b][0]));
	      else if ( is_gap(list_res1[a][0]) || is_gap (list_res2[b][0]))score+=CL->gep*list_res1[a][1]*list_res2[b][1];
	      else score+=(CL->M)[list_res1[a][0]-'a'][list_res2[b][0]-'a']*list_res1[a][1]*list_res2[b][1];
	    }
	
	
	score=score*SCORE_K/(n_used_res1*n_used_res2);
	vfree (used_res1);
	vfree (used_res2);
	
	return score;
    }


int fast_get_dp_cost_quadruplet ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
	{
	/*WARNING: WORKS ONLY WITH List to Extend*/
	  /*That function does a quruple extension beween two columns by pooling the residues together*/
	  
	  int a,b, c;
	  int n_gap1=0;
	  int n_gap2=0;
	  int score=0;
	  int s1, rs1, r1, t_r, t_s,t_w, q_r, q_s, q_w, s2, rs2, r2;
	  int **buf_pos, buf_ns, *buf_list, buf_col; 

	  static int **hasch1;
	  static int **hasch2;
	  
	  static int **n_hasch1;
	  static int **n_hasch2;
	  
	  static int **is_in_col1;
	  static int **is_in_col2;
	  

	  if (ns2>ns1)
	    {
	      buf_pos=pos1;
	      buf_ns=ns1;
	      buf_list=list1;
	      buf_col=col1;
	      
	      pos1=pos2;
	      ns1=ns2;
	      list1=list2;
	      col1=col2;
	      
	      pos2=buf_pos;
	      ns2=buf_ns;
	      list2=buf_list;
	      col2=buf_col;
	    }
	  
	CL=index_res_constraint_list ( CL, WE);
	if ( !hasch1)
	    {
	    
	    hasch1=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
            hasch2=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
	    n_hasch1=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	    n_hasch2=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
	    is_in_col1=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
	    is_in_col2=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
	    }
	
	for ( a=0; a< ns1; a++)
	    {
		rs1= list1[a];
		s1=A->order[rs1][0];
		r1=pos1[rs1][col1];
		
		if (r1<0)n_gap1++;		
		else
		   {	
		   is_in_col1[s1][r1]=1;
		   for (b=1; b< CL->residue_index[s1][r1][0]; b+=3)
		           {
			   t_s=CL->residue_index[s1][r1][b];
			   t_r=CL->residue_index[s1][r1][b+1];
			   t_w=CL->residue_index[s1][r1][b+2];
			   for ( c=1; c< CL->residue_index[t_s][t_r][0]; c+=3)
			     {
			       q_s=CL->residue_index[t_s][t_r][c];
			       q_r=CL->residue_index[t_s][t_r][c+1];
			       q_w=CL->residue_index[t_s][t_r][c+2];
			       hasch1[q_s][q_r]+=MIN(q_w, t_w);
			       n_hasch1[q_s][q_r]++;
			     }
			   }
		   }
	    }
	
	for ( a=0; a< ns2; a++)
	    {
		rs2=list2[a];
		s2=A->order[rs2][0];
		r2=pos2[rs2][col2];
	
		if (r2<0)n_gap2++;
		else
		   {
		   is_in_col2[s2][r2]=1;
		   for (b=1; b< CL->residue_index[s2][r2][0]; b+=3)
		           {
			   t_s=CL->residue_index[s2][r2][b];
			   t_r=CL->residue_index[s2][r2][b+1];
			   t_w=CL->residue_index[s2][r2][b+2];
			   for ( c=1; c< CL->residue_index[t_s][t_r][0]; c+=3)
			     {
			       q_s=CL->residue_index[t_s][t_r][c];
			       q_r=CL->residue_index[t_s][t_r][c+1];
			       q_w=CL->residue_index[t_s][t_r][c+2];
			       hasch2[q_s][q_r]+=MIN(t_w, q_w);
			       n_hasch2[q_s][q_r]++;
			     }
			   }
		   }
	    }
	

	for ( a=0; a< ns2; a++)
	  {
	    rs2=list2[a];
	    s2=A->order[rs2][0];
	    r2=pos1[rs2][col2];
	    
	    if (r2<0);
	    else
	      {
		for (b=1; b< CL->residue_index[s2][r2][0]; b+=3)
		  {
		    t_s=CL->residue_index[s2][r2][b];
		    t_r=CL->residue_index[s2][r2][b+1];

		    for ( c=1; c< CL->residue_index[t_s][t_r][0]; c+=3)
		      {
			q_s=CL->residue_index[t_s][t_r][c];
			q_r=CL->residue_index[t_s][t_r][c+1];
			if ( hasch2[q_s][q_r] && hasch1[q_s][q_r]&& !(is_in_col1[q_s][q_r] || is_in_col2[q_s][q_r]))
			  {
			    score+=MIN(hasch2[q_s][q_r]*(n_hasch1[q_s][q_r]),hasch1[q_s][q_r]*(n_hasch2[q_s][q_r]));
			  }
			else if ( hasch2[q_s][q_r] && is_in_col1[q_s][q_r])
			  {
			    score+=hasch2[q_s][q_r]*(n_hasch1[q_s][q_r]+1);
			  }
			else if (hasch1[q_s][q_r] && is_in_col2[q_s][q_r])
			  {
			    score+=hasch1[q_s][q_r]*(n_hasch2[q_s][q_r]+1);
			  }
			hasch2[q_s][q_r]=0;
			n_hasch2[q_s][q_r]=0;
		      }
		  }
		hasch2[s2][r2]=0;
		is_in_col2[s2][r2]=0;
	      }
	  }
	
	
	for ( a=0; a< ns1; a++)
	  {
	    rs1= list1[a];
	    s1=A->order[rs1][0];
	    r1=pos1[rs1][col1];
	    
	    if (r1<0);
	    else
	      {
		is_in_col1[s1][r1]=0;
		hasch1[s1][r1]=0;
		for (b=1; b< CL->residue_index[s1][r1][0]; b+=3)
		  {
		    t_s=CL->residue_index[s1][r1][b];
		    t_r=CL->residue_index[s1][r1][b+1];
		    for ( c=1; c< CL->residue_index[t_s][t_r][0]; c+=3)
		      {
			q_s=CL->residue_index[t_s][t_r][c];
			q_r=CL->residue_index[t_s][t_r][c+1];
			hasch1[q_s][q_r]=0;
			n_hasch1[q_s][q_r]=0;
		      }
		  }
	      }
	  }
	

	score=(score*SCORE_K)/((ns1-n_gap1)*(ns2-n_gap2));
	score=(CL->normalise)?((score*CL->normalise)/CL->max_ext_value):score;

	return score;
	}

	    
int fast_get_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
	{
	/*WARNING: WORKS ONLY WITH List to Extend*/
	  
	  
	  int a,b;
	  int n_gap1=0;
	  int n_gap2=0;
	  int score=0;
	  int s1, rs1, r1, t_r, t_s, s2, rs2, r2;
	  static int **hasch1;
	  static int **hasch2;
	  
	  static int **n_hasch1;
	  static int **n_hasch2;
	  
	  static int **is_in_col1;
	  static int **is_in_col2;
	    
	CL=index_res_constraint_list ( CL, WE);
	if ( !hasch1)
	    {
	    
	    hasch1=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
            hasch2=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
	    n_hasch1=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	    n_hasch2=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
	    is_in_col1=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
	    is_in_col2=declare_int( (CL->S)->nseq, (CL->S)->max_len+1);
	    }
	
	for ( a=0; a< ns1; a++)
	    {
		rs1= list1[a];
		s1=A->order[rs1][0];
		r1=pos1[rs1][col1];
		
		if (r1<0)n_gap1++;		
		else
		   {	
		   is_in_col1[s1][r1]=1;
		   for (b=1; b< CL->residue_index[s1][r1][0]; b+=3)
		           {
			   t_s=CL->residue_index[s1][r1][b];
			   t_r=CL->residue_index[s1][r1][b+1];			   
			   hasch1[t_s][t_r]+=CL->residue_index[s1][r1][b+2];
			   n_hasch1[t_s][t_r]++;
			   }
		   }
	    }
	
	for ( a=0; a< ns2; a++)
	    {
		rs2=list2[a];
		s2=A->order[rs2][0];
		r2=pos2[rs2][col2];
	
		if (r2<0)n_gap2++;
		else
		   {
		   is_in_col2[s2][r2]=1;
		   for (b=1; b< CL->residue_index[s2][r2][0]; b+=3)
		           {
			   t_s=CL->residue_index[s2][r2][b];
			   t_r=CL->residue_index[s2][r2][b+1];

			   hasch2[t_s][t_r]+=CL->residue_index[s2][r2][b+2];
			   n_hasch2[t_s][t_r]++;
			   }
		   }
	    }

	if ( ns2<ns1)
	    {
		for ( a=0; a< ns2; a++)
		    {
		    rs2=list2[a];
		    s2=A->order[rs2][0];
		    r2=pos1[rs2][col2];
		    
		    if (r2<0);
		    else
		        {
			for (b=1; b< CL->residue_index[s2][r2][0]; b+=3)
			    {
			    t_s=CL->residue_index[s2][r2][b];
			    t_r=CL->residue_index[s2][r2][b+1];
			    
			    if ( hasch2[t_s][t_r] && hasch1[t_s][t_r]&& !(is_in_col1[t_s][t_r] || is_in_col2[t_s][t_r]))
			        {
				score+=MIN(hasch2[t_s][t_r]*(n_hasch1[t_s][t_r]),hasch1[t_s][t_r]*(n_hasch2[t_s][t_r]));
				}
			    else if ( hasch2[t_s][t_r] && is_in_col1[t_s][t_r])
			        {
				score+=hasch2[t_s][t_r]*(n_hasch1[t_s][t_r]+1);
				}
			    else if (hasch1[t_s][t_r] && is_in_col2[t_s][t_r])
			        {
				score+=hasch1[t_s][t_r]*(n_hasch2[t_s][t_r]+1);
				}
			    hasch2[t_s][t_r]=0;
			    n_hasch2[t_s][t_r]=0;
			    }
			hasch2[s2][r2]=0;
			is_in_col2[s2][r2]=0;
			}
		    }

	
		for ( a=0; a< ns1; a++)
		    {
		    rs1= list1[a];
		    s1=A->order[rs1][0];
		    r1=pos1[rs1][col1];
		    
		    if (r1<0);
		    else
		        {
			is_in_col1[s1][r1]=0;
			hasch1[s1][r1]=0;
			for (b=1; b< CL->residue_index[s1][r1][0]; b+=3)
			    {
			    t_s=CL->residue_index[s1][r1][b];
			    t_r=CL->residue_index[s1][r1][b+1];
			    
			    hasch1[t_s][t_r]=0;
			    n_hasch1[t_s][t_r]=0;
			    }
			}
		    }
	    }
	else
	   {
		for ( a=0; a< ns1; a++)
		    {
		    rs1=list1[a];
		    s1=A->order[rs1][0];
		    r1=pos1[rs1][col1];
		    
		    if (r1<0);
		    else
		        {
			for (b=1; b< CL->residue_index[s1][r1][0]; b+=3)
			    {
			    t_s=CL->residue_index[s1][r1][b];
			    t_r=CL->residue_index[s1][r1][b+1];
			    
			    if ( hasch1[t_s][t_r] && hasch2[t_s][t_r]&& !(is_in_col2[t_s][t_r] || is_in_col1[t_s][t_r]))
			        {
				score+=MIN(hasch1[t_s][t_r]*(n_hasch2[t_s][t_r]),hasch2[t_s][t_r]*(n_hasch1[t_s][t_r]));
			        }
			    else if ( hasch1[t_s][t_r] && is_in_col2[t_s][t_r])
			        {
				score+=hasch1[t_s][t_r]*(n_hasch2[t_s][t_r]+1);
			        }
			    else if (hasch2[t_s][t_r] && is_in_col1[t_s][t_r])
			        {
				score+=hasch2[t_s][t_r]*(n_hasch1[t_s][t_r]+1);
			        }
			    hasch1[t_s][t_r]=0;
			    n_hasch1[t_s][t_r]=0;
			    }
			hasch1[s1][r1]=0;
			is_in_col1[s1][r1]=0;
		        }
		    }

	
		for ( a=0; a< ns2; a++)
		    {
		    rs2= list2[a];
		    s2=A->order[rs2][0];
		    r2=pos1[rs2][col2];
		    
		    if (r2<0);
		    else
		        {
			is_in_col2[s2][r2]=0;
			hasch1[s2][r2]=0;
			for (b=1; b< CL->residue_index[s2][r2][0]; b+=3)
			    {
			    t_s=CL->residue_index[s2][r2][b];
			    t_r=CL->residue_index[s2][r2][b+1];
			    
			    hasch2[t_s][t_r]=0;
			    n_hasch2[t_s][t_r]=0;
			    }
		        }
		    }
	    }
	score=(score*SCORE_K)/((ns1-n_gap1)*(ns2-n_gap2));
	score=(CL->normalise)?((score*CL->normalise)/CL->max_ext_value):score;

	return score;
	}

int slow_get_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
	{
	int a, b, s1, s2,r1, r2;
	static int *entry;

	int score=0, s;
	int gap_gap=0;
	int gap_res=0;
	int res_res=0;
	int rs1, rs2;
	
	if ( !CL->evaluate_residue_pair){fprintf ( stderr, "\nWARNING: CL->evaluate_residue_pair Not set\nSet to: extend_residue_pair\n");CL->evaluate_residue_pair= extend_residue_pair;}
	
	if ( entry==NULL) entry=calloc (LIST_N_FIELDS , sizeof (int));
	
	
	for ( a=0; a< ns1; a++)
		{
		for ( b=0; b<ns2; b++)
			{
			
			s1 =list1[a];
			rs1=A->order[s1][0];
			r1 =pos1[s1][col1];
			
			s2 =list2[b];
			rs2=A->order[s2][0];
			r2 =pos2[s2][col2];
			
			
			/*
			if ( ns1+ns2>2)
			  {
			    fprintf ( stderr, "\nWARNING: TEST in slow_get_dp_cost");
			    fprintf ( stderr, "\n[%d %d][%d %d]->%d", rs1, r1, rs2, r2,(CL->evaluate_residue_pair)(CL, rs1, r1, rs2, r2));
			  }
			*/
	        	if ( rs1>rs2)
			   {			    
			   SWAP (rs1, rs2);
			   SWAP (r1, r2);
			   }
			
			if (r1==0 && r2==0)gap_gap++;			     
			else if ( r1<0 || r2<0) gap_res++;			
			else 
			    {					
			    res_res++;
			    if ((s=(CL->evaluate_residue_pair)(CL, rs1, r1, rs2, r2))!=UNDEFINED) score+=s;				    
			    else 
			      {
				
				return UNDEFINED;
			      }
			    }
		
			}
		}
	
	score=(res_res==0)?0:( (score*SCORE_K)/res_res);

	return score;	
	} 
int sw_get_dp_cost     ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
	{
	  int a, b;
	  int s1,r1,rs1;
	  int s2,r2,rs2;
	  
	  
	  

	  for ( a=0; a< ns1; a++)
	    {	
	      s1 =list1[a];
	      rs1=A->order[s1][0];
	      r1 =pos1[s1][col1];
	      if ( r1<=0)continue;
	      for ( b=0; b< ns2; b++)
		{
		  
		  
		  s2 =list2[b];
		  rs2=A->order[s2][0];
		  r2 =pos2[s2][col2];
		  
		  if (r2<=0)continue;
		  
		  
		  if (sw_pair_is_defined (CL, rs1, r1, rs2, r2)==UNDEFINED)return UNDEFINED;	       
		}
	    }	  
			
	  return slow_get_dp_cost ( A, pos1, ns1, list1, col1, pos2, ns2, list2,col2, CL);
	  	  
	}




  


int get_domain_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2,Constraint_list *CL , int scale , int gop, int gep)

	{
	int a, b, s1, s2,r1, r2;
	static int *entry;
	int **r;
	int score=0;
	int gap_gap=0;
	int gap_res=0;
	int res_res=0;
	int rs1, rs2;
	int flag_list_is_aa_sub_mat=0;
	int p;

/*Needs to be cleanned After Usage*/
	


	if ( entry==NULL) entry=calloc (LIST_N_FIELDS , sizeof (int));
	
	for (a=0; a< ns1; a++)
		{
		s1=list1[a];
		rs1=A->order[s1][0];
		for ( b=0; b<ns2; b++)
			{
			s2 =list2[b];
			rs2=A->order[s2][0];
			
			entry[SEQ1]=rs1;
			entry[SEQ2]=rs2;
			r1=entry[R1]=pos1[s1][col1];
			r2=entry[R2]=pos2[s2][col2];
			
			if ( !flag_list_is_aa_sub_mat)
			    {
			    if ( r1==r2 && rs1==rs2)
			       {
				 
				 return UNDEFINED;
			       }
			    else if (r1==0 && r2==0)
			       {			    
			       gap_gap++;		
			       }
			    else if ( r1<=0 || r2<=0)
			       {
			       gap_res++;
			       }
			    else if ((r=main_search_in_list_constraint ( entry,&p,4,CL))!=NULL)
				{				
				res_res++; 
				
				if (r[0][WE]!=UNDEFINED) 
				    {
				    score+=(r[0][WE]*SCORE_K)+scale;
				    }
				else 
				  {
				    fprintf ( stderr, "**");
				    return UNDEFINED;
				  }
				}
			    }			    			  
			}
		}
	return score;
	score=((res_res+gap_res)==0)?0:score/(res_res+gap_res);
	return score;	
	} 

/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR ANALYSING AL OR MATRIX                                              */
/*                                                                                           */
/*********************************************************************************************/

int aln2n_res ( Alignment *A, int start, int end)
    {
    int a, b;
    int score=0;

    for ( a=start; a<end; a++)for ( b=0; b< A->nseq; b++)score+=!is_gap(A->seq_al[b][a]);
    return score;
    }

float get_gop_scaling_factor ( int **matrix,float id, int l1, int l2)
    {
	return id* get_avg_matrix_mm(matrix, AA_ALPHABET);
    }

float get_avg_matrix_mm ( int **matrix, char *alphabet)
    {
	int a, b;
	float naa;
	float gop;
	int l;

	
	l=MIN(20,strlen (alphabet));
	for (naa=0, gop=0,a=0; a<l; a++)
	    for ( b=0; b<l; b++)
	        {
		    if ( a!=b)
		        {
			gop+=matrix[alphabet[a]-'a'][alphabet[b]-'a'];
			naa++;
			}
		}
	
	gop=gop/naa;
	return gop;
    }
float get_avg_matrix_match ( int **matrix, char *alphabet)
    {
	int a;
	float gop;
	int l;

	

	
	l=MIN(20,strlen (alphabet));
	for (gop=0,a=0; a<l; a++)
	  gop+=matrix[alphabet[a]-'a'][alphabet[a]-'a'];
	
	gop=gop/l;
	return gop;
    }	
		      
float get_avg_matrix_diff ( int **matrix1,int **matrix2, char *alphabet)
    {
	int a, b;
	float naa;
	float gop;
	int l,v1;

	

	
	l=MIN(20,strlen (alphabet));
	for (naa=0, gop=0,a=0; a<l; a++)
	  for (b=0; b<l; b++)
	    {
	      v1=matrix1[alphabet[a]-'a'][alphabet[a]-'a']-matrix2[alphabet[b]-'a'][alphabet[b]-'a'];
	      gop+=v1*v1;
	      naa++;
	    }
	
	gop=gop/l;
	return gop;
    }	

float* set_aa_frequencies ()
   {
     
     float *frequency;
     /*frequencies tqken from psw*/
     
     frequency=vcalloc (100, sizeof (float));
     frequency ['x'-'a']=0.0013;
     frequency ['a'-'a']=0.0076;
     frequency ['c'-'a']=0.0176;
     frequency ['d'-'a']=0.0529;
     frequency ['e'-'a']=0.0628;
     frequency ['f'-'a']=0.0401;
     frequency ['g'-'a']=0.0695;
     frequency ['h'-'a']=0.0224;
     frequency ['i'-'a']=0.0561;
     frequency ['k'-'a']=0.0584;
     frequency ['l'-'a']=0.0922;
     frequency ['m'-'a']=0.0236;
     frequency ['n'-'a']=0.0448;
     frequency ['p'-'a']=0.0500;
     frequency ['q'-'a']=0.0403;
     frequency ['r'-'a']=0.0523;
     frequency ['s'-'a']=0.0715;
     frequency ['t'-'a']=0.0581;
     frequency ['v'-'a']=0.0652;
     frequency ['w'-'a']=0.0128;
     frequency ['y'-'a']=0.0321; 
   
     return frequency;
   }

float measure_matrix_pos_avg (int **matrix,char *alphabet)
{
	float naa=0, tot=0;
	int a, b;
	
	
	for ( tot=0,a=0; a< 20; a++)
	   for ( b=0; b<20; b++)
	     {
	       if (matrix[alphabet[a]-'a'][alphabet[b]-'a']>=0){naa++;tot+=matrix[alphabet[a]-'a'][alphabet[b]-'a'];}

	     }
	return tot/naa;
}
  
float measure_matrix_enthropy (int **matrix,char *alphabet)
   {
     
     int a, b;
     double s, p, q, h=0, tq=0;
     float lambda;
     float *frequency;
     /*frequencies tqken from psw*/
     
     frequency=set_aa_frequencies ();

     
     lambda=compute_lambda(matrix,alphabet);
     fprintf ( stderr, "\nLambda=%f", (float)lambda);
	   
     for ( a=0; a< 20; a++)
       for ( b=0; b<=a; b++)
	 {
	   s=matrix[alphabet[a]-'a'][alphabet[b]-'a'];
	   
	   
	   p=frequency[alphabet[a]-'a']*frequency[alphabet[b]-'a'];
	   
	   if ( p==0)continue;
	   
	   q=exp(lambda*s+log(p));
	   
	   tq+=q;
	   h+=q*log(q/p)*log(2);
	  
	 }
	   
     fprintf ( stderr,"\ntq=%f\n", (float)tq);
   
     return (float) h;
    }	
float compute_lambda (int **matrix,char *alphabet)
{

  int a, b;
  double lambda, best_lambda, delta, best_delta, p, tq,s;
  static float *frequency;
  
  if ( !frequency)frequency=set_aa_frequencies ();

  for ( lambda=0.001; lambda<1; lambda+=0.005)
     {
       tq=0;
       for ( a=0; a< 20; a++)
	 for ( b=0; b<20; b++)
	   {	     
	     p=frequency[alphabet[a]-'a']*frequency[alphabet[b]-'a'];
	     s=matrix[alphabet[a]-'a'][alphabet[b]-'a'];
	     tq+=exp(lambda*s+log(p));
	   }
       delta=fabs(1-tq);
       if (lambda==0.001)
	 {
	   best_delta=delta;
	   best_lambda=lambda;
	 }
       else
	 {
	   if (delta<best_delta)
	     {
	       best_delta=delta;
	       best_lambda=lambda;
	     }
	 }   
       
       fprintf ( stderr, "\n%f %f ", lambda, tq);
       if ( tq>1)break;
     }
   fprintf ( stderr, "\nRESULT: %f %f ", best_lambda, best_delta);
   return (float) best_lambda;
}



float evaluate_random_match (char  *mat, int n, int len,char *alp)
{
  int **matrix;
  matrix=read_matrice ( mat); 
  fprintf ( stderr, "Matrix=%15s ", mat);
  return evaluate_random_match2 (matrix, n,len,alp);
  
}

float evaluate_random_match2 (int **matrix, int n, int len,char *alp)
{
  int a, b, c, d, c1, c2, tot;
  static int *list;
  static float *freq;
  float score_random=0;
  float score_id=0;
  float score_good=0;
  float tot_len=0;
  float tot_try=0;


  if ( !list)
    {
      srand(time(NULL));
      freq=set_aa_frequencies ();
      list=vcalloc ( 10000, sizeof (char));
    }
  
  for (tot=0,c=0,a=0;a<20; a++)
    {
      b=freq[alp[a]-'a']*1000;
      tot+=b;
      for (d=0; d<b; d++, c++)
	{
	  list[c]=alp[a];
	}
    }
  
  
  for (a=0; a< len*n; a++)
    {
      c1=rand()%tot;
      c2=rand()%tot;
      score_random+=matrix[list[c1]-'a'][list[c2]-'a'];
      score_id+=matrix[list[c1]-'a'][list[c1]-'a'];
    }
  while (tot_len< len*n)
    {
      tot_try++;
      c1=rand()%tot;
      c2=rand()%tot;
      if ( matrix[list[c1]-'a'][list[c2]-'a']>=0){score_good+=matrix[list[c1]-'a'][list[c2]-'a']; tot_len++;}
    }
      

  score_random=score_random/tot_len;
  score_id=score_id/tot_len;
  score_good=score_good/tot_len;
  
  fprintf ( stderr, "Random=%8.3f Id=%8.3f Good=%8.3f [%7.2f]\n",score_random, score_id, score_good, tot_len/tot_try);
  
  return score_random;
}
float compare_two_mat (char  *mat1,char*mat2, int n, int len,char *alp)
{
  int **matrix1, **matrix2;
  
 evaluate_random_match (mat1, n, len,alp);
 evaluate_random_match (mat2, n, len,alp);
 matrix1=read_matrice ( mat1);
 matrix2=read_matrice ( mat2);
 matrix1=rescale_matrix(matrix1, 10, alp);
 matrix2=rescale_matrix(matrix2, 10, alp);
 compare_two_mat_array(matrix1,matrix2, n, len,alp);
 return 0;
} 


int ** rescale_two_mat (char  *mat1,char*mat2, int n, int len,char *alp)
{
  float lambda;
  int **matrix1, **matrix2;

  lambda=measure_lambda2 (mat1, mat2, n, len, alp)*10;

  fprintf ( stderr, "\nLambda=%.2f", lambda);
  matrix2=read_matrice(mat2);
  matrix2=neg_matrix2pos_matrix(matrix2);
  matrix2=rescale_matrix( matrix2, lambda,"abcdefghiklmnpqrstvwxyz");
  
  matrix1=read_matrice(mat1);
  matrix1=neg_matrix2pos_matrix(matrix1);
  matrix1=rescale_matrix( matrix1,10,"abcdefghiklmnpqrstvwxyz");
  
  output_matrix_header ( "stdout", matrix2, alp);
  evaluate_random_match2(matrix1, 1000, 100, alp);
  evaluate_random_match2(matrix2, 1000, 100, alp);
  compare_two_mat_array(matrix1,matrix2, n, len,alp);

  return matrix2;
}  
float measure_lambda2(char  *mat1,char*mat2, int n, int len,char *alp) 
{
  int **m1, **m2;
  float f1, f2;

  m1=read_matrice (mat1);
  m2=read_matrice (mat2);

  m1=neg_matrix2pos_matrix(m1);
  m2=neg_matrix2pos_matrix(m2);
  
  f1=measure_matrix_pos_avg( m1, alp);
  f2=measure_matrix_pos_avg( m2, alp);
  
  return f1/f2;
}
  

float measure_lambda (char  *mat1,char*mat2, int n, int len,char *alp) 
{
  int c;
  int **matrix1, **matrix2, **mat;
  float a;
  float best_quality=0, quality=0, best_lambda;
  
  matrix1=read_matrice ( mat1);
  matrix2=read_matrice ( mat2);
  matrix1=rescale_matrix(matrix1, 10, alp);
  matrix2=rescale_matrix(matrix2, 10, alp);

  for (c=0, a=0.1; a< 2; a+=0.05)
    {
      fprintf ( stderr, "Lambda=%.2f\n", a);
      mat=duplicate_int (matrix2,-1,-1);
      mat=rescale_matrix(mat, a, alp);
      quality=compare_two_mat_array(matrix1,mat, n, len,alp);
      quality=MAX((-quality),quality); 
      
      if (c==0 || (best_quality>quality))
	{
	  c=1;
	  fprintf ( stderr, "*");
	  best_quality=quality;
	  best_lambda=a;
	}
      
      
      evaluate_random_match2(mat, 1000, 100, alp);
      evaluate_random_match2(matrix1, 1000, 100, alp); 
      free_int (mat, -1);
    }
  
  return best_lambda;
  
}

float compare_two_mat_array (int  **matrix1,int **matrix2, int n, int len,char *alp)
{
  int a, b, c, d, c1, c2, tot;
  static int *list;
  static float *freq;
  float delta_random=0;
  float delta2_random=0;

  float delta_id=0;
  float delta2_id=0;

  float delta_good=0;
  float delta2_good=0;

  float delta;

  float tot_len=0;
  float tot_try=0;

 

  if ( !list)
    {
      srand(time(NULL));
      freq=set_aa_frequencies ();
      list=vcalloc ( 10000, sizeof (char));
    }
  
  for (tot=0,c=0,a=0;a<20; a++)
    {
      b=freq[alp[a]-'a']*1000;
      tot+=b;
      for (d=0; d<b; d++, c++)
	{
	  list[c]=alp[a];
	}
    }
  
  
  for (a=0; a< len*n; a++)
    {
      c1=rand()%tot;
      c2=rand()%tot;
      delta=matrix1[list[c1]-'a'][list[c2]-'a']-matrix2[list[c1]-'a'][list[c2]-'a'];
      delta_random+=delta;
      delta2_random+=MAX(delta,(-delta));
      
      delta=matrix1[list[c1]-'a'][list[c1]-'a']-matrix2[list[c1]-'a'][list[c1]-'a'];
      delta_id+=delta;
      delta2_id+=MAX(delta,(-delta));
    }
  while (tot_len< len*n)
    {
      tot_try++;
      c1=rand()%tot;
      c2=rand()%tot;
      if ( matrix1[list[c1]-'a'][list[c2]-'a']>=0 || matrix2[list[c1]-'a'][list[c2]-'a'] )
	{
	  delta=matrix1[list[c1]-'a'][list[c2]-'a']-matrix2[list[c1]-'a'][list[c2]-'a']; 
	  delta_good+=delta;
	  delta2_good+=MAX(delta,(-delta));
	  tot_len++;
	}
    }
      

  delta_random=delta_random/tot_len;
  delta2_random=delta2_random/tot_len;

  
  delta_id=delta_id/tot_len;
  delta2_id=delta2_id/tot_len;

  delta_good=delta_good/tot_len;
  delta2_good=delta2_good/tot_len;
  
    
  fprintf ( stderr, "\tRand=%8.3f %8.3f\n\tId  =%8.3f %8.3f\n\tGood=%8.3f %8.3f\n",delta_random, delta2_random, delta_id,delta2_id, delta_good,delta2_good);
  
  return delta_good;
}



int ** rescale_matrix ( int **matrix, float lambda, char *alp)
{
  int a, b;


  for ( a=0; a< 20; a++)
    for ( b=0; b< 20; b++)
      {
	matrix[alp[a]-'a'][alp[b]-'a']=	matrix[alp[a]-'a'][alp[b]-'a']*lambda;
      }
  return matrix;
}
void output_matrix_header ( char *name, int **matrix, char *alp)
{
  int a, b;
  FILE *fp;
  char *nalp;
  int l;
  nalp=vcalloc ( 1000, sizeof (char));
  

  
  sprintf ( nalp, "ABCDEFGHIKLMNPQRSTVWXYZ");
  l=strlen (nalp);

  fp=vfopen ( name, "w");
  fprintf (fp, "\nint []={\n");
  
  for (a=0; a<l; a++)
    {
      for ( b=0; b<=a; b++)
	fprintf ( fp, "%3d, ",matrix[nalp[a]-'A'][nalp[b]-'A']);
      fprintf ( fp, "\n");
    }
  fprintf (fp, "};\n");
  
  vfclose (fp);
}
  
		  
