#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

int gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
	{
/*******************************************************************************/
/*                NEEDLEMAN AND WUNSCH (GOTOH)                                 */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
	

/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	    int TG_MODE;
	    int l_gop, l_gep;
	    int gop, gep;
	    int maximise;
/*VARIANLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int a, b, i, j;

	int *cc;	
	int *dd,*ddg;
	int   e, eg;

	int lenal[2], len;
	int t, c,s, ch;
	int sub;
	int fop;
	int score=0;
        int **pos0;
	static char **al;
	char **aln;
	int ala, alb,LEN;
	char *buffer;
	char *char_buf;
/*trace back variables       */
	FILE       *long_trace=NULL;
	TRACE_TYPE *buf_trace=NULL;
	static TRACE_TYPE **trace;
	TRACE_TYPE k;
	TRACE_TYPE *tr;
	int long_trace_flag=0;
	int dim;
/********Prepare penalties*******/
	gop=CL->gop*SCORE_K;
	gep=CL->gep*SCORE_K;
	TG_MODE=CL->TG_MODE;
	maximise=CL->maximise;
	
	
/********************************/	
/*CLEAN UP AFTER USE*/
	if ( A==NULL)
	   {
	   free_int (trace,-1);
	   trace=NULL;
	   free_char (al,-1);
	   al=NULL;
	   return 0;
	   }		

/*DO MEMORY ALLOCATION FOR DP*/

	lenal[0]=strlen (A->seq_al[l_s[0][0]]);
	lenal[1]=strlen (A->seq_al[l_s[1][0]]);
	len= MAX(lenal[0],lenal[1])+1;
	buf_trace=vcalloc ( len, sizeof (TRACE_TYPE));	
	buffer=vcalloc ( 2*len, sizeof (char));	
        al=declare_char (2, 2*len);  
	
	char_buf= vcalloc (2*len, sizeof (char));	
	

	dd = vcalloc (len, sizeof (int));
	

	cc = vcalloc (len, sizeof (int));
	ddg=vcalloc (len, sizeof (int));
	

	
	if ( len>=MAX_LEN_FOR_DP)
	    {
	    long_trace_flag=1;
	    long_trace=vtmpfile();
	    }
	else
	    {
	   
	    dim=(trace==NULL)?0:read_size_int ( trace[-1],0);	   
	    trace    =realloc_int ( trace,dim,dim,MAX(0,len-dim), MAX(0,len-dim));
	    }
	
/*END OF MEMORY ALLOCATION*/
	
	
		/*
		0(s)   +(dd)
		  \      |
		   \     |
		    \    |
		     \   |
		      \  |
		       \ |
		        \|
		-(e)----O
		*/ 
		       
	pos0=aln2pos_simple ( A,-1, ns, l_s);


	cc[0]=0;		
	tr=(long_trace_flag)?buf_trace:trace[0];
	tr[0]=(TRACE_TYPE)1;
	for ( j=1; j<=lenal[1]; j++)tr[j]=(TRACE_TYPE)-1;
	if (long_trace_flag)fwrite (buf_trace, sizeof ( TRACE_TYPE),lenal[1]+1, long_trace);
	
	
	t=(TG_MODE==0)?gop:0;
	

	for (cc[0]=0,j=1; j<=lenal[1]; j++)
	    {
	    
	    l_gop=(TG_MODE==0)?gop:0;
	    l_gep=(TG_MODE==2)?0:gep;

	    cc[j]=t=t+l_gep;
	    dd[j]=  t+  gop;
	    }

	t=(TG_MODE==0)?gop:0;	
	
	for (i=1; i<=lenal[0];i++)
			{			
			tr=(long_trace_flag)?buf_trace:trace[i];
			s=cc[0];

			l_gop=(TG_MODE==0)?gop:0;
			l_gep=(TG_MODE==2)?0:gep;
			
			
			
			cc[0]=c=t=t+l_gep;
			e=t+  gop;
			tr[0]=(TRACE_TYPE)1;

			

			for (eg=0,j=1; j<=lenal[1];j++)
				{				   

				sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);	
				
				/*get the best Insertion*/
				l_gop=(i==lenal[0] || i==1 )?((TG_MODE==0)?gop:0):gop;
				l_gep=(i==lenal[0] || i==1)?((TG_MODE==2)?0:gep):gep;
			

				if ( a_better_than_b ( e,c+l_gop, maximise))eg++;
				else eg=1;	
				e=best_of_a_b (e, c+l_gop, maximise)+l_gep;
				
				/*Get the best deletion*/
				l_gop=(j==lenal[1] || j==1)?((TG_MODE==0)?gop:0):gop;
				l_gep=(j==lenal[1] || j==1)?((TG_MODE==2)?0:gep):gep;
				

				if ( a_better_than_b ( dd[j], cc[j]+l_gop, maximise))ddg[j]++;
				else ddg[j]=1;
				dd[j]=best_of_a_b( dd[j], cc[j]+l_gop,maximise)+l_gep;
				


				c=best_int(3,maximise,&fop, e, s+sub,dd[j]);
				/*Chose Substitution for tie breaking*/
				if ( fop==0 && (s+sub)==e)fop=1;
				else if ( fop==2 && (s+sub)==dd[j])fop=1;
				/*Chose Deletion for tie breaking*/
				else if ( fop==2 && e==dd[j])fop=1;

				fop-=1;
				s=cc[j];
				cc[j]=c;	

	
				if ( fop<0)
					{tr[j]=(TRACE_TYPE)fop*eg;
					}
				else if ( fop>0)
				        {tr[j]=(TRACE_TYPE)fop*ddg[j];
					}
				else if (fop==0)
					{tr[j]=(TRACE_TYPE)0;	
					}					
				fop= -2;
				}
			if (long_trace_flag)
			    {
			    fwrite ( buf_trace, sizeof (TRACE_TYPE), lenal[1]+1, long_trace);
			    }
			}
	
		
	score=c;
	
        i=lenal[0];
	j=lenal[1];
	ala=alb=0;
	

	while (i>=0 && j>=0 && ((i+j)!=0))
			{
			if ( i==0)
				k=-1;
			else if ( j==0)
				k=1;
			else if ( j==0 && i==0)
				k=1;	
			else
			        {
				if (long_trace_flag)
				   {
				   fseek ( long_trace, sizeof (TRACE_TYPE)*((lenal[1]+1)*(i)+j),SEEK_SET);
				   fread ( &k, sizeof (TRACE_TYPE), 1, long_trace);
				   }
				else
				   {
				   
				   k=trace[i][j];
				   }
				}
				
				
			if (k==0)
				{
				
				al[0][ala++]=1;
				al[1][alb++]=1;
				i--;
				j--;
				}		
			else if (k>0)
				{
				
				for ( a=0; a< k; a++)
					{
					al[0][ala++]=1;
					al[1][alb++]=0;
					i--;
					}
				}
			else if (k<0)
				{
				
				for ( a=0; a>k; a--)
					{
					al[0][ala++]=0;
					al[1][alb++]=1;
					j--;
					}
				}
			}
      
	LEN=ala;	
	c=LEN-1;  
	
	

	invert_list_char ( al[0], LEN);
	invert_list_char ( al[1], LEN);	
	if ( A->declared_len<=LEN)A=realloc_aln2  ( A,A->max_n_seq, 2*LEN);	
	aln=A->seq_al;

	for ( c=0; c< 2; c++)
	    {
	    for ( a=0; a< ns[c]; a++) 
		{		
		ch=0;
		for ( b=0; b< LEN; b++)
		    {		   
		    if (al[c][b]==1)
			char_buf[b]=aln[l_s[c][a]][ch++];
		    else
			char_buf[b]='-';
		   }
		char_buf[b]='\0';
		sprintf (aln[l_s[c][a]],"%s", char_buf);
	        }
	     }
	
	
	A->len_aln=LEN;
	A->nseq=ns[0]+ns[1];
	

	free ( cc);
	free (dd);		
	free (ddg);
	free (buffer);
	free (char_buf); 
	free (buf_trace);
	free_char ( al, -1);
	if ( long_trace_flag)fclose (long_trace);	



	return score;
	}
     

