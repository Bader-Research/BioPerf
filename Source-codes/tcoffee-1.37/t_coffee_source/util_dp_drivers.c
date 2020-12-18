#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

/******************************************************************/
/*                   MAIN DRIVER                                  */
/*                                                                */
/*                                                                */
/******************************************************************/




Constraint_list *seq2list     ( Constraint_list *CL, char *command, char *seq, char *weight)
    {
	char mode[100];
	static Alignment *A;
	static Constraint_list *PW_CL;
	static Constraint_list *PW_CL_BUF;
	int n, all_pdb, a;
	char *buf;
	static char *mat;
	static char *group_mat;
	static char *dp_mode;
	int ktup;
	static char *weight_buf;
	static char *seq_buf;
	
	

	if ( !weight_buf)weight_buf=vcalloc ( 100, sizeof (int));
	if ( seq_buf) free(seq_buf);
	seq_buf=vcalloc ( strlen (seq)+1, sizeof (char));
	
	if ( !mat) mat=vcalloc ( STRING, sizeof (char));
	if ( !group_mat) group_mat=vcalloc ( STRING, sizeof (char));
	if ( !dp_mode) dp_mode=vcalloc ( STRING, sizeof (char));

	sscanf ( command, "%s", mode);


/*STRUCTURE OR SEQUENCES*/

	all_pdb=1;
	buf=vcalloc ( strlen(seq)+1, sizeof (char));
	sprintf ( buf, "%s", seq);
	n=atoi(strtok (buf, SEPARATORS));

	for ( a=0; a< n; a++)
	  {
	    all_pdb=all_pdb*is_pdb ( (CL->S)->file[atoi(strtok (NULL, SEPARATORS))]);
	  }
	vfree (buf);


/*DNA OR PROTEIN*/	
	if (  strm ( (CL->S)->type, "PROTEIN"))
	   {
	       strcpy(mat, "blosum62mt");
	       strcpy(group_mat, "vasiliky");
	       ktup=2;
	   }
	else if (  strm ( (CL->S)->type, "DNA"))
	   {
	      strcpy(group_mat, "idmat"); 
	      strcpy(mat, "idmat");
	      ktup=5;
	   }  
/*FUNCTION CHOSER*/
/*Proteins*/	

	    
	if ( strm2 ( mode,"fast_pair", "ifast_pair"))
	    {
	      if (! PW_CL)
		   {
		   PW_CL=vcalloc ( 1, sizeof (Constraint_list));
		   PW_CL->pw_parameters_set=1;
		   PW_CL->maximise=1;
		   PW_CL->TG_MODE=1;
		   PW_CL->S=CL->S;
	
		   PW_CL->use_fragments=0;
		   if ( !PW_CL->use_fragments)PW_CL->diagonal_threshold=0;
		   else PW_CL->diagonal_threshold=6;
		   
		   sprintf (PW_CL->dp_mode, "fasta_pair_wise");
		   PW_CL->S=CL->S;
		   PW_CL->DO_S=CL->DO_S;	
		   PW_CL->ktup=ktup;
		   sprintf (PW_CL->matrix_for_aa_group, group_mat);
		   if ( strm ( mode, "fast_pair"))
		       {
		       PW_CL->M=read_matrice (mat);
		       PW_CL->get_dp_cost=slow_get_dp_cost;
		       PW_CL->evaluate_residue_pair=evaluate_matrix_score;
			 
		       PW_CL->extend_jit=0;
		       PW_CL->gop= get_avg_matrix_mm (PW_CL->M, AA_ALPHABET)*10;
		       PW_CL->gep=-1;
		       }
		   else if ( strm ("ifast_pair", mode))
		       {
		       PW_CL->extend_jit=1;
		       PW_CL=CL;	
		       PW_CL->get_dp_cost=slow_get_dp_cost;
		       PW_CL->evaluate_residue_pair=residue_pair_extended_list;
		       }
		   }		   
	    sprintf ( seq_buf, "%s", seq);	    
	    A=fast_pair (A,seq_buf, command+strlen (mode),CL, PW_CL);
	    if ( PW_CL==CL){PW_CL=PW_CL_BUF;}
	    
	    return aln2constraint_list (A, CL,weight);
	    }
	else if ( strm2 ( mode,"diag_fast_pair","idiag_fast_pair"))
	    {
	      if (! PW_CL)
		   {
		   PW_CL=vcalloc ( 1, sizeof (Constraint_list));
		   PW_CL->pw_parameters_set=1;
		   PW_CL->maximise=1;
		   PW_CL->TG_MODE=1;
		   PW_CL->S=CL->S;
	
		   PW_CL->use_fragments=1;
		   PW_CL->diagonal_threshold=3;
		   
		   sprintf (PW_CL->dp_mode, "fasta_pair_wise");
		   PW_CL->S=CL->S;
		   PW_CL->DO_S=CL->DO_S;	
		   PW_CL->ktup=1;
		   sprintf (PW_CL->matrix_for_aa_group, group_mat);
		   PW_CL->M=read_matrice (mat);
		   PW_CL->get_dp_cost=slow_get_dp_cost;
		   PW_CL->evaluate_residue_pair=evaluate_matrix_score;
		   
		   PW_CL->extend_jit=0;
		   PW_CL->gop= get_avg_matrix_mm (PW_CL->M, AA_ALPHABET)*10;
		   PW_CL->gep=-1;		   
		   }		   
	    sprintf ( seq_buf, "%s", seq);	    
	    A=fast_pair (A,seq_buf, command+strlen (mode),CL, PW_CL);
	    if ( PW_CL==CL){PW_CL=PW_CL_BUF;}
	    
	    return aln2constraint_list (A, CL,weight);
	    }
	else if ( strm2 ( mode,"slow_pair","islow_pair" ))
	    {
	      if (! PW_CL)
		   {
		   PW_CL=vcalloc ( 1, sizeof (Constraint_list));
		   PW_CL->pw_parameters_set=1;
		   PW_CL->maximise=1;
		   PW_CL->TG_MODE=1;
		   PW_CL->S=CL->S;
		   PW_CL->use_fragments=0;
		   sprintf (PW_CL->dp_mode, "myers_miller_pair_wise");
		   PW_CL->S=CL->S;
		   PW_CL->DO_S=CL->DO_S;	
		   PW_CL->ktup=ktup;
		   sprintf (PW_CL->matrix_for_aa_group, group_mat);
		   
		   if ( strm ("slow_pair", mode))
		       {
		       PW_CL->M=read_matrice (mat);
		       PW_CL->get_dp_cost=slow_get_dp_cost;
		       PW_CL->evaluate_residue_pair=evaluate_matrix_score;
			 
		       PW_CL->extend_jit=0;
		       PW_CL->gop= get_avg_matrix_mm (PW_CL->M, AA_ALPHABET)*10;
		       PW_CL->gep=-1;
		       }
		   else if ( strm ("islow_pair", mode))
		       {
		       PW_CL->extend_jit=1;
		       PW_CL=CL;	
		       PW_CL->get_dp_cost=slow_get_dp_cost;
		       PW_CL->evaluate_residue_pair=residue_pair_extended_list;
		       }
		   }		   
	    sprintf ( seq_buf, "%s", seq);	    
	    A=fast_pair (A,seq_buf, command+strlen (mode),CL, PW_CL);
	    if ( PW_CL==CL){PW_CL=PW_CL_BUF;}
	    return aln2constraint_list (A, CL,weight);
	    }
/*CDNA*/
	else if ( strm ( mode, "cdna_cfast_pair"))
	    {	
	      if (! PW_CL)
		   {
		   PW_CL=vcalloc ( 1, sizeof (Constraint_list));
		   PW_CL->pw_parameters_set=1;
		   PW_CL->maximise=1;
		   PW_CL->TG_MODE=1;
		   PW_CL->S=CL->S;
		   PW_CL->use_fragments=0;
		   sprintf (PW_CL->dp_mode, "cfasta_cdna_pair_wise");
		   PW_CL->S=CL->S;
		   PW_CL->DO_S=CL->DO_S;	
		   PW_CL->M=read_matrice (strcpy ( mat, "blosum62mt"));
		   PW_CL->extend_jit=0;
		   PW_CL->gop= get_avg_matrix_mm (PW_CL->M, AA_ALPHABET)*10;
		   PW_CL->gep=-1;
		   PW_CL->f_gop=CL->f_gop;
		   PW_CL->f_gep=CL->f_gep;
		   PW_CL->get_dp_cost=get_dp_cost;
		   PW_CL->evaluate_residue_pair=evaluate_cdna_matrix_score;
		   PW_CL->ktup=1;
		   PW_CL->get_dp_cost (NULL,NULL,0,NULL,0,NULL,0,NULL,0,NULL);
		   }
	      sprintf ( seq_buf, "%s", seq);	      
	      A=fast_pair (A,seq_buf, command+strlen (mode),CL, PW_CL);
	      if ( PW_CL==CL){PW_CL=PW_CL_BUF;}
	      return aln2constraint_list (A, CL,"cdna");	      
	    }
	else if ( strm ( mode, "cdna_fast_pair") ||  strncmp (mode,"cdna_fast_pair",14)==0)
	    {	
	      if (! PW_CL)
		   {
		   PW_CL=vcalloc ( 1, sizeof (Constraint_list));
		   PW_CL->pw_parameters_set=1;
		   PW_CL->maximise=1;
		   PW_CL->TG_MODE=1;
		   PW_CL->S=CL->S;
		   PW_CL->use_fragments=0;
		   sprintf (PW_CL->dp_mode, "fasta_cdna_pair_wise");
		   PW_CL->S=CL->S;
		   PW_CL->DO_S=CL->DO_S;	
		   PW_CL->M=read_matrice (strcpy ( mat, "blosum62mt"));
		   PW_CL->extend_jit=0;
		   PW_CL->gop= get_avg_matrix_mm (PW_CL->M, AA_ALPHABET)*10;
		   PW_CL->gep=-1;

		   if ( !strm ( mode, "cdna_fast_pair"))
		     {
		       sscanf (mode+14, "%d_%d", &PW_CL->f_gop, &PW_CL->f_gep);
		       PW_CL->f_gop*=(-1);
		       PW_CL->f_gep*=(-1);
		     }
		   else 
		     {
		       PW_CL->f_gop=2*PW_CL->gop;
		       PW_CL->f_gep=2*PW_CL->gep;
		     }
		   
		   PW_CL->get_dp_cost=get_dp_cost;
		   PW_CL->evaluate_residue_pair=evaluate_cdna_matrix_score;
		   PW_CL->ktup=1;
		   PW_CL->get_dp_cost (NULL,NULL,0,NULL,0,NULL,0,NULL,0,NULL);
		   }
	      sprintf ( seq_buf, "%s", seq);	      
	      A=fast_pair (A,seq_buf, command+strlen (mode),CL, PW_CL);
	      
	     
	      if ( PW_CL==CL){PW_CL=PW_CL_BUF;}
	      return aln2constraint_list (A, CL,"cdna");	      
	    }
/*MULTIPLE ALIGNMENTS*/
	else if ( strm ( mode, "prrp_aln"))
	    {
		A=compute_prrp_aln (A,CL);
		return aln2constraint_list (A, CL,weight);
	    }
	else if ( strm ( mode, "clustalw_aln"))
	    {
		A=compute_clustalw_aln (A,CL);
		return aln2constraint_list (A, CL,weight);
	    }


/*STRUCTURAL METHODS*/
	else if ( strm3 ( mode, "sap_pair", "align_pdb_pair","align_pdb_pair_2" ) && !all_pdb)
	    {
	      sprintf (mode, "fast_pair");
	      CL=seq2list     ( CL,mode, seq, weight);
	      sprintf (mode, "lalign_id_pair");
	      return seq2list     ( CL,mode, seq, weight);
	    }
	      
	else if ( strm ( mode, "sap_pair"))
	    {
		
		return sap_pair (seq, CL);
	    }
	else if ( strm ( mode, "align_pdb_pair"))
	    {
		return align_pdb_pair ( seq, CL);
	    }
	else if ( strm ( mode, "align_pdb_pair_2"))
	    {
		return align_pdb_pair_2 ( seq, CL);
	    }
	else if ( strm ( mode, "custom1_align_pdb_pair"))
	    {
		return custom1_align_pdb_pair ( seq, CL);
	    }
	else if ( strm ( mode, "custom2_align_pdb_pair"))
	    {
		return custom2_align_pdb_pair ( seq, CL);
	    }
	else if ( strm ( mode, "custom3_align_pdb_pair"))
	    {
		return custom3_align_pdb_pair ( seq, CL);
	    }
	else if ( strm ( mode, "custom4_align_pdb_pair"))
	    {
		return custom4_align_pdb_pair ( seq, CL);
	    }
	else if ( strm ( mode, "custom5_align_pdb_pair"))
	    {
		return custom5_align_pdb_pair ( seq, CL);
	    }
	else if ( strm ( mode, "custom6_align_pdb_pair"))
	    {
		return custom6_align_pdb_pair ( seq, CL);
	    }
	else if ( strm ( mode, "custom7_align_pdb_pair"))
	    {
		return custom7_align_pdb_pair ( seq, CL);
	    }
	else if ( strm ( mode, "custom8_align_pdb_pair"))
	    {
		return custom8_align_pdb_pair ( seq, CL);
	    }
	
	else if ( strm ( mode, "custom9_align_pdb_pair"))
	    {
		return custom9_align_pdb_pair ( seq, CL);
	    }
	
	else if ( strm ( mode, "custom10_align_pdb_pair"))
	    {
		return custom10_align_pdb_pair ( seq, CL);
	    }
	else
	    {
		fprintf ( CL->local_stderr, "\nWARNING: THE FUNCTION %s DOES NOT EXIST", mode);
		return CL;
	    }
	    
	}
/******************************************************************/
/*                   MULTIPLE ALIGNMENTS                          */
/*                                                                */
/*                                                                */
/******************************************************************/
Alignment * compute_prrp_aln (Alignment *A, Constraint_list *CL)
   {
       char *tmpseq=NULL;
       char *tmpaln=NULL;
       char command[10000];
       Sequence *S;

       tmpseq=vtmpnam(tmpseq);
       tmpaln=vtmpnam(tmpaln);
       

       A=seq2aln ( CL->S, A, 1);
       output_gotoh_seq (tmpseq, A);
       sprintf ( command, "prrp -E/dev/null -o%s -F9 %s >/dev/null", tmpaln, tmpseq);
       system (command);
       if (!check_file_exists(tmpaln)){return NULL;}
       S=get_fasta_sequence (tmpaln, NULL);

       S->contains_gap=0;
       A=seq2aln(S, A,0);
       free_sequence (S, S->nseq);
       
       return A;
   }
Alignment * aln2clustalw_aln (Alignment *B, Constraint_list *CL)
{
     char *tmpseq=NULL;
     char *tmpaln=NULL;
     char command[10000];
     Alignment *A=NULL;
     FILE *fp;
     int a;

     

     A=copy_aln (B,A); 
     tmpseq=vcalloc (100, sizeof (char));
     tmpaln=vcalloc (100, sizeof (char));

     sprintf ( tmpseq, "tmp%d", rand()%10000);
     sprintf ( tmpaln, "tmp%d", rand()%10000);
     
     fp=vfopen (  tmpseq, "w");
     for ( a=0; a< A->nseq; a++)
       {
	 ungap ( A->seq_al[a]);
	 fprintf ( fp, ">%s\n%s\n", A->name[a], A->seq_al[a]);
       }
     vfclose ( fp);
     
     sprintf ( command, "clustalw %sinfile=%s %soutorder=input %soutfile=%s %s", CWF,tmpseq, CWF, tmpaln, CWF, TO_NULL_DEVICE);
    
     system (command);
     if (!check_file_exists(tmpaln))return NULL;
     
    
     
     A->nseq=0;
     A=main_read_aln(tmpaln,A);
     
     for ( a=0; a< B->nseq; a++)
       sprintf (B->seq_al[a], "%s", A->seq_al[a]);
     B->len_aln=A->len_aln;
     
     
     A=free_aln (A);
     
     remove( tmpseq);
     remove (tmpaln);
     free ( tmpseq);
     free ( tmpaln);
     return B;
   }  
Alignment * compute_clustalw_aln (Alignment *A, Constraint_list *CL)
   {
       char *tmpseq=NULL;
       char *tmpaln=NULL;
       char command[10000];
     

       tmpseq=vcalloc (100, sizeof (char));
       tmpaln=vcalloc (100, sizeof (char));

       sprintf ( tmpseq, "tmp%d", rand()%10000);
       sprintf ( tmpaln, "tmp%d", rand()%10000);
       

       A=seq2aln ( CL->S, A, 1);
       output_fasta_seq (tmpseq, A);
       
       sprintf ( command, "clustalw %sinfile=%s %soutfile=%s %s",CWF, tmpseq,CWF, tmpaln,TO_NULL_DEVICE);
       
       system (command);
       if (!check_file_exists(tmpaln))return NULL;
       A->nseq=0;
       A=main_read_aln(tmpaln,A);
       
       remove( tmpseq);
       remove (tmpaln);
       free ( tmpseq);
       free ( tmpaln);
       return A;
   }  


/******************************************************************/
/*                  DNA                                           */
/*                                                                */
/*                                                                */
/******************************************************************/


/******************************************************************/
/*                   STRUCTURES                                   */
/*                                                                */
/*                                                                */
/******************************************************************/


	


Constraint_list * align_pdb_pair_2 (char *seq, Constraint_list *CL)
        {
	    char *tmp_name=NULL;
	    int s1, s2;
	    
	    
	    static char *command;
	    static char *program;
	    
	    tmp_name=vtmpnam ( NULL);
	    
	    if ( !program)program=vcalloc ( LONG_STRING, sizeof (char));
	    if ( !command)command=vcalloc ( LONG_STRING, sizeof (char));
	    
#ifndef     ALIGN_PDB_4_TCOFFEE
	    if ( getenv ( "ALIGN_PDB_4_TCOFFEE")==NULL)crash ("ALIGN_PDB_4_TCOFFEE IS NOT DEFINED");
	    else  sprintf ( program, "%s", (getenv ( "ALIGN_PDB_4_TCOFFEE")));
#else
	    if ( getenv ( "ALIGN_4_TCOFFEE")==NULL)sprintf (program, "%s", ALIGN_PDB_4_TCOFFEE);
	    else sprintf ( program, "%s", (getenv ( "ALIGN_PDB_4_TCOFFEE")));
#endif

	    atoi(strtok (seq,SEPARATORS));
	    s1=atoi(strtok (NULL,SEPARATORS));
	    s2=atoi(strtok (NULL,SEPARATORS));

	    sprintf ( command , "%s -in P%s P%s -gapopen=-40 -max_delta=2.5 -gapext=0 -scale=0 -hasch_mode=hasch_ca_trace_bubble -maximum_distance=10 -output pdb_constraint_list -outfile stdout> %s",program, (CL->S)->file[s1], (CL->S)->file[s2], tmp_name);	    
	    
	    system  ( command);
	    CL=read_constraint_list_file(CL, tmp_name);
	   

	    remove ( tmp_name);
	    free ( tmp_name);

	    return CL;
	}

Constraint_list * align_pdb_pair (char *seq, Constraint_list *CL)
        {
	    char *tmp_name=NULL;
	    int s1, s2;
	    
	    
	    static char *command;
	    static char *program;
	    
	    tmp_name=vtmpnam ( NULL);
	    
	    if ( !program)program=vcalloc ( LONG_STRING, sizeof (char));
	    if ( !command)command=vcalloc ( LONG_STRING, sizeof (char));
	    
#ifndef     ALIGN_PDB_4_TCOFFEE
	    if ( getenv ( "ALIGN_PDB_4_TCOFFEE")==NULL)crash ("ALIGN_PDB_4_TCOFFEE IS NOT DEFINED");
	    else  sprintf ( program, "%s", (getenv ( "ALIGN_PDB_4_TCOFFEE")));
#else
	    if ( getenv ( "ALIGN_4_TCOFFEE")==NULL)sprintf (program, "%s", ALIGN_PDB_4_TCOFFEE);
	    else sprintf ( program, "%s", (getenv ( "ALIGN_PDB_4_TCOFFEE")));
#endif

	    atoi(strtok (seq,SEPARATORS));
	    s1=atoi(strtok (NULL,SEPARATORS));
	    s2=atoi(strtok (NULL,SEPARATORS));

	    sprintf ( command , "%s  -in P%s P%s -output pdb_constraint_list -outfile stdout> %s",program, (CL->S)->file[s1], (CL->S)->file[s2],tmp_name);	    
	    
	    system  ( command);
	    
	    
	    CL=read_constraint_list_file(CL, tmp_name);
	    remove ( tmp_name);
	    free ( tmp_name);
	
	    return CL;
	}

Constraint_list * custom1_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function1");
	}
Constraint_list * custom2_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function2");
	}
Constraint_list * custom3_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function3");
	}
Constraint_list * custom4_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function4");
	}
Constraint_list * custom5_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function5");
	}
Constraint_list * custom6_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function6");
	}
Constraint_list * custom7_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function7");
	}
Constraint_list * custom8_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function8");
	}
Constraint_list * custom9_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function9");
	}
Constraint_list * custom10_align_pdb_pair (char *seq, Constraint_list *CL)
        {
	  return custom_align_pdb_pair (seq, CL, "custom_pair_score_function10");
	}


Constraint_list * custom_align_pdb_pair (char *seq, Constraint_list *CL, char *mode)
        {
	    char *tmp_name=NULL;
	    int s1, s2;
	    
	    
	    static char *command;
	    static char *program;
	    
	    tmp_name=vtmpnam ( NULL);
	    
	    if ( !program)program=vcalloc ( LONG_STRING, sizeof (char));
	    if ( !command)command=vcalloc ( LONG_STRING, sizeof (char));
	    
#ifndef     ALIGN_PDB_4_TCOFFEE
	    if ( getenv ( "ALIGN_PDB_4_TCOFFEE")==NULL)crash ("ALIGN_PDB_4_TCOFFEE IS NOT DEFINED");
	    else  sprintf ( program, "%s", (getenv ( "ALIGN_PDB_4_TCOFFEE")));
#else
	    if ( getenv ( "ALIGN_4_TCOFFEE")==NULL)sprintf (program, "%s", ALIGN_PDB_4_TCOFFEE);
	    else sprintf ( program, "%s", (getenv ( "ALIGN_PDB_4_TCOFFEE")));
#endif

	    atoi(strtok (seq,SEPARATORS));
	    s1=atoi(strtok (NULL,SEPARATORS));
	    s2=atoi(strtok (NULL,SEPARATORS));

	    sprintf ( command , "%s  -in P%s P%s -output pdb_constraint_list -hasch_mode=%s -outfile stdout> %s",program, (CL->S)->file[s1], (CL->S)->file[s2],mode, tmp_name);	    
	    system  ( command);
	    
	    
	    CL=read_constraint_list_file(CL, tmp_name);
	    remove ( tmp_name);
	    free ( tmp_name);
	
	    return CL;
	}
Constraint_list * sap_pair (char *seq, Constraint_list *CL)
        {
	    FILE *fp;
	    char *tmp_pdb1, *tmp_pdb2;
	    char *string1, *tmp_name, *string2, *string3, *string4, *string5;
	    float score;
	    int s1, s2;
	    int r1, r2,rs1, l, p1, p2;
	    
	    static CLIST_TYPE *entry;
	    char command[STRING];
	    char program[STRING];

	    if ( !entry)entry=vcalloc ( LIST_N_FIELDS, sizeof ( CLIST_TYPE ));
	    declare_name (string1);
	    declare_name (string2);
	    declare_name (string3);
	    declare_name (string4);
	    declare_name (string5);
	    
	  
	    atoi(strtok (seq,SEPARATORS));
	    s1=atoi(strtok (NULL,SEPARATORS));
	    s2=atoi(strtok (NULL,SEPARATORS));
	    

	    tmp_name=vtmpnam ( NULL);

#ifndef     SAP_4_TCOFFEE
	    if ( getenv ( "SAP_4_TCOFFEE")==NULL)crash ("SAP_4_TCOFFEE IS NOT DEFINED");
	    else  sprintf ( program, "%s", (getenv ( "SAP_4_TCOFFEE")));
#else
	    if ( getenv ( "SAP_4_TCOFFEE")==NULL)sprintf (program, "%s", SAP_4_TCOFFEE);
	    else sprintf ( program, "%s", (getenv ( "SAP_4_TCOFFEE")));
#endif


	    /*prepare the two files*/
	    tmp_pdb1=normalize_pdb_file((CL->S)->file[s1], vtmpnam (NULL));
	    tmp_pdb2=normalize_pdb_file((CL->S)->file[s2], vtmpnam (NULL));
	    

	    sprintf ( command , "%s %s %s > %s",program,tmp_pdb1,tmp_pdb2, tmp_name);	    
	    system  ( command);
	    
	    
	    rs1=2;
	    fp=find_token_in_file ( tmp_name, NULL, "Percent");
	    fp=find_token_in_file ( tmp_name, fp  , "Percent");
	    while ( (fgetc (fp))!='\n');
	    while ( fscanf (fp, "%s %s %s %s %s\n", string1, string3, string4, string5, string2)==5 && !strm ( string3, "Weighted") &&  !strm ( string1, "Weighted"))
	          {
		    
		      p1=atoi(string3);
		      score=atof(string4);
		      p2=atoi(string5);
		      l=strlen (string1);
		      r1=tolower(string1[l-1]);
		      r2=tolower(string2[0]);
		      
		      
		      if ( p1> (CL->S)->len[s2]){rs1=1; break;}
		      else if ( p2> (CL->S)->len[s1]){rs1=1; break;}
		      else if ( r1!=tolower ( (CL->S)->seq[s2][p1-1])){rs1=1; break;}
		      else if ( r2!=tolower ( (CL->S)->seq[s1][p2-1])){rs1=1; break;}
		  }
	    vfclose (fp);

	    fp=find_token_in_file ( tmp_name, NULL, "Percent");
	    fp=find_token_in_file ( tmp_name, fp  , "Percent");
	    while ( (fgetc (fp))!='\n');

	    
	    
	    while ( fscanf (fp, "%s %s %s %s %s\n", string1, string3, string4, string5, string2)==5 && !strm ( string3, "Weighted") && !strm ( string1, "Weighted"))
	          {
		      

		      p1=(rs1==1)?atoi(string3):atoi(string5);
		      score=atof(string4);
		      p2=(rs1==1)?atoi(string5):atoi (string3);

		     		      
		      entry[SEQ1]=s1;
		      entry[SEQ2]=s2;
		      entry[R1]=p1;
		      entry[R2]=p2;
		      entry[WE]=(int)score+1;
		      entry[CONS]=1;
		      entry[MISC]=0;
		      CL=add_entry2list(entry, CL);
		      
		  }  

	    vfclose (fp);
	    /*system ( "rm -rf data");*/
	    vfree ( tmp_name);
	    vfree ( string1);
	    vfree ( string2);

	    remove (tmp_pdb1);
	    remove (tmp_pdb2);
	    remove ("super.pdb");

	    free(tmp_pdb1);
	    free(tmp_pdb2);
	    return CL;
	}

/******************************************************************/
/*                   GENERIC PAIRWISE METHODS                     */
/*                                                                */
/*                                                                */
/******************************************************************/


Alignment * fast_pair      ( Alignment *A, char *seq, char *mode, Constraint_list *CL,Constraint_list *PW_CL)
        {
	    int s, n,a;	    
	    static int **l_s;
	    static int *ns;

	    n=atoi(strtok (seq, SEPARATORS));
	    

	    if (!A) 
		{
		A=declare_aln (CL->S);
		}

	    if ( !ns)
	        {
		ns=vcalloc ( 2, sizeof (int));
		l_s=declare_int (2,(CL->S)->nseq);
		}
	    
	    if ( n!=2){fprintf ( stderr, "\nERROR: fast_pw_aln can only handle two seq at a time [FATAL]\n");crash ("");}
	    
	    for ( a=0; a< n; a++)
	        {
		s=atoi (strtok (NULL, SEPARATORS));
		sprintf ( A->seq_al[a], "%s", (CL->S)->seq[s]);
		sprintf ( A->name[a], "%s", (CL->S)->name[s]);
		A->order[a][0]=s;
		}
	    
	    ns[0]=ns[1]=1;
	    l_s[0][0]=0;
	    l_s[1][0]=1;
	    pair_wise ( A, ns, l_s, PW_CL);

	    A->nseq=n;
	    
	    return A;
	    
	}
Alignment * align_two_sequences ( char *seq1, char *seq2, char *in_matrix, int gop, int gep, char *in_align_mode)
        {
	static Alignment *A;
	static Constraint_list *CL;
	int *ns;
	int **l_s;
	
	char       **seq_array;
	char       **name_array;
	static char *matrix;
	static char *align_mode;

	if (!matrix)matrix=vcalloc ( 100, sizeof (char));
	if (!align_mode)align_mode=vcalloc ( 100, sizeof (char));
	sprintf ( matrix, "%s", in_matrix);
	sprintf ( align_mode, "%s", in_align_mode);
	
	CL=vcalloc ( 1, sizeof (Constraint_list));
	
	CL->pw_parameters_set=1;
	CL->M=read_matrice (matrix);
	CL->matrices_list=declare_char (10, 10);


	if (strstr (in_align_mode, "cdna")) 
	  CL->evaluate_residue_pair=evaluate_cdna_matrix_score;
	else
	  CL->evaluate_residue_pair=evaluate_matrix_score;
	
	CL->get_dp_cost=get_dp_cost;
	CL->extend_jit=0;
	CL->maximise=1;
	CL->gop=gop;
	CL->gep=gep;
	CL->TG_MODE=2;
	sprintf (CL->matrix_for_aa_group, "vasiliky");
	CL->use_fragments=0;
	CL->ktup=5;
	if ( !CL->use_fragments)CL->diagonal_threshold=0;
	else CL->diagonal_threshold=6;

	sprintf (CL->dp_mode, "%s", align_mode);

	seq_array=declare_char ( 2, MAX(strlen(seq1), strlen (seq2))+1);
	sprintf (seq_array[0], "%s",seq1);
	sprintf (seq_array[1],"%s", seq2);
	ungap_array(seq_array,2);
	string_array_lower(seq_array,2);

	name_array=declare_char (2, STRING);
	sprintf ( name_array[0], "A");
	sprintf ( name_array[1], "B");
	  
	
	ns=vcalloc ( 2, sizeof(int));
	l_s=declare_int ( 2, 1);
	ns[0]=ns[1]=1;
	l_s[0][0]=0;
	l_s[1][0]=1;
	
	

	CL->S=fill_sequence_struc(2, seq_array, name_array);
	A=seq2aln(CL->S, NULL, 1);
	ungap (A->seq_al[0]);
	ungap (A->seq_al[1]);
	


	A->score_aln=pair_wise (A, ns, l_s,CL);

	free (ns);
	free_int (l_s, -1);
	free_char (name_array, -1);free_char ( seq_array,-1), free_int (CL->M, -1); free (CL);
	
	return A;
	}

NT_node ** make_tree ( Alignment *A,Constraint_list *CL,int gop, int gep,Sequence *S,  char *tree_file, char *tree_mode, int maximise)
	{
	int a, b;
	NT_node **T;
	int **distances;
	int out_nseq;
	char **out_seq_name;
	char **out_seq;
	Constraint_list *CM;
	
	out_nseq=S->nseq;
	out_seq_name=S->name;
	out_seq=S->seq;
	
	/*slow--->  use the same mode as the progressive aln
	  fast--->  use fasta_pairwise_dp with the same subst matrix as the progressive aln
	  very_fast use fasta_pairwise_dp with pam250
	*/

	
	
	if ( strm2 (tree_mode, "slow","fast"))
           {
	   distances=get_pw_distances (A,CL,gop,gep,out_seq,out_seq_name,out_nseq, tree_file, tree_mode, maximise);
	   }

	else if (  strm(tree_mode, "very_fast"))
	    {
	    CM=declare_constraint_list ( S,NULL, NULL, 0,NULL, NULL);
	    CM->L=NULL;
	    CM->S=CL->S;
	    CM->DO_S=CL->DO_S;
	    CM->local_stderr=stderr;
	    CM->M=(strm(S->type, "DNA"))?read_matrice ("blosum62mt"):read_matrice ("idmat");			    
	    distances=get_pw_distances (A,CM,-10,-1,out_seq,out_seq_name,out_nseq, tree_file, tree_mode, maximise);
	    }
	else if ( strm (tree_mode, "order"))
	   {
	   distances=declare_int ( S->nseq, S->nseq);
	   for ( a=0; a< S->nseq; a++)
	     for ( b=0; b< S->nseq; b++)
	     distances[b][a]=100;
	   }
	else 
	    {
	      fprintf ( stderr, "\n%s is an unknown mode making the Tree [FATAL]\n", tree_mode);
	      exit (1);
	    }
	
	for ( a=0; a< S->nseq; a++)
	  for ( b=0; b< S->nseq; b++)
	    if ( distances[a][b]==0)distances[a][b]=1;
	
	fprintf (CL->local_stderr , "\nMAKE NEIGHBOR JOINING DENDROGRAM\n\t[MODE=%s] [Output dendrogram file=%s]\n",tree_mode, (tree_file)?tree_file:"no");
	T=make_nj_tree (A,distances,gop,gep,out_seq,out_seq_name,out_nseq, tree_file, tree_mode);
		     
	free_int (distances, out_nseq);
	return T;
	}






int ** get_pw_distances ( Alignment *A,Constraint_list *CL,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode, int maximise)
	{
	int a, b, c,ra, rb;
	int *ns, **l_s;
	int **dist;
	int n,len;
	int max_name;
	
	CL=cache_dp_value4constraint_list("cache", CL);
	
	if ( strm2 (tree_mode, "very_fast", "fast"))
	   {
	       sprintf ( CL->dp_mode, "fasta_pair_wise");
	       CL->pw_parameters_set=1;
	       CL->extend_jit=1;
	       CL->maximise=1;
	       CL->gop=gop;
	       CL->gep=gep;
	       CL->TG_MODE=1;
	      
	       sprintf (CL->matrix_for_aa_group, "vasiliky");
	       CL->use_fragments=0;
	       	if ( !CL->use_fragments)CL->diagonal_threshold=0;
		else CL->diagonal_threshold=6;

	       CL->ktup=2;
	       CL->fasta_step=0;
	   }
	if ( strm (tree_mode, "very_fast"))
	  {
	    CL->get_dp_cost=slow_get_dp_cost;
	    CL->evaluate_residue_pair=evaluate_matrix_score;
	  }

	dist=declare_int (out_nseq, out_nseq);
	ns=vcalloc ( 2, sizeof(int));
	l_s=declare_int ( 2, out_nseq);
	
	fprintf ( (CL->local_stderr), "\nCOMPUTE PAIRWISE SIMILARITY USING %s\n", CL->dp_mode);	
	ns[0]=ns[1]=1;

	for ( max_name=0,a=0; a<  (CL->S)->nseq; a++)max_name=MAX(strlen ((CL->S)->name[a]), max_name);
	  

	for ( a=0; a< out_nseq-1; a++)
	    {
	    ra=name_is_in_list ( (CL->DO_S)->name[a], (CL->S)->name, (CL->S)->nseq, 100);
	    for ( b=a+1; b< out_nseq; b++)
	        {
		rb=name_is_in_list ((CL->DO_S)->name[b], (CL->S)->name, (CL->S)->nseq, 100); 
		l_s[0][0]=ra;
		l_s[1][0]=rb;
		A->nseq=out_nseq;
		
		
		
		dist[a][b]=dist[b][a]=pair_wise (A, ns, l_s,CL);
		
		len=strlen ( A->seq_al[ra]);
		for ( n=0,c=0; c<len; c++)n+=( A->seq_al[ra][c]!='-' && A->seq_al[rb][c]!='-')?1:0;
		
		if ( n>0)
		    dist[a][b]=dist[b][a]=(dist[a][b])/n;
		else
		    dist[a][b]=dist[b][a]=0;
		
		ungap ( A->seq_al[ra]);
		ungap ( A->seq_al[rb]);
		
		//	fprintf (CL->local_stderr, "\n\t%-*s %-*s: score=%5d", max_name,(CL->DO_S)->name[a], max_name,(CL->DO_S)->name[b],(CL->normalise)?((dist[a][b]*100)/(CL->normalise*SCORE_K)):(dist[a][b]) );
	        }
	}
	fprintf ((CL->local_stderr), "\n\n");
	free(ns);
	free_int (l_s, 2);
	CL=cache_dp_value4constraint_list("restaure", CL);

	return dist;
	}	
	



	
Alignment *recompute_local_aln (Alignment *A, Sequence *S,Constraint_list *CL, int scale, int gep)
    {
    int **coor;
    int a;
    Alignment *B;
   
    sort_constraint_list (CL, 0, CL->ne);
    coor=declare_int (A->nseq, 3);
    for ( a=0; a< A->nseq; a++)
        {
        coor[a][0]=A->order[a][0];
	coor[a][1]=A->order[a][1]+1;
	coor[a][2]=strlen(S->seq[A->order[a][0]])-coor[a][1];
	}
    B=build_progressive_nol_aln_with_seq_coor(CL,0,0,S,coor,A->nseq);
    A=copy_aln ( B, A);
 
    CL=compact_list (CL, 0,CL->ne, "shrink");
    free_Alignment(B);
    return A;
    }
    

Alignment *build_progressive_nol_aln_with_seq_coor(Constraint_list *CL,int gop, int gep,Sequence *S, int **seq_coor, int nseq)
    {

    static int ** local_coor1;
    static int ** local_coor2;
    if ( local_coor1!=NULL)free_int (local_coor1, -1);
    if ( local_coor2!=NULL)free_int (local_coor2, -1);

    local_coor1=get_nol_seq          ( CL,seq_coor, nseq, S);
    local_coor2=minimise_repeat_coor ( local_coor1, nseq, S);

    return build_progressive_aln_with_seq_coor(CL,gop, gep,S, local_coor2,nseq);
    }


Alignment *build_progressive_aln_with_seq_coor (Constraint_list*CL,int gop, int gep, Sequence *S, int **coor, int nseq)
    {
    Alignment *A=NULL;

    A=seq_coor2aln (S,NULL, coor, nseq);
   
    return build_progressive_aln ( A,CL, gop, gep);
    }

Alignment *est_progressive_aln(Alignment *A, Constraint_list *CL, int gop, int gep)
    {
    int a,n;
    int**group_list; 
    int *n_groups;
    char *seq;
    n_groups=vcalloc ( 2, sizeof (int));
    group_list=declare_int ( 2, A->nseq);
    
    n=A->nseq;

    n_groups[0]=1;
    n_groups[1]=1; 
    group_list[0][0]=0;
    group_list[0][1]=1;
    
    group_list[1][0]=1;
    fprintf ( stderr, "\n");  
    for ( a=1; a<n; a++)
        {
	sprintf ( A->seq_al[1], "%s", A->seq_al[a]);
	fprintf ( stderr, "\t[%30s]->[len=%5d]", A->name[a],strlen ( A->seq_al[0]));	
	pair_wise ( A,n_groups, group_list, CL);
	
	seq=dna_aln2cons_seq(A);
	
	sprintf ( A->seq_al[0], "%s", seq);
	vfree (seq);
	fprintf ( stderr, "\n");  
	}
   
    A->nseq=1;
    return A;
    }

void analyse_seq ( Alignment *A, int s)
   {
     int a, b, c;
     int r;

     int len=0;


     int state;
     int pstate=-1;
     float score=0;
     
     for ( a=0; a< A->len_aln; a++)
       {
	 for ( b=0, c=0; b< s; b++)
	   if ( !is_gap(A->seq_al[b][a])){c=1; break;}
	 
	 r=!is_gap(A->seq_al[s][a]);
	 
	 if      (  r &&  c) state=1;
	 else if ( !r && !c) state=2;
	 else if ( !r &&  c) state=3;
	 else if (  r && !c) state=4;
	 
	 if ( state !=pstate)
	   {
	     score+=len*len;
	     len=0;
	   }
	 len+=r;
	 pstate=state;
       }
     score=score/(float)(((A->S)->len[s]*(A->S)->len[s]));
     fprintf ( stderr, "[%.2f]", score);
     
     return;
   }
  
Alignment *build_progressive_aln(Alignment *A, Constraint_list *CL, int gop, int gep)
    {
    int a,n;
    int**group_list; 
    int *n_groups;
    char dp_mode[100];


    sprintf ( dp_mode, "%s", CL->dp_mode);
    sprintf (CL->dp_mode, "gotoh_pair_wise");

    n_groups=vcalloc ( 2, sizeof (int));
    group_list=declare_int ( 2, A->nseq);
    
    n=A->nseq;
    
    for ( a=0; a<n; a++)ungap(A->seq_al[a]);
    for ( a=1; a<n; a++)
        {
	n_groups[0]=a;
	n_groups[1]=1;
	group_list[0][a-1]=a-1;
	group_list[1][0]  =a;
	
	pair_wise ( A,n_groups, group_list, CL);
	fprintf ( stderr, "\n\t[%d]->[%d]", a,strlen ( A->seq_al[0]));
	}
    fprintf (stderr, "\n");
    free(n_groups);
    free_int ( group_list, -1);
    sprintf (CL->dp_mode, "%s",dp_mode);

    return A;
    }




void tree_aln ( NT_node LT, NT_node RT, Alignment*A, int nseq, Constraint_list *CL)
    {
    int *n_s;
    int ** l_s;
    int a, b;
    int score;
    static int n_groups_done;
    int *list;
    
    
    Alignment *SCORE, *SUB_A;
    
    if ( CL->translation==NULL)
       {
       CL->translation=vcalloc ( (CL->S)->nseq, sizeof (int));
       for ( a=0; a< (CL->DO_S)->nseq; a++)
	   CL->translation[a]=name_is_in_list ( (CL->DO_S)->name[a], (CL->S)->name, (CL->S)->nseq, 100);
       if (nseq>2)fprintf ( CL->local_stderr, "\nPROGRESSIVE_ALIGNMENT [Tree Based]\n");
       else fprintf ( CL->local_stderr, "\nPAIRWISE_ALIGNMENT [No Tree]\n");
       n_groups_done=(CL->S)->nseq;
       }
       
    
    n_s=vcalloc (2, sizeof ( int));
    l_s=declare_int ( 2, nseq);

    
    if ( nseq==2)
       {
       n_s[0]=n_s[1]=1;
       l_s[0][0]=name_is_in_list ((CL->DO_S)->name[0],(CL->S)->name, (CL->S)->nseq, 100);
       l_s[1][0]=name_is_in_list ((CL->DO_S)->name[1],(CL->S)->name, (CL->S)->nseq, 100);
       A->score_aln=score=pair_wise (A, n_s, l_s,CL);
       
       free ( n_s);
       free_int ( l_s, 2); 

       return;
       }
    else
       {
       if ( LT->leaf==1 && RT->leaf==0)
	   tree_aln ( RT->left, RT->right,A, nseq,CL);

       else if ( RT->leaf==1 && LT->leaf==0)
	   tree_aln ( LT->left, LT->right,A,nseq,CL);
 
       else if (RT->leaf==0 && LT->leaf==0)
          {
	  tree_aln ( LT->left, LT->right,A,nseq,CL);
	  tree_aln ( RT->left, RT->right,A,nseq,CL);
	  }

       if ( LT->leaf==1 && RT->leaf==1)
          {
	  
	  n_s[0]=LT->nseq;
	  for ( a=0; a< LT->nseq; a++)
	      {
	      l_s[0][a]=CL->translation[LT->lseq[a]];
	      }
	  if ( LT->nseq==1)LT->group=l_s[0][0];
	 
	  n_s[1]=RT->nseq;
	  for ( a=0; a< RT->nseq; a++)
	      {
	      l_s[1][a]=CL->translation[RT->lseq[a]];
	      }
	  if ( RT->nseq==1)RT->group=l_s[1][0];
	  

	  (RT->parent)->group=n_groups_done++;
	  //  fprintf (CL->local_stderr, "\n\tGroup %3d: [Group %3d (%3d seq)] with [Group %3d (%3d seq)]-->",(RT->parent)->group+1,RT->group+1, n_s[1],LT->group+1, n_s[0]);
	  A->score_aln=score=pair_wise (A, n_s, l_s,CL);
	  
	  list=vcalloc ( n_s[0]+n_s[1], sizeof (int));
	  for (b=0,a=0; a< n_s[0]; a++)list[b++]=l_s[0][a];
	  for (a=0; a< n_s[1]; a++)list[b++]=l_s[1][a];
	  SUB_A=copy_aln (A, NULL);
	  SUB_A=extract_sub_aln (SUB_A,n_s[0]+n_s[1],list);
	  SCORE=main_coffee_evaluate_output (SUB_A,CL, "tcoffee_fast");
	  score=SCORE->score_aln;
	  SUB_A=free_aln (SUB_A);
	  SCORE=free_aln (SCORE);
	  vfree (list);
	  
	  //	  fprintf (CL->local_stderr, "[Score=%4d][Len=%5d]",score, strlen ( A->seq_al[l_s[0][0]]));
	  if (nseq==(n_s[0]+n_s[1])) 
	      {
		//  fprintf (CL->local_stderr, "\n\n");
	      n_groups_done=0;
	       for (b=0, a=0; a< n_s[0]; a++, b++)sprintf ( A->tree_order[b],"%s", (CL->S)->name[l_s[0][a]]);
	       for (a=0; a< n_s[1]     ; a++, b++)sprintf ( A->tree_order[b], "%s",(CL->S)->name[l_s[1][a]]);
	     
	      }
	  }
       if (( LT->parent)->parent !=NULL)(LT->parent)->leaf=(RT->parent)->leaf=1;
       if ( LT->isseq==0)LT->leaf=0;
       if ( RT->isseq==0)RT->leaf=0;
       
       free ( n_s);
       free_int ( l_s, 2);
       return;
       }
    
    }

int pair_wise (Alignment *A, int*ns, int **l_s,Constraint_list *CL )
    {
	/*
	 CL->maximise
	 CL->gop;
	 CL->gep
	 CL->TG_MODE;
	*/
	int score;
	if (! CL->pw_parameters_set)
		   {
		       fprintf ( stderr, "\nERROR pw_parameters_set must be set in pair_wise [FATAL]\n" );crash("");
		   }


	if ( CL->get_dp_cost==NULL)CL->get_dp_cost=get_dp_cost;
	
	if ( CL->pair_wise) score=(CL->pair_wise)(A, ns, l_s, CL);	   

	else if (strm (CL->dp_mode,"fasta_cdna_pair_wise"))        score=fasta_cdna_pair_wise   (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"cfasta_cdna_pair_wise"))       score=cfasta_cdna_pair_wise   (A, ns, l_s, CL);
	
	else if (strm (CL->dp_mode,"gotoh_pair_wise"))            score=gotoh_pair_wise        (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"myers_miller_pair_wise"))     score=myers_miller_pair_wise (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"fasta_pair_wise"))	          score=fasta_gotoh_pair_wise  (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"cfasta_pair_wise"))	          score=cfasta_gotoh_pair_wise  (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"very_fast_pair_wise"))	  score=very_fast_gotoh_pair_wise  (A, ns, l_s, CL);

	else if (strm (CL->dp_mode,"gotoh_pair_wise_sw"))         score=gotoh_pair_wise_sw (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"cfasta_sw_pair_wise"))        score=cfasta_gotoh_pair_wise_sw (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"fasta_sw_pair_wise"))         score=fasta_gotoh_pair_wise_sw (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"gotoh_pair_wise_lalign"))     score=gotoh_pair_wise_lalign (A, ns, l_s, CL);
	
	else if (strm (CL->dp_mode,"domain_pair_wise"))           score=domain_pair_wise (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"ssec_pair_wise"))             score=ssec_pwaln_maln  (A, ns, l_s, CL);
	else if (strm (CL->dp_mode,"default"))                    score=gotoh_pair_wise  (A, ns, l_s, CL);            
	
	else
	    {
		fprintf ( stderr, "\n[%s] is an unknown mode for pair_wise[FATAL]\n", CL->dp_mode);
		crash ( "\n");
	    }
       	
	return  score;
    }

/*******************************************************************************/
/*                                                                             */
/*                                                                             */
/*	Util Functions                                                         */
/*                                                                             */
/*	                                                                       */
/*******************************************************************************/

char *build_consensus ( char *seq1, char *seq2, char *dp_mode)
        {
	Alignment *A;
	char *buf;
	int a;
	char c1, c2;
	static char *mat;

	
	if ( !mat) mat=vcalloc ( STRING, sizeof (char));


	A=align_two_sequences (seq1, seq2, strcpy(mat,"idmat"), 0, 0,dp_mode);
	buf=vcalloc ( A->len_aln+1, sizeof (char));
	
	for ( a=0; a< A->len_aln; a++)
	    {
		c1=A->seq_al[0][a];
		c2=A->seq_al[1][a];
		if (is_gap(c1) && is_gap(c2))buf[a]='-';
		else if (is_gap(c1))buf[a]=c2;
		else if (is_gap(c2))buf[a]=c1;
		else if (c1!=c2){free (buf);buf=NULL;free_aln(A);return NULL;}
		else buf[a]=c1;
	    }
	buf[a]='\0';
	free_aln (A);
	return buf;
	}
