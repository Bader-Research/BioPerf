#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

#include "dp_lib_header.h"



void print_atom ( Atom*A);


int evaluate_ca_trace_nb (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
    
     return (int)(neighborhood_match(CL, s1,r1, s2, r2, (CL->T[s1])->Chain,(CL->T[s2])->Chain ));
   }
int evaluate_ca_trace_sap2_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       
       
       
       return sap2_neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Bubble,(CL->T[s2])->Bubble );
       
     }
int evaluate_ca_trace_sap1_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       /*
	 Function documentation: start
	 
	 int evaluate_ca_trace_sap1_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2)
	 This function evaluates the cost for matching two residues:
	 
	 a1 is the cost for matching the two neighborood ( bubble type), using sap
	 a1: [0,+100], +100 is the best possible match.
	 a2 is the residue type weight:
	    min=worst substitution value
	    best=best of r1/r1, r2/r2-min

	    a2=(r1/r2 -min)/best --> a1:[0, 100]
	 
	 score=a1*a2-->[-inf, +10000];
       */
	 
	 
       
       float a1;
       
       
       a1=(int) sap1_neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Bubble,(CL->T[s2])->Bubble );
       
       return (int)a1;
       

     }
int evaluate_ca_trace_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       /*
	 Function documentation: start
	 
	 int evaluate_ca_trace_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2)
	 This function evaluates the cost for matching two residues:
	 
	 a1 is the cost for matching the two neighborood ( bubble type)
	 a1: [-inf,+100-scale], +100-scale is the best possible match.
	 
       */
	 
	 
       
       float a1;

       
            
       a1=(int) neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Bubble,(CL->T[s2])->Bubble )-((CL->T[s1])->pdb_param)->scale;
      
       return a1;
       

     }
int evaluate_ca_trace_transversal (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       return (int)(transversal_match (CL, s1, r1, s2, r2, (CL->T[s1])->Transversal,(CL->T[s2])->Transversal ));
     }

int evaluate_ca_trace_bubble_3 (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       /*This Mode evaluates :
	 
	 1-The Bubble
	 2-The Match of the transversal residues
       */

       int a1, l1;
       int a2, l2;
       int a;
       
       l1=MAX(( (CL->T[s1])->Chain )->nb[r1][0] ,((CL->T[s2])->Chain )->nb[r2][0]);
       l2=MAX(( (CL->T[s1])->Bubble)->nb[r1][0], ((CL->T[s2])->Bubble)->nb[r2][0]);
       
       a1=(int)(neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Bubble,(CL->T[s2])->Bubble ));
       a2=(int)(transversal_match (CL, s1, r1, s2, r2, (CL->T[s1])->Transversal,(CL->T[s2])->Transversal ));
       
       if ( !l1 && !l2)return 0;
       a=(a1+a2)/2;
       return a;
     }
int evaluate_ca_trace_bubble_2 (Constraint_list *CL, int s1, int r1, int s2, int r2)
     {
       /*This Mode evaluates :
	 1-The Ca neighborhood
	 2-The Bubble
       */

      
       return (int)((neighborhood_match (CL, s1, r1, s2, r2, (CL->T[s1])->Chain,(CL->T[s2])->Chain )));
     }


/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:START                                   */
/*                                                                                           */
/*********************************************************************************************/
float matrix_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)

     {
       /*
	 Function documentation: start
	 
	 float matrix_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
	 This function evaluates the matrix for matching two residues:
	
	    min=worst substitution value
	    best=best of r1/r1, r2/r2-min

	    a2=(r1/r2 -min)/best --> a1:[0, 100]
	 
	 score=a1*a2-->[-inf, +10000];
       */
	 
	 
       
       float a2;
       float m1, m2, m;
       static float min=0;
       int a, b;

       if ( !CL->M) 
	 {
	   CL->M=read_matrice ( "pam250mt");
	   min=CL->M[0][0];
	   for ( a=0; a< 26; a++)
	     for ( b=0; b< 26; b++)min=MIN(min, CL->M[a][b]);
	 }

       if ( r1<=0 || r2<=0)return 0;
       m1=CL->M[(CL->S)->seq[s1][r1-1]-'a'][(CL->S)->seq[s1][r1-1]-'a']-min;
       m2=CL->M[(CL->S)->seq[s2][r2-1]-'a'][(CL->S)->seq[s2][r2-1]-'a']-min;
       m=MAX(m1, m2);
       a2=(CL->M[(CL->S)->seq[s1][r1-1]-'a'][(CL->S)->seq[s2][r2-1]-'a']-min)/m;
     
       return a2;
     }

     
float transversal_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
      {
	int a, l1, l2;
	float score=0;
	float delta, max_delta;
	float max;
	Pdb_param*PP;
	
	PP=(CL->T[s1])->pdb_param;
	max_delta=PP->max_delta;

	l1=nbs1->nb[r1][0];
	l2=nbs2->nb[r2][0];

	if ( l1!=l2 || l1<(PP->N_ca)) return 0;
	

	max=MAX(l1, l2)*max_delta;
	for ( delta=0,a=0; a< l2 ; a++)
	       {
		   
		   delta+=max_delta-FABS((nbs1->d_nb[r1][a]-nbs2->d_nb[r2][a]));
	       }
	score=(delta*100)/max;
       

              
       return score;
      }
	
float neighborhood_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
      {
	  static float **table;
	  static int table_size;
	  int a, b, l1, l2;
	  float score=0;
	  float ins, del, sub;
	  float delta, max_delta;
	  float max;
	  Pdb_param*PP;

       
       PP=(CL->T[s1])->pdb_param;
       max_delta=PP->max_delta;

       
       if ( r1> 0 && r2 >0) {r1--; r2--;}
       else return 0;
       
       l1=nbs1->nb[r1][0];
       l2=nbs2->nb[r2][0];

       if (table_size< (MAX(l1, l2)+1))
	 {
	 table_size=MAX(l1, l2)+1;
	 if ( table)free_float (table, -1);
	 table=NULL;
	 }
       if ( !table) table=declare_float (table_size, table_size);
       
       	 
       max=MAX(l1, l2)*max_delta;
       if ( max==0)return 0;

       
       table[0][0]=0;
       for ( b=1; b<=l2; b++)
           { 
	   table[0][b]=0;
	   }
       for ( a=1; a<=l1; a++)
           {
	   table[a][0]=0;
	   for ( b=1; b<=l2 ; b++)
	       {
		   
		   delta=max_delta-FABS((nbs1->d_nb[r1][a]-nbs2->d_nb[r2][b]));
		   		   
		   del=table[a-1][b];
		   ins=table[a][b-1];
		   sub= table[a-1][b-1]+delta;

		   if ( del >= ins && del >= sub)score=del;
		   else if ( ins >= del && ins >= sub) score=ins;
		   else score=sub;		   
		   table[a][b]=score;
	       }
	   }
      

       score=((((score)*100)/max));
       
              
       return score;
      }

float sap1_neighborhood_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
      {
	/*
	  Function documentation: start
	  
	  float sap1_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
	  This function is adapted from Taylor, Orengo, Protein Structure Alignment JMB 1989, (208)1-22
	  It is the first function where 
	  score= A/(|dra-drb|+b)
	  
	  Function documentation: end
	*/
	  
	  static float **table;
	  static int table_size;
	  int a, b, l1, l2;
	  float score=0;
	  float ins, del, sub;
	  float delta;
	  float max;

	  int A=50;
	  int B=5;
	
       

       

       
       if ( r1> 0 && r2 >0) {r1--; r2--;}
       else return 0;
       
       l1=nbs1->nb[r1][0];
       l2=nbs2->nb[r2][0];

       if (table_size< (MAX(l1, l2)+1))
	 {
	 table_size=MAX(l1, l2)+1;
	 if ( table)free_float (table, -1);
	 table=NULL;
	 }
       if ( !table) table=declare_float (table_size, table_size);
       
       	 
       max=MAX(l1, l2)*(A/B);
       if ( max==0)return 0;

       
       table[0][0]=0;
       for ( b=1; b<=l2; b++)
           { 
	   table[0][b]=0;
	   }
       for ( a=1; a<=l1; a++)
           {
	   table[a][0]=0;
	   for ( b=1; b<=l2 ; b++)
	       {
	
		   delta=A/(FABS((nbs1->d_nb[r1][a]-nbs2->d_nb[r2][b]))+B);
	
		   del=table[a-1][b];
		   ins=table[a][b-1];
		   sub= table[a-1][b-1]+delta;

		   if ( del >= ins && del >= sub)score=del;
		   else if ( ins >= del && ins >= sub) score=ins;
		   else score=sub;		   
		   table[a][b]=score;
	       }
	   }
      

       score=((score*100))/(max);
       
              
       return score;
      }

float sap2_neighborhood_match (Constraint_list *CL, int s1, int r1, int s2, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
      {
	/*
	  Function documentation: start
	  
	  float sap1_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2)
	  This function is adapted from Taylor, Orengo, Protein Structure Alignment JMB 1989, (208)1-22
	  It is the first function where 
	  score= A/(|dra-drb|+b)
	  
	  Function documentation: end
	*/
	  
	  static float **table;
	  static int table_size;
	  int a, b, l1, l2;
	  float score=0;
	  float ins, del, sub;
	  float delta;
	  float max;

	  Amino_acid **pep1;
	  Amino_acid **pep2;
	  static Atom *vX_1, *vY_1, *vZ_1;
	  static Atom *vX_2, *vY_2, *vZ_2;
	  static Atom *ca1, *ca2;
	  float val;
 
	  int A=50;
	  int B=2;
	  
       


       if ( r1> 0 && r2 >0) {r1--; r2--;}
       else return 0;

       /*Make up the referencial*/
       pep1=(CL->T[s1])->peptide_chain;
       pep2=(CL->T[s2])->peptide_chain;

       /*Get Referencial for CA1*/ 
       if ( (pep1[r1])->C)vX_1 =diff_atom(pep1[r1]->C,pep1[r1]->CA, vX_1);
       if ( (pep1[r1])->N)vY_1 =diff_atom(pep1[r1]->N,pep1[r1]->CA, vY_1);
       if ( (pep1[r1])->CB)vZ_1=diff_atom(pep1[r1]->CB,(pep1[r1])->CA,vZ_1);
       else vZ_1=add_atom (vX_1, vY_1, vZ_1); 

       
      
      

       /*Get Referencial for CA2*/ 
       if ( (pep2[r2])->C)vX_2 =diff_atom((pep2[r2])->C,(pep2[r2])->CA, vX_2);
       if ( (pep2[r2])->N)vY_2 =diff_atom((pep2[r2])->N,(pep2[r2])->CA, vY_2);
       if ( (pep2[r2])->CB)vZ_2=diff_atom((pep2[r2])->CB,(pep2[r2])->CA, vZ_2);
       else vZ_2=add_atom (vX_2, vY_2, vZ_2); 
       
      
      

       /*END OF GETTING REFERENCIAL*/ 
	 
       /*Test
       if ( r1>1 && r2>1)
	 {
	 fprintf (stdout,"\n\t*******");
	 
	 fprintf (stdout, "RESIDUE %d %c", r1, (CL->S)->seq[s1][r1]);
	 if ( (pep1[r1])->CA)fprintf (stdout,"\n\tCA ");print_atom (pep1[r1]->CA );
	 if ( (pep1[r1])->C)fprintf (stdout,"\n\tC  ");print_atom (pep1[r1]->C );
	 if ( (pep1[r1])->N)fprintf (stdout,"\n\tN  ");print_atom (pep1[r1]->N );
	 if ( (pep1[r1])->CB)fprintf (stdout,"\n\tCB ");print_atom (pep1[r1]->CB );
	 fprintf (stdout,"\n\t*******");
	 fprintf (stdout,"\n\tvX ");print_atom ( vX_1);
	 fprintf (stdout,"\n\tvY ");print_atom ( vY_1);
	 fprintf (stdout,"\n\tvZ ");print_atom ( vZ_1);

	 ca1= copy_atom ((pep1[r1-1])->CA, ca1);
	 ca1 =diff_atom(ca1,(pep1[r1])->CA, ca1);
	 fprintf (stdout,"\n\tca ");print_atom ( ca1);
	 fprintf ( stdout, "\n\tSQ1=%d ", (int)square_atom(ca1));
	 ca1=reframe_atom(vX_1, vY_1, vZ_1, ca1, ca1);
	 fprintf ( stdout, "\n\tSQ2=%d ", (int)square_atom(ca1));
	 fprintf (stdout,"\n\tca ");print_atom ( ca1);
	 fprintf (stdout,"\n\n");
	 }
       */

       l1=nbs1->nb[r1][0];
       l2=nbs2->nb[r2][0];

       if (table_size< (MAX(l1, l2)+1))
	 {
	 table_size=MAX(l1, l2)+1;
	 if ( table)free_float (table, -1);
	 table=NULL;
	 }
       if ( !table) table=declare_float (table_size, table_size);
       
       	 
       max=MAX(l1, l2)*(A/B);
      
       if ( max==0)return 0;
       
       
       table[0][0]=0;
       for ( b=1; b<=l2; b++)
           { 
	   table[0][b]=0;
	   }

       for ( a=1; a<=l1; a++)
           {
	   ca1=copy_atom ((CL->T[s1])->structure[nbs1->nb[r1][a]], ca1);
	   ca1=diff_atom(ca1,(pep1[r1])->CA, ca1);
	   ca1=reframe_atom(vX_1, vY_1, vZ_1, ca1, ca1);
	   
	   table[a][0]=0;
	   for ( b=1; b<=l2 ; b++)
	       {
		  ca2  =copy_atom((CL->T[s2])->structure[nbs2->nb[r2][b]], ca2);
		  ca2  =diff_atom(ca2,(pep2[r2])->CA, ca2);
		  ca2  =reframe_atom(vX_2, vY_2, vZ_2, ca2, ca2);
		  
		  ca2=diff_atom(ca2,ca1,ca2);
		  val=square_atom (ca2);
		  
		  val=(float)sqrt ((double)val);
		  
		  delta=A/(val+B);
		  

		  del=table[a-1][b];
		  ins=table[a][b-1];
		  sub= table[a-1][b-1]+delta;

		  if ( del >= ins && del >= sub)score=del;
		  else if ( ins >= del && ins >= sub) score=ins;
		  else score=sub;		   
		  table[a][b]=score;
	       }
	   }
      

       score=(((score*100))/(max)-50);
       
              
       return score;
      }
/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:START                                   */
/*                                                                                           */
/*********************************************************************************************/

Alignment * analyse_pdb ( Alignment *A, Alignment *ST)
      {
	  int a, b;
	  float ***residues;

	  float tot, tot_n=0;
	  float m1, m2;
	  float tot_m1=0, tot_m2=0, tot_m3=0;
	  float seq_tot=0,seq_tot_n=0, seq_tot_m1, seq_tot_m2,seq_tot_m3=0;
	  float col_tot_n=0, col_tot_m2=0;
	  Pdb_param *PP;
	  Constraint_list *CL;
	  int **pos;
	  int p;
	  

	  CL=A->CL;
	  
	  for ( a=0; a< (A->S)->nseq; a++)
	      if ( CL->T[a]){PP=(CL->T[a])->pdb_param;break;}
	      
	  	 
	  residues=analyse_pdb_residues ( A, A->CL,PP);
	  pos=aln2pos_simple (A, A->nseq);

	  if (strpbrk( PP->comparison_io, "h*"))
	     {
		fprintf (stdout, "%s, %s (%s)\n%s\n",PROGRAM,VERSION,DATE, AUTHOR);      
	     }

	  
	  if (strpbrk (PP->comparison_io,"s*"))fprintf ( stdout, "\nINDIVIDUAL RESULTS:");
	  tot=tot_n=0;
	  for ( a=0; a< A->nseq; a++)
	  {

	      seq_tot=seq_tot_n=seq_tot_m1=seq_tot_m2=seq_tot_m3=0;
	      for ( p=0; p< A->len_aln; p++)
	          {
		  
		      if ( !is_gap(A->seq_al[a][p]) && (CL->T[A->order[a][0]]))
			         {
				  seq_tot++;
				  tot++;
				 }
		      
		      if ( pos[a][p]<=0)continue; 
		      b=pos[a][p]-1;
		      
		      if ( !(CL->T[A->order[a][0]]) || is_gap (A->seq_al[a][p]) ||pos[a][p]<=0 || residues[a][b][0]==0)
			  {
			      m2=NO_COLOR_RESIDUE;
			  }
		      else
			  {
			  tot_n++;
			  seq_tot_n++;

			  m1=residues[a] [b][1]/residues[a][b][0];
			  m2=(residues[a][b][2]*100)/residues[a][b][0];
			  seq_tot_m1+=m1;
			  seq_tot_m2+=m2;
			  tot_m1+=m1;
			  tot_m2+=m2;
			  if (m2>PP->similarity_threshold){tot_m3++;seq_tot_m3++;}
			  }
		      (ST)->seq_al[a][p]=(m2==NO_COLOR_RESIDUE)?m2:(MIN(m2/10, 9));	
		  }
	      if ( seq_tot_n==0)continue;   
	      else if (strpbrk( PP->comparison_io, "s*"))
	         {
		   		
		     fprintf (stdout, "\n\t%-6s [Len=%d Res.] (%d)", (A->S)->name[A->order[a][0]],  (A->S)->len[A->order[a][0]],(int)seq_tot_n );
		     if (strpbrk( PP->comparison_io,"0*#"))fprintf (stdout, "\n\t\t%-6s IM0: %6.2f %% Aligned Res.",(A->S)->name[A->order[a][0]], (seq_tot_n*100)/seq_tot);
		     if (strpbrk( PP->comparison_io,"1*#") )fprintf (stdout, "\n\t\t%-6s IM1: %6.2f Angs.",(A->S)->name[A->order[a][0]], seq_tot_m1/seq_tot_n);
		     if (strpbrk( PP->comparison_io,"2*#"))fprintf (stdout, "\n\t\t%-6s IM2: %6.2f %% of the %6.2f Angs. Ngb.",(A->S)->name[A->order[a][0]], seq_tot_m2/seq_tot_n,PP->maximum_distance );
		     if (strpbrk( PP->comparison_io,"3*#"))fprintf (stdout, "\n\t\t%-6s IM3: %6.2f %% of the %d Aligned Res.\n",(A->S)->name[A->order[a][0]], (seq_tot_m3*100)/seq_tot_n, (int)seq_tot_n);
		 }
	      ST->score_seq[a]=((seq_tot_m3*100)/seq_tot_n);
	  }
	  
	  for ( a=0; a< A->len_aln; a++)
	      {
		  col_tot_n=col_tot_m2=0;
		  for ( b=0;b<A->nseq; b++)
		      {
		      if (pos[b][a]<=0)continue;
		      col_tot_n +=residues[b][pos[b][a]-1][0];
		      col_tot_m2+=residues[b][pos[b][a]-1][2];
		      }
		  if (col_tot_n==0)(ST)->seq_al[b][a]=NO_COLOR_RESIDUE;
		  else  
		     {
			 m2=(col_tot_m2*10)/col_tot_n;
			 (ST)->seq_al[b][a]=MIN(9, m2);
		     }
	      }
			     
	  if (tot_n==0){fprintf ( stderr, "\nWARNING: NO RESIDUE WITH KNOWN COORDINATES\nIMPOSSIBLE TO DO ANALYSE [FATAL][%.2f]\n",tot_n);exit(1);}
	  
	  if (strpbrk( PP->comparison_io, "g*"))
	     {
		 fprintf ( stdout, "\nGLOBAL RESULTS:");
		 if (strpbrk(PP->comparison_io,"0#*"))fprintf ( stdout, "\n\t\tGM0: %6.2f %% Aligned Res.", (tot_n*100)/tot);
		 if (strpbrk(PP->comparison_io,"1*#"))fprintf ( stdout, "\n\t\tGM1: %6.2f Angs.", tot_m1/tot_n);
		 if (strpbrk(PP->comparison_io,"2*#"))fprintf ( stdout, "\n\t\tGM2: %6.2f %% of the %6.2f Angs. Ngb.", tot_m2/tot_n,PP->maximum_distance );
		 if (strpbrk(PP->comparison_io,"3*#"))fprintf ( stdout, "\n\t\tGM3: %6.2f %% of the %d Res.\n", (tot_m3*100)/tot_n, (int)tot_n);
	     }
	  if ( strpbrk( "d*", PP->comparison_io))
	      { 
		  fprintf ( stdout, "\nDEFINITIONS:");
		  if (strpbrk(PP->comparison_io,"0#*"))fprintf ( stdout, "\nM0 Proportion of aligned residues");
		  if (strpbrk(PP->comparison_io,"1*#"))fprintf ( stdout, "\nM1 average RMSD of Ca within a %.2f Angs. Ngb.",PP->maximum_distance);
		  if (strpbrk(PP->comparison_io,"2*#"))fprintf ( stdout, "\nM2 average Proportion of Ca with RMSD <%.2f Angs. Within a %.2f Angs. Ngb.", PP->rmsd_threshold, PP->maximum_distance );
		  if (strpbrk(PP->comparison_io,"3*#"))fprintf ( stdout, "\nM3 Proportion of Residues With M2>%.2f out of %d Res.\n", PP->similarity_threshold, (int)tot_n);
	      }
	  fprintf ( stdout, "\n");
	  ST->score_aln=ST->score=(tot_m3*100)/tot_n;
	  return ST;
      }

float *** analyse_pdb_residues ( Alignment *A, Constraint_list *CL, Pdb_param *pdb_param)
     {

	 int **pos;
	 int s1, s2, rs1, rs2;
	 int col1, col2;
	 float ***distances;
	 float d1, d2;
	 int real_res1_col1;
	 int real_res1_col2;
	 int real_res2_col1;
	 int real_res2_col2;
	 Pdb_param *PP;


	 PP=pdb_param;
	 
	 distances=vcalloc ( A->nseq, sizeof (float**));
	 for ( s1=0; s1< A->nseq; s1++)
	     {
		 distances[s1]=declare_float ( A->len_aln, 4); 
		 if (!PP->distance_on_request && CL->T[s1])
		     {
		     rs1=A->order[s1][0];
		     if ( !(CL->T[rs1])->ca_dist)(CL->T[rs1])->ca_dist=measure_ca_distances(CL->T[rs1]);
		     }
	     }
	 pos=aln2pos_simple (A, A->nseq);
	
	 for ( s1=0; s1< A->nseq; s1++)
	     {
		 rs1=A->order[s1][0];
		 if ( !(CL->T[rs1]))continue;
		 for ( col1=0; col1< A->len_aln; col1++)
		 {
		     if ( pos[s1][col1]<=0)continue;
		     for ( s2=0; s2<A->nseq; s2++)
		     {
			 rs2=A->order[s2][0];
			 if ( s2==s1)continue;
			 if ( !(CL->T[rs2]))continue;
			 if (pos[s2][col1]<=0)continue;
			 for ( col2=0; col2<A->len_aln; col2++)
			 {
			     
			     if ( FABS((col2-col1))<=PP->n_excluded_nb)continue;
			     if ( pos[s1][col2]<=0) continue;
			     if ( pos[s2][col2]<=0)continue;
			     
			    
			     real_res1_col1=pos[s1][col1]-1;
			     real_res1_col2=pos[s1][col2]-1;
			     
			     real_res2_col1=pos[s2][col1]-1;
			     real_res2_col2=pos[s2][col2]-1;
				
			     if ( !PP->distance_on_request)
			        {
				    d1=(CL->T[rs1])->ca_dist[real_res1_col1][real_res1_col2];
				    d2=(CL->T[rs2])->ca_dist[real_res2_col1][real_res2_col2];
				}
			     else 
			        {
				    d1=get_atomic_distance( (CL->T[rs1])->ca[real_res1_col1], (CL->T[rs1])->ca[real_res1_col2]);
				    d2=get_atomic_distance( (CL->T[rs2])->ca[real_res2_col1], (CL->T[rs2])->ca[real_res2_col2]);
				}

			     if (d1==UNDEFINED || d2==UNDEFINED)continue;
			     if ( d1<PP->maximum_distance && d2<PP->maximum_distance)
				 {
				     if ( FABS((d1-d2))<PP->rmsd_threshold)
					 distances[s1][real_res1_col1][2]++;

				     distances[s1][real_res1_col1][1]+=FABS((d1-d2));
				     distances[s1][real_res1_col1][0]++;
				 }
			 }
		     }
		 }
	     }	 
	 
     return distances;
     }

float square_atom ( Atom *X)
{
  
  return X->x*X->x + X->y*X->y + X->z*X->z;
}
Atom* reframe_atom ( Atom *X, Atom*Y, Atom *Z, Atom *IN, Atom *R)
     {
       float new_x, new_y, new_z;
       
       if ( R==NULL)R=vcalloc ( 1, sizeof (Atom));

       
        new_x= X->x*IN->x + Y->x*IN->y +Z->x*IN->z;
	new_y= X->y*IN->x + Y->y*IN->y +Z->y*IN->z;
	new_z= X->z*IN->x + Y->z*IN->y +Z->z*IN->z;

	R->x=new_x;
	R->y=new_y;
	R->z=new_z;
       return R;
     }

Atom* add_atom ( Atom *A, Atom*B, Atom *R)
{
  if ( R==NULL)R=vcalloc ( 1, sizeof (Atom));
  
  R->x=A->x+B->x;
  R->y=A->y+B->y;
  R->z=A->z+B->z;
  
  return R;
}
Atom* diff_atom ( Atom *A, Atom*B, Atom *R)
{
  if ( R==NULL)R=vcalloc ( 1, sizeof (Atom));
  
  R->x=A->x-B->x;
  R->y=A->y-B->y;
  R->z=A->z-B->z;
  
  return R;
}

Atom * copy_atom ( Atom *A, Atom*R)
{
  if ( R==NULL)R=vcalloc ( 1, sizeof (Atom));
  R->num=A->num;
  R->res_num=A->res_num;
  R->x=A->x;
  R->y=A->y;
  R->z=A->z;
  
  sprintf( R->type, "%s", A->type);
  return R;
}
 void print_atom (Atom *A)
{
  fprintf ( stdout, "%.2f %.2f %.2f", A->x, A->y, A->z);
}
