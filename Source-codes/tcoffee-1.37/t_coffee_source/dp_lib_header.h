
int gotoh_pair_wise_lalign ( Alignment *A, int*ns, int **l_s,Constraint_list *CL);
Constraint_list * undefine_sw_aln ( Alignment *A, Constraint_list *CL);
Constraint_list * undefine_sw_pair ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int sw_pair_is_defined ( Constraint_list *CL, int s1, int r1, int s2, int r2);


int gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

Alignment * get_best_local_aln ( Alignment *IN,Constraint_list *CL,int gop, int gep, int sw_t, int sw_l, int sw_z, int greedy);
Alignment * get_best_nol_local_aln ( Alignment *IN, Constraint_list *CL, int gop, int gep,int sw_t,int sw_l, int sw_z, int mode);
double compute_penalty   (Constraint_list *CL, char *mode, int len);
double compute_scale ( Constraint_list *CL,char *mode, int len);
int evaluate_penalty (Alignment *A, Constraint_list *CL, int *scale,char *scale_mode, int *penalty, char *penalty_mode, int len_seq);
Alignment ** t_coffee_lalign   (Constraint_list *CL, int scale, int penalty,int maximise,Sequence *S, int sw_t, int sw_l, int sw_z,int *sw_n, int sw_io);
Alignment * add_seq2aln   (Constraint_list *CL, Alignment *IN,Sequence  *S);







struct Dp_Model
{
  int *diag;

  int TG_MODE;
  int F_TG_MODE;
  int gop;
  int gep;
  int f_gop;
  int f_gep;
  int nstate;
  int START;
  int END;
  
  char**model_comments;
  int **model;
  int **model_properties;
  int **bounded_model;
  int (***model_emission_function)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *);

  int LEN_I;
  int LEN_J;
  int DELTA_I;
  int DELTA_J;
  int EMISSION;
  int START_EMISSION;
  int TERM_EMISSION;
  
  int ALN_TYPE;
  Constraint_list *CL;
  /*Associated Functions*/
  
  /*To Deprecate*/
  int UM;
  
  int TYPE;
  int F0;
  int F1;
  
  
  int NON_CODING;
  int INSERTION;
  int DELETION;
  int CODING0;
  int CODING1;
  int CODING2;
  

};
typedef struct Dp_Model Dp_Model;

struct Dp_Result
{
  int *traceback;
  int len;
  int score;
  Dp_Model *Dp_model;
};
typedef struct Dp_Result Dp_Result;

Dp_Result * make_fast_generic_dp_pair_wise (Alignment *A, int*ns, int **l_s,Dp_Model *M);

Constraint_list* free_dp_model  (Dp_Model *D);
Dp_Result * free_dp_result (Dp_Result *D );
int ** evaluate_diagonals_for_two_sequences ( char *seq1, char *seq2,int maximise,Constraint_list *CL,int ktup);
int ** evaluate_diagonals ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);
int ** evaluate_segments_with_ktup  ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);
int ** evaluate_diagonals_with_ktup ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);

int ** evaluate_diagonals_with_clist ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);

int * flag_diagonals (int l1, int l2, int **sorted_diag,float T, int window);                
int * extract_N_diag (int l1, int l2, int **sorted_diag, int n_chosen_diag, int window);

int hasch_seq(char *seq1, int **hs, int **lu,int ktup, char *alph);
int fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int cfasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int very_fast_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

int make_fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag);
int cfasta_gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int fasta_gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int make_fasta_gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag);

/*pair wise aln implementations*/
int gotoh_pair_wise    (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

  
      
  
 

int cfasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int fasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
Dp_Model* initialize_dna_dp_model (Constraint_list *CL);
Dp_Result * make_fast_dp_pair_wise (Alignment *A,int*ns, int **l_s, Constraint_list *CL,Dp_Model *M);
int make_fasta_cdna_pair_wise (Alignment *B,Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag);



int ** evaluate_diagonals_cdna ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup);
int cfasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
Alignment *clean_maln ( Alignment *A, Alignment *I, int T, int n_it);
Alignment *realign_segment (int seq, int start, int len,Alignment *A, Alignment *C);
Dp_Model * initialize_seg2prf_model(int left_tg_mode, int right_tg_mode, Constraint_list *CL);

int get_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

Dp_Model * initialize_sseq_model(int left_tg_mode, int right_tg_mode, Constraint_list *CL);
int ssec_pwaln_maln (Alignment *A, int *ns, int **ls, Constraint_list *CL);


int get_turn_gep_cost       (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_term_gep_cost  (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_alpha_gep_cost      (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_alpha_start_gep_cost(Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_alpha_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_beta_gep_cost       (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_term_gep_cost  (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_alpha_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_ssec_no_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
/*pair wise aln implementations*/
int myers_miller_pair_wise (Alignment *A, int *ns, int **l_s,Constraint_list *CL);
int diff (Alignment *A, int *ns, int **ls, int s1, int M,int s2, int N , int tb, int te, Constraint_list *CL, int **pos);
int evaluate_est_order (Sequence *S, char *concat, Constraint_list *CL, int ktuple);
 
Constraint_list *prepare_cl_for_moca ( Constraint_list *CL);
Alignment ** moca_aln ( Constraint_list *CL);
Alignment * extract_domain ( Constraint_list *CL);
Alignment * interactive_domain_extraction ( Constraint_list *CL);
int print_moca_interactive_choices ();

Alignment * approximate_domain ( int min_start, int max_start, int step_start,int min_len, int max_len, int step_len, int *best_start, int *best_len, int *best_score, Constraint_list *CL);
 
int measure_domain_length ( Constraint_list *CL,Alignment *IN, int start, int min_size, int max_size, int step);
Alignment *extract_domain_with_coordinates ( Alignment *RESULT,int start, int len, Constraint_list *CL);
int get_starting_point ( Constraint_list *CL);

Alignment * find_domain_coordinates (Constraint_list *CL, int *start, int *len);
Alignment * extend_domain ( Constraint_list *CL, int *start, int *len, int dstart, int dlen);
Alignment * modify_domain ( Constraint_list *CL, Alignment *IN, int *start, int *len, int dstart, int dlen);

int * analyse_sequence ( Constraint_list *CL);
/******************************************************************/
/*                   MAIN DRIVER                                  */
/*                                                                */
/*                                                                */
/******************************************************************/




Constraint_list *seq2list     ( Constraint_list *CL, char *command, char *seq, char *weight);

/******************************************************************/
/*                   MULTIPLE ALIGNMENTS                          */
/*                                                                */
/*                                                                */
/******************************************************************/
Alignment * compute_prrp_aln (Alignment *A, Constraint_list *CL);
Alignment * compute_clustalw_aln (Alignment *A, Constraint_list *CL);
Alignment * aln2clustalw_aln (Alignment *A, Constraint_list *CL);
/******************************************************************/
/*                  DNA                                           */
/*                                                                */
/*                                                                */
/******************************************************************/
Constraint_list * align_coding_nucleotides (char *seq, char *method, char *weight, char *mem_mode, Constraint_list *CL);
/******************************************************************/
/*                   STRUCTURES                                   */
/*                                                                */
/*                                                                */
/******************************************************************/


	
Constraint_list * align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * align_pdb_pair_2 (char *seq, Constraint_list *CL);

Constraint_list * custom1_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom2_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom3_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom4_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom5_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom6_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom7_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom8_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom9_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom10_align_pdb_pair (char *seq, Constraint_list *CL);
Constraint_list * custom_align_pdb_pair (char *seq, Constraint_list *CL, char *mode);



Constraint_list * sap_pair (char *seq, Constraint_list *CL);
 
/******************************************************************/
/*                   GENERIC PAIRWISE METHODS                     */
/*                                                                */
/*                                                                */
/******************************************************************/

Alignment * fast_pair      ( Alignment *A, char *seq, char *mode, Constraint_list *CL,Constraint_list *PW_CL);
Alignment * align_two_sequences ( char *seq1, char *seq2, char *matrix, int gop, int gep, char *align_mode);
NT_node ** make_tree ( Alignment *A,Constraint_list *CL,int gop, int gep,Sequence *S,  char *tree_file, char *tree_mode, int maximise);
int ** get_pw_distances ( Alignment *A,Constraint_list *CL,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode, int maximise);
Alignment *build_progressive_nol_aln_with_seq_coor(Constraint_list *CL,int gop, int gep,Sequence *S, int **seq_coor, int nseq);
Alignment *build_progressive_aln_with_seq_coor (Constraint_list*CL,int gop, int gep, Sequence *S, int **coor, int nseq);
Alignment *build_progressive_aln(Alignment *A, Constraint_list *CL, int gop, int gep);
Alignment *est_progressive_aln(Alignment *A, Constraint_list *CL, int gop, int gep);
void analyse_seq ( Alignment *A, int s);

void tree_aln ( NT_node LT, NT_node RT, Alignment*A, int nseq, Constraint_list *CL);
int pair_wise (Alignment *A, int*ns, int **l_s,Constraint_list *CL );
char *build_consensus ( char *seq1, char *seq2, char *dp_mode);

int domain_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL );
Alignment *domain_match_list2aln ( Alignment *A,int *ns,int **l_s,int **ml, int nseq, int len);
Alignment * domain_seq2domain (Constraint_list *CL,int scale,int gop,int gep,Alignment *SEQ_DOMAIN, Alignment *TARGET);


int custom_pair_score_function1  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function2  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function3  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function4  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function5  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function6  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function7  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function8  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function9  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function10 (Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_ca_trace_sap2_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_ca_trace_nb (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_sap1_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_3D_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_transversal (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble_2 (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble_3 (Constraint_list *CL, int s1, int s2, int r1, int r2);


/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:START                                   */
/*                                                                                           */
/*********************************************************************************************/
float matrix_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float transversal_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2,Struct_nb *nbs1, Struct_nb *nbs2);
float sap1_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float sap2_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);


/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:END                                     */
/*                                                                                           */
/*********************************************************************************************/
Alignment * analyse_pdb ( Alignment *A, Alignment *ST);
float *** analyse_pdb_residues ( Alignment *A, Constraint_list *CL, Pdb_param *pdb_param);

float square_atom ( Atom *X);
Atom* reframe_atom ( Atom *X, Atom*Y, Atom *Z, Atom *IN, Atom *R);
Atom* add_atom ( Atom *A, Atom*B, Atom *R);
Atom* diff_atom ( Atom *A, Atom*B, Atom *R);
Atom * copy_atom ( Atom *A, Atom*R);
int evaluate_aln_gibbs   ( Alignment *A, Constraint_list *CL);
int evaluate_moca_domain ( Alignment *A, Constraint_list *CL);
int moca_residue_pair_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int moca_evaluate_matrix_score      ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int moca_slow_get_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int **cache_cl_with_moca_domain (Alignment *A, Constraint_list *CL);
Alignment *make_moca_nol_aln ( Alignment *A, Constraint_list *CL);
/*********************************************************************************************/
/*                                                                                           */
/*         DOMAIN Z SCORE EVALUATION                                                         */
/*                                                                                           */
/*********************************************************************************************/

int evaluate_domain_aln_z_score (Alignment *A, int start, int end,Constraint_list *CL, char *alphabet);
int evaluate_domain_aln  ( Alignment *A, int start, int end,Constraint_list *CL);



int unpack_seq_residues ( int *s1, int *r1, int *s2, int *r2, int **packed_seq_lu);
Alignment * unpack_seq_aln ( Alignment *A,Constraint_list *C);
typedef struct 
{
  int N_COMPONENT;
  double *double_logB_alpha;
  double *exponant_list;
  double **ALPHA;
  double *DM_Q;
  double *alpha_tot;
  int n_aa;
  int tot_n;
}
Mixture;

double int_logB (int *i, int n);
double float_logB (float *i, int n);
double double_logB (double *x, int n);
double *** make_lup_table ( Mixture *D);
double  double_logB2(int j, double *n,Mixture *D);
double compute_exponant ( double *n, int j, Mixture *D);

double *compute_matrix_p ( double *n,int Nseq);
double* compute_dirichlet_p ( double *n,int Nseq);
void precompute_log_B ( double *table,Mixture *D);
double compute_X (double *n,int i,Mixture *D);
Mixture *read_dirichlet ( char *name);
int dirichlet_code( char aa);

/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR EVALUATING THE CONSISTENCY BETWEEN ALN AND CL                       */
/*                                                                                           */
/*********************************************************************************************/
Alignment * main_coffee_evaluate_output          ( Alignment *IN,Constraint_list *CL, const char *mode );
Alignment * fast_coffee_evaluate_output          ( Alignment *IN,Constraint_list *CL);
Alignment * slow_coffee_evaluate_output          ( Alignment *IN,Constraint_list *CL);
Alignment * non_extended_t_coffee_evaluate_output( Alignment *IN,Constraint_list *CL);
Alignment * heuristic_coffee_evaluate_output     ( Alignment *IN,Constraint_list *CL);
/*Old Function: To deprecate*/
Alignment * coffee_evaluate_output ( Alignment *IN,Constraint_list *CL);

/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR GETING THE COST : (Sequences) ->evaluate_residue_pair               */
/*                                                                                           */
/*********************************************************************************************/
int evaluate_cdna_matrix_score               ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_matrix_score                    ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_non_extended_list           ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list               ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list_quadruplet    ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list_mixt          ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int extend_residue_pair                      ( Constraint_list *CL, int s1, int r1, int s2, int r2);

/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR GETTING THE COST :  CL->get_dp_cost                                 */
/*                                                                                           */
/*********************************************************************************************/
int get_cdna_best_frame_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost                 ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost_quadruplet      ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
      

int very_fast_get_dp_cost       ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int fast_get_dp_cost            ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int fast_get_dp_cost_quadruplet ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int slow_get_dp_cost            ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int sw_get_dp_cost              ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_domain_dp_cost          ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2,Constraint_list *CL , int scale , int gop, int gep);

/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR ANALYSING AL OR MATRIX                                              */
/*                                                                                           */
/*********************************************************************************************/

int aln2n_res ( Alignment *A, int start, int end);
float get_gop_scaling_factor ( int **matrix,float id, int l1, int l2);
float get_avg_matrix_mm ( int **matrix, char *alphabet);
float get_avg_matrix_match ( int **matrix, char *alphabet);
float get_avg_matrix_diff ( int **matrix1,int **matrix2, char *alphabet);
float measure_matrix_enthropy (int **matrix,char *alphabet);
float measure_matrix_pos_avg (int **matrix,char *alphabet);
float evaluate_random_match (char *matrix, int n, int len,char *alp);
float compare_two_mat (char  *mat1,char*mat2, int n, int len,char *alp);
float compare_two_mat_array (int  **matrix1,int **matrix2, int n, int len,char *alp);
int ** rescale_two_mat (char  *mat1,char*mat2, int n, int len,char *alp);
int ** rescale_matrix ( int **mat, float lambda, char *alp);
void output_matrix_header ( char *name, int **matrix, char *alp);
float evaluate_random_match2 (int **matrix, int n, int len,char *alp);
float measure_lambda2(char  *mat1,char*mat2, int n, int len,char *alp);
float measure_lambda(char  *mat1,char*mat2, int n, int len,char *alp);
int ** show_pair(int istart, int iend, int jstart, int jend, int *seqlen_array, char **seq_array, int dna_ktup, int dna_window, int dna_wind_gap, int dna_signif,int prot_ktup, int prot_window,int prot_wind_gap,int prot_signif, int nseqs,int dnaflag, int max_aa, int max_aln_length  );
