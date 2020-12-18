#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <signal.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"
#include "t_coffee.h"

static char *get_dali_defaults(char *buf);


#define is_a_seq_file(file) (!is_matrix(file) && !is_matrix(file+1) && !is_method (file) && !is_method (file+1) &&(check_file_exists(file) || check_file_exists(file+1))) 
main ( int argc, char *argv[])
	{
	
	int a, b, c;
	
	Sequence *S;
	Fname *F=NULL;
	FILE *fp_parameters;
	Constraint_list *CL;
	Alignment  *A=NULL, *EA=NULL;
	Alignment  **aln_list;
	Sequence *SEQ_LIST;
	FILE *OUT;	

	
	NT_node **T=NULL;
	int tot_node;
	char *pc;
/*Parameters*/
	int garbage;
	int quiet;
	char *parameters;
	char *t_coffee_defaults;
	int t_coffee_defaults_flag;
	char *dali_defaults;
	char **nargv;

	FILE *le=NULL;
	char *se_name;
	char *run_name;
	
	char *mem_mode;
	
	int do_extend;
	char *extend_mode;
	int  nseq_for_quadruplet;
	char **seq_name_for_quadruplet;
	
	int do_compact;
	int **do_list;
	int n_do;

	char *compact_mode;
	int do_clean;
	char *clean_mode;
	int do_self;
	int do_normalise;
	


	int n_list;
	char **list_file;
	char *out_lib;
	char *outseqweight;
	int n_seq_to_align;
	char **seq_to_align;
	int cosmetic_penalty,gop,f_gop, gep,f_gep, nomatch;
	
	char *tree_file;
	char *use_tree;
	char *tree_mode;

	int quicktree;
	char *out_aln;
	char **tot_out_aln;
	int maximise;
	char **out_aln_format;
	int  n_out_aln_format;
	char *infile;
	char *matrix;
	char *dp_mode;
	int tg_mode;
	int   ktup;
	int   fasta_step;
	int   diag_mode;
	char *sim_matrix;

	char *type;
	char *outorder;
	char *output_res_num;
	char *residue_case;
	int extra_cpu;

	char *weight;
	char *seq_weight;
	int do_align;
	char *evaluate_mode;
	
	int get_type;
	/*Post_processing*/
	int clean_aln;
	int clean_threshold;
	int clean_iteration;
	char *clean_evaluate_mode;
	/*Profile Alignment*/
	
	char **profile;
	int n_profile;
	
	char *profile1;
	char *profile2;

	/*Domain Parameters*/
	int do_domain;
	int domain_start;
	int domain_len;
	int domain_scale;
	int domain_interactive;
	
	int do_score;
	int do_convert;
	int maxnseq;
	int maxlen;
	
	

argv=standard_initialisation (argv, &argc);
/*PARAMETER PROTOTYPE:    READ PARAMETER FILE     */
	       declare_name (parameters);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-parameters"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "R_F"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "get bottom parameters" ,\
			    /*Parameter*/ &parameters          ,\
			    /*Def 1*/     "NULL"             ,\
			    /*Def 2*/     "stdin"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );	     
	       declare_name (dali_defaults);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-dali_defaults"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "get top parameters" ,\
			    /*Parameter*/ &dali_defaults          ,\
			    /*Def 1*/     "NULL"             ,\
			    /*Def 2*/     "HARD_CODED"       ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
	       declare_name (t_coffee_defaults);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-t_coffee_defaults"        ,\
			    /*Flag*/      &t_coffee_defaults_flag     ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "get top parameters" ,\
			    /*Parameter*/ &t_coffee_defaults          ,\
			    /*Def 1*/     "NULL"             ,\
			    /*Def 2*/     "NULL"       ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

/*Parameters > Prompt>t_coffee_defaults>dali_defaults*/ 

if (t_coffee_defaults && t_coffee_defaults[0] )
   {
     nargv=vcalloc (argc+1, sizeof (char*));
     nargv[0]=argv[0];
     nargv[1]=vcalloc ( VERY_LONG_STRING, sizeof(char));  
     a=0;
     fp_parameters=vfopen (t_coffee_defaults, "r");
     while ((c=fgetc (fp_parameters))!=EOF)nargv[1][a++]=c;
     nargv[1][a]='\0';
   
     for (a=1; a<argc; a++)
       {
	 nargv[a+1]=argv[a];	 
       }
     argc++;
     argv=nargv;
   } 
else if ( t_coffee_defaults_flag)
  {
    declare_name(t_coffee_defaults);
    sprintf ( t_coffee_defaults, "%s/.t_coffee_defaults",getenv ( "HOME") );
    if (  getenv ( "TCOFFEE_DEFAULTS")==NULL && !check_file_exists ( t_coffee_defaults));
    else if (  getenv ( "TCOFFEE_DEFAULTS"))
      {
	sprintf ( t_coffee_defaults,"%s",  getenv ( "TCOFFEE_DEFAULTS"));       
      }
    
    if (   check_file_exists ( t_coffee_defaults))
      {
	nargv=vcalloc (argc+1, sizeof (char*));
	nargv[0]=argv[0];
	nargv[1]=vcalloc ( VERY_LONG_STRING, sizeof(char));  
	a=0;
	fp_parameters=vfopen (t_coffee_defaults, "r");
	while ((c=fgetc (fp_parameters))!=EOF)nargv[1][a++]=c;
	nargv[1][a]='\0';
	
	for (a=1; a<argc; a++)
	  {
	    nargv[a+1]=argv[a];	 
	  }
	argc++;
	argv=nargv;
      }
  }
    
				
    
  

if ( dali_defaults && dali_defaults[0] )
   {
     nargv=vcalloc (argc+1, sizeof (char*));
     nargv[0]=argv[0];
     nargv[1]=vcalloc ( VERY_LONG_STRING, sizeof(char));  
     if ( !strm(dali_defaults, "HARD_CODED"))
       {
	 a=0;
	 fp_parameters=vfopen (dali_defaults, "r");
	 while ((c=fgetc (fp_parameters))!=EOF)nargv[1][a++]=c;
	 nargv[1][a]='\0';
	 vfclose (fp_parameters);
       }
     else
       {
	 nargv[1]=get_dali_defaults(nargv[1]);
       }
     for (a=1; a<argc; a++)
       {
	 nargv[a+1]=argv[a];
       }
     argc++;
     argv=nargv;
   }
	 
       
if ( parameters && parameters[0])
   {
   argv[argc]=vcalloc ( VERY_LONG_STRING, sizeof(char));  
   a=0;
   fp_parameters=vfopen (parameters, "r");
   while ((c=fgetc (fp_parameters))!=EOF)argv[argc][a++]=c;
   vfclose ( fp_parameters);
   argv[argc][a]='\0';
   fprintf ( stderr, "\n%s", argv[argc]);
   argc++;  
   }

argv=break_list ( argv, &argc, "=;, \n");
/*PARAMETER PROTOTYPE:    DO SCORE      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-score"        ,\
			    /*Flag*/      &do_score        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "forces the program to make the alignment" ,\
			    /*Parameter*/ &do_score          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );
/*PARAMETER PROTOTYPE:    DO FORMAT      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-convert"        ,\
			    /*Flag*/      &do_convert        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "forces the program to make the alignment" ,\
			    /*Parameter*/ &do_convert          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );	     


/*PARAMETER PROTOTYPE*/
     declare_name (se_name);
     get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-quiet"      ,\
			    /*Flag*/      &quiet        ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &se_name      ,\
			    /*Def 1*/     (do_convert || do_score)?"/dev/null":"stderr"      ,\
			    /*Def 2*/     "/dev/null"   ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
     le=vfopen ( se_name, "w");
     fprintf ( le, "\nPROGRAM: %s (%s)\n",PROGRAM,VERSION);

/*PARAMETER PROTOTYPE: RUN NAME*/
               declare_name (run_name);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-run_name"   ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &run_name     ,\
			    /*Def 1*/     "NULL"        ,\
			    /*Def 2*/     ""            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE: MEM MODE*/
	       declare_name(mem_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-mem_mode"   ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "mem OR disk" ,\
			    /*Parameter*/ &mem_mode     ,\
			    /*Def 1*/     "mem"         ,\
			    /*Def 2*/     ""            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE: EXTEND  */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-extend"     ,\
			    /*Flag*/      &do_extend    ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Do The extention On the Fly"          ,\
			    /*Parameter*/ &do_extend    ,\
			    /*Def 1*/     "1"           ,\
			    /*Def 2*/     "1"           ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE: EXTEND  */
	       declare_name (extend_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-extend_mode"     ,\
			    /*Flag*/      &garbage    ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "Chose the type of extention"          ,\
			    /*Parameter*/ &extend_mode    ,\
			    /*Def 1*/     "triplet"           ,\
			    /*Def 2*/     ""           ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE: SEQUENCES TO EXTEND */
	seq_name_for_quadruplet=declare_char ( 200, STRING);
	nseq_for_quadruplet=get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-seq_name_for_quadruplet"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  200           ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ seq_name_for_quadruplet    ,\
			    /*Def 1*/     "",\
			    /*Def 2*/     ""       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE: COMPACT */
	       declare_name (compact_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-compact"    ,\
			    /*Flag*/      &do_compact   ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &compact_mode ,\
			    /*Def 1*/     "default"      ,\
			    /*Def 2*/     "default"      ,\
			    /*Min_value*/ "0"           ,\
			    /*Max Value*/ "1"           \
		   );


/*PARAMETER PROTOTYPE:        CLEAN*/
	       declare_name ( clean_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-clean"      ,\
			    /*Flag*/      &do_clean     ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &clean_mode   ,\
			    /*Def 1*/     "no"          ,\
			    /*Def 2*/     "shadow"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );

/*PARAMETER PROTOTYPE:        DO SELF */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-do_self"    ,\
			    /*Flag*/      &do_self      ,\
			    /*TYPE*/      "FL"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  0             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &do_self      ,\
			    /*Def 1*/     "0"           ,\
			    /*Def 2*/     "1"           ,\
			    /*Min_value*/ "0"           ,\
			    /*Max Value*/ "1"           \
		   );

/*PARAMETER PROTOTYPE:        DO NORMALISE */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-do_normalise"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &do_normalise ,\
			    /*Def 1*/     "1000"           ,\
			    /*Def 2*/     "1000"       ,\
			    /*Min_value*/ "0"           ,\
			    /*Max Value*/ "10000"           \
		   );
/*PARAMETER PROTOTYPE:        IN */
	list_file=declare_char ( 200, STRING);
	n_list=get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-in"         ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  200           ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ list_file     ,\
			    /*Def 1*/    "Mlalign_id_pair Mfast_pair",\
			    /*Def 2*/     "stdin"       ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );


/*PARAMETER PROTOTYPE:    OUT_LIB     */
	declare_name (out_lib);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-out_lib"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/&out_lib       ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    OUT_LIB     */
	declare_name (outseqweight);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-outseqweight"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/&outseqweight       ,\
			    /*Def 1*/    "no"      ,\
			    /*Def 2*/    "default"      ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    SEQ TO ALIGN     */
	       seq_to_align=declare_char ( 200, STRING);
	       n_seq_to_align=get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-seq_to_align",\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"          ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  200              ,\
			    /*DOC*/       "Indicates the name of the sequences to analyse in Fasta Format",\
			    /*Parameter*/ seq_to_align   ,\
			    /*Def 1*/    "NULL"          ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    COSMETIC PENALTY     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-cosmetic_penalty" ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "A Penalty That only affects the non stable portions of the alignmnent\nSet to -10"          ,\
			    /*Parameter*/ &cosmetic_penalty          ,\
			    /*Def 1*/    "-50"            ,\
			    /*Def 2*/    "-50"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );	       
/*PARAMETER PROTOTYPE:    GAPOPEN     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-gapopen"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &gop          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    GAPEXT     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-gapext"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &gep          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
/*PARAMETER PROTOTYPE:    F_GAPOPEN     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-fgapopen"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &f_gop          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    F_GAPEXT     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-fgapext"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &f_gep          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
/*PARAMETER PROTOTYPE:    NEW_TREE     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-nomatch"     ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &nomatch          ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "0"             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );	       
	       declare_name ( tree_file);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-newtree"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/&tree_file     ,\
			    /*Def 1*/    "default"      ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
/*PARAMETER PROTOTYPE:    USETREE     */
	       declare_name ( use_tree);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-usetree",\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "R_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &use_tree     ,\
			    /*Def 1*/    "NULL"         ,\
			    /*Def 2*/    "NULL"         ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );

/*PARAMETER PROTOTYPE:    OUT_LIB     */
	       declare_name ( tree_mode);
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-tree_mode"  ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "S"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &tree_mode    ,\
			    /*Def 1*/    "fast"         ,\
			    /*Def 2*/    "slow"         ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
/*PARAMETER PROTOTYPE:    OUT_LIB     */
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
 			    /*output*/    &le           ,\
	 		    /*Name*/      "-quicktree"  ,\
		 	    /*Flag*/      &quicktree    ,\
			    /*TYPE*/      "FL"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  0             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &quicktree    ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "1"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );
	       if ( quicktree)sprintf ( tree_mode, "very_fast");
/*PARAMETER PROTOTYPE:    OUTFILE     */
	       declare_name ( out_aln);
	       tot_out_aln=declare_char (200, STRING);
	        get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-outfile"    ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "W_F"         ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &out_aln      ,\
			    /*Def 1*/    "default"      ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    MAXIMISE     */		
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-maximise"   ,\
			    /*Flag*/      &maximise     ,\
			    /*TYPE*/      "FL"          ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  0             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &maximise     ,\
			    /*Def 1*/    "1"            ,\
			    /*Def 2*/    "1"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"         \
		   );
	       
/*PARAMETER PROTOTYPE:    OUTPUT_FORMAT    */
	       out_aln_format=declare_char ( 200, STRING);
	       n_out_aln_format=get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-output"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  200              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ out_aln_format,\
			    /*Def 1*/    "clustalw"           ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (infile);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-infile"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &infile        ,\
			    /*Def 1*/    ""              ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	       
/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (matrix);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-matrix"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &matrix        ,\
			    /*Def 1*/    ""              ,\
			    /*Def 2*/    "blosum62mt"    ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	    
/*PARAMETER PROTOTYPE:    TG_MODE    */
	       
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-tg_mode"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "0: Penalise Term gap with gapopen and gapext\n1: gapopen only\n2: No penalty\n",\
			    /*Parameter*/ &tg_mode        ,\
			    /*Def 1*/    "1",\
			    /*Def 2*/    "0",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	   

/*PARAMETER PROTOTYPE:    DP_MODE    */
	       declare_name (dp_mode);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-dp_mode"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &dp_mode        ,\
			    /*Def 1*/    "cfasta_pair_wise",\
			    /*Def 2*/    "cfasta_pair_wise",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	       
/*PARAMETER PROTOTYPE:    KTUP    */
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-ktuple"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &ktup          ,\
			    /*Def 1*/    "1",\
			    /*Def 2*/    "1",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	       
/*PARAMETER PROTOTYPE:    FASTA_STEP    */
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-ndiag"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &fasta_step          ,\
			    /*Def 1*/    "0",\
			    /*Def 2*/    "10",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	
/*PARAMETER PROTOTYPE:    diag_mode    */
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-diag_mode"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "D"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "0: Use the whole Diag\n1: Use the best match\n"           ,\
			    /*Parameter*/ &diag_mode          ,\
			    /*Def 1*/    "0",\
			    /*Def 2*/    "1",
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	
/*PARAMETER PROTOTYPE:    SIM_MATRIX    */
	       declare_name (sim_matrix);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-sim_matrix"      ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &sim_matrix        ,\
			    /*Def 1*/    "vasiliky",\
			    /*Def 2*/    "idmat",\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	       

/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (type);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-type"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "protein or dna"           ,\
			    /*Parameter*/ &type          ,\
			    /*Def 1*/    ""              ,\
			    /*Def 2*/    ""              ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );	
/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (outorder);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-outorder"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "input or aligned"           ,\
			    /*Parameter*/ &outorder       ,\
			    /*Def 1*/    "aligned"          ,\
			    /*Def 2*/    "input"        ,\
			    /*Min_value*/ "any"          ,\
			    /*Max Value*/ "any"           \
		   );
/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (output_res_num);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-seqnos"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &output_res_num ,\
			    /*Def 1*/    "off"            ,\
			    /*Def 2*/    "on"             ,\
			    /*Min_value*/ "any"           ,\
			    /*Max Value*/ "any"            \
		   );
/*PARAMETER PROTOTYPE:    INFILE    */
	       declare_name (residue_case);
	       get_cl_param(\
			    /*argc*/      argc           ,\
			    /*argv*/      argv           ,\
			    /*output*/    &le            ,\
			    /*Name*/      "-case"        ,\
			    /*Flag*/      &garbage       ,\
			    /*TYPE*/      "S"            ,\
			    /*OPTIONAL?*/ OPTIONAL       ,\
			    /*MAX Nval*/  1              ,\
			    /*DOC*/       "ND"           ,\
			    /*Parameter*/ &residue_case         ,\
			    /*Def 1*/    "upper"            ,\
			    /*Def 2*/    "lower"             ,\
			    /*Min_value*/ "any"           ,\
			    /*Max Value*/ "any"            \
		   );

/*PARAMETER PROTOTYPE:    CPU     */
	       
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-cpu"        ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &extra_cpu    ,\
			    /*Def 1*/    "0"            ,\
			    /*Def 2*/    "0"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    MAXNSEQ     */
	       
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-maxnseq"        ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &maxnseq    ,\
			    /*Def 1*/    "-1"            ,\
			    /*Def 2*/    "-1"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );

/*PARAMETER PROTOTYPE:    MAXLEN     */
	       
	       get_cl_param(\
			    /*argc*/      argc          ,\
			    /*argv*/      argv          ,\
			    /*output*/    &le           ,\
			    /*Name*/      "-maxlen"        ,\
			    /*Flag*/      &garbage      ,\
			    /*TYPE*/      "D"           ,\
			    /*OPTIONAL?*/ OPTIONAL      ,\
			    /*MAX Nval*/  1             ,\
			    /*DOC*/       "ND"          ,\
			    /*Parameter*/ &maxlen    ,\
			    /*Def 1*/    "-1"            ,\
			    /*Def 2*/    "-1"            ,\
			    /*Min_value*/ "any"         ,\
			    /*Max Value*/ "any"          \
		   );


/*PARAMETER PROTOTYPE:    WEIGHT      */
	       declare_name ( weight);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-weight"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "sim OR  sim_matrix" ,\
			    /*Parameter*/ &weight          ,\
			    /*Def 1*/    "sim"             ,\
			    /*Def 2*/    "sim"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );	  /*PARAMETER PROTOTYPE:    WEIGHT      */
	       declare_name ( seq_weight);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-seq_weight"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "t_coffee" ,\
			    /*Parameter*/ &seq_weight          ,\
			    /*Def 1*/    "t_coffee"             ,\
			    /*Def 2*/    "no_seq_weight"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );	    

/*PARAMETER PROTOTYPE:    DO ALIGN      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-align"        ,\
			    /*Flag*/      &do_align        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  0                ,\
			    /*DOC*/       "forces the program to make the alignment" ,\
			    /*Parameter*/ &do_align          ,\
			    /*Def 1*/     "1"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );	   
/*PARAMETER PROTOTYPE:    DO DOMAIN      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-mocca"        ,\
			    /*Flag*/      &do_domain        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  0                ,\
			    /*DOC*/       "forces the program to extract domains" ,\
			    /*Parameter*/ &do_domain          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );  
/*PARAMETER PROTOTYPE:    Domain Param      */
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-start"        ,\
			    /*Flag*/      &domain_start        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "start of the domain" ,\
			    /*Parameter*/ &domain_start          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );  	     
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-len"        ,\
			    /*Flag*/      &domain_len        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "length of the domain" ,\
			    /*Parameter*/ &domain_len          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-scale"        ,\
			    /*Flag*/      &domain_scale        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "scale factor for match value" ,\
			    /*Parameter*/ &domain_scale          ,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-mocca_interactive"        ,\
			    /*Flag*/      &domain_interactive        ,\
			    /*TYPE*/      "FL"             ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  0                ,\
			    /*DOC*/       "make the domain in an interactive manneer" ,\
			    /*Parameter*/ &domain_interactive,\
			    /*Def 1*/     "0"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );

 /*PARAMETER PROTOTYPE:    WEIGHT      */
	       declare_name (evaluate_mode);
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-evaluate_mode"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Indictes the mode used to produce the color output" ,\
			    /*Parameter*/ &evaluate_mode          ,\
			    /*Def 1*/    "t_coffee_fast"             ,\
			    /*Def 2*/    "dali"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );	     
	       get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-get_type"        ,\
			    /*Flag*/      &get_type        ,\
			    /*TYPE*/      "FL"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "forces t_coffee top get the type of the sequences" ,\
			    /*Parameter*/ &get_type          ,\
			    /*Def 1*/    "0"             ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );	     

get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-clean_aln"        ,\
			    /*Flag*/      &clean_aln        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Forces weak portion of aln to be realigned" ,\
			    /*Parameter*/ &clean_aln          ,\
			    /*Def 1*/    "1"             ,\
			    /*Def 2*/    "1"             ,\
			    /*Min_value*/ "0"            ,\
			    /*Max Value*/ "1"             \
		   );
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-clean_threshold"        ,\
			    /*Flag*/      &clean_threshold        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Threshold for aln cleaning" ,\
			    /*Parameter*/ &clean_threshold          ,\
			    /*Def 1*/     "1"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-clean_iteration"        ,\
			    /*Flag*/      &clean_iteration        ,\
			    /*TYPE*/      "D"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Number of cleaning rounds" ,\
			    /*Parameter*/ &clean_iteration          ,\
			    /*Def 1*/     "1"             ,\
			    /*Def 2*/     "1"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
declare_name (clean_evaluate_mode);
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-clean_evaluate_mode"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Mode used to score residues (see evaluate_mode)" ,\
			    /*Parameter*/ &clean_evaluate_mode          ,\
			    /*Def 1*/    "t_coffee_fast"             ,\
			    /*Def 2*/    "t_coffee_fast"             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );	    
profile=declare_char ( 200, STRING);
n_profile=get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-profile"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  200                ,\
			    /*DOC*/       "Input a profile" ,\
			    /*Parameter*/ profile          ,\
			    /*Def 1*/    ""             ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );	    
declare_name (profile1);
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-profile1"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Input a profile" ,\
			    /*Parameter*/ &profile1          ,\
			    /*Def 1*/    ""             ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
declare_name (profile2);
get_cl_param(\
			    /*argc*/      argc             ,\
			    /*argv*/      argv             ,\
			    /*output*/    &le              ,\
			    /*Name*/      "-profile2"        ,\
			    /*Flag*/      &garbage        ,\
			    /*TYPE*/      "S"              ,\
			    /*OPTIONAL?*/ OPTIONAL         ,\
			    /*MAX Nval*/  1                ,\
			    /*DOC*/       "Input a profile" ,\
			    /*Parameter*/ &profile2          ,\
			    /*Def 1*/    ""             ,\
			    /*Def 2*/    ""             ,\
			    /*Min_value*/ "any"            ,\
			    /*Max Value*/ "any"             \
		   );
 get_cl_param( argc, argv,&le, NULL,NULL,NULL,0,0,NULL);	     
	       
	     
/*PRE-PROCESS THE DATA FILE*/
	       if (  infile && infile[0] && !do_score)
		   {
		   
		   sprintf ( list_file[n_list++], "%s",infile);		   
		   }
	       else  if (  infile && infile[0] && do_score)
	           {
		   
		       if ( infile[0]=='A' ||infile[0]=='S')
			   {
			   infile[0]='S';
			   sprintf ( list_file[n_list++], "%s",infile);
			   }
		       else sprintf ( list_file[n_list++], "S%s",infile);
		       
		   }
	       
	       if ( profile1 && profile1[0])
		 {
		   sprintf ( list_file[n_list++], "R%s",profile1);
		 }
	       if ( profile2 && profile2[0])
		 {
		   sprintf ( list_file[n_list++], "R%s",profile2);
		 }
	       for ( a=0; a< n_profile; a++)
		 {
		   sprintf ( list_file[n_list++], "R%s",profile[a]);
		 }
	       

	       a=1;
	       while ( a< argc && argv[a][0]!='-')
		   {
		   if ( !do_score)sprintf ( list_file[n_list++], "%s", argv[a]); 
		   else sprintf ( list_file[n_list++], "S%s", (argv[a][0]=='A')?argv[a]+1:argv[a]);
  
		   if (!is_method (argv[a]))sprintf (infile, "%s", (is_in_set ( argv[a][0], "ALSXPR"))?argv[a]+1:argv[a]);
		   a++;
		   }	       
	       
	       if (n_list==0 || argc<=1)exit(1);
	       else if ( run_name)F=parse_fname(run_name);
	       else if ( argv[1][0]!='-')
		 {
		   if (check_file_exists(argv[1]))F=parse_fname(argv[1]);
		   else if ( check_file_exists(argv[1]+1))F=parse_fname(argv[1]+1);
		 }
	       else 
	          {		      
		  for ( a=0; a< n_list; a++)
		      {
			  
		
			  if (!is_method(list_file[a])) 
			     {
			     if ( check_file_exists( list_file[a])){F=parse_fname(list_file[a]);break;}			 
			     else if ( is_in_set ( list_file[a][0], "ASLXPR") && check_file_exists( list_file[a]+1)){F=parse_fname(list_file[a]+1);break;}
			     }
		      }
		
		  }
	      

	       /*FATAL: NO SEQUENCES*/
	       if (!F)
		 {
		   fprintf ( stderr, "\nFATAL: No Sequence 1 [%s]", argv[0]);
		   exit (1);
		 }

	       if (!run_name)F->path[0]='\0';
	       if (  matrix && matrix[0]){sprintf ( list_file[n_list++], "X%s",matrix);}
	       
	       
	       identify_list_format      (list_file, n_list);
	      
	       
	       fprintf (le, "\nINPUT FILES\n");
	       for ( a=0; a< n_list; a++)
		   {
		     fprintf (le, "\tType=%c", list_file[a][0]);
		     fprintf (le, " File=%s", list_file[a]+1); 
		     if ( list_file[a][0]=='A' || list_file[a][0]=='S' || list_file[a][0]=='P'|| list_file[a][0]=='R' ) fprintf (le, " Format=%s\n", identify_seq_format ( list_file[a]+1));
		     else fprintf (le, "\n");
		   }
	       

/*SET GAP PENALTIES POSITIVE OR NEGATIVE*/
	       if (  maximise && (gop>0 || gep>0)){gop=(gop>0)?-gop:gop; gep=(gep>0)?-gep:gep;}
	       if ( !maximise && (gop<0 || gep<0)){gop=(gop<0)?-gop:gop; gep=(gep<0)?-gep:gep;}
	       if (  maximise && (f_gop>0 || f_gep>0)){f_gop=(f_gop>0)?-f_gop:f_gop; f_gep=(f_gep>0)?-f_gep:f_gep;}
	       if ( !maximise && (f_gop<0 || f_gep<0)){f_gop=(f_gop<0)?-f_gop:f_gop; f_gep=(f_gep<0)?-f_gep:f_gep;}
	       



		            
/*CONVERT, ALIGN OR SCORE: CHOSE THE RIGHT VERB*/
	       /*Set the Hierarchy of the verbs*/
	      
	       
	       do_list=vcalloc ( 100, sizeof (int*));
	       n_do=0;
	       do_list[n_do++]=&do_convert;
	       do_list[n_do++]=&do_score;
	       do_list[n_do++]=&do_domain;
	       do_list[n_do++]=&do_align;

	       for ( a=0; a< n_do; a++)
		 {
		 if ( do_list[a][0])
		   {
		   for ( b=0; b< n_do; b++)if ( b!=a)do_list[b][0]=0;
		   break;
		   }
		 }

	       
	     
/*SET THE DEFAULT NAMES*/
	       if ( do_convert && infile && infile[0])
	           {
		       if (strm (out_lib, "default"))
		          {
			  
			      n_list=0;
			      sprintf ( list_file[n_list++], "S%s", infile); 
			      sprintf ( out_lib, "no");
			  }		       
		       else 
		          {
			  sprintf ( list_file[n_list++], "S%s", infile); 
			  }
		       
		       

		       if ( strm (tree_file, "default"))sprintf ( tree_file, "no");		   
		   }

	       if (  do_score)
		   {
		   sprintf ( out_lib, "no");
		   sprintf ( tree_file, "no");
		   }
 
     
	       if ( F && strm ( tree_file, "default"))sprintf ( tree_file ,"%s%s.dnd",F->path     ,F->name);

	       if ( F && strm ( out_aln  , "default"))
	          {
		  for (a=0; a< n_out_aln_format; a++)
		      sprintf ( tot_out_aln[a]   ,"%s%s.%s"      ,F->path,F->name, (strm (out_aln_format[a], "gcg")?"msf":(strm(out_aln_format[a],"clustalw")?"aln":out_aln_format[a])));
		  }
	       else if (  n_out_aln_format==1)
		   sprintf ( tot_out_aln[0], "%s", out_aln);
	       else
	          {
		  for (a=0; a< n_out_aln_format; a++)
		      sprintf ( tot_out_aln[a]   ,"%s%s.%s", F->path  ,out_aln, (strm (out_aln_format[a], "gcg")?"msf":(strm(out_aln_format[a],"clustalw")?"aln":out_aln_format[a])));    
		  }
	       
	       if ( F && strm ( out_lib  , "default"))sprintf ( out_lib   ,"%s%s.tc_lib",F->path     , F->name);
	       if ( type && type[0])
	          {
		      if (strm2 (type,"Protein", "protein"))sprintf ( type, "PROTEIN");
		      if (strm2 (type,"DNA", "dna"))sprintf ( type, "DNA");
		      
		  }
	       
	       if (   !use_tree && check_file_exists (tree_file))remove (tree_file);
	       else if ( !use_tree || (use_tree && strm (use_tree, "default")));
	       else sprintf ( tree_file, "%s", use_tree);
	       


/*******************************************************************************************************/	     
/*                                                                                                     */	      
/*                           Input Sequences and Library                                               */
/*                                                                                                     */	      
/*******************************************************************************************************/
	       
/*START*/
	       /*1 READ THE SEQUENCES*/
	      
	       S=read_seq_in_n_list (list_file, n_list, type);
	       le=display_sequences_names ( S, le);
	       if ( maxnseq!=-1 && S->nseq> maxnseq)
	          {
		      fprintf ( stderr, "\nTOO MANY SEQUENCES [NSEQ=%d][MAX=%d][FATAL]\n", S->nseq,maxnseq);
		      crash ("");
		  }
	       
	       if (maxlen!=-1 && S->max_len>maxlen)
		  {
		      fprintf ( stderr, "\nSEQUENCES TOO LONG [Longuest=%d][MAX=%d][FATAL]\n", S->max_len,maxlen);
		      crash ("");
		  }

	       if ( get_type)
		 {
		   S=get_sequence_type (S);
		   fprintf ( stdout , "%s\n", S->type);
		   free_sequence(S, S->nseq);
		   return 1;
		 }

	       free_sequence (S, S->nseq);
	       
	       /*2 PREPARE THE CONSTRAINT LIST*/
	       
	       CL=read_n_constraint_list (list_file,n_list,NULL, mem_mode,weight,type, le);   
	       if ( CL->M)clean_aln=0;
	       

	       /*If the List is empty*/
	       if ( CL->ne==0 && !CL->M &&!(do_convert && infile[0]))
		 {
		   fprintf ( stderr, "\n******************ERROR*****************************************\n");
		   
		   fprintf ( stderr, "\nYou have not provided any method or enough Sequences[FATAL]");
		   fprintf ( stderr, "\nIf you have used the '-in' Flag, ADD the methods you wish to use:");
		   fprintf ( stderr, "\n\t-in <your sequences> Mlalign_id_pair Mfast_pair\n");
		   fprintf ( stderr, "\nAnd make sure you provide at least TWO sequences\n");		   
		   fprintf ( stderr, "\n*****************************************************************\n");
		   exit(EXIT_FAILURE);
		 }

	       CL->normalise=do_normalise;
 
	       if ( type && type[0])sprintf ( (CL->S)->type, "%s", type); 
	       CL->extend_jit=(do_extend>0)?1:0;
	       
	       CL->extend_threshold=(do_extend==1)?0:do_extend;
	       CL->do_self=do_self;
	       sprintf (CL->extend_clean_mode,   "%s", clean_mode);
	       sprintf (CL->extend_compact_mode, "%s", compact_mode);
	       if ( CL->extend_jit && CL->extend_threshold !=0)filter_list (CL,0, CL->ne, CL->extend_threshold); 	       
	       CL->pw_parameters_set=1;



	       CL->nomatch=nomatch;

	       /*Gep and Gop*/
	       if ( !gep &&  !gop && CL->M)
		   {
		    CL->gop=get_avg_matrix_mm ( CL->M, (strm3((CL->S)->type,"PROTEIN", "Protein", "protein")?AA_ALPHABET:"gcuta"))*10;
		    CL->gep=CL->gop/10;
		    fprintf ( CL->local_stderr, "\nAUTOMATIC PENALTIES: gapopen=%d gapext=%d", CL->gop, CL->gep);
		   }
	       else if ( !CL->M && cosmetic_penalty && !gep && !gop)
	           {
		    CL->gep=0;
		    CL->gop=cosmetic_penalty;
		   }
	       else
	           {
		     CL->gep=gep;
		     CL->gop=gop;
		   }
	       
	       /*Frame Penalties*/
	       CL->f_gep=f_gep;
	       CL->f_gop=f_gop;

	       
	       CL->maximise=maximise;
	       if (strm((CL->S)->type,"DNA") )
		   CL->ktup=MAX(2,ktup);
	       else
		   CL->ktup=ktup;
	       
	       CL->use_fragments=diag_mode;
	       CL->fasta_step=fasta_step;

	       sprintf ( CL->matrix_for_aa_group, "%s", sim_matrix);
	       sprintf ( CL->dp_mode, "%s", dp_mode);
	       CL->TG_MODE=tg_mode;
	       	       
	       if (n_seq_to_align) 
	          {
		  if (n_seq_to_align==1)
		    {
		      SEQ_LIST=read_sequences (seq_to_align[0]);
		    }
		  else
		    {
		     SEQ_LIST=declare_sequence ( 1, 1,n_seq_to_align);
		     for ( a=0; a< n_seq_to_align; a++)sprintf ( SEQ_LIST->name[a], "%s", seq_to_align[a]);
		    }
		  CL->DO_S=duplicate_sequence(CL->S);
		  CL->DO_S=reorder_seq(CL->DO_S, SEQ_LIST->name,SEQ_LIST->nseq);
		  }
	      else
	          {
		  CL->DO_S=CL->S;
		  }
	      fprintf (le, "\n\n\tLibrary Total Size: [%d]\n", CL->ne); 
	      
	      if (out_lib[0]!='\0' && !strm (out_lib, "no"))
	         {
		    fprintf (CL->local_stderr, "\tLibrary    file [%s] with %d pairs\n",out_lib, CL->ne);
		    OUT=save_constraint_list ( CL, 0, CL->ne, out_lib, NULL, "ascii",CL->DO_S);
		    vfclose (OUT);
		 }
	      
	      
	      /*WEIGHT CONSTRAINT LIST*/

	      if ( !do_convert)
		{
		  CL=weight_constraint_list(CL, seq_weight);
		  if (output_seq_weights (CL->W, outseqweight))fprintf (CL->local_stderr, "\tSeq Weight     file [%s]",outseqweight);
		  le=display_weights(CL->W, le);
		}
	      /*Prepare quadruplets*/
	      CL->nseq_for_quadruplet=nseq_for_quadruplet;
	      for (a=0; a< CL->nseq_for_quadruplet; a++)
		{
		  if ( (b=name_is_in_list (seq_name_for_quadruplet[a],(CL->S)->name,(CL->S)->nseq, 100))!=-1)CL->seq_for_quadruplet[b]=1;
		  else fprintf ( stderr, "\nWarning: Sequence %s is not in the set and cannot be used for quqdruplet extension\n",seq_name_for_quadruplet[a]); 
		}
	      
		    
	      /*CONSTRAINT LIST READY*/
	      /*Chose the right Mode for comparing residues*/
	      if ( strm ( extend_mode, "triplet") && !CL->M)
		{
		  CL->evaluate_residue_pair=residue_pair_extended_list;
		  CL->get_dp_cost      =get_dp_cost;
		}
	    else if ( strm ( extend_mode, "slow_triplet") && !CL->M)
		{
		  CL->evaluate_residue_pair=residue_pair_extended_list;
		  CL->get_dp_cost      =slow_get_dp_cost;
		}
	      else if (  strm ( extend_mode, "mixt") && !CL->M)
		{
		  CL->evaluate_residue_pair=residue_pair_extended_list_mixt;
		  CL->get_dp_cost=slow_get_dp_cost;
		}
	      else if (  strm ( extend_mode, "quadruplet") && !CL->M)
		{
		  CL->evaluate_residue_pair=residue_pair_extended_list_quadruplet;
		  CL->get_dp_cost      =get_dp_cost_quadruplet;
		}
	      else if (  strm ( extend_mode, "matrix") || CL->M)
		{
		  CL->evaluate_residue_pair=evaluate_matrix_score;
		  CL->get_dp_cost=slow_get_dp_cost;
		  CL->normalise=1;
		}
	      else 
		{
		  fprintf ( stderr, "\nERROR: %s is an unknown extend_mode[FATAL]\n", extend_mode);
		  crash ("");
		}
	      

/*******************************************************************************************************/	     
/*                                                                                                     */	      
/*                           Prepare The Alignment                                                     */
/*                                                                                                     */	      
/*******************************************************************************************************/

	      if ( do_align)
		   {
		   A=seq2aln  ((CL->S),NULL,1);
		   ungap_array(A->seq_al,A->nseq);		   
		   
		   /*Chose the right Mode for evaluating Columns*/
	

		   
		   if ( (CL->DO_S)->nseq>2)
		      {	      
		      pc=tree_file;
		      if ( strm (tree_file, "default") || !check_file_exists (tree_file))
			  T=make_tree ( A,CL,gop, gep,(CL->DO_S),pc, tree_mode, maximise);	      
		      else if ( strm (tree_file, "no"))
			  T=make_tree ( A,CL,gop, gep,(CL->DO_S),NULL, tree_mode, maximise);
		      else
			T=read_tree (pc,&tot_node,(CL->DO_S)->nseq,  (CL->DO_S)->name);
		      
		      tree_aln ((T[3][0])->left,(T[3][0])->right,A,(CL->DO_S)->nseq, CL);
		      }
		   else
		      {
		      tree_aln (NULL,NULL,A, (CL->DO_S)->nseq,CL);
		      }
		   A->nseq=(CL->S)->nseq;
		   }
	       else if ( (do_score || do_convert))
	           {
		   
		   A=(infile && infile[0])?main_read_aln ( infile, declare_aln(CL->S)):NULL;
		   
		   if (A)
		     {
		       if ( CL->DO_S==CL->S){CL->DO_S=aln2seq (A);}
		       A->S=CL->S;		  
		       A->nseq=(CL->DO_S)->nseq;
		     }
		   }
	      else if (do_domain)
		   {
		     CL->moca=vcalloc ( 1, sizeof ( Moca));
		     if (strm ( "cfasta_pair_wise", dp_mode))sprintf (CL->dp_mode, "%s","domain_pair_wise"); 
		     (CL->moca)->moca_start=domain_start;
		     (CL->moca)->moca_len  =domain_len;
		     (CL->moca)->moca_scale=(domain_scale==0)?-(CL->normalise/20):domain_scale;
		     (CL->moca)->moca_interactive=domain_interactive;
		     


		     if (!cosmetic_penalty && !gep && !gop)
		       {
			 CL->gop=-200;
			 CL->gep=-100;
		       }

		     CL=prepare_cl_for_moca (CL);
		     aln_list=moca_aln (CL);
		     free_int ( CL->packed_seq_lu, -1);
		     CL->packed_seq_lu=NULL;

		     a=0;
		     while ( aln_list[a])
		       {			 
		       for ( b=0; b< n_out_aln_format; b++)
			  {			    
			    output_format_aln (out_aln_format[b],aln_list[a],EA=fast_coffee_evaluate_output(aln_list[a], CL), tot_out_aln[b]);
			    fprintf (le, "\tOutput  file [%15s] in [%10s] format\n", tot_out_aln[b],out_aln_format[b]);
			  }
		       a++;
		       }
		   exit (EXIT_SUCCESS);
		   }
	    
/*******************************************************************************************************/	     
/*                                                                                                     */	      
/*                           PREPARE THE ALIGNMENT FOR OUTPUT                                          */
/*                                                                                                     */	      
/*******************************************************************************************************/
	      if (A)
	          {
		      for ( a=0; a< A->nseq; a++)
		         {
			   for ( b=0; b< A->len_aln ; b++)
			     if ( A->seq_al[a][b]=='O' || A->seq_al[a][b]=='o')A->seq_al[a][b]='-';
			 }
		    
		     
		      A=reorder_aln ( A, (CL->DO_S)->name,(CL->DO_S)->nseq);      
		      
		      if ( strm(outorder, "aligned") && T)A=reorder_aln ( A,A->tree_order,A->nseq);
		      else A=reorder_aln ( A, (CL->DO_S)->name,(CL->DO_S)->nseq);  
		      
		     
		      A->output_res_num=strm3 ( output_res_num, "on", "On", "ON");
		      A->residue_case=strm3 (residue_case, "lower", "Lower", "LOWER");
		      
		     
		      if ( clean_aln)
			{
			  EA=main_coffee_evaluate_output(A, CL,clean_evaluate_mode);			  
			  A=clean_maln(A, EA,clean_threshold,clean_iteration);
			  free_aln (EA);
			  A=ungap_aln(A);
			}
		      EA=main_coffee_evaluate_output(A, CL, evaluate_mode);
		     
		      for ( a=0; a< n_out_aln_format; a++)
			  output_format_aln (out_aln_format[a],A,EA, tot_out_aln[a]);
		      
		      if (!strm2(out_aln, "stdout", "stderr") && le==stderr && !do_convert)output_format_aln ("t_coffee_aln",A,NULL,"stdout");
		     
		      fprintf (le, "\n\nOUTPUT RESULTS\n\tDendrogram file [%s]\n", tree_file);
		      for ( a=0; a< n_out_aln_format; a++)
			fprintf ( le, "\tOutput  file [%15s] in [%10s] format\n", tot_out_aln[a],out_aln_format[a]);
		  }
	
	      
		      
	fprintf (le, "\n\n"); 
		      
		  
	free_char (list_file, -1); 	
	free_Alignment (A); 		
	free_Alignment (EA); 		
	S=free_constraint_list (CL);
	free_sequence (S, S->nseq);
	remove ( "core");
	exit (EXIT_SUCCESS);
	}

/*Specialized set of Parameters*/
char *get_dali_defaults(char *buf)
   {
     
     buf=strcat (buf,"-cosmetic_penalty=-50 ");
     buf=strcat (buf,"-tree_mode=slow ");
     buf=strcat (buf,"-output clustalw,score_ascii ");
     buf=strcat (buf,"-evaluate_mode=non_extended_t_coffee ");
     buf=strcat (buf,"-clean_aln 0 ");
     return buf;
   }
