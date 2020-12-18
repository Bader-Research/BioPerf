#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

/*********************************************************************/
/*                                                                   */
/*                         PRODUCE IN LIST                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

	 
Constraint_list *produce_list ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode)	
    {
    int b;
    int n_aln;



    static char **aln_command;
    static char **seq_command;
    static char *seq;
    static char *aln;
    static char *list;
    static char *main_list;
    static char *aln_mode;
    static char *out_mode;


    int file_exists;


    int n_elements_in, n_new_elements;
    Constraint_list *BCL;
    FILE *local_stderr;

    /*OUT_MODE:
              A->  alignment provided in a file
	      L->  list      provided in a file
	      aln-> alignment computed and provided in a file
	      list->list      computed and provided in a file
    */
    
    if ( !aln_mode)aln_mode=vcalloc ( STRING, sizeof (char));
    if ( !out_mode)out_mode=vcalloc ( STRING, sizeof (char));
    
    if ( CL==NULL)CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?vtmpfile():NULL, NULL);

    local_stderr=(CL->local_stderr!=NULL)?CL->local_stderr:stderr;
 
    CL->local_stderr=vfopen("/dev/null", "w");

    if ( aln_command==NULL)
        {
	    
	    seq=vtmpnam (seq);
	    aln=vtmpnam (aln);
	    list=vtmpnam (list);
	    main_list=vtmpnam (main_list);

	    aln_command=declare_char (1,-1);
	    seq_command=declare_char (1,-1);
	}
    
   

    n_aln=parse_method ( method,S,&aln_command,&seq_command, seq, aln, aln_mode, out_mode);
    
    n_elements_in=CL->ne;    
    /* Cf. parse method for possible out_mode flags*/

    for ( b=0; b< n_aln; b++)
        { 
	file_exists=1;
	
	if ( !strm4 (out_mode,"A","L","fA", "fL"))
	   {	       
	    local_stderr=output_completion ( local_stderr, b, n_aln,1);
	    seq_list2fasta_file( S,  seq_command[b], seq);
	    system ( aln_command[b]);
	    remove (seq);
	    file_exists=evaluate_sys_call_io ( aln, aln_command[b], "");	       
	   }

	
	if ( !file_exists)
	   {
	   fprintf ( stderr, "\nFAILED TO EXECUTE:\n\t%s\n", aln_command[b]); 
	   }
	else
	   {
	   
	   if (!strm2 (out_mode, "A", "L"))CL->cpu=-1;
	   else CL->cpu=1;
	   
	   if (strnm(out_mode, "aln",3) ||strm(out_mode, "A"))
	       CL=read_constraint_list (CL, aln,"aln","disk",weight);
	   else if (strm2(out_mode, "lib","L"))
	       CL=read_constraint_list (CL, aln,"lib","disk",weight);
	   else if ( strm (out_mode, "fL"))
	       {
	       local_stderr=output_completion ( local_stderr, b, n_aln,1);
	       
	       BCL= seq2list(CL,aln_command[b],seq_command[b], weight);
	       if ( !BCL) 
	           {
		   fprintf ( stderr, "\nFAILED TO EXECUTE:\n\t%s\n", aln_command[b]); 
	           }
	       else 
	           {
		       CL=BCL;
		   }
	       }
	   if ( !strm4 (out_mode, "A", "L", "fA", "fL"))remove(aln);
	   }	    
        }

    n_new_elements=CL->ne - n_elements_in;
    compact_list (CL, n_elements_in,n_new_elements, "best");
    compact_list (CL, 0, CL->ne, "default");
    CL->local_stderr=local_stderr;
   
    return CL;
    }


int parse_method ( char *fname, Sequence *S, char ***aln_c, char ***seq_c, char *seq, char *aln, char *aln_mode, char *out_mode)
    {
    int a, b, c, d, x, y;


    char *** param_list;
    int * n_param;
    int ** combination_table;
    int n_combinations=1;

    char buf [LONG_STRING];
    char executable [LONG_STRING];
    char aln_c_buf[VERY_LONG_STRING];
    char seq_c_buf[VERY_LONG_STRING];
    int tot_n_param;
    int preset_method=0;
    char *perl_pg_param[3];
    int n_perl_pg_param=3; 
    int do_mirror;

    /*A method can be:
      1- a pre computed alignment out_mode=A
      2- a precomputed Library    out_mode=L
      3- a method producing an alignment out_mode=aln
      4- a method producing an alignment out_mode=list
      5- a function producing an alignment out_mode=faln
      6- a function producing a library    out_mode=flist
    */

    if ( fname[0]=='A')
	{
	sprintf ( out_mode, "A");
	sprintf ( aln, "%s", fname+1);
	return 1;
	}
    else if ( fname[0]=='L')
	{
	sprintf ( out_mode, "L");
	sprintf ( aln, "%s", fname+1);
	return 1;
	}
    else if ( fname[0]=='M' && is_in_pre_set_method_list (fname+1))
	{	
	preset_method=1;
	fname++;
	}
    else if ( is_in_pre_set_method_list (fname))
        {
	preset_method=1;	
	}
	  
   
    

    n_param   =vcalloc (4, sizeof (int));
    param_list=vcalloc (4, sizeof (char **));
    
    sprintf ( executable, "EXECUTABLE");
    perl_pg_param[0]=executable;

    sprintf ( aln_mode, "ALN_MODE");
    perl_pg_param[1]=aln_mode;
    
    sprintf ( out_mode, "OUT_MODE");
    perl_pg_param[2]=out_mode;
    
    
    for (tot_n_param=0; tot_n_param<n_perl_pg_param; tot_n_param++)
        {
	param_list[tot_n_param]=get_parameter (perl_pg_param[tot_n_param],&n_param[tot_n_param],fname);
	sprintf (perl_pg_param[tot_n_param] , "%s", param_list[tot_n_param][1]);
	n_param[tot_n_param]--;
	}   
  
   while ((param_list[tot_n_param]=get_parameter ("PARAM",&n_param[tot_n_param],NULL))!=NULL)
	{
	n_param[tot_n_param]-=2; 
	tot_n_param++;
	param_list=vrealloc ( param_list, sizeof(char**)*(tot_n_param+1));
	n_param=   vrealloc ( n_param   , sizeof(int)   *(tot_n_param+1));
	}    
   
    combination_table=make_recursive_combination_table (tot_n_param, n_param, &n_combinations, NULL,0);

    for (c=0, a=0; a< n_combinations; a++)
        {
	

	if ( strcmp (aln_mode, "multiple")==0)
	       {
	       seq_c[0]=realloc_char(seq_c[0], read_size_char(seq_c[0][-1],0),0, 1,0);
	       aln_c[0]=realloc_char(aln_c[0], read_size_char(aln_c[0][-1],0),0, 1,0);
	       
	       sprintf (seq_c_buf, "%d",S->nseq);
	       
	       for (d=0; d< S->nseq; d++)
	           {
		   sprintf ( buf," %d",d);
		   strcat (seq_c_buf, buf);
		   }
	       seq_c[0][c]=vcalloc (strlen(seq_c_buf)+1, sizeof (char));
	       sprintf ( seq_c[0][c],"%s",seq_c_buf);

	       make_aln_command (aln_c_buf, fname, seq, aln);
	       for ( b=n_perl_pg_param; b< tot_n_param; b++)
	            {
		    sprintf ( buf, " %s %s ", (strm (param_list[b][1],"no_name"))?"":param_list[b][1], param_list[b][combination_table[a][b]+2]);
		    strcat ( aln_c_buf, buf);
		    } 
	       if (strrchr(aln_c_buf, '>')==NULL)strcat (aln_c_buf,TO_NULL_DEVICE);
			
	       aln_c[0][c]=vcalloc ( strlen (aln_c_buf)+1, sizeof (char));			 
	       sprintf (aln_c[0][c], "%s", aln_c_buf);
	       c++;	       	       
	       }
	       
	else if ( strm2 (aln_mode, "pairwise", "m_pairwise"))
	       {
	       do_mirror=strm (aln_mode, "m_pairwise");
	       for (x=0; x< S->nseq-(1-do_mirror); x++)
		   for ( y=((do_mirror)?0:(x+1)); y< S->nseq; y++)
		       for ( d=0; d<n_combinations && x!=y ; d++, c++) 
		         { 
			 seq_c[0]=realloc_char(seq_c[0], read_size_char(seq_c[0][-1],0),0, 1,0);
			 aln_c[0]=realloc_char(aln_c[0], read_size_char(aln_c[0][-1],0),0, 1,0);
			 			 
			 sprintf (seq_c_buf, "2 %d %d",x,y);
			 seq_c[0][c]=vcalloc (strlen(seq_c_buf)+1, sizeof (char));
			 sprintf ( seq_c[0][c],"%s",seq_c_buf);
			 			 			 
			 make_aln_command (aln_c_buf, fname, seq, aln);
			 for ( b=n_perl_pg_param; b< tot_n_param; b++)
			     {
			     sprintf ( buf, " %s %s ",  (strm (param_list[b][1],"no_name"))?"":param_list[b][1], param_list[b][combination_table[d][b]+2]);
			     strcat ( aln_c_buf, buf);
			     } 
			 if (strrchr(aln_c_buf, '>')==NULL)strcat (aln_c_buf, TO_NULL_DEVICE);
			 aln_c[0][c]=vcalloc ( strlen (aln_c_buf)+1, sizeof (char));			 
			 sprintf (aln_c[0][c], "%s", aln_c_buf);
			 }
	       }
	else if ( strcmp (aln_mode, "s_pairwise")==0)
	       {
	       for (x=0; x< S->nseq; x++)
		   for ( y=x; y< S->nseq; y++)
		       for ( d=0; d<n_combinations; d++, c++) 
		         {
			 seq_c[0]=realloc_char(seq_c[0], read_size_char(seq_c[0][-1],0),0, 1,0);
			 aln_c[0]=realloc_char(aln_c[0], read_size_char(aln_c[0][-1],0),0, 1,0);
			 			 
			 sprintf (seq_c_buf, "2 %d %d",x,y);
			 seq_c[0][c]=vcalloc (strlen(seq_c_buf)+1, sizeof (char));
			 sprintf ( seq_c[0][c],"%s",seq_c_buf);
			 			 			 			 
			 make_aln_command (aln_c_buf, fname, seq, aln); 
			 for ( b=n_perl_pg_param; b< tot_n_param; b++)
			     {
			     sprintf ( buf, " %s %s ", (strm (param_list[b][1],"no_name"))?"":param_list[b][1], param_list[b][combination_table[d][b]+2]);
			     strcat ( aln_c_buf, buf);
			     } 
			 if (strrchr(aln_c_buf, '>')==NULL)strcat (aln_c_buf, TO_NULL_DEVICE);
			 
			 aln_c[0][c]=vcalloc ( strlen (aln_c_buf)+1, sizeof (char));			 
			 sprintf (aln_c[0][c], "%s", aln_c_buf);			 
			 }
	       }
	}
        
    for ( a=0; a< tot_n_param; a++)free_char ( param_list[a],-1);
    free ( param_list);
    free ( n_param);

    free_int ( combination_table,-1);
    
    if ( preset_method==1)remove (fname);
   
    return c;
    }
    
int is_in_pre_set_method_list (char *fname)
    {
    FILE *fp;
    static char *fname2;
    char clustalw[LONG_STRING];
    char lalign  [LONG_STRING];



    
#ifndef     CLUSTALW_4_TCOFFEE
    if ( strncmp ( fname, "clustalw", 8)!=0);
    else if ( getenv ( "CLUSTALW_4_TCOFFEE")==NULL)crash ("CLUSTALW_4_TCOFFEE IS NOT DEFINED");
    else  sprintf ( clustalw, "%s", (getenv ( "CLUSTALW_4_TCOFFEE")));
#else
    
     if ( getenv ( "CLUSTALW_4_TCOFFEE")==NULL)sprintf ( clustalw, "%s", CLUSTALW_4_TCOFFEE);
     else  sprintf ( clustalw, "%s", (getenv ( "CLUSTALW_4_TCOFFEE")));
#endif


#ifndef     LALIGN_4_TCOFFEE
    if ( strncmp ( fname, "lalign",6)!=0);
    else if ( getenv ( "LALIGN_4_TCOFFEE")==NULL)crash ("LALIGN_4_TCOFFEE IS NOT DEFINED");
    else  sprintf ( lalign, "%s", (getenv ( "LALIGN_4_TCOFFEE")));
#else
    
     if ( getenv ( "LALIGN_4_TCOFFEE")==NULL)sprintf (lalign, "%s", LALIGN_4_TCOFFEE);
     else sprintf ( lalign, "%s", (getenv ( "LALIGN_4_TCOFFEE")));
#endif

   if (fname2==NULL)fname2=vcalloc (1000, sizeof (char));

    sprintf ( fname2, "%s",fname);
    fname=vtmpnam ( fname);
    if (strm (fname2, "?"))fprintf ( stderr, "\nAVAILABLE METHODS:\n\n");
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tclustalw_pair     \t[%5s INSTALLED]\n"      ,(pg_is_installed(clustalw))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tclustalw_pair_perl\t[%5s INSTALLED]\n" ,(pg_is_installed("clustalw.perl"))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tclustalw_aln      \t[%5s INSTALLED]\n"       ,(pg_is_installed(clustalw))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tpraline_perl_perl \t[%5s INSTALLED]\n"  ,(pg_is_installed("praline.perl"))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tpraline_aln_perl  \t[%5s INSTALLED]\n"   ,(pg_is_installed("praline.perl"))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tprrp_aln_perl     \t[%5s INSTALLED]\n"      ,(pg_is_installed("prrp.perl"))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tprrp_pair_perl    \t[%5s INSTALLED]\n"      ,(pg_is_installed("prrp.perl"))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tdialign2_aln_perl \t[%5s INSTALLED]\n"  ,(pg_is_installed("dialign2.perl"))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tdialign2_pair_perl\t[%5s INSTALLED]\n"  ,(pg_is_installed("dialign2.perl"))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tlalign_rs_pair    \t[%5s INSTALLED]\n"     ,(pg_is_installed(lalign))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tlalign_id_pair    \t[%5s INSTALLED]\n"     ,(pg_is_installed(lalign))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tlalign_id_m_pair  \t[%5s INSTALLED]\n"     ,(pg_is_installed(lalign))?"":"NOT" );
    if (strm (fname2, "?"))fprintf ( stderr, "METHOD:\tlalign_s_pair     \t[%5s INSTALLED]\n"      ,(pg_is_installed(lalign))?"":"NOT" );
    if (strm (fname2, "?"))exit (1);

    /*IMPORTANT METHODS MUST START WITH A LOWERCASE*/
    /*1 INTERNAL METHODS*/
   
    if ( strm ( fname2, "fast_pair"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S fast_pair\n");
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
    else if ( strm ( fname2, "diag_fast_pair"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S diag_fast_pair\n");
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
     else if ( strm ( fname2, "slow_pair"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S slow_pair\n");
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
    /*STRUCTURE SUPER-IMPOSITION*/
    else if ( strm ( fname2, "align_pdb_pair"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S align_pdb_pair\n");
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
    else if ( strm ( fname2, "align_pdb_pair_2"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S align_pdb_pair_2\n");
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
     else if ( strm4 ( fname2, "custom1_align_pdb_pair", "custom2_align_pdb_pair","custom3_align_pdb_pair","custom4_align_pdb_pair"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n", fname2);
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
    else if ( strm4 ( fname2, "custom5_align_pdb_pair", "custom6_align_pdb_pair","custom7_align_pdb_pair","custom8_align_pdb_pair"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n", fname2);
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
    else if ( strm2 ( fname2, "custom9_align_pdb_pair", "custom10_align_pdb_pair"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n", fname2);
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
    else if ( strm ( fname2, "sap_pair"))
           {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S sap_pair\n");
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   }
    
     else if ( strm ( fname2, "cdna_fast_pair") || strncmp ( fname2,"cdna_fast_pair",14)==0 )
	   {
	     /*cdna_fast_pair(fgop_fgep) fgop and fgep are turned into neg val later*/
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n", fname2);
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   } 
    else if ( strm ( fname2, "cdna_cfast_pair"))
	   {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S cdna_cfast_pair\n");
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	   } 
    
    else if ( strm ( fname2, "ifast_pair"))
      {
	    fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S ifast_pair\n");
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
      }
     
    /*2 EXTERNAL METHODS*/
    
    else if ( strm (fname2, "clustalw_pair")  )
	{
	  
	    
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n",clustalw );
	   fprintf ( fp, "ALN_MODE   S pairwise\n");
	   fprintf ( fp, "OUT_MODE   S aln\n");
	   fprintf ( fp, "IN_FLAG    S %sINFILE=\n", CWF);
	   fprintf ( fp, "OUT_FLAG   S %sOUTFILE=\n",CWF);
	   fprintf ( fp, "DEFAULT    S %sOUTORDER=INPUT %sNEWTREE=core %salign\n",CWF,CWF,CWF);
	   fclose (fp); 	   
	   return 1;
	}
    
     else if ( strm ( fname2, "clustalw_aln"))
        {
	    fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S clustalw_aln\n");
	   fprintf ( fp, "ALN_MODE   S multiple\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	}
    else if ( strm ( fname2, "prrp_aln"))
        {
	    fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S prrp_aln\n");
	   fprintf ( fp, "ALN_MODE   S multiple\n");
	   fprintf ( fp, "OUT_MODE   S fL\n");
	   fprintf ( fp, "IN_FLAG    S no_name\n");
	   fprintf ( fp, "OUT_FLAG   S no_name\n");
	   fclose (fp); 	   
	   return 1;
	}
   
    else  if ( strcmp (fname2, "praline_pair_perl")==0)
	{
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S praline.perl\n");
	   fprintf ( fp, "ALN_MODE S pairwise\n");
	   fprintf ( fp, "OUT_MODE S aln\n");
	   fprintf ( fp, "IN_FLAG   S -in&\n");
	   fprintf ( fp, "OUT_FLAG  S -out&\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   
	    return 1;
	}
    
    else if ( strcmp (fname2, "praline_aln_perl")==0)
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S praline.perl\n");
	   fprintf ( fp, "ALN_MODE S multiple\n");
	   fprintf ( fp, "OUT_MODE S aln\n"); 
	   fprintf ( fp, "IN_FLAG   S -in&\n");
	   fprintf ( fp, "OUT_FLAG  S -out&\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   
	    return 1;    
	}
    else if ( strcmp (fname2, "prrp_aln_perl")==0)
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S prrp.perl\n");
	   fprintf ( fp, "ALN_MODE S multiple\n");
	   fprintf ( fp, "OUT_MODE S aln\n"); 
	   fprintf ( fp, "IN_FLAG   S -in&\n");
	   fprintf ( fp, "OUT_FLAG  S -out&\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   
	    return 1;    
	}  
    
    else if ( strcmp (fname2, "prrp_pair_perl")==0)
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S prrp.perl\n");
	   fprintf ( fp, "ALN_MODE S pairwise\n");
	   fprintf ( fp, "OUT_MODE S aln\n"); 
	   fprintf ( fp, "IN_FLAG   S -in&\n");
	   fprintf ( fp, "OUT_FLAG  S -out&\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   
	    return 1;    
	} 
    else if ( strcmp (fname2, "dialign2_aln_perl")==0)
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S dialign2.perl\n");
	   fprintf ( fp, "ALN_MODE S multiple\n");
	   fprintf ( fp, "OUT_MODE S aln\n");
	   fprintf ( fp, "IN_FLAG   S -in&\n");
	   fprintf ( fp, "OUT_FLAG  S -out&\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   
	    return 1;    
	}
     else if ( strcmp (fname2, "dialign2_pair_perl")==0)
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S dialign2.perl\n");
	   fprintf ( fp, "ALN_MODE S pairwise\n");
	   fprintf ( fp, "OUT_MODE S aln\n");
	   fprintf ( fp, "IN_FLAG   S -in&\n");
	   fprintf ( fp, "OUT_FLAG  S -out&\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   
	    return 1;    
	}
    else if ( strm2 (fname2, "lalign_pair", "lalign_rs_pair"))
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n",lalign);
	   fprintf ( fp, "ALN_MODE S pairwise\n");
	   fprintf ( fp, "OUT_MODE S lib\n");
	   fprintf ( fp, "IN_FLAG   S no_name\n");
	   fprintf ( fp, "OUT_FLAG  S no_name\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   fclose (fp);
	   return 1;    
	}  
    else if ( strm (fname2,"lalign_id_pair"))
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n", lalign);
	   fprintf ( fp, "ALN_MODE S pairwise\n");
	   fprintf ( fp, "OUT_MODE S lib\n");
	   fprintf ( fp, "PARAM S no_name 0\n");
	   fprintf ( fp, "PARAM S no_name 1\n");
	   fprintf ( fp, "IN_FLAG   S no_name\n");
	   fprintf ( fp, "OUT_FLAG  S no_name\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   return 1;    
	}  
    else if ( strm (fname2,"lalign_id_m_pair"))
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n", lalign);
	   fprintf ( fp, "ALN_MODE S m_pairwise\n");
	   fprintf ( fp, "OUT_MODE S lib\n");
	   fprintf ( fp, "PARAM S no_name 0\n");
	   fprintf ( fp, "PARAM S no_name 1\n");
	   fprintf ( fp, "IN_FLAG   S no_name\n");
	   fprintf ( fp, "OUT_FLAG  S no_name\n");
	   fprintf ( fp, "DEFAULT   S \n");
	   fclose (fp);
	   return 1;
	}
	 
    else if ( strcmp (fname2, "lalign_s_pair")==0)
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n",lalign);
	   fprintf ( fp, "ALN_MODE  S s_pairwise\n");
	   fprintf ( fp, "OUT_MODE  S lib\n");
	   fprintf ( fp, "PARAM S no_name -30\n");
	   fprintf ( fp, "PARAM S no_name   1\n");
	   fprintf ( fp, "IN_FLAG   S no_name\n");
	   fprintf ( fp, "OUT_FLAG  S no_name\n");
	   fprintf ( fp, "DEFAULT   S \n");
	  
	   fclose (fp);
	   return 1;    
	}
    else if ( strcmp (fname2, "lalign_rs_s_pair")==0)
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n",lalign);
	   fprintf ( fp, "ALN_MODE  S s_pairwise\n");
	   fprintf ( fp, "OUT_MODE  S lib\n");
	   fprintf ( fp, "PARAM S no_name -20\n");
	   fprintf ( fp, "PARAM S no_name   0\n");
	   fprintf ( fp, "IN_FLAG   S no_name\n");
	   fprintf ( fp, "OUT_FLAG  S no_name\n");
	   fprintf ( fp, "DEFAULT   S \n");
	  
	   fclose (fp);
	   return 1;    
	}
    else if ( strcmp (fname2, "lalign_rs_s_dna_pair")==0)
        {
	   fp=vfopen (fname, "w");
	   fprintf ( fp, "EXECUTABLE S %s\n",lalign);
	   fprintf ( fp, "ALN_MODE  S s_pairwise\n");
	   fprintf ( fp, "OUT_MODE  S lib\n");
	   fprintf ( fp, "PARAM S no_name  -40\n");
	   fprintf ( fp, "PARAM S no_name   0\n");
	   fprintf ( fp, "IN_FLAG   S no_name\n");
	   fprintf ( fp, "OUT_FLAG  S no_name\n");
	   fprintf ( fp, "DEFAULT   S \n");
	  
	   fclose (fp);
	   return 1;    
	}
    else if ( check_file_exists ( fname2)==1)
        {
	    sprintf ( fname, "%s", fname2);
	    return 0;
	}
    else
        {	
	    sprintf ( fname, "%s", fname2);
	    return 0;
        }
    }

char *make_aln_command(char *command, char *parafile, char *seq, char *aln)
    {
    int a;
    char **executable;
    char **translate_list;
    char **in_flag;
    char **out_flag;
    char **default_p;
    int n;
    int n_default_p;

    translate_list=declare_char ( 1,3);

    
    sprintf (translate_list[0], "& ");
    executable=get_parameter   ( "EXECUTABLE",&n        ,parafile);
    default_p   =get_parameter ( "DEFAULT"   ,&n_default_p,parafile);
    in_flag   =get_parameter   ( "IN_FLAG"   ,&n        ,parafile);
    out_flag  =get_parameter   ( "OUT_FLAG"  ,&n        ,parafile);
    get_parameter ( "CLOSE THE FILE",&n, NULL);

    sprintf ( command, "%s ", executable[1]);

    if ( !strm ( in_flag[1], "no_name")){string_convert (in_flag[1], 1, translate_list);strcat ( command, in_flag[1]);}
    strcat ( command, seq);
    strcat ( command, " ");

   

    for ( a=1; a< n_default_p; a++){strcat ( command," "); strcat ( command, default_p[a]);}
    strcat ( command, " ");
    if ( !strm ( out_flag[1], "no_name")){string_convert (out_flag[1], 1,translate_list);strcat ( command, out_flag[1]);}
    strcat ( command, aln);
    
    
    free_char ( executable, -1);
    free_char ( translate_list, -1);
    free_char ( in_flag, -1);
    free_char ( out_flag,-1);
    free_char ( default_p, -1);   
    
    return command;
    }
	

     
     

/*********************************************************************/
/*                                                                   */
/*                         WRITE IN LIST                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

       
    
int vread_clist ( Constraint_list *CL, int a, int b )
    {
    int x;
    
    if ( a>= CL->ne)return UNDEFINED_2;
    else if (CL->fp)
       {       
       fseek (CL->fp, a*CL->el_size*CL->entry_len+b*CL->el_size, SEEK_SET);
       fread (&x, CL->el_size, 1, CL->fp);
       return x;
       }
    else if ( CL->L)
       {
       return (CL->L)[a][b];
       }
    else if (CL->M)
       {
       return (CL->M)[a][b];
       }
    else 
       {
       return UNDEFINED_2; 
       
       }
    }
int vwrite_clist ( Constraint_list *CL, int a, int b, CLIST_TYPE x)
    {
    
    CL->seq_indexed=0;
    CL->residue_indexed=0;

    if (CL->fp)
       {
       fseek (CL->fp, a*CL->el_size*CL->entry_len+b*CL->el_size, SEEK_SET);
       fwrite(&x, CL->el_size, 1, CL->fp);
       }
    else if (!CL->M)
       {
       if (a>=CL->max_L_len)
          {
	  (CL->L)=realloc_int ( (CL->L), CL->max_L_len, CL->entry_len,CL->chunk, 0);
	  CL->max_L_len+=CL->chunk;
	  }
       (CL->L)[a][b]=x;
       }
    return x;
    }


/*********************************************************************/
/*                                                                   */
/*                        INDEXING FUNCTIONS                         */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

Constraint_list *index_res_constraint_list ( Constraint_list *CL, int field)
        {
	/*
	  function documentation: start
	  Constraint_list *index_res_constraint_list ( Constraint_list *CL, int field)
	  This function reorganises the content of the CL->L matrix, so that a single look up gives
	  the constraints associated with any residue
	  
	  1-if CL->residue_indexed=1 return
		 2-CL->residue_index[Seq X][Res Y of Seq X][0]=Z
		 Z=Number of residues matching X(Y)*3+1
		 CL->residue_index[Seq X][Res Y of Seq X][0]=Z
		 CL->residue_index[Seq X][Res Y of Seq X][c+0]=seq W
		 CL->residue_index[Seq X][Res Y of Seq X][c+1]=res V
		 CL->residue_index[Seq X][Res Y of Seq X][c+2]=weight W(V) Vs X(Y)		 
	  
	  NOTE: Works All right with a sequence against itself
	  NOTE: Any modification of CL->L should result in residue_indexed to be reset	  
	  function documentation: end
	*/
	int a, b, s1, s2, r1, r2, w;


       

	if ( CL->residue_indexed && CL->residue_field==field);
	else
	   {
	       if ( CL->residue_index)
	          {
		  for ( a=0; a< (CL->S)->nseq; a++)
		      for ( b=0; b<= (CL->S)->len[a]; b++)
		          {
			      free(CL->residue_index[a][b]);
			      CL->residue_index[a][b]=vcalloc (1, sizeof (int));
			      CL->residue_index[a][b][0]=1;
			  }    
		  }
	       else if ( !CL->residue_index)
	          {
		      CL->residue_index=vcalloc ( (CL->S)->nseq, sizeof (int**));
		      for ( a=0; a< (CL->S)->nseq; a++)
		          {
			      CL->residue_index[a]=vcalloc ( (CL->S)->len[a]+1, sizeof (int*));
			      for ( b=0; b<= (CL->S)->len[a]; b++)
				  {
				  CL->residue_index[a][b]=vcalloc (1, sizeof (int));
				  CL->residue_index[a][b][0]=1;
				  }
			  }
		  }
	       for (a=0;a<CL->ne; a++)
	          {
		  s1=vread_clist (CL, a, SEQ1);
		  s2=vread_clist (CL, a, SEQ2);
		  r1=vread_clist (CL, a, R1);
		  r2=vread_clist (CL, a, R2);
		  w=vread_clist (CL, a, field);
		  
		  CL->residue_index[s1][r1][0]+=3;
		  CL->residue_index[s1][r1]=vrealloc ( CL->residue_index[s1][r1], CL->residue_index[s1][r1][0]*sizeof (int));
		  CL->residue_index[s1][r1][CL->residue_index[s1][r1][0]-3]=s2;
		  CL->residue_index[s1][r1][CL->residue_index[s1][r1][0]-2]=r2;
		  CL->residue_index[s1][r1][CL->residue_index[s1][r1][0]-1]=w;

		  CL->residue_index[s2][r2][0]+=3;
		  CL->residue_index[s2][r2]=vrealloc ( CL->residue_index[s2][r2], CL->residue_index[s2][r2][0]*sizeof (int));
		  CL->residue_index[s2][r2][CL->residue_index[s2][r2][0]-3]=s1;
		  CL->residue_index[s2][r2][CL->residue_index[s2][r2][0]-2]=r1;
		  CL->residue_index[s2][r2][CL->residue_index[s2][r2][0]-1]=w;
		  
		  }
	   CL->residue_indexed=1;
	   CL->residue_field=field;
	   }
	return CL;
	}
      
Constraint_list *index_constraint_list ( Constraint_list *CL)
        {
	  /*
	  Function Documentation: start
	  Constraint_list *index_constraint_list ( Constraint_list *CL);
	  Indexes a constraint list
	     1-Checks the flag seq_indexed
	     2-if flag set to 0
	       CL->start_index[seq1][seq2] indicatse the first position for seq1 Vs seq2
	       CL->end_index[seq1][seq3]   indicatse the last  position for seq1 Vs seq2
	  Any modif to CL->L should cause the flag 1 to be set to 0;
	  Function Documentation: end
	*/
	int a, csA, csB, sA, sB;


	if ( CL->seq_indexed);
	else
	    {
	    
	    if ( CL->start_index!=NULL)free_int ( CL->start_index,-1);
	    CL->start_index=declare_int ( (CL->S)->nseq, (CL->S)->nseq);
	    
	    if ( CL->end_index!=NULL)free_int ( CL->end_index,-1);
	    CL->end_index=declare_int ( (CL->S)->nseq, (CL->S)->nseq);
	    
	    csA=vread_clist (CL, 0, SEQ1);
	    csB=vread_clist (CL, 0, SEQ2);

	    CL->start_index[csA][csB]=0;
	    CL->start_index[csB][csA]=0;
	    for ( a=1; a<CL->ne; a++)
	        {
		sA=vread_clist (CL, a, SEQ1);
		sB=vread_clist (CL, a, SEQ2);
		if (sA!=csA || sB!=csB)
		   {
		   CL->end_index[csA][csB]=a;
		   CL->end_index[csB][csA]=a;
		   csA=sA;
		   csB=sB;
		   CL->start_index[csA][csB]=a;
		   CL->start_index[csB][csA]=a;
		   }
		}		
	    CL->end_index[csB][csA]=CL->ne;
	    CL->end_index[csA][csB]=CL->ne;
	    CL->seq_indexed=1;
	    }
	return CL;
	}
	
/*********************************************************************/
/*                                                                   */
/*                         LIST EXTENTION                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *extend_list_pair (Constraint_list *CL,char *store_mode, int s1, int s2)
        {
	static Sequence *S;
	Constraint_list *CLout;
	/*
	  function documentation: start
	  Constraint_list *extend_list_pair (Constraint_list *CL,char *store_mode, int s1, int s2)
	  This function takes a pair of sequences s1, s2 and perrforms the extention
	  It returns the incoming list CL, with CL->L[s1][s2] now extended
	  See main documentation for store_mode
	  function documentation: end
	*/

	if ( S==NULL)S=declare_sequence ((CL->S)->min_len, (CL->S)->max_len,(CL->S)->nseq); 	
	sprintf ( S->name[0], "%s",(CL->S)->name[s1]);
	sprintf ( S->name[1],"%s",(CL->S)->name[s2]);
	S->nseq=2;
	
	CLout=extend_list (CL, store_mode, CL->extend_clean_mode, CL->extend_compact_mode,CL->do_self, S);
	return CLout;		
	}
Constraint_list *extend_list (Constraint_list *CLin, char *store_mode,char *clean_mode, char *compact_mode,int do_self, Sequence *DO_LIST)
	{
	int a, b, c, d, e, f;
	int wA, wC,w, rA, rC, miscA, miscC, misc;
	static int **posA;
	static int **posC;
	int start_ba, end_ba, start_bc, end_bc, start_ac, end_ac;
	int len;
	int lenA=0;
	int lenC=0;
	int *translation;
	Constraint_list *CLout=NULL;


	/*Do not extend if the List is a Matrix*/
	if ( !CLin->L && CLin->M)
	   {
	   CLin->extend_jit=0;
	   return CLin;
	   }
	
	translation=vcalloc ( (CLin->S)->nseq, sizeof (int));
	for ( a=0; a<(CLin->S)->nseq; a++)
	    {
	    translation[a]=name_is_in_list ((CLin->S)->name[a],DO_LIST->name, DO_LIST->nseq, 100);
	    translation[a]++;/* set translation to -1+1=0 if seq not in list*/
	    }
	
	CLout=declare_constraint_list (CLin->S, NULL,NULL,0, strm("disk", store_mode)?tmpfile():NULL, NULL);
        	
	for ( a=0; a<(CLin->S)->nseq-(1-do_self); a++)
		{	        
		fprintf (CLin->local_stderr, "\nSeq %3d: %5d", a+1,CLout->ne);
		for ( c=a+(1-do_self); c<(CLin->S)->nseq; c++)
			{
			if ( translation[a] && translation[c])
			       {			       
			       get_bounds (CLin,a, c, &start_ac, &end_ac); 						      
			       for ( d=start_ac; d<end_ac; d++)
			           {
				   for ( e=0; e< CLin->entry_len; e++)
				       vwrite_clist(CLout,CLout->ne, e, vread_clist(CLin,d, e));
			           CLout->ne++;
				   }
			       
			       for ( b=0; b<(CLin->S)->nseq; b++)
			           {
				   len=strlen ( (CLin->S)->seq[b]);
				   
				   get_bounds (CLin,b, a, &start_ba, &end_ba);  
				   posA=fill_pos_matrix (CLin,start_ba, end_ba, len, posA, &lenA,(b>a));
				   
				   if ((c!=b && a!=b) ||(do_self==1))
				       {
				   
				       get_bounds (CLin, b, c, &start_bc, &end_bc);
				       posC=fill_pos_matrix (CLin, start_bc, end_bc, len, posC, &lenC, (b>c));
				       
				       for (d=1; d<=len; d++)
				           {
					   if ( posA[d][1]==0 || posC[d][1]==0);
					   else
					       {
					       for (e=2; e<=posA[d][1]+1; e+=(CLin->entry_len-4)) 
						   for ( f=2; f<=posC[d][1]+1; f+=(CLin->entry_len-4))
						       {
						       wA   =posA[d][e+1];
						       miscA=posA[d][e+2];

						       wC   =posC[d][f+1];
						       miscC=posC[d][f+2];

						       rA=posA[d][e];
						       rC=posC[d][f];
						       
						       w   =MIN(wA,wC);
						       
						       misc=MAX(miscA, miscC);
						       
						       vwrite_clist( CLout, CLout->ne, SEQ1, a);
						       vwrite_clist( CLout, CLout->ne, SEQ2, c);
						       vwrite_clist( CLout, CLout->ne, R1  ,rA);
						       vwrite_clist( CLout, CLout->ne, R2  ,rC);
						       vwrite_clist( CLout, CLout->ne, WE  , w);
						       vwrite_clist( CLout, CLout->ne, CONS, 1);
						       vwrite_clist( CLout, CLout->ne, MISC,misc);
        					       CLout->ne++;
						       }
					       }
					   }
				       }
				   }
		
			       CLout=compact_list (CLout,0,CLout->ne,"mirror_sum");
			       CLout=clean ( clean_mode,CLout, 0, CLout->ne);
			       }
			}
		}

	
	free (translation);
	return CLout;
	}
void get_bounds (Constraint_list *CL, int s1, int s2, int *start, int *end)
	{

	CL=index_constraint_list (CL);
	
	if ( s1>s2)SWAP(s1, s2);
	
	start[0]=CL->start_index[s1][s2];
	end  [0]=CL->end_index  [s1][s2];
	}

	
int ** fill_pos_matrix (Constraint_list *CL, int beg, int end, int slen, int **pos, int *len, int mirrored)
	{
	int small_chunck;
	int a, r1,r2;
	

	small_chunck=2*CL->entry_len;

	if ( pos==NULL)
		{
		pos=declare_int (slen+1, small_chunck);
		for ( a=0; a<=slen; a++)pos[a][0]=small_chunck;
		len[0]=slen+1;
		}
	else if ( len[0]<=slen)
		{
		free_int ( pos, len[0]);
		pos=declare_int (slen+1, small_chunck);
		for ( a=0; a<=slen; a++)pos[a][0]=small_chunck;
		len[0]=slen+1;
		}
	else
		{
		for ( a=0; a<=slen; a++)pos[a][1]=0;
		}
			
			
	
	
	for ( a=beg; a<end; a++)
		{
		
		if (!mirrored)     {r1=vread_clist (CL, a, R1);r2=vread_clist (CL, a, R2);}
		else if ( mirrored){r1=vread_clist (CL, a, R2);r2=vread_clist (CL, a, R1);}

	       if ( ((pos[r1][1])+(CL->entry_len))>pos[r1][0])
			{
			pos[r1]=vrealloc (pos[r1], (pos[r1][0]+small_chunck)*sizeof (int));
			pos[r1][0]+=small_chunck;
			}
		pos[r1][pos[r1][1]+2]=r2;
		pos[r1][pos[r1][1]+3]=vread_clist(CL,a,WE);
		pos[r1][pos[r1][1]+4]=vread_clist(CL,a,MISC);
		pos[r1][1]+=(CL->entry_len-4);		
		}
	return pos;
	}
Constraint_list * evaluate_constraint_list_reference ( Constraint_list *CL)
        {
	    static CLIST_TYPE *entry;
	    
	    int a, b, c, s1, s2, r1, r2, w;
	    int ***max_res;
	    

	    CL->max_value=CL->max_ext_value=0;
	    max_res=vcalloc ( (CL->S)->nseq, sizeof (int**));
	    for ( a=0; a< (CL->S)->nseq; a++)
	        {
		    
		    max_res[a]=vcalloc ( strlen ((CL->S)->seq[a])+1, sizeof (int*));
		    for ( b=0; b<=(CL->S)->len[a]; b++)
			max_res[a][b]=vcalloc ( (CL->S)->nseq+1, sizeof (int));
		}

	     for ( a=0; a< CL->ne; a++)
	         {
		     
		     entry=extract_entry ( entry, a, CL);
		     s1=entry[SEQ1];
		     s2=entry[SEQ2];
		     r1=entry[R1];
		     r2=entry[R2];
		     w= entry[WE];
		     
		     if ( w==UNDEFINED || ( (CL->moca) && (CL->moca)->forbiden_residues  && ((CL->moca)->forbiden_residues[s1][r1]==UNDEFINED || (CL->moca)->forbiden_residues[s2][r2]==UNDEFINED)));
		     else
		       {
			 max_res[s1][r1][s2]+=w;
			 max_res[s2][r2][s1]+=w;
			 CL->max_value=MAX(w, CL->max_value);
		       }
		 }
	    
	     for ( a=0; a< (CL->S)->nseq; a++)
		 for ( b=1; b<=(CL->S)->len[a]; b++)
		     {
		     for ( c=0; c< (CL->S)->nseq; c++)
		         {
			     max_res[a][b][(CL->S)->nseq]+= max_res[a][b][c];
			 }
		     CL->max_ext_value=MAX(max_res[a][b][c],CL->max_ext_value);
		     }           
	     
	     for ( a=0; a< (CL->S)->nseq; a++)
	         {
		     for ( b=0; b<=(CL->S)->len[a]; b++)
			 vfree ( max_res[a][b]);
		     vfree (max_res[a]);
		 }
	     CL->max_ext_value=MAX(1,CL->max_ext_value);
	     vfree ( max_res);
	     return CL;
	}
			       
/*********************************************************************/
/*                                                                   */
/*                         ENTRY MANIPULATION                        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list * add_list_entry2list (Constraint_list *CL, int n_para, ...)
	{
	int a;
	int *entry;
	int field, val;
	va_list ap;

	if (n_para>LIST_N_FIELDS)
	   {
	       crash ("Too Many Fields in List [FATAL/add_list_entry2list]");
	   }
	
	va_start (ap,n_para);
	entry=vcalloc (CL->entry_len, sizeof (int));

	for ( a=0; a<n_para; a++)
	    {
	    field=va_arg(ap, int);
	    val  =va_arg(ap, CLIST_TYPE);
	    entry[field]=val;
	    }
	va_end (ap);
	add_entry2list(entry, CL);
	free(entry);
	return CL;
	}

Constraint_list *add_entry2list ( CLIST_TYPE *entry, Constraint_list *CL)
	{
	return insert_entry2list (entry, CL->ne++, CL);
	}
Constraint_list* insert_entry2list(CLIST_TYPE * entry, int pos, Constraint_list *CL)
        {
	int a;
        for ( a=0; a< CL->entry_len; a++)
	    vwrite_clist ( CL,pos, a,entry[a]);
	return CL;
	}
CLIST_TYPE* extract_entry(CLIST_TYPE * entry, int pos, Constraint_list *CL)
        {
	int a;
	
	if ( entry==NULL)entry=vcalloc ( CL->entry_len, CL->el_size);
	
	for (a=0; a< CL->entry_len; a++)entry[a]=vread_clist(CL, pos, a);
	return entry;
	}

	
/*********************************************************************/
/*                                                                   */
/*                         SEARCH IN LIST (ARRAY AND FILE)           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
FILE * compare_list (FILE *OUT, Constraint_list *CL1,Constraint_list *CL2)
	{
	int a;
	float nw_score=0;
	float w_score=0;
	int **l;

	CLIST_TYPE  *entry=NULL;
	int p;

	entry=vcalloc ( CL1->entry_len, CL1->el_size);
	for ( a=0; a<CL1->ne; a++)
		{
		entry=extract_entry (entry,a,CL1);
		if ((l=main_search_in_list_constraint (entry,&p,4,CL2))!=NULL)
			{
			vwrite_clist ( CL2, p,MISC, 1);
			vwrite_clist ( CL1, a,MISC, 1);
			nw_score++;
			w_score+=l[0][WE];
			}
		}
	fprintf ( OUT, "%-15s:%d pairs (Evaluated matrix), %d pairs in the other (%s)\n", CL2->list_name, CL2->ne, CL1->ne, CL1->list_name);
        fprintf ( OUT, "%-15s:%d pairs\n", CL1->list_name, CL1->ne);
        fprintf ( OUT, "Acurracy=%.2f%%\n", (nw_score/CL1->ne)*100);
        fprintf ( OUT, "Sensitiv=%.2f%%\n\n", (nw_score/CL2->ne)*100);
	return OUT;
	}


CLIST_TYPE **main_search_in_list_constraint ( int *key,int *p,int k_len,Constraint_list *CL)
	{

	CLIST_TYPE **l=NULL;
	int start, end;
		
	CL=index_constraint_list (CL);
	
	start=CL->start_index[key[SEQ1]][key[SEQ2]];
	end  =CL->end_index  [key[SEQ1]][key[SEQ2]];
		
	if ( CL->fp)
		    {
		    fseek(CL->fp, (long)start*CL->el_size*CL->entry_len, SEEK_SET);
		    l=(int **)search_in_list_file (key,p,k_len,CL->fp, end-start,CL->el_size, CL->entry_len);  		   
		    }
	else if ( CL->L)
		    {
		    l=(int **)search_in_list_array (key,p,k_len,(void **)((CL->L)+start), end-start,CL->el_size, CL->entry_len);  
		    }		
	
	return l;
	}
CLIST_TYPE return_max_constraint_list ( Constraint_list *CL, int field)
        {
	CLIST_TYPE max=0;
	int a;
	for ( a=0; a< CL->ne; a++)max=MAX( vread_clist(CL,a,field), max);
	return max;
	}
	
 /*********************************************************************/
/*                                                                   */
/*                                                                   */
/*      LIST SORTING                                                 */
/*                                                                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *sort_constraint_list_inv (Constraint_list *CL, int start, int len)
	{
	CL=sort_constraint_list   (CL, start,len);
	

	CL=invert_constraint_list (CL, start,len);
	if ( start+len==CL->ne)
	    {	   
	    while (vread_clist(CL,CL->ne-1, SEQ1)==-1)CL->ne--;
	    }
	

	return CL;
	}

Constraint_list *invert_constraint_list (Constraint_list *CL, int start,int len)
        {
	int a, b, c;
	CLIST_TYPE tp;
	
	
	for ( a=start, b=start+len-1; a<=b; a++, b--)
	    {
	    for (c=0; c< CL->entry_len; c++)
	        {
		tp=vread_clist(CL, a, c);
		vwrite_clist(CL,a, c, vread_clist(CL, b, c));
		vwrite_clist(CL,b, c, tp);
		}
	    }
	return CL;
	}
	
Constraint_list * sort_constraint_list(Constraint_list *CL, int start, int len)
        {

	CL=sort_constraint_list_on_n_fields (CL, start, len, 0, CL->entry_len);

	return CL;
	}
Constraint_list * sort_constraint_list_on_n_fields (Constraint_list *CL, int start, int len, int first_field, int n_fields)
	{

	if (CL->fp)
	   {	   
	   rewind( CL->fp);
	   fseek      ( CL->fp, start*CL->el_size*CL->entry_len , SEEK_SET);
	   hsort_list_file ( CL->fp,        len, CL->el_size, CL->entry_len,first_field,n_fields);
	   }
	else if ( CL->L)
	   {
	  
	   hsort_list_array ((void**)(CL->L)+start, len, CL->el_size, CL->entry_len,first_field,n_fields);
	   }
	return CL;
	}

/*********************************************************************/
/*                                                                   */
/*                         LIST PARSING                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
Constraint_list* read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode, char *type,FILE *local_stderr)
    {
    int a, b;
    Sequence *S;
    Constraint_list *CL=NULL;






    if ((S=read_seq_in_n_list (fname, n_list,type))==NULL)
	{
	fprintf ( stderr, "\nNO SEQUENCE WAS SPECIFIED[FATAL]\n");
        exit(1);
	}
   
   /*CHECK IF THERE IS A MATRIX AND GET RID OF OTHER METHODS*/
    for (b=0, a=0; a< n_list; a++)if (is_matrix(fname[a]) ||is_matrix(fname[a]+1) )b=a+1;

    if ( b)
       {
	if ( b==1);
	else sprintf ( fname[0], "%s", fname[b-1]);
	n_list=1;

        }   
    

    
    CL=declare_constraint_list ( S,NULL, NULL, 0,(strm(mem_mode, "disk"))?tmpfile():NULL, NULL);    
    CL->local_stderr=local_stderr;
    fprintf ( CL->local_stderr,"\nREAD/MAKE LIBRARIES:[%d]\n",n_list );
	   
    CL=read_constraint_list (CL, fname[0], in_mode, mem_mode,weight_mode);
    compact_list (CL, 0, CL->ne, "default");
    for ( a=1; a< n_list; a++)
        {	
	CL=read_constraint_list (CL, fname[a], in_mode, mem_mode,weight_mode);	
	compact_list (CL, 0, CL->ne, "default");
	}
    CL->local_stderr=local_stderr;

    

    CL=evaluate_constraint_list_reference (CL);

    return CL;
    }
Constraint_list* read_constraint_list(Constraint_list *CL,char *in_fname,char *in_mode, char *mem_mode,char *weight_mode)
        {
	FILE *fp;
	char *name;
	Constraint_list *CLtp=NULL;
	static char *read_mode;
	char *fname;

	fname=in_fname;	
	if ( !read_mode)read_mode=vcalloc ( STRING, sizeof (char));
		
	if ( in_mode)sprintf (read_mode, "%s", in_mode);     
	else if ( fname[0]=='A'){sprintf ( read_mode, "aln");fname++;}
	else if ( fname[0]=='L'){sprintf ( read_mode, "lib");fname++;}
	else if ( fname[0]=='M'){sprintf ( read_mode, "method");fname++;}
	else if ( fname[0]=='S'){sprintf ( read_mode, "sequence");return CL;}
	else if ( fname[0]=='P'){sprintf ( read_mode, "pdb")     ;return CL;}
	else if ( fname[0]=='R'){sprintf ( read_mode, "profile") ;return CL;}
	else if ( fname[0]=='X'){sprintf ( read_mode, "matrix");++fname;}    

	else
	         {
		     fprintf ( stderr, "\nERROR: The descriptor %s could not be identified as a file or a method.[FATAL]\nIf it is a method file please indicate it with M%s\n", fname, fname);
		     exit (1);
		 }

	fprintf (CL->local_stderr, "\n\t%s [%s]\n", fname, read_mode);
	
	
	if (strm(read_mode, "binary"))
	      {
	      name=vcalloc ( L_tmpnam+1, sizeof (char));
	      CLtp=declare_constraint_list ( (CL->S),NULL, NULL, 0,NULL, NULL);
	      CLtp->fp=vfopen ( fname, "rb");
	      fseek( CLtp->fp, 0, SEEK_END);
	      CLtp->ne=(int)ftell(CLtp->fp)/(CLtp->entry_len*CLtp->el_size);	
	      name=vtmpnam(name);
	      fp=bin_file2constraint_list(CL, NULL, name);
	      vfclose (fp);
	      vfclose (CLtp->fp);
	      CL=read_constraint_list_file ( CL,name);
	      remove (name); free (name);
	      free_constraint_list ( CLtp);	 
	      }
	   else if ( strm2 (read_mode,"ascii","lib"))
              {	      
	      if ( CL->cpu!=-1)increase_ref_time(read_cpu_in_list(fname));	     
	      CL=read_constraint_list_file(CL, fname);	      
	      }
	   else if (strm(read_mode, "method"))
	      {

	      CL=produce_list ( CL, CL->S, fname,weight_mode,mem_mode);
	      }
	   else if (strm(read_mode, "matrix"))
	      {
	      CL->L=NULL;
	      CL->extend_jit=0;
	      CL->M=read_matrice ( fname);
	      }
	   else if (strm (read_mode, "aln"))
	      {
		CL=aln_file2constraint_list ( fname, CL,weight_mode);
	      }
	   else 
	      {
	      if ( CL->cpu!=-1)increase_ref_time(read_cpu_in_list(fname));
	      CL=read_constraint_list_file(CL, fname);
	      }		
	return CL;
	}


Sequence * read_seq_in_n_list(char **fname, int n, char *type)
        {
	int nseq=0;
	int a, b;
	Alignment *A;
	char **sequences=NULL;
	char **seq_name=NULL;
	Sequence *S=NULL;
	Sequence *S1;
	char mode;
	char *lname;

	/*THE TYPE OF EACH FILE MUST BE INDICATED*/



	if ( n==0)
	   {
	   fprintf ( stderr, "\nNO IN FILE [FATAL]\n");
	   exit (1);
	   }
	else
	   {	   
	   for ( a=0; a< n ; a++)
	       {	       
		   
	       if (is_in_set (fname[a][0], "ALSMXPR"))
		       {
		       mode=fname[a][0];
		       lname=fname[a]+1;
		       }
	       else
	               {
			   mode=identify_format (&fname[a]);
			   lname=fname[a];
		       }
	       if ( mode=='A')
	          {
		  A=declare_aln (NULL);
		  A=main_read_aln (lname, A);
		  S1=aln2seq(A);
		 
		  if ((S=merge_seq ( S1, S))==NULL){fprintf ( stderr, "\nSequence Error in %s [FATAL]\n",lname); exit(1);} 
		  free_aln (A);
		  free_sequence (S1, S1->nseq);
		  }
	       else if ( mode=='R')
	          {
		  A=declare_aln (NULL);
		  A=main_read_aln (lname, A);
		  
		  S=add_prf2seq(A,S);
		  free_aln (A);
		  
		  }
	       else if ( mode=='P')
	          {		  
		  A=declare_aln (NULL);
		  A=main_read_aln (lname, A);
		  
		  S1=get_pdb_sequence (lname);
		  if ((S=merge_seq ( S1, S))==NULL){fprintf ( stderr, "\nSequence Error in %s [FATAL]\n",lname); exit(1);} 
		  free_sequence (S1, S1->nseq);

		  }
	       else if ( mode=='M');
	       else if ( mode=='X');
	       else if ( mode=='S')
	          {
		  /*1 Try with my routines (read t_coffee and MSF)*/ 
		  
		  if ( (A=main_read_aln ( lname, NULL))!=NULL)  
		     {
		     
		     S1=aln2seq(A);
		     free_aln(A);
		     }		  
		  else
		     { 
		     S1=main_read_seq (lname);
		     }
		  
		  for ( b=0; b< S1->nseq; b++)ungap(S1->seq[b]);
		  if ((S=merge_seq ( S1, S))==NULL){fprintf ( stderr, "\nSequence Error in %s [FATAL]\n",lname); exit(1);} 
		  
		  free_sequence (S1, S1->nseq);
		  }
	       else if (mode=='L')
                  {
		  
                  read_seq_in_list (lname,&nseq,&sequences,&seq_name);             
                  S1=fill_sequence_struc ( nseq, sequences, seq_name);
		  for ( b=0; b< S1->nseq; b++)sprintf ( S1->file[b], "%s", lname);
		  nseq=0;free_char (sequences, -1); free_char ( seq_name, -1);
		  sequences=NULL;
		  seq_name=NULL;
		  
                  if ((S=merge_seq( S1, S))==NULL){fprintf ( stderr, "\nSequence Error in %s [FATAL]\n",lname); exit(1);} 
		  free_sequence(S1, S1->nseq);
                  }

	       else 
	          {
		      fprintf ( stderr, "\nERROR: %s is neither a file nor a method [FATAL]\n", lname);
		      exit (1);
		  }
	       }
	   
	   if ( type && type[0])sprintf ( S->type, "%s", type);
	   else S=get_sequence_type (S);
	   
	   if ( strm (S->type, "PROTEIN_DNA"))
	      {
		  for ( a=0; a< S->nseq; a++)
		      {
			  if (strm ( get_string_type ( S->seq[a]), "DNA"));
			  else if ( strm ( get_string_type ( S->seq[a]), "PROTEIN"))
			      {
				  S->seq[a]=thread_aa_seq_on_dna_seq (S->seq[a]);
				  S->len[a]=strlen (S->seq[a]);
				  S->max_len=MAX(S->max_len, S->len[a]);
			      }
		      }
	      }
			
	   
	   return S;
	   }


	return NULL;
	}

int read_cpu_in_list ( char *fname)
        {
	FILE *fp;
	int c;
	int cpu=0;

	fp=vfopen ( fname, "r");
	while ( (c=fgetc(fp))!='#');
	while ( (c=fgetc(fp))!='C' && c!=EOF);
	if ( c=='C')fscanf( fp, "PU %d\n", &cpu);
	vfclose ( fp);
	return cpu;
	}


Constraint_list * read_constraint_list_file(Constraint_list *CL, char *fname)
        {
	int a, b, c,e,n,z;
	int seq_len, sn;
	int s1, s2;
	FILE *fp;
	static char *name;
	char *sequence;
	static char *mat;
	static char *dp_mode;
	int max_nseq=0;
	static int *sn_list;
	static int line=2;
	int list_nseq;
	static CLIST_TYPE *entry;
	Sequence *S;
	int seq_1_to_n=0;

	int **cache;
	Alignment *B;
	char *buf;
	int r, d;
	int lline;

	
	lline=measure_longest_line_in_file (fname)+1;

	if ( !mat) mat=vcalloc (STRING, sizeof (char));
	if ( !dp_mode) dp_mode=vcalloc (STRING, sizeof (char));
	fp=vfopen (fname, "r");
	while((c=fgetc(fp))!='#')if ( c=='\n')max_nseq++;
	vfclose (fp);
	
	buf=vcalloc (lline, sizeof (char));	
	sequence=vcalloc (lline, sizeof (char));
	if ( !name)name=vcalloc ( 100, sizeof (char));
	if ( !entry)entry=vcalloc ( CL->entry_len, CL->el_size);
	if ( !sn_list)sn_list=vcalloc (max_nseq, sizeof (int));
	else 
	  {
	    sn_list=vrealloc (sn_list, max_nseq*sizeof (int));
	  }
	S=CL->S;
	cache=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	for ( a=0; a< S->nseq; a++)for ( b=0; b<=(CL->S)->max_len; b++)cache[a][b]=b;

	seq_1_to_n=((fp=find_token_in_file (fname, NULL, "SEQ_1_TO_N"))!=NULL);
	vfclose (fp);
	
	fp=vfopen(fname,"r");
	
	if ( sn_list==NULL)sn_list=vcalloc (max_nseq, sizeof (int));
	
	fscanf ( fp, "%d\n", &list_nseq);
	for ( a=0; a<list_nseq; a++)
		{
		fscanf ( fp, "%s %d %s\n", name, &seq_len, sequence);
		line++;
	
		lower_string (sequence);
		
		if ((sn=name_is_in_list (name,S->name, S->nseq, 100))==-1)
			{
			fprintf (stderr, "\nERROR: Attempt to read sequence %s that is not in the sequence list [FATAL]", name);
			fprintf (stderr, "\nThe sequence name %s could contain special character that cause this problem, please rename it\n", name);
			
			exit (EXIT_FAILURE);
			}
		else
			{
			if (!strm (sequence,S->seq[sn]))
			       {
			       sprintf ( buf, "%s", sequence);
			       B=align_two_sequences (S->seq[sn],buf, strcpy ( mat,"idmat"), 0, 0, strcpy (dp_mode,"fasta_pair_wise"));
			       if ( S->len[sn]!=B->len_aln)
			          {
				  sprintf ( B->name[0], "%s", S->name[sn]);
				  sprintf ( B->name[1], "%s", name);
				  fprintf (stderr, "\nERROR: IMPOSSIBLE TO RECONCILIATE SEQ %s in file %s\nCHECK THAT TWO != SEQUENCES DO NOT HAVE THE SAME NAME\n",S->name[a], fname);
				  print_aln (B);
				  exit (EXIT_FAILURE);
			          }
			       else 
			          {
				      fprintf (stderr, "\nSequence Reconciliation [%s]", S->name[sn]);
				      for (d=0,e=0; e< B->len_aln; e++)
				          {
					      r=!is_gap(B->seq_al[1][e]);
					      d+=r;
					      if ( r)cache[sn][d]=e+1;
					  }
				  }
			       free_aln (B);
			       }
			sn_list[a]=sn;
			}
		
		}
	
	
	while ( (c=fgetc(fp))!=EOF)
		{
		
		ungetc(c, fp);
		if ( c=='#')
			{
			fscanf ( fp, "#%d %d\n", &s1, &s2);line++;
		

			/*Check If the sequence numbering is legal*/
			if ( seq_1_to_n){s1--; s2--;}
			
			if (s1<0 || s2 <0)
			  {
			    fprintf (stderr, "ERROR: Wrong Sequence Numbering in %s [FATAL:%s]\n",fname, PROGRAM); 
			    exit (EXIT_FAILURE);
			  }


			

			s1=sn_list[s1];
			s2=sn_list[s2];
			
			while (isdigit((c=fgetc(fp))))
				{
				for ( z=0; z<  CL->entry_len; z++)entry[z]=0;
				ungetc(c, fp);				
				n=0;
				entry[n++]=s1;
				entry[n++]=s2;
				while ( (c=fgetc(fp))!='\n')
					{
					if ( isspace (c));
					else 
						{
						ungetc(c, fp);
						fscanf ( fp, "%d", &entry[n]);
						n++;
						}
					
					if ( n>CL->entry_len)
						{
						fprintf ( stderr, "\nWARNING1:PARSING ERROR IN %s AT LINE %d: C=%c n=%d\n", fname,line, c,n);
						for ( e=2; e<LIST_N_FIELDS; e++)
							fprintf ( stderr, "%d ", entry[e]);
					
						exit (EXIT_FAILURE);
						}
					}
				if (c=='\n')line++;
				 
				if ( n<=CONS)entry[CONS]=1;
				
				  
				/*Check The legality of the entry*/
				if ( n>0 && n<3)
					{
					fprintf ( stderr, "\nWARNING2:PARSING ERROR IN %s AT LINE %d: C=%c\n", fname,line-1, c);
					for ( e=2; e<LIST_N_FIELDS; e++)
						fprintf ( stderr, "%d ",entry[e]);
					
					exit (EXIT_FAILURE);
					}
				
				
				if ( entry[R1]<=0 || entry[R1]>(CL->S)->len[s1])
				  {
				   fprintf ( stderr, "\nERROR: Seq1: %d, Seq2 %d, Res1 %d, Res2 %d\n", entry[SEQ1], entry[SEQ2], entry[R1], entry[R2]);
				   fprintf ( stderr, "\nERROR: Library %s, line %d, Field 1: Bad residue numbering (%d)[FATAL:%s]\n", fname, line-1,entry[R1], PROGRAM);
				   exit (EXIT_FAILURE);
				  }
				else if (entry[R2]<=0 || entry[R2]>(CL->S)->len[s2])
				  {
				    fprintf ( stderr, "\nERROR: Seq1: %d, Seq2 %d, Res1 %d, Res2 %d\n", entry[SEQ1], entry[SEQ2], entry[R1], entry[R2]);
				    fprintf ( stderr, "\nERROR: Library %s, line %d, Field 2: Bad residue numbering (%d)[FATAL:%s]\n", fname, line-1, entry[R2],PROGRAM);
				    exit (EXIT_FAILURE);
				  }
				fscanf ( fp, "\n");
				if ( (entry[SEQ1]>entry[SEQ2])|| (entry[SEQ1]==entry[SEQ2] && entry[R1]>entry[R2]))
				   {
				   SWAP(entry[SEQ1],entry[SEQ2]);
				   SWAP(entry[R1], entry[R2]);
				   }
				
				entry[R1]=cache[entry[SEQ1]][entry[R1]];
				entry[R2]=cache[entry[SEQ2]][entry[R2]];
				for ( z=0; z< CL->entry_len; z++)vwrite_clist( CL,CL->ne, z, entry[z]);
				CL->ne++;
				}
			 ungetc ( c, fp);
			 
			 }
		else if (c=='!' || c=='C')
			{
			while ((c=fgetc(fp))!='\n');
			line++;
			}
		else
			{
			fprintf ( stderr, "\nPARSING ERROR IN %s AT LINE %d: %c", fname,line,c); 	 
			exit (1);
			}
		}
	vfree(buf);
	vfree(sequence);
	vfclose (fp);	 
	return CL;
	} 
	
int read_seq_in_list ( char *fname,  int *nseq, char ***sequences, char ***seq_name)
	{
	int a,c;
	int seq_len, sn;

	FILE *fp;
	char name[1000];
	char *sequence;	
	static int max_nseq;
	static int *sn_list;
	int list_nseq;
	int lline;

	lline=measure_longest_line_in_file (fname);
	sequence=vcalloc (lline, sizeof (char)+1);
	fp=vfopen (fname, "r");
	while((c=fgetc(fp))!='#')if ( c=='\n')max_nseq++;
	vfclose (fp);
	
	
	
	if ( seq_name[0]==NULL)
		{
		seq_name[0]= declare_char (max_nseq,0);
		sequences[0]=declare_char (max_nseq,0);		
		}
	if ( sn_list==NULL)sn_list=vcalloc ( max_nseq, sizeof (int));
	else sn_list=vrealloc (sn_list, max_nseq*sizeof (int));

	fp=vfopen (fname,"r");
	
	
	if (fscanf ( fp, "%d\n", &list_nseq)!=1)return 0;
	

	

	for ( a=0; a<list_nseq; a++)
		{
		fscanf ( fp, "%s %d %s\n", name, &seq_len, sequence);
		lower_string (sequence);
		
		if ((sn=name_is_in_list (name, seq_name[0], nseq[0], 100))==-1)
			{
			seq_name[0][nseq[0]]=vcalloc (strlen (name)+1, sizeof (char));
			sprintf (seq_name[0][nseq[0]], "%s", name);
			sequences[0][nseq[0]]=vcalloc (strlen (sequence)+1, sizeof (char));
			sprintf (sequences[0][nseq[0]], "%s", sequence);
			sn_list[a]=nseq[0];
			nseq[0]++;
			}
		else
			{
			sn_list[a]=sn;
			}
		}
	vfclose (fp);
	free (sequence);
	return 1;
	}


/*********************************************************************/
/*                                                                   */
/*                         LIST OUTPUT                               */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	


FILE * save_constraint_list ( Constraint_list *CL,int start, int len, char *fname, FILE *fp,char *mode, Sequence *S)
        {
	int a, b;
	static int* translation;

	
	
	
	if ( fp==NULL)
	   {
	   if ( translation!=NULL)free(translation);
	   translation=vcalloc ( (CL->S)->nseq+1, sizeof (int));
	   for ( b=0,a=0; a< (CL->S)->nseq; a++)
	       {
	       if ( name_is_in_list((CL->S)->name[a],S->name,S->nseq, 100)==-1)
	          {
		  (CL->S)->len[a]=-1;
		  translation [a]=-1;
		  }
	       else 
	          {
		  translation[a]=b++;
		  }
	       }
	   
	   }
	if (strm2(mode, "lib","ascii"))
	   {
	   if ( fp==NULL)fp=vfopen ( fname, "w");
	   fp=save_list_header (fp,CL);
	   fp=save_constraint_list_ascii(fp, CL, 0, CL->ne, translation);
	   }
	else if (strm(mode, "binary"))
	   {
	   if ( fp==NULL)fp=vfopen ( fname, "wb");
	   fp=save_constraint_list_bin  (fp, CL, 0, CL->ne, translation);
	   }
	else
	    {
	    fprintf (stderr,"\nUNKOWN MODE FOR OUTPUT: %s [FATAL]\n",mode);
	    crash ("");
	    }
	return fp;
	}
	    	    
FILE * save_list_header ( FILE *OUT,Constraint_list *CL)
	{	
	int a;
	int nseq=0;
	
	for ( a=0; a<(CL->S)->nseq; a++)nseq+=((CL->S)->len[a]!=-1);

	fprintf ( OUT, "%d\n",nseq);
	for ( a=0; a<(CL->S)->nseq; a++)
		if ((CL->S)->len[a]!=-1) fprintf ( OUT, "%s %d %s\n", (CL->S)->name[a], (CL->S)->len[a],(CL->S)->seq[a]);
	
	return OUT;			
	}	   
FILE * save_constraint_list_ascii ( FILE *OUT,Constraint_list *CL, int start,int len, int *translation)
	{	
	int a, b, s1, s2;
        CLIST_TYPE x1, x2;

	if (len==start && CL->cpu!=-1)
	    {
	    fprintf (OUT, "CPU %d\n",get_time());
	    return OUT;
	    }
	else
	    {
	    
	    s1=translation[vread_clist(CL,start,SEQ1)];
	    s2=translation[vread_clist(CL,start,SEQ2)];
	   
	    
	    if ( s1!=-1 && s2!=-1)fprintf ( OUT, "#%d %d\n", s1+1, s2+1);
	    for ( a=start; a<(len+start); a++)
		   {
		   x1=translation[vread_clist(CL,a,SEQ1)];
		   x2=translation[vread_clist(CL,a,SEQ2)];
		   if ( x1==-1 || x2==-1);
		   else 
		       {
		       if ( x1!=s1 || x2!=s2)
			  {
			  s1=x1;
			  s2=x2;
			  fprintf ( OUT, "#%d %d\n", s1+1, s2+1);
			  }
		       for ( b=2; b<CL->entry_len; b++) fprintf ( OUT, "%5d ", vread_clist(CL, a, b));
		       fprintf (OUT, "\n");
		       }
		   }
	    }
	if ( CL->cpu)fprintf (OUT, "CPU %d\n",get_time());
	fprintf (OUT, "! SEQ_1_TO_N\n");
	return OUT;			
	}
FILE * save_constraint_list_bin ( FILE *OUT,Constraint_list *CL, int start,int len, int *translation)
	{	
	int a, b;

	CLIST_TYPE x1, x2;
	

	if (len==start && CL->cpu!=-1)
	    {
	    
	    return OUT;
	    }
	else
	    {
	    for ( a=start; a<(len+start); a++)
		   {
		   x1=translation[vread_clist(CL,a,SEQ1)];
		   x2=translation[vread_clist(CL,a,SEQ2)];
		   if ( x1==-1 || x2==-1);
		   else 
		       {
		       for ( b=2; b<CL->entry_len; b++)
		           {
			   x1=vread_clist(CL,a,b);
			   fwrite (&x1, CL->el_size, 1, OUT);
			   }
		       }
		   }
	    }
	return OUT;			
	}

/*********************************************************************/
/*                                                                   */
/*                         LIST CONVERTION                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/		

Constraint_list *aln_file2constraint_list (char *alname, Constraint_list *CL,char *weight_mode)
        {
	Alignment *A;
	

	A=declare_aln (CL->S);

	main_read_aln ( alname, A);
	CL=aln2constraint_list (A, CL, weight_mode);
	free_aln (A);
	return CL;
	}

Constraint_list *aln2constraint_list      (Alignment *A, Constraint_list *CL,char *weight_mode)
	{
	Constraint_list *CLB=NULL;
	int a, b, c,nres1, nres2;
	int **weight=NULL;
	int ***pw_weight=NULL;
	int ref_weight;
	int *col=NULL;
	int s1, s2;
	int do_swap;
	int fixed_nres1, fixed_nres2;
	int do_pdb=0;
	int pdb_weight;

	
	CLB=(Constraint_list *)A->CL;
	if ( !(strm ( weight_mode, "pdb")))
	  {
	    pw_weight=vcalloc ( A->nseq, sizeof (int**));
	    for ( a=0; a< A->nseq; a++)
	      {
		pw_weight[a]=vcalloc ( A->nseq, sizeof (int*));
		for ( b=0; b< A->nseq; b++)
		  {
		    pw_weight[a][b]=vcalloc ( A->len_aln, sizeof (int));
		  }
	      }
	  }
	
	if ( weight_mode==NULL || strcmp (weight_mode, "no")==0 || is_number (weight_mode))
	    {
	    if (is_number (weight_mode))ref_weight=atoi(weight_mode);
	    else ref_weight=1;
	    	    	    
	    for ( a=0; a< A->nseq; a++)
	      for ( b=0; b< A->nseq; b++)
		for ( c=0; c< A->len_aln; c++)
		  pw_weight[a][b][c]=ref_weight;
	    }
	else if ( strncmp ( weight_mode, "sim", 3)==0)
	    {
	     weight=get_sim_aln_array(A, weight_mode+3);
	     for ( a=0; a< A->nseq; a++)
	      for ( b=0; b< A->nseq; b++)
		for ( c=0; c< A->len_aln; c++)
		  pw_weight[a][b][c]=weight[a][b];
	     
	     free_int (weight, -1);
	    }
	else if ( strncmp (weight_mode, "winsim", 6)==0)
	    {
	      pw_weight=get_winsim_aln_array(A,weight_mode+6,pw_weight);
	    }
	
	else if (  strncmp ( weight_mode, "cdna", 4)==0)
	  {
	 
	   weight=get_sim_aln_array(A, weight_mode);
	   col=vcalloc ( A->len_aln+1, sizeof (int));
	   if (A->cdna_cache)
	       for ( a=0; a<=A->len_aln; a++)col[a]=A->cdna_cache[0][a];
	   else
	     for ( a=0; a<=A->len_aln; a++)col[a]=1;
	  
	   for ( a=0; a< A->nseq; a++)
	      for ( b=0; b< A->nseq; b++)
		for ( c=0; c< A->len_aln; c++)
		  pw_weight[a][b][c]=weight[a][b]*col[c];
	   
	    
	     free_int (weight, -1);
	     free (col);

	  }
	else if ( strncmp ( weight_mode, "col_sim", 7)==0)
	    {
	    col =get_aln_col_weight (A, weight_mode+7);
	    weight=get_sim_aln_array(A, weight_mode+7);
	    for ( a=0; a< A->nseq; a++)
	      for ( b=0; b< A->nseq; b++)
		for ( c=0; c< A->len_aln; c++)
		  pw_weight[a][b][c]=weight[a][b]*col[c];
	    
	    
	    free_int (weight, -1);
	    free (col);
	    }
	else if ( strm ( weight_mode, "pdb"))
	   {
	       if ( !CLB || !(CLB->T))
	          {
		      fprintf ( stderr, "\nCould not find the PDB structure: [FATAL]\n");
		      crash ("");
		  }
	       else do_pdb=1;
	   }
	else
	  {
	    fprintf ( stderr, "\nERROR: Weight Mode %s is unknown [FATAL:%s]", weight_mode, PROGRAM);
	    crash ("");
	  }
	      
	

	for ( a=0; a<A->nseq-1; a++)
		
		{
		for ( do_swap=0,b=a+1; b< A->nseq; b++)
			{	
			s1=name_is_in_list (A->name[a], (CL->S)->name, (CL->S)->nseq, 100);
			s2=name_is_in_list (A->name[b], (CL->S)->name, (CL->S)->nseq, 100);
			
			if ( s1>s2){do_swap=1;SWAP(s1, s2);}
			else do_swap=0;

			if ( s1!=-1 && s2!=-1)
			    {
			    for (nres1=A->order[a][1], nres2=A->order[b][1], c=0; c< A->len_aln; c++)
				{
				nres1+=!is_gap(A->seq_al[a][c]);
				nres2+=!is_gap(A->seq_al[b][c]);
				
				if ( do_pdb)
				  {
				     
				      pdb_weight=MAX(0,(CLB->evaluate_residue_pair)(CLB,0, nres1,1,nres2));
				  }


				if ( !is_gap(A->seq_al[a][c]) && !is_gap(A->seq_al[b][c]) && A->seq_al[b][c]!=UNDEFINED_RESIDUE && A->seq_al[a][c]!=UNDEFINED_RESIDUE && !(do_pdb && pdb_weight==0))
					{
					fixed_nres1=(!A->seq_cache)?nres1:A->seq_cache[(do_swap)?s2:s1][nres1];
					fixed_nres2=(!A->seq_cache)?nres2:A->seq_cache[(do_swap)?s1:s2][nres2];
					

					if ( do_swap){SWAP(fixed_nres1, fixed_nres2);}
					vwrite_clist (CL,CL->ne, SEQ1, s1);
					vwrite_clist (CL,CL->ne, SEQ2, s2);
					vwrite_clist (CL,CL->ne, R1, fixed_nres1);
					vwrite_clist (CL,CL->ne, R2, fixed_nres2);
					
					if (pw_weight)
					  {
					    vwrite_clist (CL,CL->ne, WE, pw_weight[a][b][c] );
					  }
					else if ( do_pdb) 
					  {
					    vwrite_clist (CL,CL->ne, WE,pdb_weight );
					  }
					vwrite_clist (CL,CL->ne, CONS,1);
					vwrite_clist (CL,CL->ne, MISC,0);
					if (do_swap){SWAP(fixed_nres1, fixed_nres2);}
				
					CL->ne++;
					}
					
				}
			    }
			}
		}
	
	if ( pw_weight)
	  {
	    for ( a=0; a< A->nseq; a++)
	      {
		for ( b=0; b< A->nseq; b++)
		  free (pw_weight[a][b]);
		free(pw_weight[a]);
	      }
	    free (pw_weight);
	  }

	return CL;
       } 
double **list2mat (Constraint_list *CLin,int s1,int s2, double *min, double *max)
        {
	double ** mat;
	int a, r1, r2;
	int min_def=0;
	Constraint_list *CL;
	static Sequence *S;

	
	int row, column;
	if ( S==NULL)S=declare_sequence ((CLin->S)->min_len, (CLin->S)->max_len,(CLin->S)->nseq); 	
	sprintf ( S->name[0], "%s",(CLin->S)->name[s1]);
	sprintf ( S->name[1],"%s",(CLin->S)->name[s2]);
        S->nseq=2;
	
        row   =(CLin->S)->len[s1];
	column=(CLin->S)->len[s2];
	
	if ( CLin->extend_jit)	
	    CL=extend_list(CLin,"mem",CLin->extend_clean_mode, CLin->extend_compact_mode, CLin->do_self, S);
	else
	    CL=CLin;


	min[0]=max[0];
	mat=declare_double ( row, column);

	for ( a=0; a<CL->ne; a++)
	    {
	    r1=vread_clist(CL,a,R1)-1;
	    r2=vread_clist(CL,a,R2)-1;
	    if ( vread_clist(CL,a,SEQ1)==s1 &&vread_clist(CL,a,SEQ2)==s2)
		{
		mat[r1][r2]=(double)vread_clist(CL,a,WE);
		if (min_def==0)
		   {
		   min_def=1;
		   min[0]=mat[r1][r2];
		   max[0]=mat[r1][r2];
		   }
		else
		   {
		   min[0]=(min[0]<mat[r1][r2])?min[0]:mat[r1][r2];
		   max[0]=(max[0]>mat[r1][r2])?max[0]:mat[r1][r2];
		   }
		}
	    else if (vread_clist(CL,a,SEQ2)==s1 &&vread_clist(CL,a,SEQ1)==s2) 
		{
		mat[r2][r1]=(double)vread_clist(CL,a,WE);
			if (min_def==0)
		   {
		   min_def=1;
		   min[0]=mat[r2][r1];
		   max[0]=mat[r2][r1];
		   }
		else
		   {
		   min[0]=(min[0]<mat[r2][r1])?min[0]:mat[r2][r1];
		   max[0]=(max[0]>mat[r2][r1])?max[0]:mat[r2][r1];
		   }
		}
	    }
	return mat;
	}

Constraint_list * constraint_list2bin_file(Constraint_list *clist)
        {
	int a,b;
	
	clist->fp=tmpfile();
	for ( a=0; a< clist->ne; a++)
	    for ( b=0; b<clist->entry_len; b++)
	        {
		fwrite (&clist->L[a][b],clist->el_size, 1,clist->fp);
		}
	return clist;
	}

FILE * bin_file2constraint_list ( Constraint_list *CL, FILE *fp, char *name)
        {
	int a, b, s1, s2;
	CLIST_TYPE *entry;
	
	if ( fp==NULL)fp=vfopen ( name, "w");
	entry=vcalloc ( CL->entry_len, CL->el_size);
	fprintf ( fp, "%d\n", (CL->S)->nseq);
	for ( a=0; a< (CL->S)->nseq; a++)fprintf (fp, "%s %d %s\n", (CL->S)->name[a], (CL->S)->len[a], (CL->S)->seq[a]);
	

	rewind ( CL->fp);
	fread(entry, CL->el_size, CL->entry_len, CL->fp);
	s1=entry[SEQ1];
	s2=entry[SEQ2];
	fprintf (fp, "#%d %d\n", s1, s2);
	for ( b=2; b< CL->entry_len; b++)fprintf (fp, "%5d ",entry[b]);
	fprintf (fp, "\n");
	for ( a=1; a< (CL->ne); a++)
	    {
	    fread(entry, CL->el_size, CL->entry_len, CL->fp);
	    if ( entry[SEQ1]!=s1 || entry[SEQ2]!=s2)
	       {
	       s1=entry[SEQ1];
	       s2=entry[SEQ2];
	       fprintf (fp, "#%d %d\n", s1, s2);
	       }
	    for ( b=2; b< CL->entry_len; b++)fprintf (fp, "%5d ",entry[b]);
	    fprintf (fp, "\n");
	    }
	fprintf (fp, "CPU %d\n",get_time());
	
	return fp;
	}
int **list2residue_total_weight ( Constraint_list *CL)
        {
	  /*Returns 
	    tot_weight[nseq][maxlen]
	    where each residue is associated with the total of its weights in CL
	    ####IMPORTANT
            
	          -the numbering of the residues  goes from 1 to L:
	          -the numbering of the sequences goes from 0 to N-1:
	  */

	int **tot_weight;
	int s1, s2, r1, r2, w, a;


	tot_weight=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	for ( a=0; a<CL->ne; a++)
	    {
	    r1=vread_clist(CL,a,R1)-1;
	    r2=vread_clist(CL,a,R2)-1;
	    s1=vread_clist(CL,a,SEQ1);
	    s2=vread_clist(CL,a,SEQ2);
	    w=vread_clist(CL,a,WE);
	    
	    tot_weight[s1][r1]+=w;
	    tot_weight[s2][r2]+=w;
	    }
	return tot_weight;
	}

int **list2residue_total_extended_weight ( Constraint_list *CL)
        {
	  /*Returns 
	    tot_extended_weight[nseq][maxlen]
	    where each residue is associated with the total of its weights in CL
	    ####IMPORTANT
            
	          -the numbering of the residues  goes from 1 to L:
	          -the numbering of the sequences goes from 0 to N-1:
	  */

	int **tot_extended_weight;
	int s1, s2, r1, r2, w;


	tot_extended_weight=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	
	for ( s1=0; s1< (CL->S)->nseq-1; s1++)
	  for ( s2=s1+1; s2< (CL->S)->nseq; s2++)
	      for (r1=1; r1<=(CL->S)->len[s1]; r1++)
		for (r2=1; r2<=(CL->S)->len[s2]; r2++)
		  {
		    w=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
		    tot_extended_weight[s1][r1]+=w;
		    tot_extended_weight[s2][r2]+=w;
		  }
	return tot_extended_weight;
	}
int **list2residue_partial_extended_weight ( Constraint_list *CL)
        {
	  /*Returns 
	    tot_extended_weight[nseq][maxlen]
	    where each residue is associated with the total of its weights in CL
	    ####IMPORTANT
            
	          -the numbering of the residues  goes from 1 to L:
	          -the numbering of the sequences goes from 0 to N-1:
	  */

	int **tot_extended_weight;
	int s1, s2, r1, r2, w1, w2, a;


	tot_extended_weight=declare_int ( (CL->S)->nseq, (CL->S)->max_len+1);
	for ( a=0; a<CL->ne; a++)
	    {
	    r1=vread_clist(CL,a,R1);
	    r2=vread_clist(CL,a,R2);
	    s1=vread_clist(CL,a,SEQ1);
	    s2=vread_clist(CL,a,SEQ2);
	    w1=(CL->evaluate_residue_pair)( CL, s1, r1, s2, r2);
	    w2=(CL->evaluate_residue_pair)( CL, s2, r2, s1, r1);
	    if ( w1!=w2)fprintf ( stderr, "*");

	    tot_extended_weight[s1][r1]+=w1;
	    tot_extended_weight[s2][r2]+=w2;
	    }
	return tot_extended_weight;
	}
  

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              clean functions                                            */
/*                                                                                         */
/*                                                                                         */
/*                                                                                         */
/*******************************************************************************************/
Constraint_list *clean ( char *clean_mode,Constraint_list *CL,int start, int len)
	{
	
	if ( strm ( clean_mode, "shadow"))   CL=clean_shadow (CL,start,len);	
	else if ( strm5( clean_mode, "","NO","no","No","default"));
 	else	fprintf ( CL->local_stderr, "\nWARNING: The %s CLEANING MODE DOES NOT EXIST\n", clean_mode);
	
	return CL;
	}


Constraint_list * clean_shadow ( Constraint_list *CL, int start, int len)
	{
	int s1, s2, r1, a, b, end;
	int max, min;
		
	s1=vread_clist (CL, start, SEQ1);
	s2=vread_clist (CL, start, SEQ2);
	r1=vread_clist (CL, start, R1);
	
	
	for ( a=start; a<(start+len);)
		{
		
		max=min=vread_clist (CL, a, WE);
		while ( a<CL->ne && vread_clist (CL, a, SEQ1)==s1 && vread_clist (CL, a, SEQ2)==s2 && vread_clist (CL, a, R1)==r1)
			{
			max=(vread_clist (CL, a, WE)>max)?vread_clist (CL, a, WE):max;
			min=(vread_clist (CL, a, WE)<min)?vread_clist (CL, a, WE):min;
			a++;
			}
		end=a;
		
		if ((end-start)>1)
			{
			for ( b=start; b<end; b++)
				if ( vread_clist (CL, b, WE)<max)vwrite_clist (CL, b, SEQ1,-1);
			}
		start=end;
		if ( start<CL->ne)
			{
			s1=vread_clist (CL, start, SEQ1);
			s2=vread_clist (CL, start, SEQ2);
			r1=vread_clist (CL, start, R1);
			}
		}
	CL=sort_constraint_list_inv (CL, start, (CL->ne-start));
	CL=sort_constraint_list (CL,start,(CL->ne-start)  );	

	return CL;
	}
/*********************************************************************/
/*                                                                   */
/*                         LIST FUNCTIONS                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
Constraint_list *modify_weight( Constraint_list *CL,int start, int end,  char *modify_mode)
        {
	int a;
	CLIST_TYPE x;
	
	if ( strm(modify_mode, "default"))return CL;
	for ( a=start; a<end; a++)
	    {
	    x=vread_clist(CL, a, WE);
	    
	    if (strm2 (modify_mode,"sc_eq_cons", "we_eq_cons"))
	        if(x!=UNDEFINED)
		    vwrite_clist(CL, a, WE,  vread_clist(CL, a, CONS));
	    
	    if (strm2(modify_mode,"sc_eq_wePcons","sc_eq_consPwe"))
	        if(x!=UNDEFINED)
		     vwrite_clist(CL, a, WE, vread_clist(CL, a, CONS)*x);		
	    }
	return CL;
	}
	
Constraint_list *compact_list (Constraint_list *CL, int start, int len, char *compact_mode)
	{
	int a;
	int r1, r2, rr1, rr2, s1, rs1, s2, rs2, ra;
	CLIST_TYPE x;
	int debug_compact=0;

	if (debug_compact)fprintf ( stderr, "\n[In: %d %s", CL->ne, compact_mode);

	if ( len==0  || strm3(compact_mode, "no", "No", "NO"))return CL;
	else if ( strm2(compact_mode,"mirror","mirror_sum"));
	else if ( strm4(compact_mode, "default","shrink","shrink_best","shrink_worst"))
	        {
		
		for ( a=start; a<(start+len) ; a++)
		    {
		    
		    if ( vread_clist(CL, a, SEQ1)> vread_clist(CL, a, SEQ2) ||\
		       ( vread_clist(CL, a, SEQ1)==vread_clist(CL, a, SEQ2) &&\
			 vread_clist(CL, a, R1)  > vread_clist(CL, a, R2)     ))
			
			
		        {
			s1=vread_clist(CL, a, SEQ1);
			s2=vread_clist(CL, a, SEQ2);
			r1=vread_clist(CL, a, R1);
			r2=vread_clist(CL, a, R2);
			vwrite_clist(CL, a, SEQ1,s2);
			vwrite_clist(CL, a, SEQ2,s1);
			vwrite_clist(CL, a, R1,r2);
			vwrite_clist(CL, a, R2,r1);
			}
		    }
		}
	

	sort_constraint_list ( CL, start, len);
	
	rs1=vread_clist(CL, start, SEQ1);	
	rs2=vread_clist(CL, start, SEQ2);
	rr1=vread_clist(CL, start, R1);
	rr2=vread_clist(CL, start, R2);
	ra=start;

	

	if ( (rs1==rs2) && (rr1==rr2))vwrite_clist(CL, start, SEQ1,-1);		
	for ( a=start+1; a<(start+len); a++)
		{
		s1=vread_clist(CL, a, SEQ1);
		s2=vread_clist(CL, a, SEQ2);
		r1=vread_clist(CL, a, R1);
		r2=vread_clist(CL, a, R2);
		
		if ( (s1==s2) && (r1==r2))vwrite_clist(CL, a, SEQ1, -1);
		else if ( s1==rs1 && s2==rs2 && r1==rr1 && r2==rr2)
			{
			x=vread_clist(CL, ra, WE);
			if (strm ( compact_mode, "shrink"));
			else if ( strm2 ( compact_mode, "default","mirror_sum"))    
			    vwrite_clist(CL, ra, WE, vread_clist(CL, a, WE)+x);
			else if (strm2 ( compact_mode,"best", "shrink_best"))
			    vwrite_clist(CL, ra, WE,MAX(vread_clist(CL, a, WE), vread_clist(CL, a, WE)));
			else if (strm2 ( compact_mode, "worst","shrink_worst"))
			    vwrite_clist(CL, ra, WE,MIN(vread_clist(CL, a, WE), vread_clist(CL, a, WE)));
		
			if (  strm(compact_mode, "shrink"));
			else
			    {
			    vwrite_clist(CL, ra, CONS, vread_clist(CL, ra, CONS)+ vread_clist(CL, a, CONS));
			    vwrite_clist(CL, ra, MISC, vread_clist(CL, ra, MISC)+ vread_clist(CL, a, MISC));
			    }
			vwrite_clist(CL,a, SEQ1, -1);
			
			}
		else
			{
			rs1=s1;
			rs2=s2;
			rr1=r1;
			rr2=r2;
			ra=a;
			}
		}
	
	
	sort_constraint_list_inv(CL,0,CL->ne);
	
	sort_constraint_list    (CL,0,CL->ne);
	
	
	if ( strm3 (compact_mode, "consPwe", "wePcons","cons"))
		{
		for ( a=start; a<(start+len); a++)
			{
			if ( strm2(compact_mode,"consPwe", "wePcons"))
			    vwrite_clist(CL, a, WE, vread_clist(CL,a,WE)* vread_clist(CL,a,CONS));
			else if  (strm (compact_mode, "cons"))
			     vwrite_clist(CL, a, WE, vread_clist(CL,a,CONS));
			}
		}
	if (debug_compact)fprintf ( stderr, "....OUT: %d]\n", CL->ne);	
	return CL;
	}


Constraint_list *rescale_list_simple (Constraint_list *CL,int start, int len,int new_min, int new_max)
	{
	int a, min, max;
	double x;
	/*Rescales between 0 and max2
	  Any value above max1 is set to max1 first and then to max2
	*/



	min=max=vread_clist ( CL,start, WE);
	
	for ( a=start; a<(start+len); a++)
		{
		x=(double)vread_clist ( CL,a, WE);
		if ( x>max)max=(int)x;
		if ( x< min)min=(int)x;
		}
	
	fprintf ( CL->local_stderr, "\n[%d-%d]=>[%d-%d]", min, max, new_min, new_max);

	for ( a=start; a<(start+len); a++)
	    {
	      
	    x=vread_clist(CL,a, WE);
	   
	    if ((max-min)==0)x=100;
	    else x=(((x-min)/(max-min))*new_max)+new_min;
	    
	    vwrite_clist(CL, a, WE,(CLIST_TYPE) x);
	    }
	return CL;
	}
Constraint_list *rescale_list (Constraint_list *CL,int start, int len,int max1, int max2)
	{
	int a, min_val, max_val;
	CLIST_TYPE x;
	
	/*Rescales between 0 and max2
	  Any value above max1 is set to max1 first and then to max2
	*/


	min_val=0;
	max_val=max1;

	
	
	for ( a=start; a<len; a++)
		{
		if (vread_clist ( CL,a, WE)>max1)vwrite_clist(CL, a, WE, max1);
		}
	
	for ( a=start; a<len; a++)
	    {
	    x=vread_clist(CL,a, WE);
	    vwrite_clist(CL, a, WE, (((x-min_val)*max2)/(max_val-min_val)));
	    }
	return CL;
	}


Constraint_list* filter_list (Constraint_list *CL, int start, int len,int T)
	{
	int a;
	int field;
	
	if (T==0)return CL;
	
	field=WE;
	if (T>0)field=WE;
	else if ( T<0)
		{
		field=CONS;
		T=-T;
		}
	
	
	for ( a=start; a<len; a++)
	    if (vread_clist(CL, a, field)<=T)vwrite_clist(CL,a,SEQ1,-1);
	
	CL=sort_constraint_list_inv (CL, 0, CL->ne);
	CL=sort_constraint_list     (CL, 0, CL->ne);
	return CL;
	}





Constraint_list *undefine_list (Constraint_list *CL)
      {
      int a, b;
      int undefined_flag;
      
      for ( a=0;a<CL->ne; a++)
          {
	  for ( b=0, undefined_flag=0; b< LIST_N_FIELDS; b++)
	      {

	      if ( vread_clist (CL, a, b)==UNDEFINED)undefined_flag=1;	  
	      if ( undefined_flag)
	          { 
		  for ( b=0; b< LIST_N_FIELDS; b++)
		      if ( b!=SEQ1 && b!=SEQ2 && b!=R1 && b!=R2)
			  {
			  vwrite_clist(CL, a, b, UNDEFINED);
			  }
		  }
	      }
	  }
      return CL;
      }

/*********************************************************************/
/*                                                                   */
/*                        DEBUG CONSTRAINT_LIST                       */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	  
void check_seq_pair_in_list(Constraint_list *CLin,int seq1, int seq2)
      {
      int a, s1, s2, r1, r2;
	   
      for ( a=0; a< CLin->ne; a++)
          {
	  s1=vread_clist(CLin,a,SEQ1);
	  s2=vread_clist(CLin,a,SEQ2);
	  if ( s1==seq1 && s2==seq2)
	     {
	     r1=vread_clist(CLin,a,R1);
	     r2=vread_clist(CLin,a,R2);
	     fprintf ( stderr, "\n[%d][%d %d] [%d %d]",a,s1, r1, s2, r2);
	     }
	  }
      }

void print_CL_mem(Constraint_list *CL, char *function)
     {
      fprintf ( stderr,"%s\n", function);
     if ( CL->fp==NULL && CL->L==NULL) fprintf ( stderr, "\n\tNOTHING");
     if ( CL->fp)fprintf ( stderr, "\n\tFILE SET");
     if ( CL->L)fprintf ( stderr, "\n\tMEM SET\n");
     }

int constraint_list_is_sorted ( Constraint_list *CL)
     {
     int a,b, x1, x2;
     for ( a=0; a< CL->ne-1; a++)
         {
	 for ( b=0; b< CL->entry_len; b++)
	     {
	     x1=vread_clist( CL, a, b);
	     x2=vread_clist( CL, a+1,b);
	     if ( x1<x2)break;
	     else if ( x1>x2)
		 {
		 fprintf ( stderr, "\n[%d][%d]=>%d\n[%d][%d]=>%d\n\n",a, b, x1, a+1, b, x2);
		 return 0;
		 }
	     
	     }
	 }
     return 1;
     }
/*********************************************************************/
/*                                                                   */
/*                        WEIGHT CONSTRAINT_LIST                     */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

Constraint_list *weight_constraint_list(Constraint_list * CL, char *seq_weight)

    {
      Weights *W;
      
      if ( CL->ne==0)return CL;
      else if ( strm(seq_weight, "t_coffee"))      W=compute_t_coffee_weight(CL);
      else if (check_file_exists (seq_weight))W=read_seq_weight ((CL->S)->name, (CL->S)->nseq, seq_weight);
      else return CL;
      
      CL=re_weight_constraint_list (CL,W);
      CL->W=W;
      
      return CL;
      
      
    }
      


Weights* compute_t_coffee_weight(Constraint_list * CL)
    {
      int a, b, c;
      float p, d;
      Weights *W;
      Alignment *B;
      int nseq;
      

      


      if (!CL->L)return NULL;
      
      nseq=(CL->S)->nseq;
      W=declare_weights(nseq);
      sprintf ( W->mode, "t_coffee");
      for ( a=0; a< nseq; a++)
	 {
	   sprintf ( W->seq_name[a], "%s", (CL->S)->name[a]);
	   W->SEQ_W[a]=1;
	 }
            
      for (a=0; a< (CL->S)->nseq; a++)
	for ( b=a; b< (CL->S)->nseq; b++)
	  {
	    if ( b==a){d=1;}
	    else
	      {

		B=align_two_sequences ((CL->S)->seq[a], (CL->S)->seq[b],"pam250mt", -10, -1, "fasta_pair_wise");
		for (d=0,c=0; c<B->len_aln; c++)d+=(B->seq_al[0][c]==B->seq_al[1][c]);
		d=d/B->len_aln;
		free_aln(B);
	      }
	    p=pow(d,3);
	    
	    W->SEQ_W[a]+=p;
	    W->SEQ_W[b]+=p;
	    
	  }
      for ( p=0,b=0; b< (CL->S)->nseq; b++)
	   {
	     W->SEQ_W[b]=2/W->SEQ_W[b];
	     p+=W->SEQ_W[b];
	   }
      for ( b=0; b< (CL->S)->nseq; b++)
	{
	  W->SEQ_W[b]=W->SEQ_W[b]*((float)W->nseq/p);
	}
      
      return W;
    }
      
Constraint_list *re_weight_constraint_list(Constraint_list * CL,Weights *W)
    {
      int a;
      float w;
      float *weight;
      int sA, sB;



      weight=W->SEQ_W;

      if (!CL->L)return CL;
 
      
      
      for ( a=0; a< CL->ne; a++)
	   {
	     sA=CL->L[a][SEQ1];
	     sB=CL->L[a][SEQ2];
	     
	     w=MIN(weight[sA], weight[sB]);
	     
	     CL->L[a][WE]*=w;
	   } 
	 CL=evaluate_constraint_list_reference (CL);
	 return CL;
    }
            
