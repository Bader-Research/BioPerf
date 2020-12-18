#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h" 
#include "define_header.h"



#define FATAL fatal:seq_reformat
int SeqGCGCheckSum(char *seq, int len);

/**************************************************************************************************/
/*****************************    FORMAT GUESSING     ******************************************/
/**************************************************************************************************/
Sequence_data_struc *read_data_structure ( char *in_format, char *in_file,	Action_data_struc  *RAD) 
         
        {
	Sequence_data_struc *D;
	char **seq_name=NULL, **sequences=NULL;
	int nseq=0;

	
	D=vcalloc ( 1, sizeof (Sequence_data_struc));
	
	if (!in_file[0])return NULL; 
	if (!in_format[0])
	  {
	    in_format=identify_seq_format(in_file);
	  }
	if ( !in_format[0] && is_lib(in_file))sprintf ( in_format, "tc_lib");
	
	if (!in_format[0])return NULL;
	

	D->A=declare_Alignment(NULL);		      
	if ( RAD->keep_case)(D->A)->residue_case=2;
	
	D->rm_gap=RAD->rm_gap;
	sprintf ( D->format, "%s", in_format);
	sprintf ( D->file, "%s", in_file);
	    
	if ( strm2(in_format,"saga_aln","clustal_aln"))
		{
		read_aln (in_file, D->A);
		D->S=aln2seq(D->A);
		
		}
	else if ( strm( in_format,"number_aln"))
		{		
		read_number_aln (in_file, D->A);
		D->S=aln2seq(D->A);
		}
	
	else if ( strm( in_format,"gotoh_aln"))
		{		
		read_gotoh_aln (in_file, D->A);
		D->S=aln2seq(D->A);
		}
	
	else if ( strm ( in_format, "msf_aln"))
		{
		read_msf_aln (in_file, D->A);
		D->S=aln2seq(D->A);
		}
	else if ( strm ( in_format, "amps_aln"))
		{
		read_amps_aln (in_file, D->A);
		D->S=aln2seq(D->A);
		}
	else if ( strm (in_format, "number_fasta"))
		{		
		D->S=get_fasta_sequence_num (in_file, NULL);
		(D->S)->contains_gap=0;
		D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}	
	else if ( strm2 (in_format, "fasta_aln", "fasta_seq"))
		{
		
		D->S=get_fasta_sequence (in_file, NULL);
		if ( strcmp (in_format, "fasta_aln")==0)(D->S)->contains_gap=0;
		D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}	
	else if ( strm (in_format, "pdb") || strm (in_format, "pdb_struc"))
		{
		    D->S=get_pdb_sequence (in_file);
		    D->A=seq2aln(D->S, D->A,RAD->rm_gap);
		}
	else if ( strm2(in_format, "pir_seq", "pir_aln"))
		{
		D->S=get_pir_sequence ( in_file,NULL );
		seq2aln(D->S, D->A, RAD->rm_gap);
		}
        else if ( strm(in_format, "gor_seq") )
		{
		D->S=get_gor_sequence ( in_file,NULL );
		seq2aln(D->S, D->A, RAD->rm_gap);
		}
	else if ( strm2 ( in_format, "dali_aln", "dali_seq"))
		{
		D->S=get_sequence_dali ( in_file);
		seq2aln(D->S, D->A, RAD->rm_gap);
		}
	else if ( strm (in_format, "barton_list_tc"))
		{
		get_barton_list_tc_seq ( in_file);
		}
	else if ( strm (in_format, "amps_sd_scores"))
		{
		D->W=get_amps_sd_scores ( in_file);
		}
	
	else if ( strm ( in_format, "pima_aln"))
		{
		D->S=get_pima_sequence ( in_file);
		seq2aln (D->S, D->A, RAD->rm_gap);
		}
	else if ( strm( in_format, "gor_struc"))
	        {
		D->S=get_struc_gor ( in_file);
		seq2aln(D->S, D->A, RAD->rm_gap);
		}
	else if ( strm( in_format, "dialign_aln"))
		{
		D->S=get_dialign_sequence ( in_file);
		seq2aln (D->S, D->A, RAD->rm_gap);
		}
	else if ( strm( in_format, "tc_lib") ||  strm( in_format, "mocca_lib") ||  strm( in_format, "lib"))
	        {
		  read_seq_in_list (in_file,&nseq,&sequences,&seq_name); 
		  D->S=fill_sequence_struc ( nseq, sequences, seq_name);
		  seq2aln (D->S, D->A, RAD->rm_gap);
		  free_char (sequences,-1);
		  free_char (seq_name, -1);
		}
	else if ( strm( in_format,"swissprot_seq"))
	        {
		  D->S=get_swissprot_sequence ( in_file,NULL);
		  seq2aln (D->S, D->A, RAD->rm_gap);
		}
	else
	        {
		
		return NULL; 
		}
	
	return D;
	}
Sequence  * main_read_seq ( char *name)
       {
       char *format=NULL;
       Sequence *S=NULL;
       Alignment *A=NULL;
       int a, b;

       
       format=identify_seq_format (name);
      
      

       if (format &&strm(format, "fasta_seq"))        S= get_fasta_sequence ( name, NULL);
       else if (format &&strm(format, "pir_seq"))     S= get_pir_sequence ( name, NULL);
       else if (format &&strm(format,"swissprot_seq"))S= get_swissprot_sequence (name, NULL); 
       else if (format && strstr (format, "aln")) 
	 {
	   A=main_read_aln ( name, NULL);
	   S=aln2seq(A);
	   ungap_seq(S);
	   free_aln(A);
	 }
       else
	  {
	  /*Use The ClustalW routine*/
	    S=read_sequences (name);
	  }
       for ( b=0,a=0; a< S->nseq; a++)
	 {
	   b+=( S->seq[a][0]!='\0');
	 }
       S->nseq=b;
       
       for ( a=0; a<S->nseq; a++)sprintf ( S->file[a], "%s", name);
       vfree(format);
       
       ungap_seq(S);
       
       return S;

       }
Alignment * main_read_aln ( char *name, Alignment *A)
       {
       int a;

       char *format=NULL;
       Sequence *S=NULL;
       Sequence *IN_SEQ;


      
       if (!A)A=declare_aln(NULL);
       format=identify_seq_format (name);
       IN_SEQ=A->S;
       
              
       if      (format && strm(format, "saga_aln" ))read_aln ( name, A);
       else if (format &&strm(format, "msf_aln"  ))read_msf_aln ( name, A);
       else if (format &&(strm(format, "fasta_aln")))
                {
		S=get_fasta_sequence ( name, NULL);
		
		S->contains_gap=0;
		seq2aln (S, A, 0);		
		}
       else if (format &&strm(format, "pir_aln"))
                {
		S=get_pir_sequence ( name, NULL);
		S->contains_gap=0;
		seq2aln (S, A, 0);
		} 
       else if (format && strm(format, "fasta_seq") && A)
	   {
	   S=get_fasta_sequence ( name, NULL);
	   
	   for ( a=1; a<S->nseq; a++)if ( strlen (S->seq[a-1])!=strlen (S->seq[a])){free_sequence (S, S->nseq); free_aln (A); return NULL;}
	   S->contains_gap=0;
	   seq2aln (S, A, 0);
	   }
       else if (format && strm(format, "pir_seq") && A)
	   {
	   S=get_pir_sequence ( name, NULL);
	  
	   for ( a=1; a<S->nseq; a++)if ( strlen (S->seq[a-1])!=strlen (S->seq[a])){free_sequence (S, S->nseq); free_aln (A); return NULL;}
	   S->contains_gap=0;
	   seq2aln (S, A, 0);
	   }
       else
          {
	      free_aln(A);
	      vfree(format);
	      return NULL;	  
	  }
       if ( check_list_for_dup( A->name, A->nseq))
          {
	      fprintf ( stderr, "\nERROR: %s is duplicated in file %s[FATAL]\n", check_list_for_dup( A->name, A->nseq), A->file[0]);
	      exit (1);
	  }
       if (IN_SEQ)A->S=IN_SEQ;
       else if (!A->S){A->S=aln2seq(A);}
       
       A->S=ungap_seq(A->S);
       A=fix_aln_seq(A, A->S);     
       compress_aln (A);
       for ( a=0; a< A->nseq; a++) sprintf ( A->file[a], "%s", name);
       vfree(format);
       return A;
       }

      
char * identify_aln_format ( char *file)
       {
	/*This function identify known sequence and alignmnent formats*/
	 return identify_seq_format (file);
       }
char * identify_seq_format ( char *file)
       {
       char *format=NULL;
       /*This function identify known sequence and alignmnent formats*/
       
       if ( format==NULL)format=vcalloc ( 100, sizeof (char));
       else format[0]='\0';
       
       
       if ( !check_file_exists(file));
       else if ( format_is_msf      (file))sprintf ( format, "msf_aln");
       else if ( format_is_fasta_seq(file))sprintf ( format, "fasta_seq");
       else if ( format_is_fasta_aln(file))sprintf ( format, "fasta_aln");       
       else if ( format_is_pir_aln  (file))sprintf ( format, "pir_aln");
       else if ( format_is_pir_seq  (file))sprintf ( format, "pir_seq");
       else if ( format_is_oligo    (file))sprintf ( format, "oligo_aln");
       else if ( format_is_swissprot     (file))sprintf ( format, "swissprot_seq");
       else if ( format_is_saga     (file))sprintf ( format, "saga_aln");
       else if ( is_pdb(file))sprintf ( format, "pdb_struc");  
      
      

       return format;
       }
char **identify_list_format ( char **list, int n)
       {
	   int a;
	   char *name;
	   char *string;
	   char mode;
	
	   declare_name (name);
	   for ( a=0; a< n; a++)
	       {
		
	       sprintf (name, "%s", list[a]);
	       string=list[a];
	       if ((mode=identify_format ( &string))!='?')
		   {
		       sprintf ( name, "%s", string);
		       sprintf ( list[a], "%c%s", mode,name);
		   }
	       else
	           {
		       fprintf ( stderr, "\nERROR: %s not recognised [FATAL]", name);
		   }
	     
	       }
	
	   vfree(name);
	   return list;
       }
	       
	   
char identify_format (char **fname)
       {
	   char mode='?';
	   mode=fname[0][0];
	   
	   if ((is_in_set (mode, "ALMSPR") && check_file_exists(fname[0]+1)) ||(mode=='X' && is_matrix ( fname[0]+1)) ||(mode=='M' && is_method(fname[0]+1)) )
	     {
	       
	       fname[0]++;
	     }
	   else
	       {
		   
		      if      (is_method(fname[0]))mode='M';
                      else if (is_lib(fname[0]))mode='L';
		      else if (is_pdb(fname[0]))mode='P';
		      else if (is_seq(fname[0]))mode='S';
		      else if (is_aln(fname[0]))mode='A';
		      else if (is_matrix(fname[0]))mode='X';
		      
		      else mode='?';
		  }
	   return mode;
       }

int is_pdb ( char *name)
       {
       FILE *fp;


       

       if ((fp=find_token_in_file (name, NULL, "HEADER"))!=NULL)
           {vfclose (fp);

	   }
	else return 0;
	

        if ((fp=find_token_in_file (name, NULL, "ATOM")))
           {vfclose (fp);
	   
	   }
	else return 0;
      
	if ((fp=find_token_in_file (name, NULL, "SEQRES")))
           {vfclose (fp);
	   
	   }
	else return 0;
       
	return 1;	
       }
int is_seq ( char *name)
       {
	 char *format;
	 if ( !check_file_exists(name))return 0;
	 
	 format= identify_seq_format(name);
	 if(!format || format[0]=='\0'){vfree (format);return 0;}
	 else if (strstr(format, "seq")){vfree (format);return 1;}
	 else return 0;
       }
int is_aln ( char *name)
       {
       char *format;	 
       if ( !check_file_exists       (name))return 0;   
	
       format= identify_seq_format(name);
       if ( !format || format[0]=='\0'){vfree (format);return 0;}
       else if (strstr(format, "aln")){vfree (format); return 1;}
       else return 0;
       }   

int is_matrix (char *name)
       {
       int **m;
       
       if ((m=read_matrice (name))!=NULL){free_int (m, -1); return 1;}
       return 0;       	      
       }
int is_lib (char *name)
       {
       int **list;
       int n_blocks;
       int a, b;
       FILE *fp;

       if ( !check_file_exists(name))return 0;
       list=get_file_block_pattern (name,&n_blocks, 100);
       
       
       if ( n_blocks>1) {free_int (list, -1);return 0;}
       if ( list[0][1]!=1) {free_int (list, -1);return 0;}
       fp=vfopen ( name, "r");
       if (!fscanf ( fp, "%d", &a)){vfclose(fp);free_int (list, -1);return 0;}
       vfclose ( fp);
       for ( b=0; b<a; b++)
           {
	       if (list[0][b+2]!=3){free_int (list, -1);return 0;}
	   }
       

       if ( list[0][b+2]!=2){free_int (list, -1);return 0;}
       for ( b=a+3; b<list[0][0] && b< 100; b++)
	   if ( !(list[0][b]==LIST_N_FIELDS-2 || list[0][b]==2)){free_int(list, -1);return 0;}
       free_int (list, -1);
       return 1;
       }
int is_method ( char *file)
    {
	char new_file[200];
	
	sprintf ( new_file, "%s", file);
	if ( is_in_pre_set_method_list(new_file)) 
	    {
		remove ( new_file);
		return 1;
	    }
	else

	    return 0;
    }

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              SEQUENCE FORMAT IDENTIFIERS                                */
/*                                                                                         */
/***************************************************************************************** */

int format_is_oligo(char *file)
    { 
    char buf[1000];
    FILE *fp;
    
    fp=vfopen ( file, "r");
    fscanf (fp , "%s", buf);
    vfclose ( fp);

    if ( strm (buf, "ALPHABET"))return 1;
    return 0;
    }
int format_is_msf ( char *file)
    {
    char buf[1000];
    FILE *fp;

   
    
    if ( (fp=find_token_in_file_nlines (file,NULL,"MSF:", 30))!=NULL){vfclose (fp);return 1;}
    else
      {
       return 0;
      }
    
    fp=vfopen ( file, "r");
    fscanf (fp , "%s", buf);
    vfclose ( fp);

    if ( strm (buf, "MSF:"))return 1;
    return 0;
    }


int format_is_fasta_aln ( char *file)

    {
      if ( format_is_fasta(file) && !format_is_fasta_seq(file))return 1;
      else return 0;
    }

  
int format_is_fasta_seq ( char *file)
    {
      int a, l1, l2,l;
      Sequence *S;


    if ( format_is_fasta (file))
      {	
	S=get_fasta_sequence (file, NULL);
	l=strlen ( S->seq[0]);
	for ( a=0; a< S->nseq; a++)if(strlen(S->seq[a])!=l){free_sequence (S, S->nseq);return 1;}
	for ( a=0; a< S->nseq; a++)
	  {
	    l1=strlen ( S->seq[a]);
	    ungap (S->seq[a]);
	    l2=strlen ( S->seq[a]);
	    if ( l1!=l2)
	      {
		free_sequence (S, S->nseq);
		return 0;
	      }
	  }
	free_sequence (S, S->nseq);
	return 1;
      }
    else
      {
	return 0;
      }
    }

int format_is_fasta ( char *file)
    {
    int c;
    FILE *fp;
    
    
    fp=vfopen ( file, "r");
    while (isspace(c=fgetc(fp)));
    
    if ( c=='>')
        {	
	while ( (c=fgetc(fp))!='\n');
	while ( (c=fgetc(fp))!='>' && c!='*' && c!=EOF);
	vfclose (fp);
	if ( c!='*')return 1;
	}
    vfclose (fp);
    return 0;
    }
int format_is_pir_aln ( char *file)

    {
      if ( format_is_pir(file) && !format_is_pir_seq(file))return 1;
      else return 0;
    }

int format_is_pir_seq ( char *file)
    {
      int a, l1, l2;
      Sequence *S;


    if ( format_is_fasta (file))
      {
	S=get_pir_sequence (file, NULL);
	for ( a=0; a< S->nseq; a++)
	  {
	    l1=strlen ( S->seq[a]);
	    ungap (S->seq[a]);
	    l2=strlen ( S->seq[a]);
	    if ( l1!=l2)
	      {
		free_sequence (S, S->nseq);
		return 0;
	      }
	  }
	return 1;
      }
    else
      {
	return 0;
      }
    }
    
int format_is_pir ( char *file)
    {

    FILE *fp;
    int c;
    
    fp=vfopen ( file, "r");
    while (isspace(c=fgetc(fp)));
    
    if ( c=='>')
        {
	while ( (c=fgetc(fp))!='\n');
	while ( (c=fgetc(fp))!='>' && c!='*' && c!=EOF);
	vfclose (fp);
	if ( c=='*')return 1;
	}
    vfclose ( fp);
    return 0;
    }      

int format_is_saga ( char *file) 
    {
    FILE *fp;
    int **list;
    int n_blocks;
    int n_seq;
    int a, b;
    
    if ( (fp=find_token_in_file (file, NULL, "SAGA"))){vfclose (fp); return 1;}
    else if  ((fp=find_token_in_file (file, NULL, "CLUSTAL"))){vfclose (fp); return 1;}
    else if  ((fp=find_token_in_file (file, NULL, "T-COFFEE"))){vfclose (fp); return 1;}
    else if  ((fp=find_token_in_file (file, NULL, "SAGA_FORMAT"))){vfclose (fp); return 1;}
    else if  ((fp=find_token_in_file (file, NULL, "GARP"))){vfclose (fp); return 1;}
    else if  ((fp=find_token_in_file (file, NULL, "INTERLEAVED"))){vfclose (fp); return 1;}
    
      else 
       {
	   list=get_file_block_pattern (file,&n_blocks,100); 
	   if (n_blocks<=2){free_int (list, -1);return 0;}
	   else 
	       {		  
	       n_seq=list[1][0];
	       for ( a=1; a< n_blocks-1; a++)
	           {
		       if ( list[a][0]!=n_seq){free_int (list, -1);return 0;}
		       else
		       {
			   for ( b=1; b<=list[a][0]; b++)
			       if ( list[a][b]!=2){free_int (list, -1);return 0;}
		       }
		   }
	       }
	   return 1;
       }
    
    return 0;
    }


int format_is_swissprot (char *name)
    {
      FILE *fp;
      
      if ( !check_file_exists(name))return 0;
	 
	 
   
    
      if (   (fp=find_token_in_file_nlines (name,NULL,"\nID",10))!=NULL\
	   &&(fp=find_token_in_file (name,NULL,"\nSQ"))!=NULL)
	{
	  
	  vfclose (fp);return 1;
	}
      else
	{
	  return 0;
	}
    } 

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT STUFF                                              */
/*                                                                                         */
/***************************************************************************************** */		    
int output_format_aln ( char *format, Alignment *inA, Alignment *inEA,char *name)
        {
	
	Sequence_data_struc *D1=NULL;
	Sequence_data_struc *D2=NULL;
	Alignment *A;
	Alignment *EA;

	
	

	A =copy_aln (inA, NULL);
	EA=copy_aln (inEA,NULL);
		
	
	A =expand_aln(A);
	EA=expand_number_aln(EA);
		
	  
        D1=vcalloc ( 1, sizeof (Sequence_data_struc));
	D1->A=A;
	
	if (EA)
	   {
	   D2=vcalloc ( 1, sizeof (Sequence_data_struc));
	   D2->A=EA;
	   }
       main_output ( D1, NULL,D2, format, name);
       return 1;
       
       }



void main_output  (Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char *out_format, char *out_file)

	{  
	FILE *fp;
	int value;

	if ( D1->rm_gap)ungap_aln ((D1->A));
	
	if ( strm (out_format, ""))return;


	else if      ( strncmp (out_format, "score",5)==0)
		{
	
		if ( !DST) 
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL]\n");	
		   exit(1);
		   }
		while ( out_format[0]!='_' && out_format[0]!='\0' )out_format++;
				
		(DST->A)=aln2number (DST->A);
			
		if     ( strm ( out_format, "_html"  ))output_reliability_html  ( D1->A,  DST->A, out_file);
		else if( strm ( out_format, "_ps"    ))output_reliability_ps    ( D1->A,  DST->A, out_file);
		else if( strm ( out_format, "_pdf"   ))output_reliability_pdf   ( D1->A,  DST->A, out_file);	
		else if( strm ( out_format, "_ascii" ))output_reliability_ascii ( D1->A,  DST->A, out_file);	
		
		}
       
       else if      ( strncmp (out_format, "color",5)==0)
	 {

	   
		if ( !DST) 
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL]\n");	
		   exit(1);
		   }
		while ( out_format[0]!='_' && out_format[0]!='\0' )out_format++;
		

		(DST->A)=aln2number (DST->A);
			
		if     ( strm ( out_format, "_html"  ))output_color_html  ( D1->A,  DST->A, out_file);
		else if( strm ( out_format, "_ps"    ))output_color_ps    ( D1->A,  DST->A, out_file);
		else if( strm ( out_format, "_pdf"   ))output_color_pdf   ( D1->A,  DST->A, out_file);	
		
		
		}
	else if ( strm4  ( out_format, "tc_aln","t_coffee_aln", "t_coffee", "tcoffee"))
	        {
		vfclose (output_aln ( D1->A, vfopen (out_file, "w")));
		}
	else if ( strm  ( out_format, "analyse_pdb"))
	        {
		 if ( !DST) 
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL]\n");	
		   exit(1);
		   }
		analyse_pdb ( D1->A,DST->A);
		(DST->A)=aln2number (DST->A);
		output_reliability_ps    ( D1->A,  DST->A, out_file);
		}
	else if ( strm4 ( out_format, "lower0", "lower1", "lower2", "lower3") || strm4(out_format, "lower4", "lower5", "lower6", "lower7") || strm4 (out_format,"lower8", "lower9", "align_pdb", "malign_pdb") )
	       {
	       if ( !DST) 
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL]\n");	
		   exit(1);
		   }
	        

	       
	       (DST->A)=aln2number (DST->A);
	       if ( strm (out_format, "align_pdb"))value=0;
	       else if (  strm (out_format, "malign_pdb"))value=5;
	       else value=atoi(out_format+5);

	       D1->A=filter_aln_lower_upper (D1->A, DST->A,value);
	       output_clustal_aln ( out_file, D1->A);
	       }
	else if ( strm4 ( out_format, "filter0", "filter1", "filter2", "filter3"))
	       {
	       if ( !DST) 
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL]\n");	
		   exit(1);
		   }
	        
	       (DST->A)=aln2number (DST->A);
	       
	       D1->A=filter_aln (D1->A, DST->A, atoi(out_format+6));
	       output_clustal_aln ( out_file, D1->A);
	       }
	
	else if ( strm3 ( out_format, "phylip_aln", "phylip", "phy"))
		{
		output_phylip_aln ( out_file, D1->A);
		}
	else if ( strm ( out_format, "mocca_aln"))
	        {
		output_mocca_aln ( out_file, D1->A, DST->A);
		}       
	else if ( strm ( out_format, "saga_pw_sd_weights") )
		{
		output_pw_weights4saga ((D1->W),(D1->W)->PW_SD, out_file);
		}
	else if ( strm ( out_format, "saga_aln"))
		{
		output_saga_aln (out_file, D1->A);
		}
	else if (strm5 ( out_format, "aln","clustal_aln", "clustalw","clustal", "clustalw_aln"))
		{		  
		output_clustal_aln (out_file, D1->A);
		}
	else if ( strm2 ( out_format, "lalign_aln","lalign"))
	        {
		output_lalign_aln (out_file, D1->A);
		}
	else if ( strm2 ( out_format, "fasta_aln","fasta" ))
		{
		output_fasta_aln( out_file, D1->A);
		}
	else if ( strm ( out_format, "est_prf" ))
		{
		output_est_prf( out_file, D1->A);
		}
	else if ( strm ( out_format, "clean_est_fasta_seq" ))
		{
		
		D1->A=clean_est(D1->A);
		output_fasta_seq(out_file, D1->A);
		
		}
	
	else if ( strm3 ( out_format, "msf_aln", "gcg", "msf"))
		{
		output_msf_aln( out_file, D1->A);
		}
	else if ( strm ( out_format, "rnalign"))
		{
		output_rnalign (out_file, D1->A, DST->S);
		}
	else if ( strm ( out_format, "fasta_seq"))
		{
		output_fasta_seq (out_file,D1->A);
		}
	else if ( strm ( out_format, "gotoh_seq"))
		{
		output_gotoh_seq (out_file,D1->A);
		}
	else if ( strm (out_format, "fasta_seq1"))
		{
		output_fasta_seq1 (out_file, D1->A);
		}
	else if ( strm2 (out_format, "pir_aln", "pir"))
		{
		output_pir_aln (out_file, D1->A);
		}
	else if ( strm (out_format, "pir_seq"))
		{
		output_pir_seq (out_file, D1->A);
		}
        else if ( strm (out_format, "gor_seq"))
		{
		output_gor_seq (out_file, D1->A);
		}
	else if ( strm (out_format, "pir_seq1"))
		{
		output_pir_seq1 (out_file, D1->A);
		}
	else if ( strm (out_format, "pw_lib_saga_aln"))
		{
		output_pw_lib_saga_aln (out_file, D1->A);
		}
	else if ( strm (out_format, "lib"))
		{
		output_lib (out_file, D1->A);
		}
	else if ( strm (out_format, "pdb_constraint_list"))
	        {
		output_constraints (out_file, "pdb",D1->A);
		}
	else if ( strm (out_format, "constraint_list"))
	        {
		output_constraints (out_file,"sim", D1->A);
		}
	else if ( strm (out_format, "cache_id"))
		{
		cache_id (D1->A);
		output_saga_aln (out_file, D1->A);
		}	
        else if ( strm (out_format, "compress_aln"))
		{
                compress_aln (D1->A);
		output_saga_aln (out_file, D1->A);
		} 
	else if (strm (out_format, "n_seq"))
		{
		fp=vfopen ( out_file, "w");
		fprintf ( fp, "%d\n", (D1->A)->nseq);
                vfclose (fp);
		}
	else if ( strm ( out_format, "thread_dna_on_prot_aln"))
	        {
		D1->A=thread_dnaseq_on_prot_aln (D1->S, D2->A);
		output_saga_aln ( out_file, D1->A);
		}
	else if ( strm ( out_format, "tdna_fasta_seq1"))
	        {
		D1->A=translate_dna_aln (D1->A,0);
		output_fasta_seq1 (out_file, D1->A);
		}
	else if ( strm ( out_format, "tdna_aln"))
	        {		
		D1->A=translate_dna_aln (D1->A,0);
		output_saga_aln ( out_file, D1->A);
		}
	else if ( strm ( out_format, "cdna_fasta_seq1"))
	        {		
		D1->A= gene2prot(D1->A);
		output_fasta_seq1 ( out_file, D1->A);
		}
	else if ( strm ( out_format, "mutate_cdna_aln"))
	        {
		    D1->A= mutate_cdna_aln ( D1->A);
		    output_clustal_aln ( out_file, D1->A);
		}
	else if ( strm ( out_format, "tdna_sp_aln"))
	        { 
	        if ( !DST) 
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL]\n");	
		   exit(1);
		   }	        
	       (DST->A)=aln2number (DST->A);
		D1->A=translate_splice_dna_aln (D1->A, DST->A);
		output_saga_aln ( out_file, D1->A);
		}
	else if (out_format && out_format[0] && (strcmp ( out_format,"rna_graph_fasta")==0))
		{
		sprintf ( (D1->A)->seq_al[0], "%s",(DST->S)->seq[0]);
		(D1->A)->nseq=0;
		output_fasta_seq (out_file, DST->A);
		}
	else if (strm ( out_format, "freq_mat"))
	        {
		output_freq_mat (out_file, D1->A);  
		}
	else if (strm ( out_format, "maln_pval"))
	        {
		output_maln_pval ( out_file, D1->A);
		}
	else if ( strm ( out_format, "model_aln"))
	        {
		output_model_aln ( out_file, D1->A);
		}
	else if (strncmp (out_format, "mult",4)==0)
	        {
		output_mult_fasta_seq ( out_file, D1->A, atoi(out_format+4));
		}
	else if (strm (out_format, "len"))
	        {
		output_statistics (out_file, D1->A, "nr");
		}
	else if ( strm (out_format, "name"))
	        {
		output_statistics (out_file, D1->A, "n");
		}
	else if ( strncmp (out_format, "statistics",10)==0)
	        {
		output_statistics (out_file, D1->A,out_format+10);
		}       
	else 
	        {
		    fprintf ( stderr, "\n%s is an UNKNOWN FORMAT[FATAL]\n",out_format); 
		    crash ("");
		    
		}
	
	}
int is_in_format_list ( char *name)
	{
	if ( strcmp ( name, "saga_aln")==0)return 1;
	if ( strcmp ( name, "number_aln")==0)return 1;
	if ( strcmp ( name, "clustal_aln")==0)return 1;	
	if ( strcmp ( name, "fasta_aln")==0)return 1;
	if ( strcmp ( name, "number_fasta")==0)return 1;
	if ( strcmp ( name, "fasta_seq")==0)return 1;
	if ( strcmp ( name, "pdb")==0)return 1;
	if ( strcmp ( name, "msf_aln")==0)return 1;
	if ( strcmp ( name, "dali_aln")==0)return 1;
	if ( strcmp ( name, "dali_seq")==0)return 1;
	if ( strcmp ( name, "barton_list_tc")==0)return 1;
	if ( strcmp ( name, "est_prf")==0)return 1;
	
	if ( strcmp ( name, "gotoh_aln")==0)return 1;
	if ( strcmp ( name, "amps_aln")==0)return 1;
	if ( strcmp ( name, "pir_aln")==0)return 1;
	if ( strcmp ( name, "pir_seq")==0)return 1;
	if ( strcmp ( name, "est_fasta")==0)return 1;
	if ( strcmp ( name, "amps_sd_scores")==0)return 1;
	if ( strcmp ( name, "pima_aln")==0)return 1;
	if ( strcmp ( name, "dialign_aln")==0)return 1;
	if ( strcmp ( name, "gor_seq")==0)return 1;
	if ( strcmp ( name, "gor_struc")==0)return 1;
	return 0;
	}
int is_struc_in_format_list ( char *name)
	{
	if ( strcmp ( name, "rna_number")==0)return 1;
	if ( strcmp ( name, "fasta_seq")==0)return 1;
	return 0;
	}
int is_out_format_list ( char *name)
	{

	if ( strcmp ( name, "saga_aln")==0)return 1;
	if ( strcmp ( name, "clustal_aln")==0)return 1;
	if ( strcmp ( name, "fasta_aln")==0)return 1;

	if ( strcmp ( name, "fasta_seq")==0)return 1;
	if ( strcmp ( name, "fasta_seq1")==0)return 1;
	if ( strcmp ( name, "msf_aln")==0)return 1;
	if ( strcmp ( name, "phylip_aln")==0)return 1;
	if ( strcmp ( name, "rnalign")==0)return 1;
	if ( strcmp ( name, "pw_lib_saga_aln")==0)return 1;
	if ( strcmp ( name, "lib")==0)return 1;
	if ( strcmp ( name, "pir_aln")==0)return 1;
	if ( strcmp ( name, "pir_seq")==0)return 1;
	if ( strcmp ( name, "pir_seq1")==0)return 1;
	if ( strcmp ( name, "gotoh_seq")==0)return 1;
	if ( strcmp ( name, "saga_pw_sd_weights")==0)return 1;
	if ( strcmp ( name, "constraint_list")==0)return 1;
	if ( strcmp ( name, "pdb_constraint_list")==0)return 1;
	if ( strcmp ( name, "cache_id")==0)return 1;
	if ( strcmp ( name, "gor_seq")==0)return 1;
	if ( strcmp ( name, "n_seq")==0)return 1;
	if ( strcmp ( name, "thread_dna_on_prot_aln")==0)return 1;
	if ( strcmp ( name, "tdna_aln")==0)return 1;
	if ( strcmp ( name, "cdna_fasta_seq1")==0)return 1;
	if ( strcmp ( name, "tdna_fasta_seq1")==0)return 1;
	
	if ( strcmp ( name, "filter0")==0)return 1;
	if ( strcmp ( name, "filter1")==0)return 1;
	if ( strcmp ( name, "filter2")==0)return 1;
	if ( strcmp ( name, "filter3")==0)return 1;
	
	if ( strcmp ( name, "lower0")==0)return 1;
	if ( strcmp ( name, "lower1")==0)return 1;
	if ( strcmp ( name, "lower2")==0)return 1;
	if ( strcmp ( name, "lower3")==0)return 1;
	if ( strcmp ( name, "lower4")==0)return 1;
	if ( strcmp ( name, "lower5")==0)return 1;
	if ( strcmp ( name, "lower6")==0)return 1;
	if ( strcmp ( name, "lower7")==0)return 1;
	if ( strcmp ( name, "lower8")==0)return 1;
	if ( strcmp ( name, "lower9")==0)return 1;
	if ( strcmp ( name, "align_pdb")==0)return 1;
	if ( strcmp ( name, "malign_pdb")==0)return 1;
	if ( strm ( name, "maln_pval"))return 1;
	if ( strm ( name, "freq_mat"))return 1;

	return 0;
	}
	
int is_struc_out_format_list ( char *name)
	{
	if ( strcmp ( name, "rna_graph_fasta")==0)return 1;
	return 0;
	}	

/**************************************************************************************************/
/*************************************REFORMAT UTIL*************************************************/
/**************************************************************************************************/





void read_check ( Alignment *A, char *check_file)
	{
	FILE *fp;
	int a;
	int l;
	static char *buf;
	
	if (buf==NULL) buf=calloc ( 1001, sizeof (char));
	fp=vfopen ( check_file, "r");
	
	for ( a=0; a< A->nseq; a++)
		{
		fscanf ( fp, "%s", A->name[a]);
		fgets ( buf, 1000, fp);
		}
	vfclose ( fp);
	fp=vfopen ( check_file, "r");
	
	for ( a=0; a< A->nseq; a++)
		{
		fgets ( A->comment[a], 1000, fp);
		l=strlen ( A->comment[a]);
		A->comment[a][l-1]='\0';
		}
	
	vfclose ( fp);
	}
	
/**************************************************************************************************/
/*************************************REFORMAT IN**************************************************/
/**************************************************************************************************/
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT WEIGHTS                                            */
/*                                                                                         */
/***************************************************************************************** */
	
Weights* get_amps_sd_scores ( char *fname)
	{
	FILE *fp;
	char *buf;
	char *buf2;
	int nseq;
	Weights *W;
	int a, b,e;
	int c;
	float array[20];
	
	buf=calloc ( 1001, sizeof (char));
	buf2=calloc ( 1001, sizeof (char));
	
	fp=vfopen ( fname, "r");
	set_fp_id ( fp, "Index");
	buf=fgets ( buf, 1000, fp);
	fscanf ( fp, "%s", buf2);	
	
	nseq=0;
	while ( isalnum(buf2[0]) && !isalpha(buf2[0]))
		{
		nseq++;
		buf=fgets ( buf, 1000, fp);
		fscanf ( fp, "%s", buf2);
		}
	vfclose ( fp);
	
	W=declare_weights (nseq);
	
	fp=vfopen ( fname, "r");
	set_fp_id ( fp, "Index");
	buf=fgets ( buf, 1000, fp);
	fscanf ( fp, "%s", buf2);	
	
	a=0;
	while ( isalnum(buf2[0]) && !isalpha(buf2[0]))
		{
		fp=set_fp_after_char (fp, '>');
		fscanf ( fp, "%s",W->seq_name[a]);
		buf=fgets ( buf, 1000, fp);
		fscanf ( fp, "%s", buf2);
		a++;
		}
	buf=fgets ( buf, 1000, fp);
	c=1;
	while ( c!=0)
		{
		for ( e=0; e< 16; e++)
			{
			c=fscanf ( fp, "%f", &array[e]);
			}
		fscanf ( fp, "\n");
		if ( c!=0)
			{
			
			a=(int)array[0]-1;
			b=(int)array[1]-1;
			W->PW_ID[b][a]=W->PW_ID[a][b]=array[9];
			W->PW_SD[b][a]=W->PW_SD[a][b]=array[14];
			}
		
		}
	vfclose ( fp);
	sprintf ( W->comments, "SD WEIGHTS GENERATED WITH THE PROGRAM AMPS IN PAIRWISE MODE");
	free ( buf);
	return W;
	}

Weights *read_seq_weight (char **name, int nseq, char* seq_weight)
       {
       int a, p;
       Weights *W;
       float w;
       
       FILE *fp;
       char line[LONG_STRING];
       char sname[MAXNAMES];
       
       
       /*Read sequence weights:
	* comment
	name1 weight1
	.....


	NOTE:
	weights must be between 0 and 1;
	
	sequences not in S do not get any weight
	sequences in S but not in file get a weight of 1
       */
       W=declare_weights(nseq);
       for ( a=0; a< nseq; a++)
	 {
	   sprintf ( W->seq_name[a], "%s", name[a]);
	   W->SEQ_W[a]=1;
	 }
       sprintf ( W->mode, "%s", seq_weight);
       fp=vfopen (seq_weight, "r");


       while ( fgets( line,LONG_STRING-1, fp))
	 {
	   if ( line[0]=='*' ||line[0]=='#' || isblanc(line));
	   else
	     {
	       if (sscanf(line, "%s %f", sname, &w)!=2)continue;
	       if ( (p=name_is_in_list ( sname, W->seq_name, nseq, MAXNAMES-1))!=-1)
		 {
		   W->SEQ_W[p]=w;
		 }
	     }
	 }
       vfclose (fp);
       return W;
       }
       
       
  
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT MISC                                               */
/*                                                                                         */
/***************************************************************************************** */


	

void get_barton_list_tc_seq ( char *in_file)
	{
	FILE *fp, *fp_make, *fp_length, *fp_long;
	FILE *fp_small[9];
	
	static char *buf;
	int len_buf=10000;
	char name[100];
	
	char pwd[100];
	int a,c,nseq;
	int k=0;
	int *length;
	int longest=0;
	
	c=0;
	length=calloc ( 1000, sizeof(int));
	if ( buf==NULL)buf=calloc ( len_buf, sizeof (char));
	fp=vfopen (in_file, "r");
	fp_long=vfopen ( "barton_seq_list_large", "w");
	fp_make=vfopen ( "make_dir", "w");
	fp_length=vfopen ( "barton_length", "w");
	for ( a=0; a< 9; a++)
		{
		sprintf ( name, "barton_nseq%d",a);
		fp_small[a]=vfopen ( name, "w");
		}
	get_pwd (pwd);
	
	
	while ( c!=EOF)
		{a=0;
		while ( (c=fgetc(fp))!='#');
		while ( (c=fgetc(fp))=='#');
		ungetc ( c, fp);
		while ( (c=fgetc(fp))!='#')buf[a++]=c;
		buf[a]='\0';
		
		sprintf ( name, "%s", buf);
	
		while ( (c=fgetc(fp))=='#');
		
		if ( c!=EOF)
			{
			a=0;
			while ( (c=fgetc(fp))!='#' && c!=EOF)
				{
				buf[a++]=c;
				if (a==len_buf)
					{
					len_buf+=10000;
					buf=realloc ( buf, len_buf*sizeof (char));
					}
				} 
			buf[a]='\0';
			if (c!=EOF)
				{
				
				nseq=process_barton_entry ( buf,name);
				length[nseq]++;
				longest=(longest<nseq)?nseq:longest;
				
				if ( nseq<=8) fprintf ( fp_small[nseq], "%s.pep\n", name);
				else fprintf ( fp_long, "%s.pep\n",name);
				fprintf ( fp_make, "mkdir %s\nmv %s.pep %s\nmv %s.check %s\n", name, name, name, name, name);
				k++;
				}
			
				
			}
		}
	
	vfclose (fp);
	vfclose (fp_long);
	for ( a=0; a< 9; a++)vfclose (fp_small[a]);
	vfclose (fp_make);
	for ( a=0; a<= longest; a++)fprintf ( fp_length, "%d: %d\n", a, length[a]);
	vfclose ( fp_length);
	
	}
	
int process_barton_entry (char *buf, char *name)			
    {	    
    Alignment *LA;
    Sequence *LS;
    int a,c;
    static char *buf2;
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0;
    int l;
    char fname[100];
    char com_name[100];
    int rm_gap=1;

    sprintf ( fname, "%s.pep", name);
    sprintf ( com_name, "%s.check",name);
    
    if ( buf2==NULL)buf2=calloc ( 10000, sizeof (char));
    a=0;		
    while (buf[a]!='\0')
	 	{
		 if (buf[a]=='>')
			{
			a=get_string_line (a,2, buf, buf2); 
			while ((c=buf[a++])!='*')
				if (isalnum (c)|| c=='.' || c=='-')
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			 clen=0;
			}
		if ( buf[a]!='\0')a++;
		}
    
    
    LS=declare_sequence (  min_len_seq,  max_len_seq,  nseq); 
    LS->nseq=nseq;
    
    
    for (a=0, current=0; current< nseq; current++) 
    	{
    	a=get_string_line ( a, 1, buf, buf2);
    	sscanf ( buf2, ">P1;%s", LS->name[current]);
    	a=get_string_line ( a, 1, buf, buf2);
    	l=strlen ( buf2);
	buf2[l-1]='\0';
    	sprintf ( LS->comment[current], buf2);
    	
    	p=0;
    	while ( (c=buf[a++])!='*')
    		{
    		if (isalpha (c))
			LS->seq[current][p++]=tolower (c);
		else if ( isgraph(c))
			LS->seq[current][p++]=(c);
		}	
    	a++;
    	}
    
    LA=declare_Alignment(LS);
    seq2aln ( LS, LA,rm_gap);
    output_fasta_seq (fname,LA);
    output_pir_check (com_name,LA->nseq, LA->comment);
    free_Alignment ( LA);
    free_sequence ( LS, nseq);   
    
    return nseq;
    }

	
	

Sequence *read_rna_struc_number (Alignment *A,char *fname)
	{
	FILE *fp;
	int a;
	char x,y;
	float f;
	Sequence *SA;
	int first, last;
	
	SA=declare_sequence ( A->len_aln, A->len_aln, 1);
	
	
	SA->len[0]=A->len[0];
	for ( a=0; a< SA->len[0]; a++)
		SA->seq[0][a]='.';
	SA->ST=declare_rna_structure_num (SA);
	
	fp=vfopen ( fname, "r");
	fscanf ( fp, "%c\n%d\n",&x, &(SA->ST)->tot_list);
	for ( a=0; a<(SA->ST)->tot_list; a++)
		{
		fscanf ( fp, "%d %d %d %c %c %f\n", &(SA->ST)->list[a][0],&(SA->ST)->list[a][1],&(SA->ST)->list[a][2], &x, &y, &f);
		(SA->ST)->list[a][0]--;
		(SA->ST)->list[a][1]--;
		(SA->ST)->list[a][2]--;
		if ( a==0)
			{
			(SA->ST)->stem[0][0]=(SA->ST)->list[a][0];
			(SA->ST)->stem[0][1]=a;
			}
		else if ( (SA->ST)->stem[(SA->ST)->tot_stem][0]==(SA->ST)->list[a][0]);
		else if ( (SA->ST)->stem[(SA->ST)->tot_stem][0]!=(SA->ST)->list[a][0])
			{
			(SA->ST)->stem[(SA->ST)->tot_stem][2]=a-1;
			(SA->ST)->tot_stem++;
			(SA->ST)->stem[(SA->ST)->tot_stem][0]=(SA->ST)->list[a][0];
			(SA->ST)->stem[(SA->ST)->tot_stem][1]=a;
			}
			
		SA->seq[0][(SA->ST)->list[a][1]]='-';
		SA->seq[0][(SA->ST)->list[a][2]]='-';
		}
	(SA->ST)->stem[(SA->ST)->tot_stem][2]=a-1;	
	(SA->ST)->tot_stem++;
	for ( a=0; a< (SA->ST)->tot_stem; a++)
     		{
     	
     		first=(SA->ST)->stem[a][1];
     		last=(SA->ST)->stem[a][2];
     		SA->seq[0][(SA->ST)->list[first][1]]='>';
     		SA->seq[0][(SA->ST)->list[first][2]]='<';
     		SA->seq[0][(SA->ST)->list[last][1]]='>';
     		SA->seq[0][(SA->ST)->list[last][2]]='<';	
     		}
	
	return SA;	
	}
		  
Structure * declare_rna_structure_num (Sequence *SA)
	{
	Structure *ST;
	ST=calloc ( 1, sizeof ( Structure));
	ST->list=declare_int ( SA->len[0], 3);
	ST->stem=declare_int ( SA->len[0], 3);
	return ST;
	}
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT SEQ                                                */
/*                                                                                         */
/***************************************************************************************** */
Sequence* get_pdb_sequence   ( char *fname)
     {
	 char *tp_name;
	 char *command;
	 Sequence *S;

	 command=vcalloc ( LONG_STRING, sizeof (char));
	 tp_name=vtmpnam (NULL);
	 
	 sprintf ( command, "extract_from_pdb -chain FIRST -infile %s -mode fasta > %s", fname, tp_name);
	 system ( command);
	

	 S=get_fasta_sequence ( tp_name, NULL);
	 
	 S->nseq=1;
	
	 sprintf ( S->file[0], "%s", fname);
	 S->max_len=S->min_len=S->len[0];
	 remove ( tp_name);
	 free ( command);
	 return S;
     }

char * get_pdb_file   ( char *fname)
     {
	 char *file;
	 int a, c;
	 FILE *fp;
	 

	 a=0;
	 file=vcalloc ( sizeof (char),count_n_char_in_file ( fname)+1);
	 fp=vfopen ( fname, "r");
	 while ( (c=fgetc(fp))!=EOF)file[a++]=c;
	 file[a]='\0'; 
	 return file;
     }
	 
Sequence* get_struc_gor ( char *fname)
    {
    int nseq, min_len, max_len;
    int a, c;
    int len;
    char name[STRING];
    

    FILE *fp;
    Sequence *S;

    min_len=max_len=-1;
    fp=vfopen ( fname, "r");
    nseq=0;
    while ( (c=fgetc(fp))!=EOF)
	    {
	    if ( c!='!');
	    else
		{
		nseq++;
		fscanf ( fp, "%s %d", name, &len);
		if (min_len==-1)min_len=max_len=len;
		else
		    {
		    min_len=(len>min_len)?min_len:len;
		    max_len=(len>max_len)?len:max_len;
		    }
		}
	    
	    }
    vfclose (fp);
   
    S=declare_sequence (  min_len,  max_len+1,nseq); 
    S->nseq=0;
    
    fp=vfopen (fname,"r");	
     while ( (c=fgetc(fp))!=EOF)
	     {
	     if ( c!='!');
	     else
	        {
		fscanf ( fp, "%s %d\n",S->name[S->nseq], &(S->len[S->nseq]));
		
		while ( (c=fgetc(fp))!='\n');
	
		for ( a=0; a<S->len[S->nseq]; a++)
		    fscanf ( fp, " %*c %c %*f %*f %*f\n",&(S->seq[S->nseq][a]));
		
		S->seq[S->nseq][a]='\0';
		while ( (c=fgetc(fp))!='!' && c!=EOF);
		ungetc (c, fp);
		S->nseq++;
		}
	     
	     }
    vfclose (fp);
    return S;		
    }
    	        
Sequence* get_sequence_dali (char *fname)
    {
    Sequence *LS;
    FILE *fp;
    int c;

    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0;
    
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  exit(1);
	 }  
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (isdigit(c))
			{
			ungetc(c, fp);
			fscanf (fp, "%s",name);
			while (!isdigit(c=fgetc(fp)) && c!=EOF)
				if (isalnum (c) || c=='.' || c=='-')
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			clen=0;
			}
		else
		    c=fgetc (fp);
		}
    vfclose (fp);
		
    LS=declare_sequence (  min_len_seq,  max_len_seq+1,nseq); 
    LS->nseq=nseq;
    
    fp=vfopen (fname,"r");
    
    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (isdigit(c))
			{
			ungetc(c, fp);
			fscanf (fp, "%s",LS->name[current]);
			p=0;
			while (!isdigit(c=fgetc(fp)) && c!=EOF)
				{
				if (isalpha (c))
				    LS->seq[current][p++]=tolower (c);
				else if ( c=='.')
				    LS->seq[current][p++]='-';
				else if ( c=='-')
				    LS->seq[current][p++]='-';
				}	    
			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);
			current++;
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);
    
    
    return LS;
    }	

Sequence* get_dialign_sequence (char *fname)
    {
    Sequence *LS;
    FILE *fp;
    int c;

    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0, l=0;
    char *buf;
    
    buf=calloc ( 1000, sizeof (char));
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  exit(1);
	 }  
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{fscanf (fp, "%s",name);
			
			buf=fgets ( buf, 1000, fp);
			while ((c=fgetc(fp))!='>' && c!=EOF && c!=' ' && c!='\t')
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			clen=0;
			}
		else
		    c=fgetc (fp);
		}
    vfclose (fp);
		
    LS=declare_sequence (  min_len_seq,  max_len_seq, nseq); 
    LS->nseq=nseq;
    
    fp=vfopen (fname,"r");
    
    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (c=='>')
			{
			fscanf (fp, "%s",LS->name[current]);
			l=strlen ( LS->name[current]);
			if ( LS->name[current][l-1]==','||LS->name[current][l-1]==',')LS->name[current][l-1]='\0';
			buf=fgets ( buf, 1000, fp);
			p=0;
			while ((c=fgetc(fp))!='>' && c!=EOF && c!=EOF && c!=' ' && c!='\t')
				if (isalpha (c))
				    LS->seq[current][p++]=tolower (c);
				else if ( isgraph(c))
				    LS->seq[current][p++]=(c);
			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);
			current++;
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);
    return LS;
    }

Sequence* get_pima_sequence (char *fname)
    {
    Sequence *LS;

    FILE *fp;
    int c;

    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0, l=0, len=0;
    char *buf, *buf2;
    char prefix[1000];
    
    sprintf (  prefix, "%s",fname);
    
    buf=strstr(prefix, "-");
    buf[0]='\0';
    len=strlen (prefix);
    	
   
    
    buf=calloc ( 1000, sizeof (char));
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  exit(1);
	 }  
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{fscanf (fp, "%s",name);
			if ( strlen(name)>=len && strncmp ( name, prefix, len)==0)
				{
				c=fgetc(fp);
				}
			else
				{
				
				buf=fgets ( buf, 1000, fp);
				while ((c=fgetc(fp))!='>' && c!=EOF)
					if (isalnum (c)|| is_gap(c))
						clen++;
				 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 	 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 	nseq++;
				clen=0;
				}
			}
		else
		    	c=fgetc (fp);
		}
    vfclose (fp);
		
    LS=declare_sequence (  min_len_seq,  max_len_seq, nseq); 
    LS->nseq=nseq;
    
    fp=vfopen (fname,"r");
    
    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (c=='>')
			{
			fscanf (fp, "%s",LS->name[current]);
			if ( strlen(LS->name[current])>=len && strncmp ( LS->name[current], prefix, len)==0)
				c=fgetc (fp);
			else
				{
				buf2=strstr (LS->name[current], ".");
				if ( buf2!=NULL) buf2[0]='\0';
				 
				l=strlen ( LS->name[current]);
				if ( LS->name[current][l-1]==','||LS->name[current][l-1]==',')LS->name[current][l-1]='\0';
				buf=fgets ( buf, 1000, fp);
				p=0;
				while ((c=fgetc(fp))!='>' && c!=EOF)
					if (isalpha (c))
					    LS->seq[current][p++]=tolower (c);
					else if ( isgraph(c))
					    LS->seq[current][p++]=(c);
				LS->seq[current][p]='\0';
				LS->len[current]=strlen ( LS->seq[current]);
				current++;
				}
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);
    return LS;
    }

Sequence* get_fasta_sequence_num (char *fname, char *comment_out)
    {
    Sequence *LS;
    char *buffer;
    FILE *fp;

    int   c;
    char *name;
    int clen=0;
    int current=0;
    int p=0;
    int max;
    int max_len_seq=0;
    int min_len_seq=0;
    int nseq=0, l=0;
 
    
    
    
    int *sub;
    
    buffer=vcalloc (1000, sizeof (char)); 
    name=vcalloc ( 100, sizeof (char));

    nseq=count_n_char_x_in_file(fname, '>');
    min_len_seq=max=count_n_char_in_file(fname);
    sub=vcalloc (max+1, sizeof (int));

    fp=vfopen (fname,"r");

    
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{
			fscanf (fp, "%s",name);
			while ((c=fgetc(fp))!='\n' && c!=EOF);
			while ((c=fgetc(fp))!='>' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 clen=0;
			}
		else
		    c=fgetc (fp);
		 
		}  

    vfclose (fp);		
    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq); 
    LS->nseq=nseq;
    
    fp=vfopen (fname,"r");
    current=0;
    c=fgetc(fp);
    while (c!=EOF)
		{
	 	if (c=='>')
			{
			
			fscanf (fp, "%s",LS->name[current]);
			l=strlen ( LS->name[current]);
			if ( LS->name[current][l-1]==','||LS->name[current][l-1]==';')LS->name[current][l-1]='\0';
			translate_name ( LS->name[current]);			
			while ((c=fgetc(fp))!='\n' && c!=EOF);
			
			p=0;
			while ((c=fgetc(fp))!='>' && c!=EOF)
			        {
				if (isalnum (c))
				    LS->seq[current][p++]=c;
				else if (is_gap(c))
				    LS->seq[current][p++]=c;				
				}

			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);

			current++;
		
			}
			
		else
		    c=fgetc ( fp);
		}
			
    
    vfclose (fp);
    

    vfree (sub);
    free (name);
    free (buffer);
    return LS;
    }
Sequence* get_fasta_sequence (char *fname, char *comment_out)
    {
    Sequence *LS;
    char *buffer;
    FILE *fp;

    int   c;
    char *name;
    int clen=0;
    int current=0;
    int p=0;
    int max;
    int max_len_seq=0;
    int min_len_seq=0;
    int nseq=0, l=0;

    
    
    
    int *sub;
    
    buffer=vcalloc (1000, sizeof (char)); 
    name=vcalloc ( 100, sizeof (char));

    nseq=count_n_char_x_in_file(fname, '>');
    min_len_seq=max=count_n_char_in_file(fname);
    sub=vcalloc (max+1, sizeof (int));

    fp=vfopen (fname,"r");

    
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{
			fscanf (fp, "%s",name);
			while ((c=fgetc(fp))!='\n' && c!=EOF);
			while ((c=fgetc(fp))!='>' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 clen=0;
			}
		else
		    c=fgetc (fp);
		 
		}  

    vfclose (fp);		
    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq); 
    LS->nseq=nseq;
    
    fp=vfopen (fname,"r");
    current=0;
    c=fgetc(fp);
    while (c!=EOF)
		{
	 	if (c=='>')
			{
			
			fscanf (fp, "%s",LS->name[current]);
			l=strlen ( LS->name[current]);
			if ( LS->name[current][l-1]==','||LS->name[current][l-1]==';')LS->name[current][l-1]='\0';
			translate_name ( LS->name[current]);			
			while ((c=fgetc(fp))!='\n' && c!=EOF);
		
			p=0;
			while ((c=fgetc(fp))!='>' && c!=EOF)
			        {
				if (isalnum (c))
				    LS->seq[current][p++]=tolower (c);
				else if (is_gap(c))
				    LS->seq[current][p++]=(c);				
				}

			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);

			current++;		
			}
			
		else
		    c=fgetc ( fp);
		}
			
  
    vfclose (fp);
    

    vfree (sub);
    free (name);
    free (buffer);
    return LS;
    }
Sequence* get_sub_fasta_sequence (char *fname, char *comment_out)
    {
    Sequence *LS;
    
    FILE *fp;

    int c;
    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max;
    int max_len_seq=0;
    int min_len_seq=0;
    int nseq=0, l=0;
    char *buf;
    
    
    
    int *sub;

    nseq=count_n_char_x_in_file(fname, '>');
    min_len_seq=max=count_n_char_in_file(fname);
    sub=vcalloc (max+1, sizeof (int));
    buf=calloc ( max+1, sizeof (char));
    fp=vfopen (fname,"r");

    
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{
			fscanf (fp, "%s",name);
			while ((c=fgetc(fp))!='\n' && c!=EOF);
			buf=fgets ( buf,max, fp);
			while ((c=fgetc(fp))!='>' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 clen=0;
			}
		else
		    c=fgetc (fp);
		 
		}  

    vfclose (fp);		
    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq); 
    LS->nseq=nseq;
    
    fp=vfopen (fname,"r");
    current=0;
    c=fgetc(fp);
    while (c!=EOF)
		{
	 	if (c=='>')
			{
			
			fscanf (fp, "%s",LS->name[current]);
			l=strlen ( LS->name[current]);
			if ( LS->name[current][l-1]==','||LS->name[current][l-1]==';')LS->name[current][l-1]='\0';
			translate_name ( LS->name[current]);
			while ((c=fgetc(fp))!='\n' && c!=EOF);
		
			p=0;
			while ((c=fgetc(fp))!='>' && c!=EOF)
			        {
				if (isalpha (c))
				    LS->seq[current][p++]=tolower (c);
				else if (is_gap(c))
				    LS->seq[current][p++]=(c);				
				}

			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);

			current++;
		
			}
			
		else
		    c=fgetc ( fp);
		}
			
    
    vfclose (fp);
    

    vfree (sub);
    return LS;
    }
Sequence* get_pir_sequence (char *fname, char *comment_out)
    {
    Sequence *LS;

    FILE *fp;
    int c;

    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=999999;
    int nseq=0, l=0;
    char *buf;
    
    buf=calloc ( 1000, sizeof (char));
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  exit(1);
	 }  
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{
			if ( (c=fgetc(fp))=='P')while ( (c=fgetc(fp))!=';');
			else ungetc ( c, fp);
			fscanf (fp, "%s",name);
			
			buf=fgets ( buf, 1000, fp);
			while ((c=fgetc(fp))!='>' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			clen=0;
			}
		else
		    c=fgetc (fp);
		}
    vfclose (fp);


   
    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq); 
    LS->nseq=nseq;
    
    fp=vfopen (fname,"r");
    
    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (c=='>')
			{
			if ( (c=fgetc(fp))=='P')while ( (c=fgetc(fp))!=';');
			else ungetc ( c, fp);

			fscanf (fp, "%s",LS->name[current]);
		
			l=strlen ( LS->name[current]);
			if ( LS->name[current][l-1]==','||LS->name[current][l-1]==',')LS->name[current][l-1]='\0';
			translate_name ( LS->name[current]);
			buf=fgets ( buf, 1000, fp);
			
			LS->comment[current]=fgets ( LS->comment[current], 1000, fp);
			LS->comment[current][strlen(LS->comment[current])-1]='\0';
			p=0;
			while ((c=fgetc(fp))!='>' && c!=EOF)
				if (isalpha (c))
				    LS->seq[current][p++]=tolower (c);
				else if ( !isspace(c) && c!='*')
				    LS->seq[current][p++]=(c);
			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);
			current++;
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);
    if (comment_out!=NULL) output_pir_check ( comment_out,LS->nseq, LS->comment);
    return LS;
    }

Sequence* get_gor_sequence (char *fname, char *comment_out)
    {
    Sequence *LS;

    FILE *fp;
    int c;

    char name[100];
    int clen=0;
    int current=0;
    int p=0;
    int max_len_seq=0;
    int min_len_seq=99999;
    int nseq=0;
    char *buf;
    
    buf=calloc ( 1000, sizeof (char));
    if ((fp=vfopen (fname,"r"))==NULL)
	 {printf ( "\nCOULDN'T OPEN %s",fname);
	  exit(1);
	 }  
    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='!')
			{
			fscanf (fp, "%s",name);
			
			buf=fgets ( buf, 1000, fp);
			while ((c=fgetc(fp))!='!' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
					clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 nseq++;
			clen=0;
			}
		else
		    c=fgetc (fp);
		}
    vfclose (fp);
		
    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq); 
    LS->nseq=nseq;
    
    fp=vfopen (fname,"r");
    
    current=0;
    c=fgetc(fp);
	while (c!=EOF)
		{
	 	if (c=='!')
			{
			
		       
			fscanf (fp, "%s",LS->name[current]);
			translate_name ( LS->name[current]);
			buf=fgets ( buf, 1000, fp);
			
			p=0;
			while ((c=fgetc(fp))!='!' && c!=EOF)
				if (isalnum (c)|| is_gap(c))
				    LS->seq[current][p++]=tolower (c);
				
			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);
			current++;
			}
		else
		    c=fgetc ( fp);
		}

    vfclose (fp);

    return LS;
    }
Sequence* get_swissprot_sequence (char *fname, char *comment_out)
    {
    Sequence *LS;
    FILE *fp;
    int c;
    char *buf;    
    int nseq=0;
    int len, max_len_seq=0, min_len_seq;
    
    if ( !check_file_exists(fname))
      {printf ( "\nCOULDN'T OPEN %s",fname);
	  exit(1);
      }  

    buf=vcalloc (LONG_STRING+1, sizeof (char));
    fp=NULL;   
    while ( (fp=find_token_in_file(fname,fp,"\nSQ")))
      {
	nseq++;
	fgets (buf, LONG_STRING, fp);
	len=0;
	while ((c=fgetc(fp))!='/')if(isalpha(c))len++;
	if ( max_len_seq==0)max_len_seq=min_len_seq=len;
	else
	  {
	    max_len_seq=MAX(len, max_len_seq);
	    min_len_seq=MIN(len, min_len_seq);
	  }
      }

    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq);     
    LS->nseq=0;
    
    fp=NULL;
    while ( (fp=find_token_in_file(fname,fp,"\nID")))
      {
	fscanf (fp, "%s", LS->name[LS->nseq]);
	fp=find_token_in_file(fname,fp,"\nSQ");
	fgets (buf, LONG_STRING, fp);
	while ((c=fgetc(fp))!='/')if (isalpha(c))LS->seq[LS->nseq][LS->len[LS->nseq]++]=tolower(c);
	LS->seq[LS->nseq][LS->len[LS->nseq]]='\0';
	LS->nseq++;
      }

   
    return LS;
    }

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT ALN                                                 */
/*                                                                                         */
/***************************************************************************************** */
void read_aln ( char *file_name, Alignment *A)
   {
    FILE *fp, *fp2;
    int * ptr_aln;
    int a,b,d;
    int c;
    char buf[VERY_LONG_STRING+1];

    int tot=0;
    int flag=0;
    char *fname;   
    int n_comment=0;

    int nseq=0;
    int max_len=0;

    
    fp=vfopen ( file_name, "r");
    
    fname=vtmpnam(NULL);
    fp2=vfopen ( fname, "w");
    while ( (c=fgetc(fp))!=EOF)
	{
	    fprintf ( fp2, "%c", c);
	}
    vfclose (fp);
    vfclose (fp2);

  
    /*1 Count The number of sequences*/ 
    fp=vfopen ( fname, "r");
    fgets ( buf, 1000, fp);
    if ( !isblanc (buf));
    while ( isblanc (buf))
    	{
        fgets ( buf, 1000, fp);    	
    	}
    while (!isblanc (buf))
    	{
    	fgets ( buf, 1000, fp);
	}
    while ( !isalnum ((c=fgetc(fp))))
    	{    	
    	ungetc(c,fp);
    	fgets ( buf, 1000, fp);    	
    	}
    
    if ( c!='\n')ungetc(c,fp);
    
    while ( isalnum ((c=fgetc(fp))))
    	{
    	ungetc(c,fp);    	
    	a=0;
    	while ( isgraph ((c=fgetc(fp))));    		    	
        nseq++;
    	fgets ( buf, 1000, fp);
    	}    
    vfclose (fp);

    /*DONE*/
    /*2 get_max_len*/
    max_len=count_n_char_in_file(fname)/nseq;
    A=realloc_alignment2( A, nseq+1, max_len+1);

    /*DONE*/
        
   
    fp=vfopen ( fname, "r");
    fgets ( buf, 1000, fp);
    if ( !isblanc (buf))sprintf (A->comment[n_comment++], "%s", buf);
    while ( isblanc (buf))
    	{
        fgets ( buf, 1000, fp);    	
    	}
    while (!isblanc (buf))
    	{
    	fgets ( buf, 1000, fp);
	sprintf ( A->comment[n_comment++], "%s", buf);
    	
	}
    while ( !isalnum ((c=fgetc(fp))))
    	{    	
    	ungetc(c,fp);
    	fgets ( buf, 1000, fp);
    	
    	}
    
    if ( c!='\n')ungetc(c,fp);
    
    while ( isalnum ((c=fgetc(fp))))
    	{
    	ungetc(c,fp);
    	
    	a=0;
    	while ( isgraph ((c=fgetc(fp))))
    		A->name[A->nseq][a++]=c;
    	A->name[A->nseq][a]='\0';
	if ( name_is_in_list (A->name[A->nseq], A->name, A->nseq, 100)!=-1){fprintf ( stderr, "\nERROR: Sequence %s Duplicated in File %s", A->name[nseq], A->file[0]);exit(1);}      
    	A->nseq++;
    	fgets ( buf, 1000, fp);
    	}
    
    vfclose (fp);
  
    
   
    if ((fp=vfopen ( fname, "r"))==NULL)
	printf ( "\nCOULDN'T READ %s", fname);
    
    ptr_aln=calloc ( A->nseq, sizeof(int));
    while ( flag==0)
	{
	while ( !isalnum(c=fgetc(fp)));
	fgets ( buf, VERY_LONG_STRING, fp);
	while ( is_alpha_line ( buf))
	  fgets ( buf, VERY_LONG_STRING, fp);
	flag=1;
	}
    while ( !isalnum(c=fgetc(fp)));
    ungetc ( c, fp);
    while ( c!=EOF)
	{
	tot=0;
	while(tot< A->nseq)
	    {
	     b=0;
	     while ( isgraph((buf[b++]=fgetc(fp))));
	     buf[b-1]='\0';
	     for ( d=0; d< A->nseq; d++)
		if ( strcmp (A->name[d], buf)==0)
		    {a=d;
		     tot++;
		    }
	     while ( !isgraph(c=fgetc(fp)));
	     ungetc ( c, fp);
	     
	     while ( (c=fgetc(fp))!='\n' && !isspace(c))
		{
		 if ( isgraph(c) || is_gap(c))
		    {
		     if ( isalpha(c))
			A->seq_al[a][ptr_aln[a]++]=(A->residue_case==2)?c:tolower(c);
		     else if (isdigit(c))
		       A->seq_al[a][ptr_aln[a]++]=c;
		     else if (is_gap(c))
		       A->seq_al[a][ptr_aln[a]++]=c;
		    }
		}
	     if (c!='\n') while ( (c=fgetc(fp))!='\n');
	     }
	 while ( !isalnum(c=getc(fp)) && c!=EOF);
	 if ( c!=EOF)
	    ungetc (c, fp);
	 }
	 
    vfclose (fp);
    
   
    for ( a=0; a< A->nseq; a++)
	{A->seq_al[a][ptr_aln[a]]='\0';
	 A->order[a][0]=a;
	 A->order[a][1]=0;
	}
    
    A->len_aln= strlen(A->seq_al[0]);  
    
    free(ptr_aln);
    remove (fname);
    
    }	
void read_number_aln ( char *file_name, Alignment *A)
   {
    FILE *fp, *fp2;
    int * ptr_aln;
    int a,b,d;
    int c;
    char buf[1001];

    int tot=0;
    int flag=0;
    char *fname;   
    int n_comment=0;

    int nseq=0;
    int max_len=0;

    
    fp=vfopen ( file_name, "r");
    
    fname=vtmpnam(NULL);
    fp2=vfopen ( fname, "w");
    while ( (c=fgetc(fp))!=EOF)
	{
	    fprintf ( fp2, "%c", c);
	}
    vfclose (fp);
    vfclose (fp2);

  
    /*1 Count The number of sequences*/ 
    fp=vfopen ( fname, "r");
    fgets ( buf, 1000, fp);
    if ( !isblanc (buf));
    while ( isblanc (buf))
    	{
        fgets ( buf, 1000, fp);    	
    	}
    while (!isblanc (buf))
    	{
    	fgets ( buf, 1000, fp);
	}
    while ( !isalnum ((c=fgetc(fp))))
    	{    	
    	ungetc(c,fp);
    	fgets ( buf, 1000, fp);    	
    	}
    
    if ( c!='\n')ungetc(c,fp);
    
    while ( isalnum ((c=fgetc(fp))))
    	{
    	ungetc(c,fp);    	
    	a=0;
    	while ( isgraph ((c=fgetc(fp))));    		    	
        nseq++;
    	fgets ( buf, 1000, fp);
    	}    
    vfclose (fp);

    /*DONE*/
    /*2 get_max_len*/
    max_len=count_n_char_in_file(fname)/nseq;
    A=realloc_alignment2( A, nseq+1, max_len+1);

    /*DONE*/
        
   
    fp=vfopen ( fname, "r");
    fgets ( buf, 1000, fp);
    if ( !isblanc (buf))sprintf (A->comment[n_comment++], "%s", buf);
    while ( isblanc (buf))
    	{
        fgets ( buf, 1000, fp);    	
    	}
    while (!isblanc (buf))
    	{
    	fgets ( buf, 1000, fp);
	sprintf ( A->comment[n_comment++], "%s", buf);
    	
	}
    while ( !isalnum ((c=fgetc(fp))))
    	{    	
    	ungetc(c,fp);
    	fgets ( buf, 1000, fp);
    	
    	}
    
    if ( c!='\n')ungetc(c,fp);
    
    while ( isalnum ((c=fgetc(fp))))
    	{
    	ungetc(c,fp);
    	
    	a=0;
    	while ( isgraph ((c=fgetc(fp))))
    		A->name[A->nseq][a++]=c;
    	A->name[A->nseq][a]='\0';
	if ( name_is_in_list (A->name[A->nseq], A->name, A->nseq, 100)!=-1){fprintf ( stderr, "\nERROR: Sequence %s Duplicated in File %s", A->name[nseq], A->file[nseq]);exit(1);}      
    	A->nseq++;
    	fgets ( buf, 1000, fp);
    	}
    
    vfclose (fp);
  
    
     
    if ((fp=vfopen ( fname, "r"))==NULL)
	printf ( "\nCOULDN'T READ %s", fname);
   
    ptr_aln=calloc ( A->nseq, sizeof(int));
    while ( flag==0)
	{
	while (  (c=fgetc(fp))!='\n');
	if ( (c=fgetc(fp))=='\n')
	    flag=1;
	}
    while ( !isalnum(c=fgetc(fp)));
    ungetc ( c, fp);
    while ( c!=EOF)
	{
	tot=0;
	while(tot< A->nseq)
	    {
	     b=0;
	     while ( !isgraph (c=fgetc(fp)) && c!=EOF);
	     if ( c!=EOF)ungetc(c, fp);
	     while ( isgraph((buf[b++]=fgetc(fp))));
	     buf[b-1]='\0';
	     for ( d=0; d< A->nseq; d++)
		if ( strcmp (A->name[d], buf)==0)
		    {a=d;
		     tot++;
		    }
	     while ( (c=fgetc(fp))!='\n')
		{
		 if ( isgraph(c) || is_gap(c))
		    {if ( isalpha(c))
			c=(A->residue_case==2)?c:tolower(c);
		   
		    if (!isspace(c))A->seq_al[a][ptr_aln[a]++]=c;
		    }
		}
	     }
	 while ( !isalnum(c=getc(fp)) && c!=EOF);
	 if ( c!=EOF)
	    ungetc (c, fp);
	 }
	 
    vfclose (fp);
    
   
    for ( a=0; a< A->nseq; a++)
	{A->seq_al[a][ptr_aln[a]]='\0';
	 A->order[a][0]=a;
	 A->order[a][1]=0;
	}
    
    A->len_aln= strlen(A->seq_al[0]);  
    
    free(ptr_aln);
    remove (fname);
    
    }		
void read_amps_aln ( char *in_file, Alignment *A)
	{
	FILE *fp;
	int a, b, c, cont=1;
	A->nseq=get_amps_seq_name ( A->name, in_file);
	
	fp=vfopen ( in_file, "r");
	fp=set_fp_id(fp, "1*");
	while ( (c=fgetc(fp))!='\n');
	b=0;
	while ( cont==1)
		{
		c=fgetc ( fp);
		c=fgetc(fp);
		if ( c=='*')
			{
			cont=0;
			for ( a=0; a<A->nseq; a++)
				A->seq_al[a][b]='\0';
			A->len_aln=b;
			}
			 
		else
		    	{
		    	ungetc (c, fp);
		  	for ( a=0; a< A->nseq; a++)
		  		{
		  		c=fgetc(fp);
		  		if ( c==' ')A->seq_al[a][b]='-';
		  		else
		  			{
		  			A->seq_al[a][b]=c;
		  			A->len[a]++;
		  			}
		  		}
		  	while ((c=fgetc(fp))!='\n');
		  	b++;
		  	}
		}
	}






int get_amps_seq_name ( char **name, char* fname)
	{
	FILE *fp;
	int nseq=0;
	
	fp=vfopen ( fname, "r");
	fp=set_fp_id ( fp, "Index");
	while ( (fgetc(fp))!='\n');
	while ( isspace(fgetc(fp)))
		{fscanf (fp, "%*d >%s", name[nseq++]);
		 while ( (fgetc(fp))!='\n');
		}
	vfclose ( fp);
	return nseq;
	}
Alignment * read_gotoh_aln ( char *fname, Alignment *A)
   {
    FILE *fp;
    int * ptr_aln;
    int a,b,d,e;


    char *buf;
    char buf2[VERY_LONG_STRING+1];
    char buf3[VERY_LONG_STRING+1];
    char buf4[VERY_LONG_STRING+1];

    int tot=0;

    int l;
    int nseq, max_len;
    
   
    if ( !check_file_exists (fname))return NULL;
    fp=vfopen ( fname, "r");

/*1 GET THE NUMBER OF SEQUENCES*/
    nseq=0;
    buf=calloc ( VERY_LONG_STRING+1, sizeof (char));    
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while (!isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while ( !isblanc ( buf) && buf!=NULL)
    	{
    	a=-1;
     	d=sscanf ( buf, "%d %s %s %s", &a, buf2, A->name[A->nseq],buf3);
    	if ( a!=-1)
    		{
		if ( name_is_in_list (A->name[A->nseq], A->name, A->nseq, 100)!=-1){fprintf ( stderr, "\nERROR: Sequence %s Duplicated in File %s", A->name[nseq], A->file[nseq]);exit(1);}		    
    		nseq++;
    		fgets(buf, VERY_LONG_STRING, fp);
    		}
    	else ( buf=NULL);
    	}
    vfclose (fp);
/*2 Get the MAX Len and Reallocate*/
    max_len=count_n_char_in_file(fname)/nseq;
    A=realloc_aln2( A, nseq+1, max_len+1);
/*3 Get The Sequences Names*/
    A->nseq=0;
    fp=vfopen ( fname, "r");
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while (!isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while ( !isblanc ( buf) && buf!=NULL)
    	{
    	a=-1;
     	d=sscanf ( buf, "%d %s %s %s", &a, buf2, A->name[A->nseq],buf3);
    	if ( a!=-1)
    		{
    		if ( d==4)sprintf (A->name[A->nseq],"%s", buf3); 	
    		A->nseq++;
    		fgets(buf, VERY_LONG_STRING, fp);
    		}
    	else ( buf=NULL);
    	}
    vfclose (fp);   

/*READ THE ALN*/     
    fp=vfopen ( fname, "r");

    buf=calloc ( VERY_LONG_STRING+1, sizeof (char));;	
    ptr_aln=calloc ( A->nseq, sizeof(int));
    
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    while (!isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    
    
    while ( isblanc (buf=fgets ( buf, VERY_LONG_STRING, fp)));
    
    while (buf!=NULL)
	{
	tot=0;
	while(tot< A->nseq)
	    {
	    
	    e=sscanf (buf, "%d %s %s %s", &e, buf2, buf3, buf4);
	    if ( e==4)sprintf( buf3, "%s", buf4);
	    
	    
	    for ( d=0; d< A->nseq; d++)
		{
		
		if ( strcmp (A->name[d], buf3)==0)
		    {a=d;
		     tot++;
		    }
		}
	     l=strlen (buf2);
	     if ( buf2[l-1]=='|')l--;
	     buf2[l]='\0';
	    
	     for (b=0; b<l; b++)
	     	{
	     	if ( isgraph (buf2[b]))
	     	 	A->seq_al[a][ptr_aln[a]++]=(A->residue_case==2)?buf2[b]:tolower (buf2[b]);
	     	 }
	     buf=fgets(buf, VERY_LONG_STRING, fp);	
	     }
	 if ( buf!=NULL)
	 	{
	 	buf=fgets(buf, VERY_LONG_STRING, fp);
	 	while ( isblanc (buf) && buf!=NULL)
	 		{
	 		buf=fgets ( buf, VERY_LONG_STRING, fp);
	 		}
	 	}
	 
	 }
	 
    vfclose (fp);
    
   
    for ( a=0; a< A->nseq; a++)
	{A->seq_al[a][ptr_aln[a]]='\0';
	}
    
    A->len_aln= strlen(A->seq_al[0]);  
    
   
    
    for ( a=0; a< A->nseq; a++)
    	{
    	for ( b=0; b< A->len_aln; b++)
    		A->len[a]+=1-is_gap(A->seq_al[a][b]);
    	}
    for ( a=0, b=0; a< A->len_aln; a++)
    	{
    	if ( !is_gap(A->seq_al[0][a]) &&!is_gap(A->seq_al[1][a]))b++;
    	}
    return A;
    }
    




void read_msf_aln ( char *fname, Alignment *A)
   {
    FILE *fp;
    int * ptr_aln;
    int a,b;
    int c;
    char buf[VERY_LONG_STRING+1];
    char *p;
    int nseq, max_len;

/*COUNT N SEQ*/
   nseq=0;
   fp=vfopen ( fname, "r");
   while ( (c=fgetc(fp))!='/');
   fgets ( buf, VERY_LONG_STRING, fp);
   while ( !is_alpha_line ( buf))
    	fgets ( buf, VERY_LONG_STRING, fp);
   while ( is_alpha_line (buf))
    	{
    	nseq++;
    	fgets(buf,VERY_LONG_STRING, fp);
    	}
   vfclose (fp);
/*COUNT MAX LEN AND REALLOCATE*/
   
   max_len=count_n_char_in_file(fname)/nseq;
   A=realloc_aln2( A, nseq+1, max_len+1);
/*GET SEQ NAMES*/
   A->nseq=0;
   fp=vfopen ( fname, "r");
   while ( (c=fgetc(fp))!='/');
   fgets ( buf, VERY_LONG_STRING, fp);
   while ( !is_alpha_line ( buf))
    	fgets ( buf, VERY_LONG_STRING, fp);
   while ( is_alpha_line (buf))
    	{
    	sscanf ( buf, "%s",A->name[A->nseq]);
	if ( name_is_in_list (A->name[A->nseq], A->name, A->nseq, 100)!=-1){fprintf ( stderr, "\nERROR: Sequence %s Duplicated in File %s", A->name[nseq], A->file[nseq]);exit(1);}   	
        A->nseq++;
        fgets(buf,VERY_LONG_STRING, fp);
    	}
   vfclose (fp);

    
/*GET THE ALIGNMENT*/     
    fp=vfopen ( fname, "r");
    
    
    ptr_aln=calloc ( A->nseq, sizeof(int));
    while ( (c=fgetc(fp))!='/');
    fgets ( buf, VERY_LONG_STRING, fp);
    while ( !is_alpha_line ( buf))
    	fgets ( buf, VERY_LONG_STRING, fp);
    
    b=0;
    
    p=buf;
    while (p!=NULL)
    	{
    	if ( is_alpha_line (buf))
    		{
    		
    		if (b==A->nseq)b=0;
    		a=0;
    		while ( isspace ( buf[a]) && buf[a]!='\0')a++;
    		while ( !isspace( buf[a]) && buf[a]!='\0')a++;
    		while ( !isgraph( buf[a]) && buf[a]!='\0')a++;
    		a--;
    		
    		while ( buf[a]!='\0')
    			{
    			if ( isalpha(buf[a]) || is_gap(buf[a]))
    				{
    				
    				c=(is_gap(buf[a]))?'-':(A->residue_case==2)?buf[a]:tolower(buf[a]);
    				
				A->seq_al[b][ptr_aln[b]++]=c;
    				}
    			a++;
    			}
    		b++;
    		}
    	p=fgets( buf, VERY_LONG_STRING, fp); 
    	}	
   
    
         
    for ( a=0; a< A->nseq; a++)
	{
	
	A->seq_al[a][ptr_aln[a]]='\0';
	A->order[a][0]=a;
	A->order[a][1]=0;	
	}
    
    A->len_aln= strlen(A->seq_al[0]);  
    
    
    for ( a=0; a< A->nseq; a++)
    	{
    	for ( b=0; b< A->len_aln; b++)
    		A->len[a]+=1-is_gap(A->seq_al[a][b]);
    	}
    for ( a=0, b=0; a< A->len_aln; a++)
    	{
    	if ( !is_gap(A->seq_al[0][a]) && !is_gap(A->seq_al[1][a]))b++;
    	}
    free (ptr_aln);
    }		

/**************************************************************************************************/
/*************************************REFORMAT OUT*************************************************/
/**************************************************************************************************/
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT MATRICES                                           */
/*                                                                                         */
/***************************************************************************************** */



int output_freq_mat ( char *outfile, Alignment *A)
    { /*
	function documentation: start
	
	int output_freq_mat ( char *outfile, Aligmnent *A)

	This function counts the number of residues in each column of an alignment (Prot)
	It outputs these values in the following format

	A | 0 0 0 1 0
	B | 1 0 0 0 1
	- | 0 1 1 0 0

	This format can be piped into:
	The routine used for computing the p-value  gmat-inf-gc-v2c
	
	function documentation: end
      */
      
    int a, b;
    int **freq_mat;
    FILE *fp;
    
    
    freq_mat=aln2count_mat (A);
            
    fp=vfopen ( outfile, "w");
    for ( b=0; b< 26; b++)
      {
	fprintf (fp, "%c |", 'A'+b);
	for ( a=0; a< A->len_aln; a++)fprintf (fp,"%d ", freq_mat[b][a]);
	fprintf (fp, "\n");
      }
    fprintf (fp, "- |");
    for ( a=0; a< A->len_aln; a++)fprintf (fp,"%d ", freq_mat[26][a]);
    
    free_int (freq_mat, -1);
    vfclose ( fp);
    return 1;
    }
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT P-Values                                           */
/*                                                                                         */
/***************************************************************************************** */	
float output_maln_pval ( char *outfile, Alignment *A)
    {
      /*
	function documentation: start
	float output_maln_pval ( char *outfile, Aligmnent *A)

	This function outputs the p-value of a multiple alignmnet as described 
	in Hertz, Stormo, Bioinformatics, 15-7/8, 563/577
	    ftp beagle.colorado.edu /pub/cosensus
	Locally
	    packages/consensus/gmat-inf-gc-v2c
	
	
	The routine used for computing the p-value is the program gmat-inf-gc-v2c
	function documentation: end
      */

  
    char *mat;
    char *result;
    FILE *fp;
    float value;
    char command[LONG_STRING];
    char string[STRING];
    mat=vtmpnam (NULL);
    result=vtmpnam (NULL);
    
    output_freq_mat (mat,A);
    sprintf ( command, "more %s | gmat-inf-gc-v2c -A abcdefghijklmnopqrstuvwxyz> %s",mat, result);
    system ( command);
    
    if ( !check_file_exists(result))return 0;
    fp=find_token_in_file ( result, NULL, "ln(p-value):");
    
    fscanf ( fp, "%s",string);
    value=atof ( string);
    vfclose ( fp);
    
    remove ( mat);
    remove ( result);
    vfree ( mat);
    vfree( result);
    
    fp=vfopen ( outfile, "w");
    fprintf ( fp, "%.6f\n", value);
    vfclose ( fp);
    
    return value;
    }
	      
    
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT WEIGHTS                                            */
/*                                                                                         */
/***************************************************************************************** */
int output_seq_weights ( Weights *W, char *wfile)
        {
	FILE*fp;
	int a;
	
	if ( W==NULL)return 0;
	
	fp=vfopen (wfile, "w");
	if ( fp==NULL)return 0;
	
	
	for ( a=0; a< W->nseq; a++)
		{
		
		  fprintf ( fp, "%s %.2f\n", W->seq_name[a],W->SEQ_W[a]);
		}
	vfclose ( fp);
	return 1;
	}  
void output_pw_weights4saga ( Weights *W, float **w_list, char *wfile)
	{
	FILE*fp;
	int a, b;
	fp=vfopen (wfile, "w");
	
	fprintf ( fp, "%s\n$\n", W->comments); 
	for ( a=0; a< W->nseq-1; a++)
		{
		for (b=a+1; b< W->nseq; b++)
			{
			fprintf ( fp, "%s %s %f\n", W->seq_name[a], W->seq_name[b],w_list[a][b]);
			}
		}
	fprintf ( fp, "$\n");
	vfclose ( fp);
	}

FILE * display_weights (Weights *W, FILE *fp)
{
  int a;
  int max_len;
  
  if ( W==NULL)
    {
      fprintf ( fp, "\n\nUN-WEIGHTED MODE: EVERY SEQUENCE WEIGHTS 1\n");
      return fp;
    }
  fprintf ( fp, "\n\nWEIGHTED MODE:%s\n\n", (W)->mode);
  for ( a=0, max_len=0; a< W->nseq; a++)max_len=MAX(max_len, strlen (W->seq_name[a]));
  for ( a=0; a< (W->nseq); a++)
    {
      fprintf ( fp, "\t%*s %.2f\n", max_len,(W)->seq_name[a],W->SEQ_W[a]);
    }
  fprintf ( fp, "\n");
  return fp;
}

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT SEQ                                                */
/*                                                                                         */
/***************************************************************************************** */
void output_statistics (char *file, Alignment *A, char *mode)
    {
      FILE *fp;
      int a, b, c, d, n;
      int maxname=0;

      
      if (!mode || !mode[0])
	mode="hnrgl";
      else if ( mode[0]=='_')
	mode++;
      for ( a=0; a<A->nseq; a++)maxname=MAX(strlen(A->name[a]), maxname);
      maxname++;
      
      
      fp=vfopen (file, "w");
      
      if (mode[0]=='h')
	{
	  b=0;
	  while ((c=mode[b++])!='\0')
	    {
	      if ( c=='n') fprintf (fp, "%-*s ",maxname,"name");
	      if ( c=='l') fprintf (fp, "%-*s ",5,"nres");
	      if ( c=='g') fprintf (fp, "%-*s ",5,"ngap");
	      if ( c=='t') fprintf (fp, "%-*s ",5,"len");
	    }
	  fprintf (fp, "\n");
	  mode++;
	}

      for (a=0; a<A->nseq; a++)
	{
	  b=0;
	  while ((c=mode[b++])!='\0')
	    {
	      if (c=='n')fprintf ( fp, "%-*s ", maxname,A->name[a]);
	      if (c=='l')
		{
		  for (n=0,d=0; d<A->len_aln; d++)n+=!is_gap(A->seq_al[a][d]);
		  fprintf ( fp, "%-5d ",n);
		}
	      if (c=='g')
		{
		  for (n=0,d=0; d<A->len_aln; d++)n+=((is_gap(A->seq_al[a][d]) && !is_gap(A->seq_al[a][d+1]))||(is_gap(A->seq_al[a][d])&& A->seq_al[a][d+1]=='\0')) ;
		  fprintf ( fp, "%-5d ",n); 
		}
	      if (c=='t')
		{
		 fprintf ( fp, "%-5d ",strlen (A->seq_al[a]));		 
		}
	      
	    } 
	  fprintf ( fp, "\n");
	}
      vfclose (fp);
    }

FILE * display_sequences_names (Sequence *S, FILE *fp)
        {
	    int a;
	    int max_len;

	    
	    if ( !S)
	       {
		   fprintf (fp,"\nERROR: NO SEQUENCE READ [FATAL]\n"); exit (1);
	       }
	    for ( a=0, max_len=0; a< S->nseq; a++)max_len=MAX(max_len, strlen (S->name[a]));
	    fprintf ( fp, "\nINPUT: %d SEQUENCES  [%s]", S->nseq,(S->type)?S->type:"Unknown type");
	    /*   for ( a=0; a< S->nseq; a++)
	        {
		    fprintf (fp, "\n\t%-*s Length= %4d", max_len,S->name[a],strlen ( S->seq[a]));
		}
		fprintf ( fp, "\n");*/
	    return fp;
	    
	}


void output_est_prf   (char *fname, Alignment *A)
        {
	int a;
	FILE *fp;

	if ( !A->P)
	  {
	    fprintf ( stderr, "\nFormat output_est_prf Impossible: No profile\n");
	    exit(1);
	  }
	

	fp=vfopen ( fname, "w");
	fprintf ( fp, "Consensus Sequence\nReconstructed with %s (%s,%s)\n",PROGRAM,AUTHOR,DATE);
	fprintf ( fp, "%4c %4c %4c %4c %15s    Consensus\n",  'A','G','C','T', "Internal Gaps");

	for ( a=0; a< A->len_aln; a++)
	  {
	    fprintf (fp, "%4d %4d %4d %4d %15d %c\n", (A->P)->count[0][a],(A->P)->count[1][a],(A->P)->count[2][a], (A->P)->count[3][a], (A->P)->count[4][a],A->seq_al[0][a]);
	  }
	return;
	}

	  
void output_gotoh_seq (char *fname, Alignment*A )
	{
	int a;
	FILE *fp;
		
	fp=vfopen ( fname, "w");
	fprintf ( fp, "%d %d\n",A->nseq, A->max_len);
	for ( a=0; a< A->nseq; a++)
		{
		ungap ( A->seq_al[a]);
		fprintf ( fp, ">%s\n", A->name[a]);
		fp=output_string_wrap ( 50,A->seq_al[a] , fp);
		fprintf ( fp, "//\n");
		}
		
	vfclose (fp);
	}	    

void output_mult_fasta_seq (char *fname, Alignment*A, int n )
	{
	int a;
	FILE *fp;
	
	fp=vfopen (fname, "w");
	ungap(A->seq_al[0]);
	for (a=0; a<n; a++)
	  {
	    fprintf (fp, ">%s_%d\n%s\n", A->name[0],a+1, A->seq_al[0]);
	  }
	vfclose (fp);
	}

void output_fasta_seq1 (char *fname, Alignment*A )
	{
	char seq_name[VERY_LONG_STRING];
	int a;
	FILE *fp;
	char *extension;
	
	for ( a=0; a< A->nseq; a++)
		{
		if ( strncmp( fname, "name",4)==0)
		  {
		    if ( (fname+4)[0]!='\0')extension=fname+5;
		    else
		      extension=NULL;
		    
		     sprintf ( seq_name,"%s.%s", A->name[a],(extension==NULL)?"seq":extension);
		  }
		else
		   sprintf ( seq_name,"%s%d.seq", fname,a);
		
		ungap ( A->seq_al[a]);
		fp=vfopen (seq_name, "w");
		fprintf (fp, ">%s\n", A->name[a]);
		fp=output_string_wrap ( 50, A->seq_al[a],fp);
		fprintf ( fp, "\n");
		vfclose (fp);
		}
	}
void output_pir_check (char *fname,int nseq, char **comment )
	{
	int a;
	FILE *fp;
	
	if ( fname==NULL)return;
	fp=vfopen ( fname, "w");
	
	for ( a=0; a< nseq; a++)fprintf (fp, "%s\n", comment[a]);
	vfclose (fp);
	}
void output_fasta_seq (char *fname, Alignment*A )
	{
	int a;
	FILE *fp;
	
	fp=vfopen ( fname, "w");
	
	for ( a=0; a< A->nseq; a++)
		{
		ungap(A->seq_al[a]);
		fprintf ( fp, ">%s\n", A->name[a]);
		fp=output_string_wrap ( 50, A->seq_al[a],fp);
		fprintf ( fp, "\n");
		}
	vfclose (fp);
	}    
void output_gor_seq (char *fname, Alignment*A )
	{
	int a;
	FILE *fp;
	
	fp=vfopen ( fname, "w");
	
	for ( a=0; a< A->nseq; a++)
		{
		ungap(A->seq_al[a]);
		fprintf ( fp, "!%s                                                               %d \n", A->name[a], strlen(A->seq_al[a]));
		upper_string ( A->seq_al[a]);
		fp=output_string_wrap ( 50, A->seq_al[a],fp);
		fprintf ( fp, "@\n");
		}
	vfclose (fp);
	}    
void output_pir_seq (char *fname, Alignment*A )
	{
	int a;
	for ( a=0; a< A->nseq; a++)ungap(A->seq_al[a]);
	output_pir_aln (fname, A);
	} 
void output_pir_seq1 (char *fname, Alignment*A )
	{
	char seq_name[VERY_LONG_STRING];
	int a;
	FILE *fp;
	char type[20];

	
	for ( a=0; a< A->nseq; a++)
		{
		if      ( strm ( get_string_type (A->seq_al[a]),"DNA"))sprintf(type, "DL");
		else if ( strm ( get_string_type (A->seq_al[a]),"PROTEIN"))sprintf(type, "P1"); 
		sprintf ( seq_name,"%s;%s_%s.seq",type, fname,A->name[a]);
		ungap ( A->seq_al[a]);
		fp=vfopen (seq_name, "w");
		fprintf (fp, ">%s\n\n", A->name[a]);
		fp=output_string_wrap ( 50, A->seq_al[a],fp);
		fprintf ( fp, "\n*\n");
		vfclose (fp);
		}
	} 
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT ALN                                                */
/*                                                                                         */
/***************************************************************************************** */
void output_mocca_aln (char *outfile, Alignment *A, Alignment *S)
    {
    FILE *fp;
    int **score;
    char **new_name_order;
    int a, maxl;

    score=declare_int (S->nseq, 2);
    new_name_order=declare_char ( S->nseq,MAXNAMES+1); 
    for ( a=0; a<A->nseq; a++)
      {
	score[a][0]=a;
	score[a][1]=S->score_seq[a];
      }
    sort_int_inv (score+1,2,1,0,S->nseq-2);
    for ( a=0; a<A->nseq; a++)
      {
	sprintf ( new_name_order[a], "%s", A->name[score[a][0]]);
      }
    A=reorder_aln (A, new_name_order, A->nseq);

    fp=vfopen (outfile, "w");
    fprintf ( fp, "MOCCA,(%s,%s, C. Notredame)\nSCORE %d\nNSEQ  %d\nLEN   %d\n",VERSION,DATE, A->score_aln, A->nseq, A->len_aln);     
    
    maxl=return_maxlen ( new_name_order, A->nseq); 
    
   
    for (a=0; a< A->nseq; a++)
      {
	fprintf (fp, "%-*s: %3d\n", maxl, A->name[a], score[a][1]);
      }
    
    fprintf ( fp, "\n");
    
    fp=output_Alignment_without_header ( A, fp);
    vfclose (fp);
    free_int  (score, -1);
    free_char (new_name_order, -1);
    return ;
    }
  
void print_aln ( Alignment *B)
    {
    output_Alignment_without_header ( B, stderr);
    }


FILE * output_aln ( Alignment *B, FILE *fp){return output_Alignment(B, fp);}
FILE * output_Alignment ( Alignment *B, FILE *fp)
    {
      fprintf ( fp, "%s, %s (%s)\n%s\nCPU   %d sec\nSCORE %d\nNSEQ  %d\nLEN   %d\n",PROGRAM,VERSION,DATE, AUTHOR,  (B->cpu+get_time())/1000, B->score_aln, B->nseq, B->len_aln);     
      return output_Alignment_without_header ( B, fp);
    }
  
FILE * output_Alignment_without_header ( Alignment *B, FILE *fp)
    {
    int a,b, c;
    int max_len=0;
    int line;	    
    int *n_residues;
    char s;

    if (fp==NULL)return fp;
    for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
    max_len=MAX(max_len+2, 16);
    line=60;
    n_residues=vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a<B->nseq; a++)n_residues[a]=B->order[a][1];
    
    
    
    
  fprintf ( fp, "\n"); 
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<=B->nseq; b++)
	     {
		
			fprintf (fp,"%-*s",max_len,B->name[b]);
			if (B->output_res_num)fprintf (fp, " %4d ", n_residues[b]+1);
			for (c=a;c<a+line && c<B->len_aln;c++)
			        {
				  if (b==B->nseq){n_residues[b]++;s=analyse_aln_column ( B, c);}
				  else 
				    {n_residues[b]+=!is_gap(B->seq_al[b][c]);
				     s=GET_CASE(B->residue_case, B->seq_al[b][c]);
				    }
				  
				  fprintf (fp,"%c",s );
			        }
			if (B->output_res_num)fprintf (fp, " %4d", n_residues[b]);
			fprintf (fp,"\n");
	     }
	   
	   fprintf (fp,"\n");
	   }
     fprintf (fp,"\n\n");
     free (n_residues);
    return fp;
    }
FILE * output_aln_score ( Alignment *B, FILE *fp){return output_Alignment_score(B, fp);}
FILE * output_Alignment_score ( Alignment *B, FILE *fp)
    {
    int a, b, c;
    static int max_len=0;
    static int line;	    
    int ch;
    
    if (fp==NULL)return fp;
    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	max_len+=4;
	line=60;
	}	
    
   sprintf (B->name[B->nseq], "CONS"); 
   fprintf ( fp, "T_COFFEE ALIGNMENT\nCPU TIME:%d sec.\n", (B->cpu+get_time())/1000);  
   fprintf ( fp, "SCORE=%d\n", B->score_aln);
   for ( a=0;a<B->nseq; a++)fprintf ( fp, "%s: %d\n", B->name[a], B->score_seq[a]);
   fprintf ( fp, "\n"); 
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	      {
	      fprintf (fp,"%-*s",max_len,B->name[b]);
	      for (c=a;c<a+line && c<B->len_aln;c++)
		{
		ch=B->seq_al[b][c];
		if (ch==NO_COLOR_RESIDUE)fprintf (fp,"-");
		else if ( ch==NO_COLOR_GAP)fprintf (fp,"*");
		else if ( ch<10 && ch>=0)fprintf (fp,"%d",ch);
		else if ( ch>10)fprintf (fp,"#");
		else if ( ch<0)fprintf  (fp,".");
		else fprintf (fp,"9");		
		}	      
	      fprintf (fp,"\n");	      
	      }
	    fprintf (fp,"\n");
	    fprintf (fp,"%-*s",max_len,B->name[b]);
	    for (c=a;c<a+line && c<B->len_aln;c++)
	      {
	      ch=B->seq_al[b][c];
	      if (ch==NO_COLOR_RESIDUE)fprintf (fp,"-");
	      else if ( ch==NO_COLOR_GAP)fprintf ( fp, "*");
	      else if ( ch<10 && ch>=0)fprintf (fp,"%d",ch);
	      else if ( ch>10)fprintf (fp,"#");
	      else if ( ch<0)fprintf (fp,".");
	      else fprintf (fp,"9");		
	      }	      
	    fprintf (fp,"\n\n\n");
	   }
    fprintf (fp,"\n\n");
    return fp;
    }
FILE * output_aln_with_res_number ( Alignment *B, FILE *fp){return  output_Alignment_with_res_number(B, fp);}
FILE * output_Alignment_with_res_number ( Alignment *B, FILE *fp)
    {
    int a, b, c;
    static int max_len=0;
    static int line;	    
    int**order;

    if (fp==NULL)return fp;
    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	max_len+=4;
	line=60;
	}	
   order=copy_int ( B->order,declare_int ( B->nseq, 2), B->nseq, 2);
    
   fprintf ( fp, "T_COFFEE ALIGNMENT\nCPU TIME:%d sec.\n", (B->cpu+get_time())/1000);     
   fprintf ( fp, "\n"); 
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {
	      fprintf (fp,"%-*s %3d %4d ",max_len,B->name[b], order[b][0], order[b][1] );
	      for (c=a;c<a+line && c<B->len_aln;c++)
		{
		order[b][1]+=1-is_gap(B->seq_al[b][c]);
		fprintf (fp,"%c",toupper(B->seq_al[b][c]) );
		}
	      fprintf (fp," %4d\n", order[b][1] );
	      }
	    fprintf (fp,"\n");
	    }
    fprintf (fp,"\n\n");

    free_int (order, -1);
    return fp;
    }

void output_constraints ( char *fname, char *weight_mode,Alignment *A)
	{
	FILE *fp;
	Constraint_list *CL;


	if ( !A->CL || strm ( weight_mode, "pdb"))
	   {
	       if (!A->S)
	          {
		      A->S=aln2seq(A);
		  }

	       CL=declare_constraint_list ( A->S, NULL, NULL, 0, NULL, NULL);
	       CL=aln2constraint_list (A,CL, weight_mode);
	   }
	else 
	   {
	       CL=(Constraint_list *)A->CL;
	   }
	
	fp=save_constraint_list ( CL, 0, CL->ne,fname, NULL, "lib",A->S);
	vfclose ( fp);
	    
	if ( (Constraint_list *)A->CL !=CL)free_constraint_list (CL);

	return;

	}
	    
void output_model_aln (char *fname, Alignment*A )
        {
	  FILE *fp;
	  int a;
	  Dp_Model *M;
	  Dp_Result *R;
	  char *string;
	  
	  if ( A->Dp_result==NULL)
	    {
	      fprintf ( stderr, "\nWARNING Could Not Output Model %s [%s]", fname, PROGRAM);
	    }
	  R=A->Dp_result;
	  M=R->Dp_model;

	  fp=vfopen ( fname, "w");
	  for (a=0; a<M->nstate; a++)
	    {
	      if (M->model_comments[a][0])fprintf ( fp, "#STATE %c: %s\n", 'a'+a, M->model_comments[a]);
	    }
	  string=vcalloc ( R->len+1, sizeof (char));
	  for (a=0; a<R->len; a++)string[a]=R->traceback[a]+'a';
	  fprintf ( fp, ">%s\n",fname);
	  fp=output_string_wrap ( 50,string, fp);
	  vfree(string);
	  fprintf ( fp, "\n");
	
	  vfclose (fp);
	  return;
	}
void output_fasta_aln (char *fname, Alignment*A )
	{
	FILE *fp;
	int a;

	fp=vfopen ( fname, "w");
	for ( a=0; a< A->nseq; a++)
		{
		fprintf ( fp, ">%s\n", A->name[a]);
		fp=output_string_wrap ( 50,A->seq_al[a] , fp);
		fprintf ( fp, "\n");
		}
	vfclose (fp);
	}
	
void output_pir_aln (char *fname, Alignment*A )
	{
	int a;
	FILE *fp;
	char type[20];
	
	

	
	
	fp=vfopen ( fname, "w");
	for ( a=0; a< A->nseq; a++)
		{
		if      ( strm ( get_string_type (A->seq_al[a]),"DNA"))sprintf(type, "DL");
		else if ( strm ( get_string_type (A->seq_al[a]),"PROTEIN"))sprintf(type, "P1");
		fprintf ( fp, ">%s;%s\n\n",type, A->name[a]);
		fp=output_string_wrap ( 50,A->seq_al[a] , fp);
		fprintf ( fp, "\n*\n");
		}
		
	vfclose (fp);
	}	    


void output_msf_aln (char *fname,Alignment *B)
        {
	int a, b, c;	
	char *seq;
	int *all_checks;
	int i,j;
	long grand_checksum;
	FILE *fp;
	int max_len;
	int line=50;
	int block=10;
	int c_block;
	char aa;

	for ( max_len=0,a=0; a< B->nseq; a++)max_len= MAX(strlen ( B->name[a]),max_len);


	max_len+=5;
	
	fp=vfopen (fname, "w");
	
	seq =vcalloc(B->len_aln,  sizeof(char));
	all_checks =vcalloc(B->nseq, sizeof(int));
	for ( i=0; i< B->nseq; i++)
	  {
	    for ( j=0; j<B->len_aln; j++)
	      {
		if ( is_gap(B->seq_al[i][j]))seq[j]='.';
		else seq[j]=B->seq_al[i][j]=toupper(B->seq_al[i][j]);
		
	      }
	    all_checks[i] = SeqGCGCheckSum(seq, (int)B->len_aln);
	  }
	grand_checksum = 0;
	for(i=0; i<B->nseq; i++) grand_checksum += all_checks[i];
	grand_checksum = grand_checksum % 10000;
	fprintf(fp,"PileUp\n\n");
	B=get_aln_type(B);
	fprintf(fp,"\n\n   MSF:%5d  Type: ",B->len_aln);
	if(strm ( (B->S)->type, "DNA"))
		fprintf(fp,"N");
	else
		fprintf(fp,"P");
	fprintf(fp,"    Check:%6ld   .. \n\n", (long)grand_checksum);
	for (i=0; i< B->nseq; i++)
	  {
	    fprintf ( fp, " Name: %s oo  Len:%5d  Check:%6ld  Weight:  %.3f\n", B->name[i], B->len_aln,(long)all_checks[i],1.00);
	  }
	fprintf(fp,"\n//\n\n");
	
	for (a=0; a<B->len_aln; a+=line)
	   {
	     fprintf ( fp,"\n\n"); 
	     for (b=0; b<B->nseq; b++)
	       {
		 fprintf (fp,"%-*s ",max_len,B->name[b]);
		 for (c_block=0,c=a;c<a+line && c<B->len_aln;c++)
		   {
		     if ( c_block==block)
			    {
			      fprintf (fp, " ");
			      c_block=0;
			    }
			c_block++;
		     aa=(is_gap(B->seq_al[b][c]))?'.': toupper(B->seq_al[b][c]);
		     fprintf (fp,"%c",aa );
		   }
		 if ( c_block==block)
			    {
			      fprintf (fp, " ");
			      c_block=0;
			    }
		 fprintf (fp,"\n");
		 
	       }
	   }
    	fprintf ( fp,"\n"); 		 
	vfclose ( fp);
	
	
	vfree(seq);
	vfree(all_checks);
	

	return;
} 
int SeqGCGCheckSum(char *seq, int len)
{
	int  i;
        long check;
        
        for( i=0, check=0; i< len; i++,seq++)
                check += ((i % 57)+1) * toupper(*seq);

        return(check % 10000);
}  
void old_output_msf_aln (char *fname,Alignment *B)
	{
	FILE *fp;
	static int *put_seq;
	int a, b, c;
	int line=50;
	char aa;
	char *buf;
	int max_len;
    	int seq_max_len;
    
    	for ( max_len=0,a=0; a< B->nseq; a++)max_len= MAX(strlen ( B->name[a]),max_len);
	for ( seq_max_len=0,a=0; a< B->nseq; a++)seq_max_len= MAX(strlen ( B->seq_al[a]),max_len);
	

	buf=vcalloc(seq_max_len+1, sizeof (int)); 
	
	if ( put_seq==NULL)
		put_seq= vcalloc ( B->nseq, sizeof (int));
	put_seq[0]=1;
	
	
	for ( b=1; b< B->nseq; b++)
		{
		sprintf ( buf, "%s", B->seq_al[b]);
		ungap(buf);
		put_seq[b]=( strlen (buf)>0)?1:0;
		}	
	
	fp=vfopen ( fname, "w");
	fprintf ( fp, "MSF: %d Type P Check: 5083 ..\n", B->len_aln);
	for ( a=0; a< B->nseq; a++)
		{
		if ( put_seq[a]==1)
			fprintf ( fp,"Name: %s\n",B->name[a]);
		}
	 fprintf ( fp, "//\n");
	 for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {
	     if ( put_seq[b]==1)
	     	{
	     	fprintf (fp,"%-*s ",max_len,B->name[b]);
	        for (c=a;c<a+line && c<B->len_aln;c++)
			{
			
			  
			    
			aa=(B->seq_al[b][c]=='-')?'.': toupper(B->seq_al[b][c]);
			fprintf (fp,"%c",aa );
			}
	      	fprintf (fp,"\n");
	      	}
	      }
	    fprintf (fp,"\n");
	    }
    	fprintf ( fp,"\n\n"); 		 
	vfclose ( fp);

	free (buf);
	free(put_seq);
	}
	
void output_saga_aln ( char *name, Alignment *B)
    {
    int a, b, c;
    FILE *fp;



    int max_len;
    int line=50;

    
    
    for ( max_len=0,a=0; a< B->nseq; a++)max_len= (strlen ( B->name[a])>max_len)?(strlen ( B->name[a])):max_len;
	    
	


    fp= vfopen ( name, "w");
    
    fprintf (fp, "\nSAGA FORMAT\nalignement  %s nseq=%d len=%d\n", name, B->nseq, B->len_aln);
       
    fprintf (fp, "\n\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {fprintf (fp,"%-*s ",max_len,B->name[b]);
	      for (c=a;c<a+line && c<B->len_aln;c++)
		{
		  fprintf (fp,"%c",(B->seq_al[b][c]) );
		}
	      fprintf (fp,"\n");
	      }
	    fprintf (fp,"\n");
	    }
    fprintf (fp,"\n\n");
    vfclose ( fp);
    }
void output_compact_aln ( char *name, Alignment *B)
    {
    int a, b, c;
    FILE *fp;
    int do_print=0;


    int max_len;
    int line=50;

    
    
    for ( max_len=0,a=0; a< B->nseq; a++)max_len= (strlen ( B->name[a])>max_len)?(strlen ( B->name[a])):max_len;
	    
	


    fp= vfopen ( name, "w");
    
    fprintf (fp, "\nSAGA FORMAT\nalignement  %s nseq=%d len=%d", name, B->nseq, B->len_aln);
    fprintf (fp, "\n\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {
	     
	     for ( do_print=0, c=a;c<a+line && c<B->len_aln;c++)
	       do_print+=1-is_gap(B->seq_al[b][c]);
	     if ( do_print>0)
	           {
		     fprintf (fp,"%-*s ",max_len,B->name[b]);
	     
	     
	     
		     for (c=a;c<a+line && c<B->len_aln;c++)
		       {
			 if ( is_gap(B->seq_al[b][c])&& B->seq_al[b][c]!='-' )fprintf (fp,"%c", '-');
			 else fprintf (fp,"%c",(B->seq_al[b][c]) );
		       }
		     fprintf (fp,"\n");
		   }
	      }
	    fprintf (fp,"\n");
	    }
    fprintf (fp,"\n\n");
    vfclose ( fp);
    }
void output_clustal_aln ( char *name, Alignment *B)
    {
    int a, b, c;
    FILE *fp;
    int max_len=0;
    int line;	
    int *n_residues;

    n_residues=vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    n_residues[a]=B->order[a][1];
	    }
    max_len=MAX(max_len+2, 16);
    line=60;
    
    

    fp= vfopen ( name, "w");
    
    fprintf (fp, "CLUSTAL FORMAT for %s %s, CPU=%.2f sec, SCORE=%d, Nseq=%d, Len=%d", PROGRAM, VERSION, (float)(B->cpu+get_time())/1000, B->score_aln, B->nseq, B->len_aln);
    fprintf (fp, "\n\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<=B->nseq; b++)
	     {
		 if (b!=B->nseq)
		    {
			fprintf (fp,"%-*s",max_len, B->name[b]);
			for (c=a;c<a+line && c<B->len_aln;c++)
			    {
				if ( is_gap(B->seq_al[b][c]))fprintf (fp,"%c", '-');
				else 
				    {
					n_residues[b]++;
					fprintf (fp, "%c", GET_CASE(B->residue_case, B->seq_al[b][c]));
			    
				    }
			    }
			if (B->output_res_num)fprintf (fp, " %d", n_residues[b]);
			fprintf (fp,"\n");
		    }
		 else if ( b==B->nseq)
		    {
		    fprintf (fp,"%-*s",max_len," ");		    
		    for (c=a;c<a+line && c<B->len_aln;c++)
			    {
			    fprintf ( fp, "%c", analyse_aln_column (B, c));
			    }
		    fprintf (fp,"\n");
		    }
	     }
	   fprintf (fp,"\n"); 
	   }
    
    fprintf (fp,"\n\n");
    free (n_residues);
    vfclose ( fp);
    }    
void output_phylip_aln ( char *name, Alignment *B)
    {
    int a, b, c;
    FILE *fp;

    int *print_name;
    static int line=50;	    
    
    
    
    print_name=calloc ( B->nseq, sizeof (int));
    fp= vfopen ( name, "w");
    
    fprintf (fp, "%d %d", B->nseq, B->len_aln);
    fprintf (fp, "\n");
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<B->nseq; b++)
	     {if ( print_name[b]==0)
	     	{fprintf (fp,"%10s",B->name[b]);
	     	 print_name[b]=1;
	     	 }
	          
	      for (c=a;c<a+line && c<B->len_aln;c++)
		{
		if ( is_gap(B->seq_al[b][c])&& B->seq_al[b][c]!='-' )fprintf (fp,"%c", '-');
		else fprintf (fp,"%c",(B->seq_al[b][c]) );
		}
	      fprintf (fp,"\n");
	      }
	    fprintf (fp,"\n");
	    }
    fprintf (fp,"\n\n");
    vfclose ( fp);
    }

void output_rnalign (char *out_file, Alignment *A, Sequence *STRUC)
    {
    int a, b;
    FILE *fp;
    char bank_file[100];
    char pep_file[100];
    char *buf;
    
    sprintf ( bank_file, "%s.mss", out_file);
    sprintf ( pep_file, "%s.one_rna", out_file);
    
   
    buf=calloc ( strlen ( A->seq_al[0]+1), sizeof (char));
    
    for ( b=0,a=0; a< strlen(A->seq_al[0]); a++) 
    	{
    	if ( is_gap(A->seq_al[0][a]))
    		buf[a]='.';
    	else
    		buf[a]=STRUC->seq[0][b++];
    	}
    buf[a]='\0';
    
    fp=vfopen ( bank_file, "w");
    
    fprintf ( fp, "ST\n");
    fp=output_string_wrap ( 50, buf, fp);
    fprintf ( fp, "\n\n");
    
    for ( a=0; a<A->nseq-1; a++)
    	{
    	fprintf ( fp, "AS %s\n ", A->name[a]);
    	fp=output_string_wrap ( 50, A->seq_al[a], fp); 
	fprintf ( fp, "\n\n");
	}
    vfclose ( fp);
    fp=vfopen ( pep_file, "w");
    fprintf ( fp, ">%s\n", A->name[A->nseq-1]); 
    fp=output_string_wrap ( 50, A->seq_al[A->nseq-1], fp);
    fprintf ( fp, "\n");
    vfclose (fp);
    }

void output_lib (char *pw_lib_saga_aln_name, Alignment *A )
    {
    Alignment *B;
    char fname[VERY_LONG_STRING];
    int a,b;
    
    B=declare_Alignment (NULL);
    
    B->nseq=2;
    
    for ( a=0; a< A->nseq-1; a++)
    	{
    	for ( b=a+1; b<A->nseq; b++)
    		{
    		sprintf ( B->seq_al[0], "%s", A->seq_al[a]);
    		sprintf ( B->name[0], "%s", A->name[a]);
    		sprintf(B->name[1], "%s", A->name[b]);
    		sprintf ( B->seq_al[1], "%s",A->seq_al[b]);
    		B->nseq=2;
    		sprintf ( fname, "%s_%s_%s.lib",pw_lib_saga_aln_name, A->name[a], A->name[b]);
    	
    		B->len_aln=strlen ( B->seq_al[0]);
    		ungap_aln (B);
    		output_clustal_aln (fname,B);  
        	}
        }
    } 
void output_pw_lib_saga_aln (char *pw_lib_saga_aln_name, Alignment *A )
    {
    Alignment *B;
    char fname[VERY_LONG_STRING];
    int a,b;
    
    B=declare_Alignment (NULL);
    
    B->nseq=2;
    
    for ( a=0; a< A->nseq-1; a++)
    	{
    	for ( b=a+1; b<A->nseq; b++)
    		{
    		sprintf ( B->seq_al[0], "%s", A->seq_al[a]);
    		sprintf ( B->name[0], "%s", A->name[a]);
    		sprintf(B->name[1], "%s", A->name[b]);
    		sprintf ( B->seq_al[1], "%s",A->seq_al[b]);
    		B->nseq=2;
    		sprintf ( fname, "%s_%s_%s.pw_lib_saga_aln",pw_lib_saga_aln_name, A->name[a], A->name[b]);
    	
    		B->len_aln=strlen ( B->seq_al[0]);
    		ungap_aln (B);
    		output_clustal_aln (fname,B);  
        	}
        }
    } 	
void output_lalign_header( char *name, Alignment *A)
    {
    FILE *fp;
    
    fp=vfopen ( name, "w");
    fprintf ( fp, " LALIGN_3D finds the best local alignments between two sequences\n");
    fprintf ( fp, " %s(%s)\n\n", VERSION, DATE);
    fprintf ( fp, " Comparison of:\n(A) %s\t%s\t-%d aa\n", (A->S)->file[A->order[0][0]],(A->S)->name[A->order[0][0]], (A->S)->len[A->order[0][0]]);
    fprintf ( fp, "(B) %s\t%s\t-%d aa\n", (A->S)->file[A->order[1][0]],(A->S)->name[A->order[1][0]], (A->S)->len[A->order[1][0]]);
    fprintf ( fp, " using structural similarity ONLY\n"); 
    
    vfclose ( fp);
    return;
    }

void output_lalign_aln   ( char *name, Alignment *B)
    {
    int a, b, c,d=0, s;
    char col;

    float tot=0;
    float id=0;

    FILE *fp;
    int max_len=0;
    int line;	
    int *n_residues;
    int res;

    static int output_header;

    if ( B==NULL){output_header=0;return;}
    else if ( output_header==0)
      {
	output_lalign_header(name, B);
	output_header=1;
      }
    
    
    if ( B->nseq>2)crash("\nlalign_aln can only output two sequences[FATAL]");
	
    n_residues=vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    n_residues[a]=B->order[a][1];
	    }
    max_len=MAX(max_len+2, 16);
    line=60;
    
    

    fp= vfopen ( name, "a");
    
    for (a=0; a< B->len_aln; a++)
      {
        if ( !is_gap(B->seq_al[0][a]) && !is_gap(B->seq_al[1][a]))
	     {
	       tot++;
	       id+=(B->seq_al[0][a]==B->seq_al[1][a]);
	     }
      }
    
    id=(id*100)/tot;
    fprintf (fp, " %.1f%% identity in %d aa overlap; score: %d\n\n", id,(int)tot, B->score_aln);
    
    
    for (a=0; a<B->len_aln; a+=line)
	   {for (b=0; b<5; b++)
	     {
	         if ( b==0 || b==4)
		   {
		     if ( b==0)s=0;
		     if ( b==4)s=1;
		     fprintf (fp,"%-*s",max_len," ");
		     for (d=0,c=a;c<a+line && c<B->len_aln;c++)
		       {
			 res=!is_gap ( B->seq_al[s][c]);
			 n_residues[s]+=res;
			 if ( (n_residues[s]%10)==0 && res && (c-a+4)<line){fprintf (fp, "%-4d", n_residues[s]);d=-3;}	
			 else
			   {
			     if ( d==0)fprintf (fp, " ");
			     else d++;
			   }
		       }
		     fprintf (fp,"\n");
		   }
		 else if (b==1 || b==3)
		    {
		      if ( b==1)s=0;
		      if ( b==3)s=1;
			fprintf (fp,"%-*s",max_len, B->name[s]);
			for (c=a;c<a+line && c<B->len_aln;c++)
			    {
				if ( is_gap(B->seq_al[s][c]))fprintf (fp,"%c", '-');
				else 
				    {
					fprintf (fp, "%c", GET_CASE(B->residue_case, B->seq_al[s][c]));
				    }
			    }
			fprintf (fp,"\n");
		    }
		 else if ( b==2)
		    {
		    fprintf (fp,"%-*s",max_len," ");		    
		    for (c=a;c<a+line && c<B->len_aln;c++)
			    {
			    col=analyse_aln_column (B, c);
			    if ( col=='*')col=':';
			    else if ( col==':')col='.';
			    else if ( col=='.')col=' ';
			    fprintf ( fp, "%c", col);
			    }
		    fprintf (fp,"\n");
		    }
	     }
	   fprintf (fp,"\n"); 
	   }
    
    fprintf (fp,"\n\n----------\n\n");
    free (n_residues);
    vfclose ( fp);
    }    	 
/****************************************************************************************************/
/*************************************UTIL *********************************************************/
/**************************************************************************************************/


/****************************************************************************************************/
/***************************                                    *************************************/
/***************************             PROCESSING 		*************************************/
/***************************                                    *************************************/
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              THREADING                                                  */
/***************************************************************************************** */

char *thread_aa_seq_on_dna_seq( char *s)
     {
	 int l, b, c;
	 char *array;
	 

	 l=strlen ( s);
	 array=vcalloc ( l*3 +1, sizeof (char));
	 for ( b=0, c=0; b< l; b++, c+=3)
	     {
		 array[c]=s[b];
		 array[c+1]='o';
		 array[c+2]='o';
	     }
	 array[c]='\0';
	 return array;
     }

Alignment *thread_dnaseq_on_prot_aln (Sequence *S, Alignment *A)
        {
	    Alignment *B=NULL;
	    int a, b, c, n, la, ls, ln, m;

	    B=copy_aln ( A, B);
	    B=realloc_aln2 ( B, B->nseq, B->len_aln*3 +1);

	    for ( n=0,a=0; a< A->nseq; a++)
	        {
		for ( m=0,b=0; b< S->nseq; b++)
		    {
		    if (strm (A->name[a], S->name[b]) )
		       {
			   m=1;
			   n++;
			   ungap ( S->seq[b]);
			   B->seq_al[a][0]='\0';
			   for (la=0, ls=0, ln=0; la< A->len_aln; la++)
			       {
				   for (c=0; c< 3; c++)
				       B->seq_al[a][ls++]=(is_gap(A->seq_al[a][la]))?'-':S->seq[b][ln++];
			       }
		       B->seq_al[a][ls]='\0';
		       }
		    }
		if ( m==0)
		       {
		       for (la=0, ls=0, ln=0; la< A->len_aln; la++)
			       {
				
				   B->seq_al[a][ls++]=A->seq_al[a][la];
				   B->seq_al[a][ls++]='-';
				   B->seq_al[a][ls++]='-';
			       }
		       }
		}
	    
	    B->len_aln=strlen ( B->seq_al[0]);
	    return B;
	}
void thread_seq_struc2aln ( Alignment *A, Sequence *ST)
	{
	int a, b, c,d;
	int len;
	
	for ( a=0; a< A->nseq; a++)
		for ( b=0; b< ST->nseq; b++)
			{
			if ( strcmp ( A->name[a], ST->name[b])==0)
				{
				ungap (ST->seq[b]);
				len=strlen(A->seq_al[a]);
				for ( c=0, d=0; c<len; c++)
					{
					if ( !is_gap(A->seq_al[a][c]))A->seq_al[a][c]=ST->seq[b][d++];
					}
				}
			}
	}
						
void cache_id ( Alignment *A)
	{
	int a, b,n;
	char r1, r2, r3;
	
	for ( a=0; a< A->len_aln; a++)
		{
		for ( b=0, n=0; b< A->nseq; b++)if ( !is_gap(A->seq_al[b][a]))n++;
		for ( b=0; b< A->nseq; b++)
			if ( !is_gap(A->seq_al[b][a]) && n==A->nseq)A->seq_al[b][a]='h';
			else if( !is_gap(A->seq_al[b][a]))A->seq_al[b][a]='x';
		}
	for ( a=0; a< A->nseq; a++)
		{
		for ( b=1; b< A->len_aln-1; b++)
			{
			r1=A->seq_al[a][b-1];
			r2=A->seq_al[a][b];
			r3=A->seq_al[a][b+1];
			if (r2=='h')
				{
				if ( (r1=='h' || r1=='b') && (r3=='h' || r3=='b'))A->seq_al[a][b]='h';
				else A->seq_al[a][b]='b';
				}
			}
		for ( b=1; b< A->len_aln-1; b++)if ( A->seq_al[a][b]=='b')A->seq_al[a][b]='x';
		}
			
	}
			
									 
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               PROCESING OF EST                                          */
/*                                                                                         */
/***************************************************************************************** */
int process_est_sequence ( Sequence *S, int *cluster_list)
	{
	char **inverted_seq;
	int T=20;
	int a, b;
	int V1, V2;
	int **sens;
	int **a_sens;
	int **best;
	int *solution;
	char buf [VERY_LONG_STRING];
	int n_clusters=0;
	int n;
	
	sens=declare_int ( S->nseq,S->nseq);
	a_sens=declare_int ( S->nseq,S->nseq);
	best=declare_int ( S->nseq,S->nseq);
	
	
	inverted_seq=calloc ( S->nseq, sizeof (char*));
	for ( a=0; a<S->nseq; a++)
		inverted_seq[a]=invert_seq ( S->seq[a]);
		
	for ( a=0; a< S->nseq-1; a++)
		{
		
		for ( b=a+1; b<S->nseq; b++)
			             {
			             
			             V1=sens[a][b]=sens[b][a]=get_best_match ( S->seq[a], S->seq[b]);
			             V2=a_sens[a][b]=a_sens[b][a]=get_best_match ( S->seq[a],inverted_seq[b]);
				     best[a][b]=best[b][a]=(V1>V2)?V1:V2;
				     }
		}
	solution=SHC ( S->nseq, a_sens, sens);
	
	
	for ( a=0; a<S->nseq; a++)cluster_list[a]=-1;
	for ( a=0; a<S->nseq; a++)
		{
		n=search_for_cluster (a, n_clusters, cluster_list, T, S->nseq, best);
		if ( n>0)n_clusters++;
		}
	fprintf ( stderr, "\nTHERE %s %d Independant Cluster(s) in your sequences",(n_clusters>1)?"are":"is",(n_clusters));
	for (a=0; a<n_clusters; a++)
		{
		fprintf (stderr, "\n");
		for ( b=0; b<S->nseq; b++)
			{
			if ( cluster_list[b]==a)fprintf ( stderr, "%s ", S->name[b]);
			}
		}
		
	for ( a=0; a<S->nseq; a++)
		{
		if ( solution[a]==-1)
			{
			S->seq[a]=inverted_seq[a];
			sprintf ( buf, "i_%s", S->name[a]);
			sprintf ( S->name[a], "%s", buf);
			}
		}
	return n_clusters;
	}

int search_for_cluster ( int seq, int cluster_number, int *cluster_list, int T, int nseq, int **S)	
	{
	int n=0,a;
	
	if (cluster_list[seq]==-1)
		{
		cluster_list[seq]=cluster_number;
		n++;
		}
	for ( a=0; a<nseq; a++)
		if ( cluster_list[a]==-1)
			{
			
			if (S[seq][a]>T)
				{
				n++;
				cluster_list[a]=cluster_number;
				n+=search_for_cluster ( a, cluster_number, cluster_list, T, nseq, S);
				}
			}
	return n;
	}	
		
int * SHC ( int nseq, int **NST, int **ST)
	{
	int a;
	int mut;
	int score, new_score;
	int N_IT=VERY_LONG_STRING;
	int *sol;
	int count;
	
	sol=calloc ( nseq, sizeof (int));
	for ( a=0; a<nseq; a++)
		sol[a]=(addrand ((unsigned long)100)>49)?1:-1;
		
	score=evaluate_sol (sol, nseq, ST, NST);
	fprintf ( stderr, "\nI_Score=%d\n", score);
	N_IT=N_IT*nseq;
	
	for ( count=0,a=0; a< N_IT && score<VERY_LONG_STRING; a++, count++)
		{
		mut=mutate_sol ( sol,nseq);
		new_score=evaluate_sol (sol, nseq, ST, NST);
		if ( new_score>score)
			{
			score=new_score;
			}
		else if ( (addrand ((unsigned long)VERY_LONG_STRING))>score)
			{
			score=new_score;
			}
		else
			sol[mut]=sol[mut]*-1;
		if ( count==VERY_LONG_STRING)
			{
			count=0;
			fprintf ( stderr, "\nScore=%d", score);
			}	
		}
	fprintf ( stderr, "\nScore=%d\n", score);
	return sol;
	}

int mutate_sol (int *sol, int nseq)
	{
	int n;
	n=addrand ((unsigned long)nseq);
	sol[n]=sol[n]*-1;
	return n;
	}
int evaluate_sol ( int *sol, int nseq, int **ST, int **NST)
	{
	static int max_score;
	int a, b, score=0;
	
	if ( max_score==0)
		{
		for ( a=0; a<nseq-1; a++)
			for ( b=a+1; b<nseq; b++)
				{
				max_score+=(ST[a][b]>NST[a][b])?ST[a][b]:NST[a][b];
				}
		}
	
	for ( a=0; a<nseq-1; a++)
		for (b=a+1; b<nseq; b++)
			if ( (sol[a]*sol[b])<0)score+=NST[a][b];
			else score+=ST[a][b];
	return (score*VERY_LONG_STRING)/max_score;
	}		
		
     
char * invert_seq ( char *seq)
	{
	int a, b;
	
	char *nseq;
	int l;
	
	
	l=strlen ( seq);
	for ( a=0; a<l; a++)
		seq[a]=tolower ( seq[a]);
	nseq=calloc ( l+1, sizeof (char));
	
	for ( a=0, b=l-1; a<l; a++, b--)
		{
		if (seq[b]=='n')nseq[a]='n';
		else if (seq[b]=='g')nseq[a]='c';
		else if (seq[b]=='c')nseq[a]='g';
		else if (seq[b]=='a')nseq[a]='t';
		else if (seq[b]=='t')nseq[a]='a';
		}
		
	nseq[l]='\0';
	return nseq;
	}
		
		
int get_best_match ( char *seq1, char *seq2)
	{
	static int **m;
	static int ml;
	int a, b;
	int **mdiag;
	int n_mdiag=0;
	int best;
	int l1, l2;
	
	
	l1=strlen ( seq1);
	l2=strlen (seq2);
	if ( m==NULL)
		{
		ml=(l1>l2)?l1:l2;
		m=declare_int (ml, ml);
		}
	else if ( (ml<l1) || (ml<l2))
		{
		free_int (m, ml);
		ml=(l1>l2)?l1:l2;
		m=declare_int (ml, ml);
		}
		
	for ( a=0; a<l1; a++)
		{
		for ( b=0; b<l2; b++)
			m[a][b]=((seq1[a]==seq2[b])|| seq1[a]=='n' ||seq2[b]=='n')?1:0;
		}
	mdiag= extract_m_diag_streches ( m, l1, l2,seq1, seq2, &n_mdiag); 
	
	for ( best=0,a=0; a<n_mdiag; a++)
		best=(mdiag[a][0]>best)?mdiag[a][0]:best;
	
	return best;
	}

int** extract_m_diag_streches ( int ** m, int l1, int l2,char *seq1, char *seq2, int *n_mdiag)
	{
	
	int b, x, y, s1, s2;
	static int **mdiag;
	int in;
	static int max_diag=VERY_LONG_STRING;
	 
	 /*
	 diag[0]=len;
	 diag[1]=x_start;
	 diag[2]=y_start;
	 diag[3]=x_end;
	 diag[4]=y_end;
	 */
	 
	if ( mdiag==NULL)
		mdiag=declare_int ( max_diag, 5);
	
	for ( s1=l1-1, s2=0;s2<l2;)
		{
		for ( in=0,x=s1, y=s2; x<l1 && y<l2; x++, y++)
			{ 
			if (m[x][y]>0)
				{
				if (in==1)
					mdiag[n_mdiag[0]][0]++;
				else 
					{
					mdiag[n_mdiag[0]][0]=1;
					mdiag[n_mdiag[0]][1]=x;
					mdiag[n_mdiag[0]][2]=y;
					in=1;
					}
				}
			else
				if (in==1)
					{
					in=0;
					mdiag[n_mdiag[0]][3]=x-1;
					mdiag[n_mdiag[0]][4]=y-1;
					if ( !is_strech ( "ta", seq1, seq2,mdiag[n_mdiag[0]][0], mdiag[n_mdiag[0]][1],mdiag[n_mdiag[0]][2]))n_mdiag[0]++;
					}
			if (n_mdiag[0]==(max_diag-1))
				{mdiag=realloc (mdiag, (max_diag+VERY_LONG_STRING)*sizeof (int*));
				for ( b=max_diag; b<max_diag+VERY_LONG_STRING; b++)mdiag[b]=calloc ( 5, sizeof (int));
				max_diag+=VERY_LONG_STRING;
				}
			}
		s2+= (s1==0)?1:0;	
		s1-= (s1==0)?0:1;
		if (in==1)
			{
			in=0;
			mdiag[n_mdiag[0]][3]=x-1;
			mdiag[n_mdiag[0]][4]=y-1;
			if ( !is_strech ( "ta", seq1, seq2,mdiag[n_mdiag[0]][0], mdiag[n_mdiag[0]][1],mdiag[n_mdiag[0]][2]))n_mdiag[0]++;
			}
		}
	
	return mdiag;
	}				
int is_strech ( char *AA, char *seq1, char *seq2, int len, int x, int y)
	{
	int n, i, j, c,a,nr;
	int T=70;
	
	n=strlen ( AA);
	for ( a=0; a<n; a++)
		{
		for (nr=0, i=x, j=y, c=0; c<len; c++, i++, j++)
			if ((seq1[i]==AA[a]) && (seq2[j]==AA[a]))nr++;
		if ( ((nr*100)/len)>T)return 1;
		}
	return 0;
	}
	
	
/************************************************************************************/
/*                                                                                  */
/*                                      DNA                                         */
/*                                                                                  */
/*                                                                                  */
/************************************************************************************/

Alignment *code_dna_aln (Alignment *A)
       {
	 int a, b,l,r;

	 for ( a=0; a< A->nseq; a++)
	   {
	     for (l=0, b=0; b< A->len_aln; b++)
	       {
		 r=A->seq_al[a][b];
		 if ( r=='-')l++;
		 else if ( r=='~')continue;
		 else if ( r=='.')l++;
		 else if ( !islower(r))A->seq_al[a][b]='4';
		 else
		   {
		     A->seq_al[a][b]=(l+3)%3+'0';
		     l++;
		   }
	       }
	   }
	 return A;
       }

  
Alignment *back_translate_dna_aln (Alignment *A)
       {
	 /*Given a set of aligned sequences
	   starts from left to right
	   1 aa->3 nuc
	   ambiguities are randomly resolved.
	   returns the corresponding amino acid alignment
	 */	 
	  int a;
	  char *seq    ;
	 
	 ungap_aln(A);
	 A=realloc_aln (A, 10000);
	 seq=vcalloc ( 10000, sizeof (char));

	 
	 for ( a=0; a< A->nseq; a++)
	   {
	   seq=back_translate_dna_seq (A->seq_al[a], seq, RANDOM);
	   sprintf ( A->seq_al[a], "%s", seq);
	   }
	 A->len_aln=A->len_aln*3;
	 compress_aln (A);
	 vfree (seq);
	 return A;
       }
char * back_translate_dna_seq ( char *in_seq,char *out_seq, int mode)     
       {
	 int a,len;

	 len=strlen(in_seq);
	 
	 if (out_seq==NULL)out_seq=vcalloc ( len*3+1, sizeof (char));
	 
	 out_seq[0]='\0';
	 for (a=0; a<len; a++)
	   {
	   strcat (out_seq,  back_translate_dna_codon (in_seq[a],mode));
	   }
	 
	 return out_seq;
       }
	   
Alignment *translate_dna_aln (Alignment *A, int frame)
       {
	 /*Given a set of aligned sequences
	   starts from left to right
	   3 nuc->1 aa
	   2nuc+1gap, 1nuc+2gap->3 gaps
	   1 stop-> 3gaps
	   returns the corresponding amino acid alignment
	 */


	 int a, b,r;

	 
	 for ( a=0; a< A->nseq; a++)
	   for (b=0; b< frame; b++)
	     A->seq_al[a][b]='-';
	 ungap_aln(A);
	 print_aln (A);

	 for ( b=0; b< A->nseq; b++)
	   for ( a=0; a< A->len_aln;)
	     {
	       
	       r=translate_dna_codon (A->seq_al[b]+a, 'z');
	       if (is_gap(r))
		 {
		   A->seq_al[b][a++]='-';
		   A->seq_al[b][a++]='-';
		   A->seq_al[b][a++]='-';
		 }
	       else if ( r=='x')
		 {
		   A->seq_al[b][a++]='o';
		   A->seq_al[b][a++]='-';
		   A->seq_al[b][a++]='-';
		 }
	       else if ( r=='z')
		 {
		   A->seq_al[b][a++]='x';
		   A->seq_al[b][a++]='-';
		   A->seq_al[b][a++]='-';
		 }
	       else
		 {
		   A->seq_al[b][a++]=r;
		   A->seq_al[b][a++]='-';
		   A->seq_al[b][a++]='-';
		 }
	     }	 
	 compress_aln (A);
	 
	 return A;
       }

Alignment *clean_gdna_aln (Alignment *A)
       {
	   int a, b, c, r1, r2,s, p, n, tn;
	   int *col;
	   static int **mat;
	   Alignment *T=NULL;
	   int **score;
	   char *buffer;
	   

	   /*Viterbi Parameters*/
	   int AL=0;        /*Allowed Transition*/
	   int F=-1000000; /*Forbiden Transition*/
	   int SPLICE_PENALTY=100;
	   int ORF1=0, ORF2=1, ORF3=2, NC=3;
	   
	   int state, pstate, best_e, best_pstate_p,best_state_p, best_pstate_v, best_state_v, v;
	   int nstate=4;
	   int **transitions;
	   int e;
	   int **v_tab_p;
	   int **v_tab;
	   int * is_dna;
	   

	   buffer=vcalloc ( 100000, sizeof (char));
	   is_dna=vcalloc ( A->nseq, sizeof (int));
	   score=declare_int ( A->nseq+1, A->len_aln);


	   if ( !mat)mat=read_matrice("pam250mt");
	   T=copy_aln (A, T);
	   col=vcalloc ( A->nseq, sizeof (int));
	   
	   for (a=0; a<= A->len_aln; a++)
	       for ( b=0; b< A->nseq; b++){A->seq_al[b][a]=tolower(A->seq_al[b][a]); A->seq_al[b][a]=(A->seq_al[b][a]=='t')?'u':A->seq_al[b][a];}

	   for ( a=0; a< A->nseq; a++)
	       {
		   sprintf ( buffer, "%s", A->seq_al[a]);
		   ungap (buffer);
		   is_dna[a]=strm ( get_string_type (buffer), "DNA");
	       }
	   

           for (a=0; a< A->len_aln-2; a++)
	       {
	       for (b=0; b< A->nseq; b++)
		       {
		       if (is_dna[b])col[b]=translate_dna_codon (A->seq_al[b]+a, 'x');
		       else col[b]=tolower ( A->seq_al[b][a]);   
		       }

	       for (n=0,tn=0,b=0; b< A->nseq; b++)
		   for ( c=b; c< A->nseq; c++   )
		       {
			   r1=col[b];
			   r2=col[c];
			   
			   if (r1=='x' || r2=='x'){score[A->nseq][a]=F;break;}
			   else if (r1=='-' && r2=='-');
			   else if (r1=='-' || r2=='-');
			   else 
			       {
			       
				   if ( is_dna[b] && is_dna[c])score[A->nseq][a]+= mat[r1-'a'][r2-'a'];
				   else score[A->nseq][a]+=mat[r1-'a'][r2-'a']* (A->nseq*A->nseq);
			       }
			   n+=( !is_gap(r1) && !is_gap(r2));
			   score[A->nseq][a]=(((tn!=0)?score[A->nseq][a]/tn:0));
		       }
	       
	       }

	   /*initialisation*/

	   transitions=declare_int ( nstate, nstate);
	   v_tab=declare_int ( A->len_aln+2, nstate       );
	   v_tab_p=declare_int ( A->len_aln+2, nstate       );

	   for (a=0; a<nstate;a++)
	       for (b=0; b<nstate;b++)
	             {transitions[a][b]=F;}

	   transitions[ORF1][ORF2]=AL;
	   transitions[ORF2][ORF3]=AL;
	   transitions[ORF3][ORF1]=AL;	   
	   
	   transitions[ORF3][NC]  =AL-SPLICE_PENALTY;
	   transitions[NC][ORF1]  =AL-SPLICE_PENALTY;


	   for ( s=0; s<A->nseq; s++)
	       {
	       for ( p=0; p<=A->len_aln; p++){for (state=0; state< nstate; state++)v_tab_p[p][state]=-1; }
	       for (p=1+2; p<= A->len_aln; p++)
	           {

		   for (state=0; state< nstate; state++)
		       {
			   
			   if ( state==NC){e=-best_e;}
			   else
			      {
				  e=score[A->nseq][(p-1)-state];
				  if ( state==0)best_e=e;
				  else best_e=MAX(e, best_e);
			      }

			   for ( pstate=0; pstate<nstate; pstate++)
		               {
				   v=e+transitions[pstate][state]+v_tab[p-1][pstate];
				   if (pstate==0 ||(v>best_pstate_v) )
				      {
				       best_pstate_v=v;
				       best_pstate_p=pstate;
				      }
			       }
			
			   v_tab[p][state]=best_pstate_v;
			   v_tab_p[p][state]=best_pstate_p;
			   if (state==0 ||best_pstate_v>best_state_v )
			      {
			       best_state_p=state; 
			       best_state_v=best_pstate_v;
			      }
		       }

		   }

	       
       
	       for (p=0; p< A->len_aln; p++)T->seq_al[s][p]='.';
	       for (p=A->len_aln; p>0; p--)
	           {
		       
		       if ( best_state_p==0)T->seq_al[s][p-1]=translate_dna_codon (A->seq_al[s]+(p-1), 'x');
		       else if ( best_state_p==1 || best_state_p==2)T->seq_al[s][p-1]='-';
		      
		      
		       
		       best_state_p=v_tab_p[p][best_state_p];
		       
		   }
	       }
	   
	   

	   vfree (col);
	   return T;
       }

Alignment *clean_cdna_aln (Alignment *A)
       {
	 /*Given an alignmnet of nucleotides
	   Returns the same alignmnent whith non coding nucleotides replaced with dots
	   
	   at each position, the emission probability is the sum of pair of the substitution of amino-acids
	 */
	 
	   int a, b, c,s, p;
	   static int **mat;
	   int   *emission;
	   float em1, em2;
	   char *buffer;
	   Alignment *B=NULL;


	   

	   /*Viterbi Parameters*/
	   int AL=0;        /*Allowed Transition*/
	   int F=-1000000; /*Forbiden Transition*/
	   int PENALTY=30;
	   int NC, C1,C2, C3, START, END;
	   int nstate=0;

	   int state,best_state, score, best_score;
	   int p_state;
	   int e;
	   int **score_tab;
	   int **state_tab;
	 
	   int **transitions;
	   int n;
	   int r1, r2, r3;

	   NC=nstate++;
	   C1=nstate++;
	   C2=nstate++;
	   C3=nstate++;
	   START=nstate++;
	   END=nstate++;

	   
	   B=copy_aln (A, B);
	   buffer=vcalloc ( 100000, sizeof (char));
	   emission=vcalloc (A->len_aln, sizeof (int));

	   if ( !mat)
	     {
	       mat=read_matrice("pam250mt");
	     }

	   /*Computation of the emission proba for the coding state*/


	   for (a=0; a< A->len_aln; a++)
	     {

	       /*First component: % occupancy of the column*/
	       em1=0;
	       for ( b=0; b< A->nseq; b++) em1+=!is_gap(translate_dna_codon (A->seq_al[b]+a, '-'));
	       em1=em1/(float)A->nseq;
	       
	       /*Second Component: % similarity within column*/
	       em2=0;
	       for (n=0,b=0; b< A->nseq-1; b++)
		 {
		   r1=translate_dna_codon (A->seq_al[b]+a, '-');
		   
		   for (c=b+1; c<A->nseq; c++)
		     {
		       r2=translate_dna_codon (A->seq_al[c]+a, '-');
		       if (is_gap(r2) || is_gap(r1));
		       else
			 {
			   n++;
			   em2+=((mat[r1-'a'][r2-'a'])>1)?1:0;
			 }
		     }
		 }
	       em2=em2/(float)((n==0)?1:n);
	       
	       
	       emission[a]=(em1*100);

	     }
	   
	 

	   /*initialisation*/

	   transitions=declare_int ( nstate, nstate);
	   score_tab=declare_int ( A->len_aln+2, nstate       );
	   state_tab=declare_int ( A->len_aln+2, nstate       );

	   for (a=0; a<nstate;a++)
	       for (b=0; b<nstate;b++)
	             {transitions[a][b]=F;}

	   
	   transitions[START][C1]=AL;
	   transitions[START][NC]=AL;
	   transitions[C3][END]=AL;
	   transitions[NC][END]=AL;
	   transitions[C1 ][C2 ]=AL;
	   transitions[C2 ][C3 ]=AL;
	   transitions[C3 ][C1 ]=AL;
	   transitions[C3 ][NC ]=AL-PENALTY;
	   transitions[NC ][C1 ]=AL-PENALTY;
	   transitions[NC][NC]=AL-PENALTY;
	   
	   
	   	   
	   for ( s=0; s< A->nseq; s++)
	     {
	     for ( p=0; p<=A->len_aln; p++){for (state=0; state< nstate; state++){score_tab[p][state]=F;state_tab[p][state]=-1;} }
	     score_tab[0][START]=0;
	     
	     for (p=1; p<= A->len_aln; p++)
	       {
		 for (state=0; state< nstate; state++)
		   {
		     if ( state==START || state==END)continue;
		     else if      ( state==NC)  e=-10;
		     else if ( state==C1)
		       {
			 e=emission[p-1];
		       }
		     else if ( state ==C2)
		       {
			 if ( p-2<0)e=F;
			 else e=emission[p-2];
		       }
		     else if ( state==C3)
		       {
			 if ( p-3<0)e=F;
			 else e=emission[p-3];
		       }
		     
		     for (p_state=0; p_state<nstate; p_state++)
		       {
			 
			 if (e==F)score=F;
			 else 
			   {
			     score=(score_tab[p-1][p_state]==F)?F:(e+transitions[p_state][state]+score_tab[p-1][p_state]);
			   }
			 
			 if(p_state==0 || score>best_score){ best_score=score;best_state=p_state;}
			 
		       }
		     
		     score_tab[p][state]=best_score;
		     state_tab[p][state]=best_state;
		     
		   }
	       }
	     
	     best_score=best_state=UNDEFINED;
	     for (state=0; state<nstate; state++)
	       {
		 if (state==START || state==END)continue;
		 e=transitions[state][END];
		 if (e==F || score_tab[p-1][state]==F)continue;
		 
		 if (best_score==UNDEFINED || score_tab[p-1][state]>best_score)
		   {
		     best_score=score_tab[p-1][state]+e; 
		     best_state=state;
		   }
		 
	       }
	     
	     for (p=A->len_aln; p>0;)
	       {
		 B->seq_al[s][p-1]=best_state+'0';
		 best_state=state_tab[p][best_state];
		 p--;
	       }
	     }

	   for ( a=0; a< A->nseq; a++)
	     for ( b=0; b< A->len_aln;)
	       {
		 s=B->seq_al[a][b];
		 if ( s==C1+'0')
		   {
		     r1=A->seq_al[a][b];
		     r2=A->seq_al[a][b+1];
		     r3=A->seq_al[a][b+2];
		     

		     if ( is_gap(r1) ||is_gap(r2) ||  is_gap(r3))
		       {
			 A->seq_al[a][b]=(is_gap(r1))?'~':'.';
			 A->seq_al[a][b+1]=(is_gap(r2))?'~':'.';
			 A->seq_al[a][b+2]=(is_gap(r3))?'~':'.';
		       }
		     b+=3;
		   }
		 else if ( s==NC+'0')
		   {
		     A->seq_al[a][b]=(is_gap(A->seq_al[a][b]))?'~':'.';
		     b++;
		   }
		 else 
		   {
		     fprintf (stderr, "\nPROBLEM: [%d %d]->%d", a, b, s-'0');
		   }
	       }
		 

	   free_aln (B);
	   free_int (transitions, -1);
	   free_int (score_tab, -1);
	   free_int (state_tab, -1);
	   vfree (emission);
	   vfree (buffer);
	   
	   return A;
       }




Alignment *translate_splice_dna_aln (Alignment *A, Alignment *ST)
       {
	   int a, b, c, r1, r2,s, p, n, tn;
	   int *col;
	   static int **mat;
	   Alignment *T=NULL;
	   int **score;
	   
	   /*Viterbi Parameters*/
	   int AL=0;        /*Allowed Transition*/
	   int F=-1000000; /*Forbiden Transition*/
	   int ORF1=0, ORF2=1, ORF3=2,SPL1=3, SPL2=4, SPL3=5, SPL4=6, NC=7;
	   int SPLICE_PENALTY;
	   int frame1, frame2, frame3, best_frame;
	   int nstate=8;
	   char r;
	   


	   int state, pstate, best_pstate_p,best_state_p, best_pstate_v, best_state_v, v;
	   
	   int **transitions;
	   int e;
	   int **v_tab_p;
	   int **v_tab;

	   score=declare_int ( A->nseq+1, A->len_aln);


	   if ( !mat)mat=read_matrice("pam250mt");
	   T=copy_aln (A, T);
	   col=vcalloc ( A->nseq, sizeof (int));
	   
	   for (a=0; a<= A->len_aln; a++)
	       for ( b=0; b< A->nseq; b++){A->seq_al[b][a]=tolower(A->seq_al[b][a]); A->seq_al[b][a]=(A->seq_al[b][a]=='t')?'u':A->seq_al[b][a];}

	   
	   

	   for (a=0; a< A->len_aln-2; a++)
	       {
	       for (b=0; b< A->nseq; b++)
		       {
		       col[b]=translate_dna_codon (A->seq_al[b]+a, 'x');
		       }
	       
	       for (n=0,tn=0,b=0; b< A->nseq-1; b++)
		   for ( c=b+1; c< A->nseq; c++, tn++   )
		       {
			   r1=col[b];
			   r2=col[c];
			   
			   if (r1=='x' || r2=='x')score[A->nseq][a]=F;
			   else if (r1=='-' && r2=='-');
			   else if (r1=='-' || r2=='-');
			   else 
			       {
				   score[A->nseq][a]+= mat[r1-'a'][r2-'a'];
				   
			       }
			   n+=( !is_gap(r1) && !is_gap(r2));
		       }   
 	       score[A->nseq][a]=(((tn!=0)?score[A->nseq][a]/tn:0));
		 
	       }
	       
	   /*initialisation*/

	   transitions=declare_int ( nstate, nstate);
	   v_tab=declare_int ( A->len_aln+2, nstate*nstate);
	   v_tab_p=declare_int ( A->len_aln+2, nstate*nstate);

	   for (a=0; a<nstate;a++)
	     for (b=0; b<nstate;b++)
	       {transitions[a][b]=F;}

	   SPLICE_PENALTY=-1000;

	   transitions[ORF1][ORF2]    =AL;
	   transitions[ORF1][SPL1]    =AL-SPLICE_PENALTY;
	   
	   transitions[ORF2][ORF3]    =AL;
	   transitions[ORF2][SPL1]    =AL-SPLICE_PENALTY;
	   
	   transitions[ORF3][ORF1]    =AL;
	   transitions[ORF3][SPL1]    =AL-SPLICE_PENALTY;
	   
	   transitions[ORF3][ORF1]    =AL;
	   transitions[ORF3][SPL1]    =AL-SPLICE_PENALTY;
	   
	   transitions[ORF3][NC]=AL-100;
	   transitions[NC][ORF1]=AL-100;


	   transitions[SPL1][SPL2]=AL;
	   transitions[SPL2][NC  ]=AL-SPLICE_PENALTY;
	   transitions[NC  ][NC  ]=AL;
	   transitions[NC  ][SPL3]=AL-SPLICE_PENALTY;
	   transitions[SPL3][SPL4]=AL;
	   transitions[SPL4][ORF1]=AL;
	   transitions[SPL4][ORF2]=AL;
	   transitions[SPL4][ORF3]=AL;
	   

	   for ( s=0; s<A->nseq; s++)
	       {
	       for ( p=0; p<=A->len_aln; p++){for (state=0; state< nstate; state++)v_tab_p[p][state]=-1; }
	       for (p=1+2; p<= A->len_aln; p++)
	           {
		    frame1=score[A->nseq][(p-1)];
		    frame2=score[A->nseq][(p-1)-1];
		    frame3=score[A->nseq][(p-1)-2];  
		    best_frame=best_int (3, 1, &a, frame1, frame2, frame3);
		    for (state=0; state< nstate; state++)
		       {
			 r=tolower (A->seq_al[s][p-1]);
			 r=(r=='u')?'t':r;
			 
			 if      (state==ORF1)e=frame1;
			 else if (state==ORF2)e=frame2;
			 else if (state==ORF3)e=frame3;
			 else if (state==SPL1)e=(r=='g')?best_frame:F;
			 else if (state==SPL2)e=(r=='t')?best_frame:F;
			 else if (state==SPL3)e=(r=='a')?best_frame:F;
			 else if (state==SPL4)e=(r=='g')?best_frame:F;
			 else if (state==NC)e=-best_frame;
			 for ( pstate=0; pstate<nstate; pstate++)
		               {
				   v=e+transitions[pstate][state]+v_tab[p-1][pstate];
				   if (pstate==0 ||(v>best_pstate_v) ){best_pstate_v=v;best_pstate_p=pstate;}
			       }
			
			   v_tab[p][state]=best_pstate_v;
			   v_tab_p[p][state]=best_pstate_p;
			   if (state==0 ||best_pstate_v>best_state_v ){best_state_p=state; best_state_v=best_pstate_v;}
		       }
		   }

	       
       
	       for (p=0; p< A->len_aln; p++)T->seq_al[s][p]='.';
	       for (p=A->len_aln; p>0; p--)
	           {
		       if ( best_state_p==0)T->seq_al[s][p-1]=toupper(translate_dna_codon (A->seq_al[s]+(p-1), 'x'));
		       else if ( best_state_p>=SPL1  && best_state_p<=SPL4)T->seq_al[s][p-1]='-';
		       best_state_p=v_tab_p[p][best_state_p];
		   }
	       }
	   
	   

	   vfree (col);
	   return T;
       }

Alignment * mutate_cdna_aln ( Alignment *A)
{
    int a, b, c, n;
    int n1, n2, r1, r2;
    int **pos, ps;
    int neutral_substitution=50;
    int random_substitution=0;
    int random_deletion=0;
    int amino_acid_deletion=0;
    int amino_acid_substitution=0;
    char nuc_list[]="agct";
    char *new_codon;

    neutral_substitution=atoi(get_env_variable ("NEUTRAL_SUBSTITUTION", 1));
    random_substitution =atoi(get_env_variable ("RANDOM_SUBSTITUTION", 1));
    random_deletion     =atoi(get_env_variable ("RANDOM_DELETION", 1));
    amino_acid_deletion =atoi(get_env_variable ("AMINO_ACID_DELETION", 1));
    amino_acid_substitution =atoi(get_env_variable ("AMINO_ACID_SUBSTITUTION", 1));
    
    
    if (A->S)free_sequence ( A->S, (A->S)->nseq);
    A->S=aln2seq(A);

    addrandinit(time (NULL));

    
    pos=aln2pos_simple ( A, A->nseq);
    
    /* 1 Apply neutral substitutions    */
    
    if ( neutral_substitution)
        {
	for (  c=0; c< neutral_substitution; c++)
	    {
	    for (  a=0; a< A->nseq; a++)
                {
		    
		    for ( b=0; b< A->len_aln; b++)
		        {
			
			if (pos[a][b]<=0)continue; 
			ps=MAX(0,pos[a][b]-(pos[a][b]-1)%3-1);


			n1=(A->S)->seq[a][pos[a][b]-1];
			r1=translate_dna_codon ( (A->S)->seq[a]+ps, 'o');
			
			n2=nuc_list[(int)addrand((unsigned long) 4)];
			(A->S)->seq[a][pos[a][b]-1]=n2;
			r2=translate_dna_codon ( (A->S)->seq[a]+ps, 'o');
			
			
			if ( r1==r2 && r1!='o')A->seq_al[a][b]=n2;
			
			else (A->S)->seq[a][pos[a][b]-1]=n1;
			}
		}
	    }
	}

    /* 2 Apply         substitutions    */
     if ( random_substitution)
        {
	for (  a=0; a< A->nseq; a++)
            {
		for ( b=0; b< A->len_aln; b++)
		    {
		    if (pos[a][b]<=0)continue; 
		    if (addrand ((unsigned long) 100)>random_substitution)continue; 
		    
		    n1=nuc_list[(int)addrand((unsigned long)4)];
		    (A->S)->seq[a][pos[a][b]-1]=n1;
		    A->seq_al[a][b]=n1;
		    }
	    }
	}
    
    /* 3 Apply amino acid substitutions */
      if ( amino_acid_substitution)
        {
	for (  a=0; a< A->nseq; a++)
            {
		for ( b=0; b< A->len_aln; b+=3)
		    {
		    if (pos[a][b]<=0)continue; 
		    if (addrand ((unsigned long) 100)>amino_acid_substitution)continue; 
		    ps=MAX(0,pos[a][b]-(pos[a][b]-1)%3-1);
		    
		    r1=translate_dna_codon ( (A->S)->seq[a]+ps, 'o');
		    new_codon=mutate_amino_acid(r1, "clustalw_col");
		    
		    for ( c=ps; c<ps+3; c++)(A->S)->seq[a][c]=new_codon[c-ps];
		    }
		for ( b=0; b< A->len_aln; b++)
		    {
		    if (pos[a][b]<=0)continue; 
		    else A->seq_al[a][b]=(A->S)->seq[a][pos[a][b]-1];
		    }
	    }
	}  
    /* 3 Apply amino acid deletions     */
     if ( amino_acid_deletion)
        {
	for (  a=0; a< A->nseq; a++)
            {
		for ( b=0; b< A->len_aln; b+=3)
		    {
		    if (pos[a][b]<=0)continue; 
		    if (addrand ((unsigned long) 1000)>amino_acid_deletion)continue; 
		    ps=MAX(0,pos[a][b]-(pos[a][b]-1)%3-1);
		    n=addrand ((unsigned long) 4)+1;
		    
		    for ( c=ps; c<ps+(3*n) && c<A->len_aln; c++)(A->S)->seq[a][c]='-';
		    }
		for ( b=0; b< A->len_aln; b++)
		    {
		    if (pos[a][b]<=0)continue; 
		    else A->seq_al[a][b]=(A->S)->seq[a][pos[a][b]-1];
		    }
	    }
	}
    /* 4 Apply amino acid insertions    */

/*FRAMESHIFT MUTATIONS*/
    /* 5 Apply nucleotide deletions*/
     if ( random_deletion)
        {
	for (  a=0; a< A->nseq; a++)
            {
		for ( b=0; b< A->len_aln; b++)
		    {
		    if (pos[a][b]<=0)continue; 
		    if (addrand ((unsigned long) 1000)>random_deletion)continue; 
		    
		    n1='-';
		    (A->S)->seq[a][pos[a][b]-1]=n1;
		    A->seq_al[a][b]=n1;
		    }
	    }
	}
    /* 6 Apply nucleotide deletions*/
    return A;

}    
    
Alignment* clean_est  ( Alignment *A)
        {
	  /*Rules are as follow:
	    Internal Gap > 30% Requences ----> -
	    Best Residue < 50% Residues  ----> 'N'
	  */
	  int a, b,c;
	  int best;
	  int tot;

	  for ( a=0; a< A->len_aln; a++)
	    {
	      
	      for (tot=0, b=0; b<4; b++)tot+=(A->P)->count[b][a];
	      best=best_int (5,1, &c, (A->P)->count[0][a],(A->P)->count[1][a],(A->P)->count[2][a],(A->P)->count[3][a],(A->P)->count[4][a]);
	      
	      if ( tot==0)
		{
		  fprintf ( stderr, "\nWARNING: POSITION WITH NO INFORMATION [clean_est:%s]", PROGRAM);
		  A->seq_al[0][a]='-';
		}
	      else if (((A->P)->count[4][a]*100)/tot >30)A->seq_al[0][a]='-';
	      else if ( (best*100)/tot<50)A->seq_al[0][a]='n';
	      
	    }
	return A;
	}
	   
    

char **make_symbols ( char *name, int *n)
    {
    char **symbol;

    symbol=declare_char ( STRING, STRING);
    
    if ( strcmp (name, "3d_ali")==0)
        {
	sprintf ( symbol[0], "gih");
	sprintf ( symbol[1], "eb");
	sprintf ( symbol[2], "x");
	sprintf ( symbol[3], "#l");
	n[0]=4;
	}
    else if ( strcmp (name, "set1")==0)
        {
	sprintf ( symbol[0], "ilvmfywhktcagH");
	sprintf ( symbol[1], "reqdnsP");
	sprintf ( symbol[2], "--");
	sprintf ( symbol[3], "#l");
	n[0]=4;
	}
    else if ( strcmp (name, "set2")==0)
        {
	n[0]=0;
	sprintf ( symbol[n[0]++], "gsacT");
	sprintf ( symbol[n[0]++], "ndtvpS");
	sprintf ( symbol[n[0]++], "ilkreqL");
	sprintf ( symbol[n[0]++], "--");
	sprintf ( symbol[n[0]++],"#l"); 
	}
    else if ( strcmp ( name, "any")==0)
        {
	sprintf ( symbol[0], "*x");
	n[0]=1;
       	}




    return symbol;
    }

char * translate_dna_seq_on3frame (  char *dna_seq, char stop, char *prot)
       {
	  int a, l;
	  char *buf;

	  l=strlen (dna_seq);
	  if ( prot==NULL)prot=vcalloc ( l+2, sizeof (char));
	   
	   buf=vcalloc (l+4, sizeof (char));
	   sprintf (buf, "%s", dna_seq);
	   lower_string ( buf);
	   for ( a=0; a< l; a++)buf[a]=(buf[a]=='t')?'u':buf[a];
	   
	   for (a=0; a< l; a++)
	       prot[a]=translate_dna_codon (buf+a, stop);
	   free (buf);
	   prot[a]='\0';

	   return prot;
       }
char * translate_dna_seq ( char *dna_seq, int frame, char stop, char *prot)
       {
	   int a, b, l;
	   char *buf;

           l=strlen (dna_seq);
	   if ( prot==NULL)prot=vcalloc ( l/3 +2, sizeof (char));
	   
	   buf=vcalloc (l+4, sizeof (char));
	   sprintf (buf, "%s", dna_seq);
	   lower_string ( buf);
	   for ( a=0; a< l; a++)buf[a]=(buf[a]=='t')?'u':buf[a];
	   
	   for ( b=0,a=0+frame; a< l; a+=3,b++)
	       prot[b]=translate_dna_codon (buf+a, stop);
	   free (buf);
	   prot[b]='\0';

	   return prot;
       }
char * back_translate_dna_codon ( char aa, int deterministic)
        {
	static char *r;
	int choice;
	
	srand(time(NULL));
	if ( r==NULL)r=vcalloc (4, sizeof (char));
	if (!is_gap(aa))aa=tolower(aa);
	  
	if (is_gap(aa))sprintf (r, "---");
	else if ( aa=='a')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if      ( choice==0)sprintf (r, "gca");
	    else if ( choice==1)sprintf (r, "gcg");
	    else if ( choice==2)sprintf (r, "gcc");
	    else if ( choice==3)sprintf (r, "gct");
	  }
	else if ( aa=='c')
	  {
	   choice=(deterministic)?0:rand()%2;
	    if      ( choice==0)sprintf (r, "tgc");
	    else if ( choice==1)sprintf (r, "tgt");
	  } 
	else if ( aa=='d')
	  {
	  choice=(deterministic)?0:rand()%2;
	  if ( choice==0)sprintf (r, "gac");
	  else if ( choice==1)sprintf (r, "gat");
	  }
	
	else if ( aa=='e')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if ( choice==0)sprintf (r, "gaa");
	    else sprintf (r, "gag");
	  }
	else if ( aa=='f')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if ( choice==0)sprintf (r, "ttc");
	    else sprintf (r, "ttt");
	  }
	else if ( aa=='g')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if  ( choice==0)     sprintf (r, "gga");
	    else if ( choice==1) sprintf (r, "ggg");
	    else if ( choice==2) sprintf (r, "ggc");
	    else if ( choice==3) sprintf (r, "ggt");
	  }	
	else if ( aa=='h')
	  {
	    choice =rand()%2;
	    if ( choice==0)sprintf (r, "cac");
	    else sprintf (r, "cat");
	  }
	else if ( aa=='i')
	  {
	    choice=(deterministic)?0:rand()%3;
	    if  ( choice==0)     sprintf (r, "ata");
	    else if ( choice==1) sprintf (r, "atc");
	    else if ( choice==2) sprintf (r, "att");
	  }	
	else if ( aa=='k')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if  ( choice==0)     sprintf (r, "aaa");
	    else if ( choice==1) sprintf (r, "aag");
	    
	  }
	else if ( aa=='l')
	  {
	    choice=(deterministic)?0:rand()%6;
	    if  ( choice==0)     sprintf (r, "cta");
	    else if ( choice==1) sprintf (r, "ctg");
	    else if ( choice==2) sprintf (r, "ctc");
	    else if ( choice==3) sprintf (r, "ctt");
	    else if ( choice==4) sprintf (r, "tta");
	    else if ( choice==5) sprintf (r, "ttg");	    
	  }	
	else if ( aa=='m')sprintf ( r, "atg");
	else if ( aa=='n')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if  ( choice==0)     sprintf (r, "aac");
	    else if ( choice==1) sprintf (r, "aat");
	  }	
	else if ( aa=='p')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if  ( choice==0)     sprintf (r, "cca");
	    else if ( choice==1) sprintf (r, "ccg");
	    else if ( choice==2) sprintf (r, "ccc");
	    else if ( choice==3) sprintf (r, "cct");
	  }	
	else if ( aa=='q')
	  {
	    choice=(deterministic)?0:rand()%2;
	    if  ( choice==0)     sprintf (r, "caa");
	    else if ( choice==1) sprintf (r, "cag");
	  }
        else if ( aa=='r')
	  {
	    choice=(deterministic)?0:rand()%6;
	    if  ( choice==0)     sprintf (r, "cga");
	    else if ( choice==1) sprintf (r, "cgg");
	    else if ( choice==2) sprintf (r, "cgc");
	    else if ( choice==3) sprintf (r, "cgt");
	    else if ( choice==4) sprintf (r, "aga");
	    else if ( choice==5) sprintf (r, "agg");
	    
	  }
	else if ( aa=='s')
	  {
	    choice=(deterministic)?0:rand()%6;
	    if  ( choice==0)     sprintf (r, "tca");
	    else if ( choice==1) sprintf (r, "tcg");
	    else if ( choice==2) sprintf (r, "tcc");
	    else if ( choice==3) sprintf (r, "tct");
	    else if ( choice==4) sprintf (r, "agt");
	    else if ( choice==5) sprintf (r, "agc");
	    
	  }
	else if ( aa=='t')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if  ( choice==0)     sprintf (r, "aca");
	    else if ( choice==1) sprintf (r, "acg");
	    else if ( choice==2) sprintf (r, "acc");
	    else if ( choice==3) sprintf (r, "act");
	  }
	else if ( aa=='v')
	  {
	    choice=(deterministic)?0:rand()%4;
	    if  ( choice==0)     sprintf (r, "gta");
	    else if ( choice==1) sprintf (r, "gtg");
	    else if ( choice==2) sprintf (r, "gtc");
	    else if ( choice==3) sprintf (r, "gtt");
	  }
	else if ( aa=='w')
	  {
	    sprintf (r, "tgg");
	  }
	else if ( aa=='y')
	  {
	     choice=(deterministic)?0:rand()%2;
	    if  ( choice==0)     sprintf (r, "tac");
	    else if ( choice==1) sprintf (r, "tat");
	  }
	else
	  {
	    sprintf (r, "nnn");
	  }
	return r;
		
	}
int translate_dna_codon ( char *sequence, char stop)
        {
	char seq[4];
	int a,b;


	if ( (b=strlen (sequence))<3)
	  {
	    for ( a=0; a<b; a++)
	      if ( !is_gap(sequence[a]))return 'x';
	  return '-';
	  }
	else 
	  {
	    seq[0]=tolower(sequence[0]);
	    seq[1]=tolower(sequence[1]); 
	    seq[2]=tolower(sequence[2]);
	    seq[3]='\0';
	    
	    seq[0]=(seq[0]=='u')?'t':seq[0];
	    seq[1]=(seq[1]=='u')?'t':seq[1];
	    seq[2]=(seq[2]=='u')?'t':seq[2];
	   
	}
	

	
	if ( is_gap(seq[0])||is_gap(seq[1]) || is_gap(seq[2]))return '-';
	else if ( strm5(seq, "gca", "gcg", "gcc", "gct","gcn"))return 'a';
	else if ( strm2(seq, "tgc","tgt"))return 'c';
	else if ( strm2(seq, "gac","gat"))return 'd';
	else if ( strm2(seq, "gaa","gag"))return 'e';
        else if ( strm2(seq, "ttc","ttt"))return 'f';
	else if ( strm5(seq, "gga","ggg","ggc", "ggt", "ggn"))return 'g';
	else if ( strm2(seq, "cac","cat"))return 'h';
	else if ( strm3(seq, "ata","atc","att"))return 'i';
	else if ( strm2(seq, "aaa","aag"))return 'k';
        else if ( strm6(seq, "cta","ctg","ctc", "ctt", "tta", "ttg"))return 'l';
	else if ( strm (seq, "ctn"))return 'l';
	else if ( strm (seq, "atg"))return 'm';
	else if ( strm2(seq, "aac","aat"))return 'n';
	else if ( strm5(seq, "cca","ccg","ccc", "cct","ccn"))return 'p';
	else if ( strm2(seq, "cag","caa"))return 'q';
	else if ( strm6(seq, "cga","cgg","cgc", "cgt","aga","agg"))return 'r';
	else if ( strm (seq, "cgn"))return 'r';
	else if ( strm6(seq, "tca","tcg","tcc", "tct","agc","agt"))return 's';
	else if ( strm (seq, "ccn"))return 's';
        else if ( strm5(seq, "aca","acg","acc", "act", "acn"))return 't';
	else if ( strm5(seq, "gta","gtg","gtc", "gtt", "gtn"))return 'v';
	else if ( strm (seq, "tgg"))return 'w';
	else if ( strm2(seq, "tac","tat"))return 'y';
	else if ( strm3(seq, "tag","taa","tga"))return stop;
	else if ( seq[0]=='n' || seq[1]=='n' || seq[2]=='n') return stop;
	else
	  {
	    fprintf ( stderr, "\n%s is an unknown codon [FATAL]",seq);
	    exit (1);
	    return 1;
	  }
	}
	   
Alignment * mutate_aln ( Alignment *A, char *r)
{
  int a, b, c, mut,type, ratio;
  char alp[30];
  int alp_size;
  Sequence *S;
  Alignment*B;
  int n_mut, tot;

  srand(time(NULL));
  if ( r[0]=='\0')ratio=0.01*RAND_MAX;
  else ratio=atof(r)*RAND_MAX;

  S=aln2seq(A);
  S=get_sequence_type(S);
  


  if ( strm(S->type, "DNA"))sprintf (alp, "AGCT");
  else if (  strm(S->type, "PROTEIN"))sprintf (alp, "ACDEFGHIKLMNPQRSTVWY");

  alp_size=strlen(alp);

  B=copy_aln (A,NULL);
  B=realloc_aln(B, B->len_aln*2+1);

  for ( a=0, b=0; a< A->len_aln; a++, b+=2)
    {
      for ( c=0; c< A->nseq; c++)
	{
	  B->seq_al[c][b]=tolower(A->seq_al[c][a]);
	  B->seq_al[c][b+1]='~';
	}      
    }

  for ( c=0; c< A->nseq; c++)B->seq_al[c][b]='\0';
  B->len_aln=A->len_aln*2;
  

 
  tot=n_mut=0;
  for (a=0; a< B->len_aln; a+=2)
    for ( b=0; b<B->nseq; b++)
      {
	if ( is_gap(B->seq_al[b][a]))continue;
	mut=((rand()%RAND_MAX)>ratio)?0:1;
	tot++;
	n_mut+=mut;

	if (mut)
	  {
	    type=rand()%2;
	    if (type==0)/*deletion*/
	      {
		B->seq_al[b][a]='.';
	      }
	    else if ( type==1)
	      {
		B->seq_al[b][a+1]=alp[rand()%alp_size];
	      }
	    else if (type==2)
	      {
		B->seq_al[b][a]=alp[rand()%alp_size];
	      }
	    
	  }
      }
  ungap_aln (B);
  
  
  free_sequence (S, S->nseq);
  free_aln (A);
  return B;
  
}

char* mutate_amino_acid ( char aa, char *mode)

     {
	 int a, b, c, d;
	 char nucleotide[]="agct";
	 char amino_acid[]="acdefghiklmnpqrstvwy";
	 static char **triplet;
	 static char **cw_col;
	 int ng_cw_col;
	 static int **amino_acid_list;
	 static int *lu;
	 char a1, a2;
	 char *mat;
	 
	 aa=tolower(aa);
	 declare_name(mat);
	 if ( !mode)sprintf (mat, "clustalw_col");
	 else sprintf (mat, "%s", mode);
	 if (!triplet)
	    {
		triplet=declare_char ( 64, 4);
		for (d=0, a=0; a< 4;a++)
		    for ( b=0; b< 4; b++)
			for ( c=0; c< 4; c++, d++)
			    {
				triplet[d][0]=nucleotide[a];
				triplet[d][1]=nucleotide[b];
				triplet[d][2]=nucleotide[c];
			    }
	    }
	 if ( !cw_col)cw_col=make_group_aa ( &ng_cw_col,mat);
	 if ( !amino_acid_list)
	    {
		amino_acid_list=declare_int ( 20, 65);
		for ( a=0; a< 20; a++)
		    for ( b=0; b< 64; b++)
		        {
			    a1=translate_dna_codon ( triplet[b], 'x');
			    a2=amino_acid[a];
			    for ( d=0; d< ng_cw_col; d++)
				if ( is_in_set ( a1, cw_col[d]) && is_in_set ( a2, cw_col[d]))
				   {
				       amino_acid_list[a][++amino_acid_list[a][0]]=b;
				   }
			}
		lu=vcalloc ( 26, sizeof (int));
		for ( a=0; a<20; a++)
		    {
			lu[amino_acid[a]-'a']=a;
		    }
		/*
		for ( a=0; a< 20; a++)
		    {
			fprintf ( stderr, "\n%c", amino_acid[a]);
			for ( b=1; b<=amino_acid_list[a][0]; b++)
			    fprintf ( stderr, "\n\t%s %c", triplet[amino_acid_list[a][b]], translate_dna_codon (triplet[amino_acid_list[a][b]], 'x'));
		    }
		*/		    
	    }
	
	 return triplet [addrand((unsigned long)amino_acid_list[lu[aa-'a']][0])+1];
     }			
				 
/**************************************************************************************************/
/********************************                      ********************************************/
/********************************    PROCESSING        ********************************************/
/********************************                      ********************************************/


	 
void modify_data  (Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char *action, Action_data_struc *RAD)
     {
       Sequence  *COOR=NULL, *NS=NULL, *OUT_S=NULL;
       char *s;
       int value, start, end, a;
       static int keep_name;
       char separator[2];
       char sep1, sep2;
       
       if (  strm(action, "seqnos"))
	 {
	  (D1->A)->output_res_num=1;
	 }  
       else if (strm(action, "keep_name"))
	 {
	   keep_name=1-keep_name;
	 }
       else if ( strncmp (action, "rm_gap",6)==0)
	 {

	   ungap_aln_n (D1->A, (action[6]=='\0')?100:atoi(action+6));
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       else if ( strcmp (action, "clean_maln")==0)
	  {
	    if ( !DST) 
		   {
		   fprintf ( stderr,"\n[You Need an evaluation File: Change the output format][FATAL]\n");	
		   exit(1);
		   }
	    (DST->A)=aln2number (DST->A);
	    D1->A=clean_maln(D1->A, DST->A, 1, 1);   
	  }
       else if ( strm (action, "extract"))
	 {
	   
	   COOR=get_pir_sequence  (RAD->coor_file, NULL);
	   D1->S=extract_sub_seq ( COOR, D1->S);
	   free_aln (D1->A);
	   D1->A=declare_Alignment(D1->S);
	   seq2aln (D1->S, D1->A, RAD->rm_gap);
	   free_sequence (COOR, COOR->nseq);
	 }
       else if ( strncmp (action, "extract_seq",11)==0)	 
	 {
	  
	   if ( action[11]==':')
	     {
	       sep1=':';
	       sep2='|';
	       action+=12;
	     }
	   else if ( action[11]=='_')
	     {
	       sep1=action[12];
	       sep2=action[13];
	       action+=14;
	     }
	   
	   separator[0]=sep1;
	   separator[1]='\0';
	   
	   
	   
	   s =strtok (action,separator);	   	   
	   s=seq_name2coor (s, &start, &end, sep2);
	   if ( strm (s, "*"))
	     {
	       OUT_S=extract_one_seq((D1->A)->name[0],start, end, D1->A, keep_name);
	       for ( a=1; a< (D1->A)->nseq; a++)
		 {
		   NS=extract_one_seq((D1->A)->name[a],start, end, D1->A, keep_name);
		   if (count_n_res_in_array(NS->seq[0], -1))
		       OUT_S=add_sequence ( NS,OUT_S, 0);
		 }
	     }
	   else	     
	     {
	       OUT_S=extract_one_seq(s,start, end, D1->A, keep_name);
	   
	       while ((s=strtok (NULL,separator))!=NULL)
		 {
		   s=seq_name2coor (s, &start, &end, sep2);
		   NS=extract_one_seq(s,start, end, D1->A, keep_name);
		   OUT_S=add_sequence ( NS,OUT_S, 0);
		 }
	     }
	   
	   D1->S=OUT_S;
	   free_aln (D1->A);
	   D1->A=declare_Alignment(D1->S);
	   seq2aln (D1->S, D1->A, RAD->rm_gap);
	 }
       else if ( strm (action, "clean_cdna"))
	 {
	   D1->A=clean_cdna_aln ( D1->A);
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       
       else if (strm ( action, "translate"))
	 {
	  D1->A=translate_dna_aln( D1->A, 0);
	  free_sequence ( D1->S, (D1->S)->nseq);
	  D1->S=aln2seq ( D1->A);
	 }
       else if ( strncmp ( action, "translate", 9)==0)
	 {
	   D1->A=translate_dna_aln( D1->A, atoi(action+9));
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq ( D1->A);
	 }
       else if (strm2 ( action, "back_translate","backtranslate"))
	 {
	  D1->A=back_translate_dna_aln( D1->A);
	  free_sequence ( D1->S, (D1->S)->nseq);
	  D1->S=aln2seq ( D1->A);
	 }
       else if (strm ( action, "code_dna_aln"))
	 {
	  D1->A=code_dna_aln( D1->A);
	  free_sequence ( D1->S, (D1->S)->nseq);
	  D1->S=aln2seq ( D1->A);
	 }
       else if ( strncmp ( action, "mutate", 6)==0)
	 {
	   D1->A=mutate_aln( D1->A,action+6);
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
       	else if ( strm ( action, "thread_dna_on_prot_aln"))
	  {
	    D1->A=thread_dnaseq_on_prot_aln (D1->S, D2->A);
	    free_sequence (D1->S,(D1->S)->nseq);
	    D1->S=aln2seq (D1->A); 
	  }
       else if ( strm ( action, "thread_struc_on_aln"))
	 {
	   thread_seq_struc2aln ( D2->A, D1->S);
	   D1->A=copy_aln(D2->A, D1->A);
	   free_sequence ( D1->S, (D1->S)->nseq);
	   D1->S=aln2seq (D1->A);
	 }
       else if (strncmp ( action, "lower",5)==0)
	 {
	   if (action[5]=='\0')value=10;
	   else value=atoi(action+5);
	   
	   D1->A=filter_aln_upper_lower (D1->A, (DST)?DST->A:NULL,value);
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A); 
	 } 
       else if (strncmp ( action, "trim",4)==0)
	 {
	   value=(action[5])? atoi(action+4):10;
	   D1->A=trimseq(D1->A,D1->S,action+4);
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A); 
	 }
       else if (strncmp ( action, "upper",5)==0)
	 {
	   if (action[5]=='\0')value=10;
	   else value=atoi(action+5);
	   
	   D1->A=filter_aln_lower_upper (D1->A, DST?DST->A:NULL,value);
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A); 
	 }
       else if ( strncmp (action, "convert",7)==0)
	 {
	   if (action[7]=='\0')value=10;
	   else value=atoi(action+7);
	   
	   
	   D1->A=filter_aln_convert (D1->A, DST?DST->A:NULL,value,RAD->n_symbol, RAD->symbol_list);
	   free_sequence (D1->S,(D1->S)->nseq);
	   D1->S=aln2seq (D1->A); 
	   /*
	   string_array_convert ( (D1->S)->seq,   (D1->S)->nseq, RAD->n_symbol, RAD->symbol_list);
	   string_array_convert ( (D1->A)->seq_al,(D1->A)->nseq, RAD->n_symbol, RAD->symbol_list);   
	   */
	 }
       else
	 {
	   fprintf ( stderr, "\nWARNING: ACTION %s UNKNOWN and IGNORED\n", action);
	 }
     }

