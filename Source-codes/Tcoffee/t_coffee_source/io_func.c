#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "matrices.h"

#define DEFAULT_COLOR -1
#define GAP_COLOR     -2
#define INK_COLOR     -3
	
Sequence * read_sequences ( char *seq_name)
	{
	Sequence *S;


	char **seq=NULL;
	char **name=NULL;
	int *len=NULL;
	int nseq=0;
	int min_len=0;
	int max_len=0;
	int a;


	get_sequence ( seq_name, &nseq, &seq, &name, &len, &min_len, &max_len); 
	
	S=declare_sequence ( min_len, max_len, nseq);
	for ( a=0; a< nseq; a++)sprintf ( S->file[a], "%s", seq_name);
	
	for ( a=0; a<S->nseq; a++)
	    {
	    S->len[a]=len[a];
	    sprintf ( S->name[a],"%s", name[a]);
	    free ( name[a]);
	    sprintf ( S->seq[a], "%s", seq[a]);
	    free ( seq[a]);
	    }
	free (seq);
	free (name);
	free (len);
	S=get_sequence_type ( S);
	return S;
	}
char *    get_string_type   (char *S)
        {
	int a, l;
	int protein=0,  dna=0, tot=0;
	char *dna_type;
	char *protein_type;
	

	
	declare_name(dna_type);sprintf (dna_type, "DNA");
	declare_name(protein_type);sprintf (protein_type, "PROTEIN");
	
	l=strlen (S);
	
	
	
	for ( a=0; a<l; a++)     
	        {
		    if ( !is_gap(S[a]))
			{
			protein+=is_aa(S[a]);
			dna+=is_dna(S[a]);
			tot++;
			}
		}

	if ( protein==0 && dna==0);
        else if ( protein==0 ){vfree(protein_type);return dna_type;}
        else if ( dna==0){vfree(dna_type);return protein_type;}
        else if ( ((dna*100)/protein)>70){vfree(protein_type);return dna_type;}
        else 
	    {
	      vfree(dna_type);
	      return protein_type;
	    }
	vfree (dna_type);
	return protein_type;
	}

Alignment* get_aln_type (Alignment *A)
        {
	  if ( A->S!=NULL && (A->S)->type[0]!='\0')
	    {
	    ;
	    }
	  else if (A->S!=NULL && (A->S)->type[0]=='\0') 
	    {
	      A->S=get_sequence_type (A->S);
	    }
	  else if (A->S==NULL) 
	      {
		A->S=aln2seq (A);
		A->S=get_sequence_type(A->S);
	      }
	  return A;
	}
  
Sequence *get_sequence_type (Sequence *S)
        {
	int a;
	int protein=0,  dna=0;
	char *type;
	
	for ( a=0; a< S->nseq; a++)
	    {
	     
	    type=get_string_type (S->seq[a]);
	    
	    if ( strm (type, "DNA"))dna++;
	    else if strm ( type, "PROTEIN")protein++;
	    vfree(type);
	    }

	if ( protein==0 && dna==0);
	else if ( protein==0 )sprintf ( S->type, "DNA");
	else if ( dna==0)sprintf ( S->type, "PROTEIN");
	else sprintf ( S->type, "PROTEIN_DNA");
	return S;
	}

void get_sequence (char *seq_file,int *NSEQ, char ***SEQ, char ***SN, int **sl, int *min, int *max)
	{
	int a,b;
	int min_len;
	int max_len;
	int nseq;

	int **SL;
	
	

	
	
	
	nseq=NSEQ[0]= readseqs ( seq_file,  SEQ, SN, &SL);
	sl[0]=calloc ( nseq, sizeof (int));
	
	 
	min_len= max_len= (SL)[0][0];
	for ( a=0; a<NSEQ[0]; a++)
		{
		sl[0][a]=SL[a][0];
		for ( b=0; b<(SL)[a][0]; b++) 
		 	(SEQ[0])[a][b]=tolower ((SEQ[0])[a][b]);
		 } 
	for ( a=1; a<NSEQ[0]; a++)
		{
		min_len= ( min_len > (SL)[a][0])?(SL)[a][0]:min_len;
		max_len= ( max_len < (SL)[a][0])?(SL)[a][0]:max_len;
		}
	min[0]=min_len;
	max[0]=max_len;
	}
int ** get_matrix   ( char *name, char *format)
       {
       if ( strm ( "blast", format))return read_blast_matrix ( name);
       else if ( strm ( "clustalw", format))return read_matrice(name);
       else
           {
	   fprintf ( stderr, "\nError:\nUnknowm Format %s for Matrix %s[FATAL]", format, name);
	   exit (1);
	   }
       return NULL;
       }
       
int ** read_blast_matrix ( char *mat_name)
        {
	FILE *fp;
	int n_aa,aa1, aa2;
	int a, b;
	int **matrix;
	int value;
	char buf[2];

	matrix=declare_int ( SIZEOF_AA_MAT, SIZEOF_AA_MAT);
	matrix[30]=vcalloc(10000, sizeof (int));
	fp=vfopen ( mat_name, "r");
	while ( (fgetc(fp))=='#')while ( (fgetc(fp))!='\n');
	while ( (fgetc(fp))!='\n');
	n_aa=strlen ( BLAST_AA_ALPHABET);
	for ( a=0; a< n_aa; a++)
	    {
	    fscanf ( fp, "%s ", buf);
	    aa1=buf[0];
	    
	    aa1=tolower ((char)aa1);	    
	    if ( aa1!=BLAST_AA_ALPHABET[a])
		{
		fprintf ( stderr, "\nParsing_error when reading blast_matrix %s:\n%c %c",mat_name, aa1,BLAST_AA_ALPHABET[a] );
		fprintf ( stderr, "\n%c ", fgetc(fp));
		exit (1);
		}
	    for ( b=0; b<n_aa; b++)
	        {
		aa2=tolower ((char) BLAST_AA_ALPHABET[b]);
		fscanf ( fp, "%d ", &value);
		if ( aa1!='*' && aa2!='*')
		    matrix[aa1-'a'][aa2-'a']=value;
		}
	    fscanf(fp, "\n");
	    }
	fclose (fp);
	return matrix;
	}


int ** read_matrice (char *mat_name_in)
	{
	int a,b,c;
	int max, min;
	char AA[]="abcdefghiklmnpqrstvwxyz";
        FILE *fp;
	int **matrice;
	int **matrix2;
	char mat_name[200];
	int *vector=NULL;

	matrice=declare_int ( SIZEOF_AA_MAT,SIZEOF_AA_MAT);
	matrix2=declare_int ( SIZEOF_AA_MAT,SIZEOF_AA_MAT);
	



	if ( strm2 (mat_name_in, "pam", "PAM"))sprintf ( mat_name, "pam250mt");
	else if (strm2 (mat_name_in, "blosum", "BLOSUM"))sprintf ( mat_name, "blosum62mt");
	else if (strm3 (mat_name_in, "id", "ID", "idmat"))sprintf ( mat_name, "idmat");
	else sprintf ( mat_name, "%s", mat_name_in);
	

	if (strm(mat_name, "pam250mt"))vector=pam250mt;
	else if (strm(mat_name, "idmat"))vector=idmat;
	else if (strm(mat_name, "est_idmat"))vector=est_idmat;
	else if (strm(mat_name, "md_350mt"))vector=md_350mt;
	else if (strm(mat_name, "md_250mt"))vector=md_250mt;
	else if (strm(mat_name, "md_120mt"))vector=md_120mt;
	else if (strm(mat_name, "md_40mt" ))vector= md_40mt;
	else if (strm(mat_name, "pam350mt" ))vector=pam350mt;
	else if (strm(mat_name, "pam160mt" ))vector=pam160mt;
	else if (strm(mat_name, "pam120mt" ))vector=pam120mt;
	
	else if (strm(mat_name, "blosum80mt" ))vector=blosum80mt;
	else if (strm(mat_name, "blosum62mt" ))vector=blosum62mt;
	else if (strm(mat_name, "blosum62mt2" ))vector=blosum62mt2;
	else if (strm(mat_name, "blosum55mt" ))vector=blosum55mt;
	else if (strm(mat_name, "blosum45mt" ))vector=blosum45mt;
	else if (strm(mat_name, "blosum40mt" ))vector=blosum40mt;
	else if (strm(mat_name, "blosum30mt" ))vector=blosum30mt;
	else if (strm(mat_name, "beta_mat" ))vector=beta_mat;
	else if (strm(mat_name, "alpha_mat" ))vector=alpha_mat;
	else if (strm(mat_name, "coil_mat" ))vector=coil_mat;
	
	else if (strm(mat_name, "rblosum80mt" ))vector=rblosum80mt;
	else if (strm(mat_name, "rblosum62mt" ))vector=rblosum62mt;
	else if (strm(mat_name, "rblosum30mt" ))vector=rblosum30mt;
	
	else if (strm(mat_name, "rpam250mt" ))vector=rpam250mt;
	else if (strm(mat_name, "rpam350mt" ))vector=rpam350mt;
	else if (strm(mat_name, "rpam160mt" ))vector=rpam160mt;
	else if (strm(mat_name, "rpam120mt" ))vector=rpam120mt;

	else if (strm(mat_name, "tmpam250mt" ))vector=tmpam250mt;
	else if (strm(mat_name, "rtmpam250mt" ))vector=rtmpam250mt;

	else if (strm(mat_name, "rbeta_mat" ))vector=rbeta_mat;
	else if (strm(mat_name, "ralpha_mat" ))vector=ralpha_mat;
	else if (strm(mat_name, "rcoil_mat" ))vector=rcoil_mat;
	if(vector)
	  {
	  for (a=0; a<23; ++a)
		{
		 for (b=0;b<=a;++b)
				{matrice[a][b]=matrice[b][a]=(vector[(a*a+a)/2+b]);
				}
		}
	
	  }
	else if ( check_file_exists(mat_name))
		{
		  fp=vfopen ( mat_name, "r");
		  while ( (c=fgetc (fp))!='$' && c!=EOF);
		  if ( c==EOF){vfclose (fp);free_int (matrice, -1);free_int (matrix2, -1);return NULL;};
	
		  fgetc(fp);
		  for ( a=0; a< 23; a++)
		    {
		      for ( b=0; b<=a; b++)
			if (fscanf ( fp, "%d,", &matrice[a][b])==0){vfclose (fp);free_int (matrice, -1);free_int (matrix2, -1);return NULL;};
		      fscanf ( fp, "\n");
		    }		
		  fclose ( fp);
		}
	else {free_int (matrice, -1);free_int (matrix2, -1);return NULL;}
	
	a=strlen (AA);
	min=max=matrice[a][b];
	for ( b=0; b<a; b++)
		for ( c=0; c<a; c++)
			{
			min=(matrice[b][c]<min)?matrice[b][c]:min;
			max=(matrice[b][c]<max)?matrice[b][c]:max;
			}
	
			
	/*
	for ( b=0; b<a; b++)
		for ( c=0; c<a; c++)
			{
			matrice[b][c]=matrice[b][c]-min;
			}
	*/ 
		
	for ( b=0; b<a; b++)
		for ( c=0; c<a; c++)
			{
			matrix2[AA[b]-'a'][AA[c]-'a']=matrice[b][c];
			}
	free_int (matrice, -1);
	return matrix2;	
	
	}
int **neg_matrix2pos_matrix ( int **matrix)
    {
      int b,c,l, min, max;
      char AA[]="abcdefghiklmnpqrstvwxyz";
      l=strlen(AA);
      min=max=matrix[AA[0]-'a'][AA[0]-'a'];
      for ( b=0; b<l; b++)
	for ( c=0; c<l; c++)
	  {
	    min=(matrix[AA[b]-'a'][AA[c]-'a']<min)?matrix[AA[b]-'a'][AA[c]-'a']:min;
	    max=(matrix[AA[b]-'a'][AA[c]-'a']<max)?matrix[AA[b]-'a'][AA[c]-'a']:max;	    	   
	  }
      if (min>0)return matrix;
      else
	 {
	   for ( b=0; b<l; b++)
	     for ( c=0; c<l; c++)
	       {
		 matrix[b][c]=matrix[b][c]-min;
	       }
	 }
      return matrix;
    }

      

/*****************************************************************/
void get_rgb_values ( int val, Color *C)
     {

       /*Colors can be modified, see the definition of 
	 COLOR_FILE in prgogrammes_define.h
       */
       
     static char **html_code;
     static float **ps_code;
     int class, n, c;
     float *r, *g, *b;
     FILE *fp;



     if ( !html_code)
       {
	 html_code=declare_char(10, 10);
	 ps_code=declare_float  (10, 3);
	 
	 n=0;
	 /*0*/
	 sprintf (html_code[n], "#6666FF");
	 ps_code[n][0]=0.4;
	 ps_code[n][1]=0.4;
	 ps_code[n][2]=1;
	 n++;
	 /*1*/
	 sprintf (html_code[n], "#00FF00");
	 ps_code[n][0]=0.6;
	 ps_code[n][1]=1;
	 ps_code[n][2]=0;
	 n++;

	 /*2*/
	 sprintf (html_code[n], "#66FF00");
	 ps_code[n][0]=0.8;
	 ps_code[n][1]=1;
	 ps_code[n][2]=0;
	 n++;
	   
	 /*3*/
	 sprintf (html_code[n], "#CCFF00");
	 ps_code[n][0]=1.0;
	 ps_code[n][1]=1.0;
	 ps_code[n][2]=0;
	 n++;
	   
	 /*4*/
	 sprintf (html_code[n], "#FFFF00");
	 ps_code[n][0]=1;
	 ps_code[n][1]=0.85;
	 ps_code[n][2]=0;
	 n++;
	   
	 /*5*/
	 sprintf (html_code[n], "#FFCC00");
	 ps_code[n][0]=1;
	 ps_code[n][1]=0.7;
	 ps_code[n][2]=0;
	 n++;
	   
	 /*6*/
	 sprintf (html_code[n], "#FF9900");
	 ps_code[n][0]=1;
	 ps_code[n][1]=0.6;
	 ps_code[n][2]=0;
	 n++;
	   
	 /*7*/
	 sprintf (html_code[n], "#FF6600");
	 ps_code[n][0]=1;
	 ps_code[n][1]=0.4;
	 ps_code[n][2]=0;
	 n++;
         
	 /*8*/
	 sprintf (html_code[n], "#FF3300");
	 ps_code[n][0]=1;
	 ps_code[n][1]=0.2;
	 ps_code[n][2]=0;
	 n++;
	 
	   
	 /*9*/
	 sprintf (html_code[n], "#FF2000");
	 ps_code[n][0]=1;
	 ps_code[n][1]=0;
	 ps_code[n][2]=0;
	 n++;
	 
         
	 
	 

	 if ( check_file_exists(COLOR_FILE))
	   {
	     fp=vfopen( COLOR_FILE, "r");
	     
	     while ((c=fgetc(fp))!='*');
	     while ((c=fgetc(fp))!='\n');
	     
	     c=0;
	     while ((c=fgetc(fp))!=EOF)
	       {
		 ungetc(c, fp);
		 if ( fscanf (fp, "%d", &class)==0)break;
		 fscanf (fp, "%s %f %f %f", html_code[class], &ps_code[class][0], &ps_code[class][1],&ps_code[class][2]);
		 while ((c=fgetc(fp))!='\n' && c!=EOF);
		 if ( c==EOF)ungetc(c, fp);
	       }
	     vfclose(fp);
	   }
       }
     


     /*Conversions*/
       if ( val==10)val--;
       
     r=&C->r;
     g=&C->g;
     b=&C->b;

     if ( val==10)val--;
     sprintf ( C->html_color_class, "value%d",val); 
     
     
     if (val<=9 && val>=0)
       {

	 sprintf ( C->html_color, "%s", html_code[val]);
	 r[0]=ps_code[val][0];
	 g[0]=ps_code[val][1];
	 b[0]=ps_code[val][2];
       }

     else if (val==DEFAULT_COLOR || val==NO_COLOR_RESIDUE || val==NO_COLOR_GAP || (val>'A' && val<'z'))
       {
	C->html_color[0]='\0';
	sprintf ( C->html_color_class, "valuedefault");
	r[0]=1.;
	g[0]=1;
	b[0]=1;
	
	} 
     else if (val==GAP_COLOR)
       {
	C->html_color[0]='\0';
	sprintf ( C->html_color_class, "valuegap");
	r[0]=1.;
	g[0]=1;
	b[0]=1; 
	} 
     else if (val==INK_COLOR )
        {
	sprintf ( C->html_color, "000000");
	sprintf ( C->html_color_class, "valueink");
	r[0]=0.;
	g[0]=0;
	b[0]=0;
	}
     return;

     if ( val==0)
        {
	sprintf ( C->html_color, "%s", html_code[val]);/*0.4 1.0 0.0*/
	r[0]=0.4;
	g[0]=0.4;
	b[0]=1;
	}
     else if ( val==1)
        {
	sprintf ( C->html_color, "#00FF00");/*0.4 1.0 0.0*/
	r[0]=0.6;
	g[0]=1;
	b[0]=0;
	}
     else if (val==2)
        {
	sprintf ( C->html_color, "#66FF00");/*0.6 1.0 0.0*/
	r[0]=0.8;
	g[0]=1;
	b[0]=0;
	}
     else if ( val==3)
        {
	sprintf ( C->html_color, "#CCFF00");/*0.8 1.0 0.0*/
	r[0]=1.0;
	g[0]=1.0;
	b[0]=0;
	}
     else if (val==4)
        {
	sprintf ( C->html_color, "#FFFF00");/*1 1.0 0.0*/    
	r[0]=1;
	g[0]=0.85;
	b[0]=0;
	}
     else if ( val==5)
        {
	sprintf ( C->html_color, "#FFCC00");/*1.0 0.85 0.0*/  	
	r[0]=1;
	g[0]=0.7;
	b[0]=0;
	}
     else if (val==6)
        {
	sprintf ( C->html_color, "#FF9900");/*1.0 0.7 0.0*/  
	r[0]=1;
	g[0]=0.6;
	b[0]=0;
	}
     else if ( val==7)
        {
	sprintf ( C->html_color, "#FF6600");/*1.0 0.6 0.0*/  
	r[0]=1;
	g[0]=0.4;
	b[0]=0;
	}
     else if (val==8)
        {
	sprintf ( C->html_color, "#FF3300 ");/*1.0 0.4 0.0*/  	
	r[0]=1;
	g[0]=0.2;
	b[0]=0;
	}   
     else if (val==9 && val<'a')
        {
	sprintf ( C->html_color, "#FF2000");/*1.0 0.2 0.0*/  	
	r[0]=1.;
	g[0]=0;
	b[0]=0;
	}
     else if (val==DEFAULT_COLOR || val==NO_COLOR_RESIDUE || val==NO_COLOR_GAP || (val>'A' && val<'z'))
        {
	C->html_color[0]='\0';
	sprintf ( C->html_color_class, "valuedefault");
	r[0]=1.;
	g[0]=1;
	b[0]=1;
	
	}
      else if (val==GAP_COLOR)
        {
	C->html_color[0]='\0';
	sprintf ( C->html_color_class, "valuegap");
	r[0]=1.;
	g[0]=1;
	b[0]=1; 
	}
      else if (val==INK_COLOR )
        {
	sprintf ( C->html_color, "000000");
	sprintf ( C->html_color_class, "valueink");
	r[0]=0.;
	g[0]=0;
	b[0]=0;
	}
     }
int output_color_format ( Alignment *B,Alignment *S,char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *))
    {
    int a, b, c;
    int max_name_len=15;
    int max_len=0;
    
    static char *buf;
    int s;
    int *n_residues;
    static FILE_format *fps;
    Color *ink;
    Color *box_c;
    Color *white;
    
    

    
    box_c=vcalloc ( 1, sizeof (Color));    
    get_rgb_values_format (DEFAULT_COLOR, (white=vcalloc ( 1, sizeof (Color))));	
    get_rgb_values_format (INK_COLOR,     (ink  =vcalloc ( 1, sizeof (Color))));

    n_residues=vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a<B->nseq; a++)n_residues[a]=B->order[a][1];

    fps=vfopen_format( name);      
    if ( buf==NULL)
	{
	buf=vcalloc (10000, sizeof (int));
	}

    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	}	
    if ( max_len>max_name_len)max_len=max_name_len;
    
   sprintf (buf, "\n%s, %s(%s)\n%s\n",PROGRAM,VERSION,DATE, AUTHOR);     
   fps=print_format_string ( buf,white, ink, fps);

   fps=print_format_string ( "\n\n",white,ink, fps);
   
   fps->line-=max_len;
   fps->line=fps->line-fps->line%3;
   
   for (a=0; a<B->len_aln; a+=fps->line)
	   {	   
	   
	   if ( (fps->n_line+(B->nseq+4))>fps->max_line_ppage && !((B->nseq+4)>fps->max_line_ppage))
		 {
		 fps=print_format_char ( fps->eop,white, ink, fps);		
		 }
	  
	   for (b=0; b<=S->nseq; b++)
	     {
	     sprintf (buf,"%-*.*s ",max_len+2, max_len,S->name[b]);
	     fps=print_format_string ( buf,white, ink, fps);
	     if(B->output_res_num)
		 {
		 sprintf (buf, " %4d ", n_residues[b]+1);
		 fps=print_format_string ( buf,white, ink, fps);
		 }
	     		 
	     for (fps->in_seq=1,c=a;c<a+fps->line && c<B->len_aln;c++)
		{
		if (b==S->nseq)
		   {
		   n_residues[b]++;
		   get_rgb_values_format (DEFAULT_COLOR,box_c);
		   s=analyse_aln_column ( B, c);
		   }
		else
		  {
		    n_residues[b]+=!is_gap(B->seq_al[b][c]);
		    s=B->seq_al[b][c]; 
		    if (!is_gap(s) && S->seq_al[b][c]!=NO_COLOR_RESIDUE )
		      {
			get_rgb_values_format ( S->seq_al[b][c], box_c);			    			
		      }
		    else
		      {
			get_rgb_values_format (GAP_COLOR, box_c);		   	
		      }
		  }
		fps=print_format_char ( s,box_c, ink,fps);
		}
	      fps->in_seq=0;

	       if(B->output_res_num)
		 {
		 sprintf (buf, " %4d ", n_residues[b]);
		 fps=print_format_string ( buf,white, ink, fps);
		 }

	      fps=print_format_char ( '\n', white, ink, fps);	      
	     
	     }
	    fps=print_format_string ( "\n\n",white, ink, fps);
	    }
    fps=print_format_string ( "\n\n\n",white, ink,fps);
    vfclose_format( fps);
    return 1;

    }

int output_reliability_format ( Alignment *B,Alignment *S,char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *))
    {
    int a, b, c,l;
    int max_name_len=15;
    int max_len=0;
    static char *buf,*buf2;
    int s;
    static FILE_format *fps;
    Color *ink;
    Color *box_c;
    Color *white;
    int *n_residues;
   

    
    box_c=vcalloc ( 1, sizeof (Color));    
    get_rgb_values_format (DEFAULT_COLOR, (white=vcalloc ( 1, sizeof (Color))));	
    get_rgb_values_format (INK_COLOR,     (ink  =vcalloc ( 1, sizeof (Color))));
    
    n_residues=vcalloc ( B->nseq+1, sizeof (int));
    for ( a=0; a<B->nseq; a++)n_residues[a]=B->order[a][1];


    fps=vfopen_format( name);      
    if ( buf==NULL)
	{
	buf=vcalloc (10000, sizeof (int));
	buf2=vcalloc (10000, sizeof (int));
	}

    if ( max_len==0)
	{
	for ( a=0; a< B->nseq; a++)
	    {if ( strlen (B->name[a])>max_len)
		max_len= strlen ( (B->name[a]));
	    }
	}	
    if ( max_len>max_name_len)max_len=max_name_len;
    
    
    
   sprintf (buf, "%s, %s(%s)\n%s\nCPU TIME:%d sec.\n",PROGRAM,VERSION,DATE, AUTHOR,  (B->cpu+get_time())/1000);     
   fps=print_format_string ( buf,white, ink, fps);
   sprintf (buf, "SCORE=%d\n*\n", S->score_aln);
   fps=print_format_string ( buf,white, ink, fps);
   
   sprintf ( buf2, " BAD AVG GOOD");
   l=strlen(buf2);
   get_rgb_values_format ( DEFAULT_COLOR, box_c);
   fps=print_format_char ( buf2[0],box_c, ink, fps);
   for ( a=1; a<l-1; a++)
        {
	get_rgb_values_format ( MIN(9,a-1), box_c);
	fps=print_format_char ( buf2[a],box_c,ink,fps);
	}
   fps=print_format_char ( buf2[a], box_c, ink, fps);
   
   
   fps=print_format_string ( "\n*\n",white,ink, fps);
   
   for ( a=0;S->score_seq && a< B->nseq; a++)
       {
       get_rgb_values_format (S->score_seq[a]/10, box_c);
       sprintf ( buf, "%-*.*s ", max_len+2,max_len, S->name[a]);
       fps=print_format_string ( buf,box_c, ink,fps);
       sprintf ( buf, ": %3d\n", S->score_seq[a]);
       fps=print_format_string ( buf,white, ink,fps);
       }
   fps=print_format_string ( "\n",white, ink,fps);


   
   fps->line-=max_len;
   fps->line=fps->line-(fps->line%3);
   
   for (a=0; a<B->len_aln; a+=fps->line)
	   {	   
	   
	   if ( (fps->n_line+(B->nseq+4))>fps->max_line_ppage && !((B->nseq+4)>fps->max_line_ppage))
		 {
		 fps=print_format_char ( fps->eop,white, ink, fps);		
		 }
	  
	   for (b=0; b<=S->nseq; b++)
	     {
	     if ( b==S->nseq) fps=print_format_string ( "\n",white, ink, fps);
	     sprintf (buf,"%-*.*s ",max_len+2,max_len,S->name[b]);
	     fps=print_format_string ( buf,white, ink, fps);
	     if(B->output_res_num)
		 {
		 sprintf (buf, " %4d ", n_residues[b]+1);
		 fps=print_format_string ( buf,white, ink, fps);
		 }
	     
	     for (fps->in_seq=1,c=a;c<a+fps->line && c<B->len_aln;c++)
		{
		if (b==S->nseq)
		   {
		   if (S->score_seq)get_rgb_values_format ( S->seq_al[b][c],box_c);
		   else get_rgb_values_format (DEFAULT_COLOR,box_c);
		   n_residues[b]++;
		   s=analyse_aln_column ( B, c);
		   }
		else
		   {
		   n_residues[b]+=!is_gap(B->seq_al[b][c]);
		   s=toupper(B->seq_al[b][c]);		   
		   if (!is_gap(s) && S->seq_al[b][c]!=NO_COLOR_RESIDUE )
			{
			get_rgb_values_format ( S->seq_al[b][c], box_c);			    			
			
			}
		   else
		        {
			get_rgb_values_format (GAP_COLOR, box_c);		   	
			
			}
		   
		   }
		fps=print_format_char ( s,box_c, ink,fps);
		}
	      fps->in_seq=0;

	      if(B->output_res_num)
		 {
		 sprintf (buf, " %4d ",n_residues[b]);
		 fps=print_format_string ( buf,white, ink, fps);
		 }

	      fps=print_format_char ( '\n', white, ink, fps);	      
	     
	     }
	    fps=print_format_string ( "\n\n",white, ink, fps);
	    }
    fps=print_format_string ( "\n\n\n",white, ink,fps);
    vfclose_format( fps);
    return 1;

    }


/*****************************************************************************/
/*                       PDF         FUNCTIONS                               */
/*                                                                           */
/*****************************************************************************/
int       output_color_pdf    ( Alignment *B,Alignment *S, char *name)
      {
      char *tmp_name;
      char command[LONG_STRING];

      
#ifndef PS2PDF 
      fprintf (stderr, "\nPDF FORMAT IS NOT SUPPORTED: INSTALL THE PROGRAM PS2PDF\n");
      exit (1);
#else
      tmp_name=vcalloc ( L_tmpnam+1, sizeof (char));
      tmp_name=vtmpnam(tmp_name);
      
      output_color_ps (B, S, tmp_name);
      sprintf ( command, "%s %s %s", PS2PDF, tmp_name, name);
      system  ( command); 
      remove  ( tmp_name); 
      free ( tmp_name);
#endif      
      
      
      return 1;
      }
int       output_reliability_pdf    ( Alignment *B,Alignment *S, char *name)
      {
      char *tmp_name;
      char command[LONG_STRING];

      

#ifndef PS2PDF 
      fprintf (stderr, "\nPDF FORMAT IS NOT SUPPORTED: INSTALL THE PROGRAM PS2PDF\n");
      exit (1);
#else
      tmp_name=vcalloc ( L_tmpnam+1, sizeof (char));
      tmp_name=vtmpnam(tmp_name);
      
      output_reliability_ps (B, S, tmp_name);
      sprintf ( command, "%s %s %s", PS2PDF, tmp_name, name);
      system  ( command); 
      remove  ( tmp_name); 
      free ( tmp_name);
#endif      
      
      
      return 1;
      }
/*****************************************************************************/
/*                       POST SCRIPT FUNCTIONS                               */
/*                                                                           */
/*****************************************************************************/
int       output_color_ps     ( Alignment *B,Alignment *S, char *name)
      {
      output_color_format (B, S, name, vfopen_ps,print_ps_string,print_ps_char,get_rgb_values_ps, vfclose_ps);
      return 1;
      }
int       output_reliability_ps     ( Alignment *B,Alignment *S, char *name)
      {
      output_reliability_format (B, S, name, vfopen_ps,print_ps_string,print_ps_char,get_rgb_values_ps, vfclose_ps);
      return 1;
      }
FILE_format *print_ps_string( char *s, Color *box, Color *ink, FILE_format *fps)
      {
      int l;
      int a;
      
      l=strlen (s);
      
      for ( a=0; a< l; a++)
          {
	  fps=print_ps_char (s[a], box, ink, fps);
	  }
      return fps;
      }
      

FILE_format * print_ps_char ( int c, Color *box, Color *ink, FILE_format *f)
       {
       
       int ch;
       int cw;

       ch=f->font+3;
       cw=f->font-2;
        
       if ( c=='(' || c==')')return f;
       else if (c!='\n' && c!=f->eop)
          {
	  fprintf(f->fp,"%d %d moveto\n", f->x,f->y);
	  fprintf(f->fp,"0 %d rlineto\n%d 0 rlineto\n0 -%d rlineto\nclosepath\n",ch,cw,ch   );
	  fprintf(f->fp,"%3.1f %3.1f %3.1f setrgbcolor\nfill\n%3.1f %3.1f %3.1f setrgbcolor\n", box->r,box->g,box->b, ink->r, ink->g, ink->b);
	  fprintf(f->fp,"%d %d moveto\n(%c) show\n", f->x+1,f->y+3, c);
	  
	  f->x+=cw;
	  }
       else 
          {
	  f->n_line++;
	  if ( f->n_line==f->max_line_ppage || c==f->eop)
	     {
	     
	     f->n_line=0;
	     f->x=f->x0;
	     f->y=f->y0;
	     fprintf(f->fp,"showpage\n");
	     f->n_pages++;
	     fprintf ( f->fp, "%c%cPage:  %d %d\n",'%', '%', f->n_pages, f->n_pages);	    
	     }
	  else
	     {
	     f->x=f->x0;
	     f->y-=ch;	  
	     }
	  }
       return f;
       }
void get_rgb_values_ps ( int val, Color *C)
     {
     get_rgb_values ( val, C);
     }



FILE_format* vfopen_ps ( char *name)
      {
      FILE_format*fps;

      fps=vcalloc ( 1, sizeof ( FILE_format));
      fps->font=9;
      fps->max_line_ppage=60;
      fps->line=60;/*N char per line*/
      fps->x0=15;
      fps->y0=750;
      fps->eop='^';
      
      fps->fp=vfopen ( name, "w");
      fprintf(fps->fp,"%%!PS-Adobe-2.0\n/Courier findfont\n%d scalefont\nsetfont\n",fps->font);
      fprintf(fps->fp, "%%%%Pages: (atend)\n");
      fprintf(fps->fp,"newpath\n"); 
      ++(fps->n_pages);
      fprintf (fps->fp, "%%%%Page:  %d %d\n", fps->n_pages, fps->n_pages);
      fprintf (fps->fp,"%d %d translate\n",fps->x0, fps->y0);
      return fps;
      }

FILE_format* vfclose_ps ( FILE_format *fps)
      {
      
      fprintf(fps->fp,"showpage\n");	  
      fprintf ( fps->fp, "%%%%Pages:  %d\n", fps->n_pages);
      fprintf(fps->fp,"%%%%EOF"); 
      fprintf(fps->fp,"%%%%\n"); 
      vfclose ( fps->fp);
      free (fps);
      return NULL;
      }
/*****************************************************************************/
/*                       HTML FUNCTIONS                               */
/*                                                                           */
/*****************************************************************************/
int       output_color_html     ( Alignment *B,Alignment *S, char *name)
      {
      output_color_format (B, S, name, vfopen_html,print_html_string,print_html_char,get_rgb_values_html, vfclose_html);
      return 1;
      }
int       output_reliability_html     ( Alignment *B,Alignment *S, char *name)
      {
      output_reliability_format (B, S, name, vfopen_html,print_html_string,print_html_char,get_rgb_values_html, vfclose_html);
      return 1;
      }
FILE_format *print_html_string( char *s, Color *box, Color *ink, FILE_format *fhtml)
      {
      int l;
      int a;
      
      l=strlen (s);
      
      for ( a=0; a< l; a++)
          {
	  fhtml=print_html_char (s[a], box, ink, fhtml);
	  }
      fhtml=print_html_char (CLOSE_HTML_SPAN,NULL,NULL,fhtml);
      return fhtml;
      }
      

FILE_format * print_html_char ( int c, Color *box, Color *ink, FILE_format *f)
       {
       char html_color[100];
       int in_span, new_color;
       char string[1000];


        if (c==CLOSE_HTML_SPAN)
	 {
	   if (f->in_html_span)fprintf ( f->fp, "</span>");
	   f->in_html_span=0;
	   return f;
	 }
	
	
       in_span=f->in_html_span;
       new_color=1-(strm (box->html_color_class, f->previous_html_color));
     
       
 
       sprintf (f->previous_html_color, "%s", box->html_color_class);
       sprintf ( html_color, "class=%s", box->html_color_class);

       
       if ( c!=' ')sprintf ( string, "%c", c);
       else sprintf ( string, "&nbsp;");
       
       if ( !in_span &&                  c!='\n' && c!=f->eop)
          {
	  fprintf ( f->fp, "<span %s>%s",html_color,string );
	  f->in_html_span=1;
	  }
       else if (in_span && !new_color && c!='\n' && c!=f->eop)
	  {
	 
	  fprintf ( f->fp, "%s",string);
	  }
       else if (in_span &&  new_color && c!='\n' && c!=f->eop)
	  {
	  fprintf ( f->fp, "</span><span %s>%s",html_color,string);
	  }	   
       else if ( c=='\n')
          {
	  if ( f->in_html_span)fprintf ( f->fp, "</span>");
	  fprintf ( f->fp, "<br>");
	  sprintf ( f->previous_html_color, "no_color_set");
	  f->in_html_span=0;
	  f->n_line++;
	  }
       
       
       
       
     
       return f;
       }
      
void get_rgb_values_html ( int val, Color *C)
     {
     get_rgb_values ( val, C);
     }

FILE_format* vfopen_html ( char *name)
      {
      FILE_format*fhtml;
      Color *color;
      int a;

      color=vcalloc ( 1, sizeof (Color));

      fhtml=vcalloc ( 1, sizeof ( FILE_format));
      fhtml->font=11;
      fhtml->max_line_ppage=100000;
      fhtml->line=70;/*N char per line*/
      fhtml->x0=15;
      fhtml->y0=800;
      fhtml->eop='^';
      sprintf ( fhtml->previous_html_color, "no_value_set");
      fhtml->fp=vfopen ( name, "w");

      fprintf(fhtml->fp,"<html>\n<style>\n");
     
      fprintf(fhtml->fp,"SPAN { font-family: courier new, courier-new, courier; font-weight: bold; font-size: %dpt;}\n", fhtml->font);
      fprintf(fhtml->fp,"SPAN { line-height:100%%}\n");
      fprintf(fhtml->fp,"SPAN {	white-space: pre}\n");
      
      for ( a=0; a< 10; a++)
          {
	  get_rgb_values_html ( a, color);    
	  if ( !strm (color->html_color, ""))fprintf (fhtml->fp, "SPAN.%s {background: %s}\n", color->html_color_class, color->html_color );
	  else fprintf (fhtml->fp, "SPAN.%s {}\n", color->html_color_class );
	  }
      get_rgb_values_html (DEFAULT_COLOR, color);
      if ( !strm (color->html_color, ""))fprintf (fhtml->fp, "SPAN.%s {background: %s}\n", color->html_color_class, color->html_color );
      else fprintf (fhtml->fp, "SPAN.%s {}\n", color->html_color_class );
      
      get_rgb_values_html (GAP_COLOR, color);
      if ( !strm (color->html_color, ""))fprintf (fhtml->fp, "SPAN.%s {background: %s}\n", color->html_color_class, color->html_color );
      else fprintf (fhtml->fp, "SPAN.%s {}\n", color->html_color_class );
       
      get_rgb_values_html (INK_COLOR, color);
      if ( !strm (color->html_color, ""))fprintf (fhtml->fp, "SPAN.%s {background: %s}\n", color->html_color_class, color->html_color );
      else fprintf (fhtml->fp, "SPAN.%s {}\n", color->html_color_class );
      
      
      
      fprintf(fhtml->fp,"</style>");
      fprintf(fhtml->fp,"<body>");
      
      return fhtml;
      }
FILE_format* vfclose_html ( FILE_format *fhtml)
      {
      if ( fhtml->in_html_span)fprintf(fhtml->fp,"</span>");
      fprintf(fhtml->fp,"</body></html>\n");	  
      vfclose ( fhtml->fp);
      free (fhtml);
      return NULL;
      }
/*****************************************************************************/
/*                       ascii FUNCTIONS                               */
/*                                                                           */
/*****************************************************************************/
int       output_reliability_ascii     ( Alignment *B,Alignment *S, char *name)
      {
      output_reliability_format (B, S, name, vfopen_ascii,print_ascii_string,print_ascii_char,get_rgb_values_ascii, vfclose_ascii);
      return 1;
      }
FILE_format *print_ascii_string( char *s, Color *box, Color *ink, FILE_format *fascii)
      {
      int l;
      int a;
      
      l=strlen (s);
      
      for ( a=0; a< l; a++)
          {
	  fascii=print_ascii_char (s[a], box, ink, fascii);
	  }
      return fascii;
      }
      

FILE_format * print_ascii_char ( int c, Color *box, Color *ink, FILE_format *f)
       {             
       if (box->ascii_value>=0 && f->in_seq)fprintf ( f->fp, "%c", box->ascii_value);
       else fprintf ( f->fp, "%c",c);
       return f;
       }

      
void get_rgb_values_ascii ( int val, Color *C)
     {
     
     if ( val==NO_COLOR_RESIDUE)C->ascii_value='-';
     else if ( val==NO_COLOR_GAP)C->ascii_value='*';
     else if ( val>9)C->ascii_value='#';     
     else if ( val>=0 && val<=9) C->ascii_value=val+'0';
     else   C->ascii_value=val;
     }

FILE_format* vfopen_ascii ( char *name)
      {
      FILE_format*fascii;

      fascii=vcalloc ( 1, sizeof ( FILE_format));
      fascii->font=11;
      fascii->max_line_ppage=100000;
      fascii->line=70;/*N char per line*/
      fascii->x0=15;
      fascii->y0=800;
      fascii->eop='^';
      fascii->fp=vfopen ( name, "w");
     
      
      return fascii;
      }
FILE_format* vfclose_ascii ( FILE_format *fascii)
      {
       vfclose ( fascii->fp);
      free (fascii);
      return NULL;
      }
