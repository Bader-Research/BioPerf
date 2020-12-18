typedef struct
    {
    char *mode;
    char *comments;
    int nseq;
    char **seq_name;
    float **PW_SD;
    float **PW_ID;
    float *SEQ_W;
    }Weights;
typedef struct
    {
    int **list;
    int tot_list;
    int **stem;
    int tot_stem;
    int n_fields;
    int nseq;
    int *len;
    int ***struc;	
    }Structure;

typedef struct
    {
    char **file;
    char **comment;
    char **seq;
    int *len;
    int max_len;
    int min_len;
    int nseq;
    int max_nseq;
    char **name;
    Structure *ST;
    struct Alignment ** Profile;
/*Constraint list*/
    struct Constraint_list *CL;
    int contains_gap;
    char *type;

    
    }Sequence;
typedef struct
    {
    int max_len;
    int alp_size;
    char *alphabet;
    int **count;
    }Profile;

struct Alignment
    {
/*Size*/
    int max_len;
    int min_len;   
    int *  len;
    int declared_len;
    int max_n_seq;
    int nseq;
    int len_aln;

/*Sequence Information*/
    char **file;
    char **comment;
    char **name;
    char **tree_order;
    char **seq_al;
    int  **order;
    Profile *P;
    Sequence *S;
    struct Dp_Result *Dp_result;
    struct Constraint_list *CL;

    int **seq_cache; /*Contains the index of the residues:
		       The sequence Numbering is relative to the sequences, and not to the alignmnent
		       
		       seq_cache[0][1]=3
		       indicates that in the aln residue (0)1 corresponds to [order[0][0]][3]
		       residues: 1...N
		       Sequences 0...M
		     */
    int **cdna_cache; /*Contains the information about wheather a nucleotide is coding or not*/
                     /*Only defined if used */


   
/*Score*/
    int *  score_seq;
    int ** pw_score_seq;
    int ** score_res;
    int score_aln;
    int score;
    	
    int cpu;
    int finished;

/*Input/Output Options*/
    int output_res_num;
    int residue_case; /*1 for lower, 0 for Upper, 2 for keeping unchanged*/

/*Must Not be copied*/
     int used;
     int num;

    };
    
typedef struct Alignment Alignment;
typedef struct
    {
    int in_seq;	
    FILE *fp;
    int font;
    int x0;
    int y0;
    int x;
    int y;
    int n_pages;
    int max_line_ppage;
    int n_line;
    int line;
    int eop;
    int in_html_span;
    char previous_html_color[100];
   
    }
FILE_format;

typedef struct
    {
    float r;
    float g;
    float b;
    char html_color[30];
    char html_color_class[30]; 
    int ascii_value;
    }
Color;


Sequence * fill_sequence_struc ( int nseq, char **sequences, char **seq_name);
Sequence * read_sequences ( char *seq_name);
Sequence * get_sequence_type (Sequence *S);
Alignment* get_aln_type (Alignment *A);
char     * get_string_type   (char *string);


void get_sequence (char *seq_file,int *NSEQ, char ***SEQ, char ***SN, int **sl, int *min, int *max);

int ** get_matrix   ( char *name, char *format);
int ** read_matrice (char *mat_name);
int **neg_matrix2pos_matrix ( int **matrix);   
int ** read_blast_matrix ( char *mat_name);

void   print_aln ( Alignment *B);

int       output_reliability_ps     ( Alignment *B,Alignment *S, char *name);
int       output_reliability_pdf    ( Alignment *B,Alignment *S, char *name);
int       output_reliability_html   ( Alignment *B,Alignment *S, char *name);
int       output_score_ps     ( Alignment *B,Alignment *S, char *name);
int       output_score_pdf    ( Alignment *B,Alignment *S, char *name);
int       output_score_html   ( Alignment *B,Alignment *S, char *name);
void      get_rgb_values(int val, Color *C);
int       output_reliability_format ( Alignment *B,Alignment *S, char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *));
int       output_score_format ( Alignment *B,Alignment *S, char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *));


FILE_format * print_ps_string      ( char *s , Color *box, Color *ink, FILE_format *f);
FILE_format * print_ps_char        ( int   c,    Color *box, Color *ink, FILE_format *f);



void get_rgb_values_ps ( int val, Color *C);
FILE_format* vfopen_ps ( char *name);
FILE_format* vfclose_ps ( FILE_format *fps);

FILE_format *print_html_string( char *s, Color *box, Color *ink, FILE_format *fhtml);
FILE_format * print_html_char ( int c, Color *box, Color *ink, FILE_format *f);
void get_rgb_values_html ( int val, Color *C);
FILE_format* vfopen_html ( char *name);
FILE_format* vfclose_html ( FILE_format *fhtml);

int       output_reliability_ascii     ( Alignment *B,Alignment *S, char *name);
FILE_format *print_ascii_string( char *s, Color *box, Color *ink, FILE_format *fascii);
FILE_format * print_ascii_char ( int c, Color *box, Color *ink, FILE_format *f);
void get_rgb_values_ascii ( int val, Color *C);

FILE_format* vfopen_ascii ( char *name);
FILE_format* vfclose_ascii ( FILE_format *fascii);
/*********************CLUSTALW.H*********************************************/
/****************************************************************************/

   /*
   Main header file for ClustalW.  Uncomment ONE of the following 4 lines
   depending on which compiler you wish to use.
   */

#define VMS 1                 /*VAX or ALPHA VMS */

/*#define MAC 1                 Think_C for MacIntosh */

/*#define MSDOS 1               Turbo C for PC's */

/*#define UNIX 1                Ultrix/Decstation, Gnu C for 
                                Sun, IRIX/SGI, OSF1/ALPHA */

/***************************************************************************/
/***************************************************************************/





#define MAXTITLES		60      /* Title length */

	
#define UNKNOWN   0
#define EMBLSWISS 1
#define PIR 	  2
#define PEARSON   3
#define GDE    	  4
#define CLUSTAL   5	/* DES */
#define MSF       6 /* DES */
#define USER      7	/* DES */

#define PAGE_LEN       22   /* Number of lines of help sent to screen */

#ifdef VMS						/* Defaults for VAX VMS */
#define DIRDELIM ']'		/* Last character before file name in full file 
							   specs */
#define SEQ_MAX_LEN		10000	/* Max Sequence Length */
#define MAXN		500		/* Max Number of Sequences */
#define FSIZE       25000   /* Work space for pairwise alignments */
#define MAXTREE		5000	/* Max Nodes for phylogenetic tree */
#define LINELENGTH    	60  /* Output line length */
#define GCG_LINELENGTH 	50  /* Output line length for GCG output */

#elif MAC
#define DIRDELIM ':'
#define SEQ_MAX_LEN		1000
#define MAXN		30
#define FSIZE       5000
#define MAXTREE		1000
#define LINELENGTH     	50
#define GCG_LINELENGTH 	50


#elif MSDOS
#define DIRDELIM '\\'
#define SEQ_MAX_LEN		1300
#define MAXN		30
#define FSIZE           5000
#define MAXTREE		1000
#define LINELENGTH     	50
#define GCG_LINELENGTH 	50

#elif UNIX
#define DIRDELIM '/'
#define SEQ_MAX_LEN		10000
#define MAXN		500
#define FSIZE           25000
#define MAXTREE		5000
#define LINELENGTH     	60
#define GCG_LINELENGTH 	50
#endif

#define NUMRES 26		/* max size of comparison matrix */

#define INPUT 0
#define ALIGNED 1

#define LEFT 1
#define RIGHT 2

#define NODE 0
#define LEAF 1

#define GAPCOL 32		/* position of gap open penalty in profile */
#define LENCOL 33		/* position of gap extension penalty in profile */

typedef struct node {		/* phylogenetic tree structure */
        struct node *left;
        struct node *right;
        struct node *parent;
        float dist;
        int  leaf;
        int order;
        char name[64];
} stree, *treeptr;

typedef struct tnode *NT_node;

typedef struct tnode{
    char *name;
    NT_node parent;
    NT_node left;
    NT_node right;
    int isseq;
    int seq;
    int nseq;
    int *lseq;
    float dist;
    float dp;
    int order;
    int leaf;
    int group;
    }Treenode;
    
int ** make_sub_tree_list ( NT_node **T, int nseq, int n_node);
void make_all_sub_tree_list ( NT_node N, int **list, int *n);
void make_one_sub_tree_list ( NT_node T, int *list);

NT_node** read_tree(char *treefile, int *nnodes,int nseq, char **seq_names);
FILE * create_tree(NT_node ptree, NT_node parent,int *numseq,int  *ntotal,int  *nnodes,NT_node **lu, FILE *fp);
NT_node declare_tree_node ();
void set_info(NT_node p, NT_node parent, int pleaf, char *pname, float pdist);
NT_node insert_tree_node(NT_node pptr);
FILE * skip_space(FILE *fd);
void create_tree_node(NT_node pptr, NT_node parent);
float calc_mean(NT_node nptr, float *maxdist, int nseq,NT_node **lu);
NT_node insert_root(NT_node p, float diff);
float calc_root_mean(NT_node root, float *maxdist, int neq, NT_node **lu);
NT_node reroot(NT_node ptree, int nseq, int ntotal, int nnodes, NT_node **lu);
/* General purpose header file - rf 12/90 */

#ifndef _H_general
#define _H_general



#define pint int			/* cast ints in printf statements as pint */
typedef int Boolean;			/* Is already defined in THINK_C */

#undef TRUE						
#undef FALSE
#define TRUE 1
#define FALSE 0

#define EOS '\0'				/* End-Of-String */
#define MAXLINE 512			/* Max. line length */


#endif /* ifndef _H_general */
