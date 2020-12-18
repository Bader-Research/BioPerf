
/*DEBUGGING*/
/*#include "mshell.h"*/
/*MEMORY MANAGEMENT*/
#define GIVE_MEMORY_BACK 0
/*OUTPUT DEFINITIONS*/
#define  NO_COLOR_RESIDUE 127
#define  NO_COLOR_GAP 126
#define  CLOSE_HTML_SPAN -1
/*TYPE DEFINITIONS*/

/*SWITCHES*/
#define GO_LEFT            -1
#define GO_RIGHT            1

#define LOCAL            1
#define GLOBAL           2

#define TRUE             1
#define FALSE            0

#define NEW              1
#define OLD              0

#define RANDOM           0
#define DETERMINISTIC    1

#define GREEDY           1
#define NON_GREEDY       0

#define OPTIONAL         1
#define NON_OPTIONAL     0
#define GV_MAXIMISE      1
#define GV_MINIMISE      0

#define ALLOWED          0
#define FORBIDEN         -99999999

/*SIZE DEFINITIONS*/
#define SIZE_OF_INT      10
#define UNDEFINED        FORBIDEN
#define UNDEFINED_INT    UNDEFINED
#define UNDEFINED_FLOAT  UNDEFINED
#define UNDEFINED_DOUBLE UNDEFINED
#define UNDEFINED_CHAR   125
#define UNDEFINED_SHORT  -125
#define UNDEFINED_2      0
#define UNDEFINED_RESIDUE '>'

#define FACTOR           1
#define MAX_N_SEQ        1
#define MAX_N_ALN        1
#define MAX_LEN_ALN      1
#define MAX_N_LIST       100

#define MAXNAMES         100
#define FILENAMELEN 	 500            /* Max. file name length */
#define MAX_N_PARAM      200
#define MAX_PARAM_LEN    200
#define MAX_LINE_LENGTH  10000

#define SHORT_STRING     10
#define STRING           100
#define LONG_STRING      1000
#define VERY_LONG_STRING 10000

#define AA_ALPHABET            "acdefghiklmnpqrstvwy-ACDEFGHIKLMNPQRSTVWY"
#define DNA_ALPHABET           "AGCTUagctu"
#define BLAST_AA_ALPHABET      "arndcqeghilkmfpstwyvbzx*"
#define SIZEOF_AA_MAT   60
#define GAP_LIST         "-.#*~"
#define SSPACE           "   "
#define SCORE_K          10
#define CLEAN_FUNCTION   NULL

#define TRACE_TYPE       int
#define MAX_LEN_FOR_DP   600
#define MATCH            1
#define UNALIGNED        2
#define GAP              3

#define MNE 3

/*CODE SHORT CUTS*/

/*1-COMMAND LINE PROCESSING*/
#define GET_COMMAND_LINE_INFO ((strncmp ( argv[1], "-h",2)==0)||(strncmp ( argv[1], "-man",4)==0)||(strncmp ( argv[1], "-",1)!=0))
#define NEXT_ARG_IS_FLAG ((argc<=(a+1)) ||(( argv[a+1][0]=='-') && !(is_number(argv[a+1]))))


/*UTIL MACROS*/
#define GET_CASE(f,c) ((f==0)?toupper(c):((f==1)?tolower(c):c))

#define SWAP(x,y) {x=x+y;y=x+y; x=y-x; y=y-2*x;}
#define SWAPP(x,y,tp) {tp=y;y=x;x=tp;}

#define MAX(x, y) (((x) >(y)) ? (x):(y))
#define MAX3(x,y,z) (MAX(MAX(x,y),z))
#define MAX4(a,b,c,d) (MAX(MAX(a,b),MAX(c,d)))
#define MIN(x, y) (((x) <(y)) ? (x):(y))
#define FABS(x) ((x<0)?(-x):(x))
#define is_defined(x) ((x==UNDEFINED)?0:1)
#define a_better_than_b(x,y,m) ((m==1)?(((x)>(y))?1:0):(((x)<(y))?1:0))

/*#define bod_a_b(x,y,m)   ((m==1)?(MAX((x),(y))):(MIN((x),(y))))
#define bo_a_b(x,y,m)    ((x==UNEFINED)?y:((y==UNDEFINED)?x:bod_a_b(y,y,m)))
#define best_of_a_b(x,y,m)   ((x==UNDEFINED && y==UNDEFINED)?(UNDEFINED):(bo_a_b(x,y,m)))
*/

#define best_of_a_b(x,y,m) ((m==1)?(MAX((x),(y))):(MIN((x),(y))))

#define strm(x,y)            ((strcmp((x),(y))==0)?1:0)
#define strnm(x,y,n)           ((strncmp((x),(y),(n))==0)?1:0)
#define strm2(a,b,c)         (strm(a,b) || strm(a,c))
#define strm3(a,b,c,d)       (strm2(a,b,c) || strm(a,d))
#define strm4(a,b,c,d,e)     (strm2(a,b,c) || strm2(a,d,e))
#define strm5(a,b,c,d,e,f)   (strm2(a,b,c) || strm3(a,d,e,f))
#define strm6(a,b,c,d,e,f,g) (strm3(a,b,c,d) || strm3(a,e,f,g))
#define declare_name(x) x=vcalloc (MAX(FILENAMELEN,L_tmpnam)+1, sizeof (char)) 
#define is_parameter(x) (x[0]=='-' && !isdigit(x[1])) 

/*Freing functions*/
#define free_2(a, b)            free(a);free(b)
#define free_1(a)               free(a)
#define free_3(a, b, c)         free_2(a,b);free_1(c)
#define free_4(a, b, c,d)       free_2(a,b);free_2(c,d)
#define free_5(a, b, c,d,e)     free_3(a,b,e);free_2(c,d)
#define free_6(a, b, c,d,e,f)   free_3(a,b,e);free_3(c,d,f)
#define free_7(a, b, c,d,e,f,g) free_3(a,b,e);free_4(c,d,f,g)
/*2-FILE PARSING*/
#define SEPARATORS "\n \t,;"
/*END 1-*/


/*WIDOWS/UNIX DISTINCTIONS*/
#if defined(_WIN32) || defined(__WIN32__) ||  defined(__WINDOWS__) || defined(__MSDOS__) || defined(__DOS__) || defined(__NT__) || defined(__WIN32__)
#define WIN32
#define TO_NULL_DEVICE " >nul"
#define    NULL_DEVICE "nul"
#define CWF "/" /*ClustalW Flag*/
#else
#define TO_NULL_DEVICE " >/dev/null 2>&1"
#define    NULL_DEVICE "/dev/null"
#define CWF "-" /*ClustaW Flag*/
#endif

/*Generic Data*/
#define MAIL "cedric.notredame@europe.com"
/*               PROGRAM PATH                  */
#define CLUSTALW_4_TCOFFEE "/home/vipin/compilecodes/TCOFFEE/T-COFFEE_distribution_Version_1.37/bin/clustalw"

#ifdef WIN32
#define LALIGN_4_TCOFFEE "t_coflal"
#else
#define LALIGN_4_TCOFFEE "/home/vipin/compilecodes/TCOFFEE/T-COFFEE_distribution_Version_1.37/bin/lalign2list"
#endif

#define PS2PDF                "ps2pdf"
#define SAP_4_TCOFFEE         "sap"
#define ALIGN_PDB_4_TCOFFEE   "align_pdb"
#define SEQ2MSA_WEIGHT        "seq2msa_weight"


/*               PARAMETER FILES               */
#define COLOR_FILE         "seq_reformat.color"
/*This file specifies the 10 colors available to seq_reformat.
If the file is not on the system, hard coded defaults will be used.
The format is as follow:

-------------------------------------------------------------------------------------------
<Your comments (as many lines as needed >
*
<number in the range 0-9> <HTML code> <R value (float)> <G value (float)> <B value (float)>
-------------------------------------------------------------------------------------------
the RGB values are used for the post-script generation, the html code is used in html documents
*/
/*               Identification    */           
#define DATE "Wed Jul 11 14:38:06 PDT 2001"
#define PROGRAM "T-COFFEE"
#define VERSION "Version_1.37"
#define AUTHOR "Notredame, Higgins, Heringa, JMB(302)pp205-217,2000"
