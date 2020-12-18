//
//  Programmer:    Arthur L. Delcher
//  File:          codon-usage.cc
//  Last Updated:  Thu Nov 14 15:25:20 EST 2002
//                
//
//  Purpose:  Read a multi-fasta file of gene transcripts and count
//  the number of occurrences of each codon (i.e., triple of nucleotides)




#include  "delcher.h"
#include  "gene.h"



static void  Parse_Command_Line
    (int argc, char * argv []);
static int  Read_String
    (FILE * fp, char * & s, int & s_size, char * & tag, int & tag_size);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   char  * string = NULL, * tag = NULL;
   int  string_size = 0;
   int  string_num = 0;
   int  tag_size = 0;
   int  bad = 0, total = 0;
   int  ct [64] = {0};
   int  i;

   Parse_Command_Line (argc, argv);

   while  (Read_String (stdin, string, string_size, tag, tag_size))
     {
      int  len;

      string_num ++;
      len = strlen (string);

      for  (i = 0;  i < len;  i += 3)
        {
         int  j;

         j = Codon_To_Subscript (string + i);
         if  (j < 0)
             bad ++;
           else
             ct [j] ++;
         total ++;
        }
     }

   for  (i = 0;  i < 64;  i ++)
     printf ("%s  %9s\n", Subscript_To_Codon (i), Commatize (ct [i]));

   printf ("\n%s  %9s\n", "bad", Commatize (bad));
   printf ("\nTotal = %s\n", Commatize (total));

   return  0;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   int  ch, errflg = FALSE;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "h")) != EOF))
     switch  (ch)
       {
        case  'h' :
          errflg = TRUE;
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   return;
  }



int  Read_String
    (FILE * fp, char * & s, int & s_size, char * & tag, int & tag_size)

//  Read next string from  fp  (assuming FASTA format) into  s [0 .. ]
//  which has  s_size  characters.  Allocate extra memory if needed
//  and adjust  s_size  accordingly.  Return  TRUE  if successful,  FALSE
//  otherwise (e.g., EOF).  Put FASTA header line into  tag [0 .. ]
//  (and adjust  tag_size  if needed).

  {
   int  ch, ct;

   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     ;

   if  (ch == EOF)
       return  FALSE;

   ct = 0;
   while  ((ch = fgetc (fp)) != EOF && ch != '\n' && isspace (ch))
     ;
   if  (ch == EOF)
       return  FALSE;
   if  (ch != '\n' && ! isspace (ch))
       ungetc (ch, fp);
   while  ((ch = fgetc (fp)) != EOF && ch != '\n')
     {
      if  (ct >= tag_size - 1)
          {
           tag_size += INCR_SIZE;
           tag = (char *) Safe_realloc (tag, tag_size);
          }
      tag [ct ++] = char (ch);
     }
   tag [ct ++] = '\0';

   ct = 0;
   while  ((ch = fgetc (fp)) != EOF && ch != '>')
     {
      if  (isspace (ch))
          continue;

      if  (ct >= s_size - 1)
          {
           s_size += INCR_SIZE;
           s = (char *) Safe_realloc (s, s_size);
          }
      s [ct ++] = char (ch);
     }
   s [ct ++] = '\0';

   if  (ch == '>')
       ungetc (ch, fp);

   return  TRUE;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [options]\n"
           "\n"
           "Read multi-fasta sequences from stdin and count the number of each\n"
           "codon in them\n"
           "\n"
           "Options:\n"
           " -h        Print this message\n"
           "\n",
           command);

   return;
  }



