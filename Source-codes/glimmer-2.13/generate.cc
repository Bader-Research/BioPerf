//
//  Programmer:  A. Delcher
//
//        File:  generate.cc
//
//  This program reads an interpolated context model and then
//  generates an artificial genome using the model to create
//  the genes and iid sequence for the intergenic regions.
//  Command-line options specify the genome length, gc composition,
//  and gene-length, intergenic length ranges.
//


#include  "delcher.h"
#include  "gene.h"


const int  MODEL_LEN = 12;
const int  ALPHABET_SIZE = 4;
const int  MAX_NAME_LEN = 256;
const double  MAX_LOG_DIFF = -46.0;    // approx 1e-20 difference in probs


#include  "icm.h"


static int  Genome_Len_g = 0;
  // Minimum length of output string.  Will generally be longer
  // to ensure that last gene and intergenic region is written
static double  Inter_Prob_g [4] = {0.0};
  // Probabilities for a, c, g, t in intergenic regions
static bool  Inter_Probs_Set = false;
  // Set true if values of  Inter_Prob_g  array are set by  -p  option
static char  * Key_Path_g = NULL;
  // Filename of file that holds key to where the genes are
static int  Lo_Gene_Len_g = 0;
static int  Hi_Gene_Len_g = 0;
  // Specify range of lengths of coding regions
static int  Lo_Inter_Len_g = 0;
static int  Hi_Inter_Len_g = 0;
  // Specify range of lengths of intergenic regions
static char  * Model_Path_g = NULL;
  // Filename of the probability model
static double  Reverse_Prob_g = 0.50;
  // Probability for inserting each gene on the reverse strand


void  Emit_Gene_Sequence
    (char * alphabet, int alpha_len, int lo, int hi,
     char * buff, int & len, bool & reversed);
void  Emit_IID_Sequence
    (char * alphabet, int alpha_len, double prob [], int lo, int hi,
     char * buff, int & len);
void  Parse_Command_Line
    (int argc, char * argv []);
void  Read_Probability_Model
    (char *);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   FILE  * key_fp;
   float  * prob;
   char  alphabet [] = "acgt";
   char  tag [1000];
   char  * buff;
   int  gene_id = 0;
   int  i, pos, max_len;

   Parse_Command_Line (argc, argv);

   if  (Genome_Len_g == 0 || Hi_Gene_Len_g == 0)
       {
        fprintf (stderr,
          "Must specify positive genome length with -l option\n"
          "and positive coding length range with -g option\n");
        exit (EXIT_FAILURE);
       }

   max_len = Genome_Len_g + Hi_Gene_Len_g + 3 + Hi_Inter_Len_g + 1;
     // 3 is for stop codon and 1 is for null termination
   buff = (char *) Safe_malloc (max_len);

   Read_Probability_Model (Model_Path_g);
   key_fp = File_Open (Key_Path_g, "w");
   
#if  0
{
 printf ("Model 0:\n");
 prob = Partial_Window_Prob_Dist ("tag", 3, MODEL [0]);
 for  (i = 0;  i < 4;  i ++)
   printf ("%c:  %8.6f\n", alphabet [i], prob [i]);

 printf ("Model 1:\n");
 prob = Partial_Window_Prob_Dist ("tag", 3, MODEL [1]);
 for  (i = 0;  i < 4;  i ++)
   printf ("%c:  %8.6f\n", alphabet [i], prob [i]);

 printf ("Model 2:\n");
 prob = Partial_Window_Prob_Dist ("tag", 3, MODEL [2]);
 for  (i = 0;  i < 4;  i ++)
   printf ("%c:  %8.6f\n", alphabet [i], prob [i]);

 exit (-3);
}
#endif

   if  (! Inter_Probs_Set)
       { // Use independent probabilities of the model
         // for intergenic regions
        
        prob = Partial_Window_Prob_Dist ("a", 1, MODEL [0]);
        for  (i = 0;  i < 4;  i ++)
          Inter_Prob_g [i] += prob [i];

        prob = Partial_Window_Prob_Dist ("a", 1, MODEL [1]);
        for  (i = 0;  i < 4;  i ++)
          Inter_Prob_g [i] += prob [i];

        prob = Partial_Window_Prob_Dist ("a", 1, MODEL [2]);
        for  (i = 0;  i < 4;  i ++)
          Inter_Prob_g [i] += prob [i];

        Inter_Prob_g [0] = Inter_Prob_g [3]
            = (Inter_Prob_g [0] + Inter_Prob_g [3]) / 6.0;
        Inter_Prob_g [1] = Inter_Prob_g [2]
            = (Inter_Prob_g [1] + Inter_Prob_g [2]) / 6.0;
        // Make a & t equiprobable, likewise for c & g
       }
   fprintf (stderr, "Independent Intergenic Probabilities:\n"
            "  a = %6.4f  c = %6.4f  g = %6.4f  t = %6.4f\n",
            Inter_Prob_g [0], Inter_Prob_g [1], Inter_Prob_g [2], Inter_Prob_g [3]);

   // Get prefix sum of probabilities
   for  (i = 1;  i < 4;  i ++)
     Inter_Prob_g [i] += Inter_Prob_g [i - 1];

   Emit_IID_Sequence
       (alphabet, 4, Inter_Prob_g, Lo_Inter_Len_g, Hi_Inter_Len_g,
        buff, pos);

   while  (pos < Genome_Len_g)
     {
      bool  reversed;
      int  new_len;

      Emit_Gene_Sequence
          (alphabet, 4, Lo_Gene_Len_g, Hi_Gene_Len_g, buff + pos,
           new_len, reversed);
      if  (reversed)
          fprintf (key_fp, "%5d  %8d  %8d  [l=%d]\n",
            ++ gene_id, pos + new_len, pos + 1 + 3, new_len - 3);
        else
          fprintf (key_fp, "%5d  %8d  %8d  [l=%d]\n",
            ++ gene_id, pos + 1, pos + new_len - 3, new_len - 3);
#if  0
if  (gene_id == 485)
{
 int  j;

 printf ("> id=%d\n", gene_id);
 for  (j = 0;  j < new_len;  j += 3)
   printf ("%c%c%c\n", buff [pos + j], buff [pos + j + 1], buff [pos + j + 2]);
}
#endif
      pos += new_len;

      Emit_IID_Sequence
          (alphabet, 4, Inter_Prob_g, Lo_Inter_Len_g, Hi_Inter_Len_g,
           buff + pos, new_len);
      pos += new_len;
     }

   sprintf (tag, "model=%s len=%d/%d genelen=%d/%d interlen=%d/%d",
            Model_Path_g, Genome_Len_g, pos, Lo_Gene_Len_g, Hi_Gene_Len_g,
            Lo_Inter_Len_g, Hi_Inter_Len_g);
   Fasta_Print (stdout, buff, tag);

   fclose (key_fp);

   return  0;
  }



void  Emit_Gene_Sequence
    (char * alphabet, int alpha_len, int lo, int hi,
     char * buff, int & len, bool & reversed)

//  Put in  buff  a sequence of characters from  alphabet [0 .. (alpha_len - 1)]
//  where each character is chosen to fit the global ICM probability
//  model  MODEL .  The length of the sequence is uniformly
//  distributed between  lo  and  hi  and put in  len .   buff  is
//  assumed to be large enough to hold the sequence.
//   reversed  is set to indicate if the gene is put in the forward or
//  reverse-complement orientation.

  {
   float  * prob;
   double  x, sum;
   int  i, j, start;

   len = 3 * ((1 + lo + lrand48 () % (1 + hi - lo)) / 3);

   assert (len % 3 == 0);
   if  (len < 3)
       {
        reversed = false;
        return;
       }
   
   strcpy (buff, "atg");   // automatic start codon

   for  (i = 4;  i <= len && i < MODEL_LEN;  i ++)
     {
      prob = Partial_Window_Prob_Dist (buff, i, MODEL [i % 3]);
      x = drand48 ();
      sum = prob [0];
      for  (j = 0;  j < alpha_len - 1 && x > sum;  j ++)
        sum += prob [j + 1];
      buff [i - 1] = alphabet [j];
      if  (i % 3 == 0
             && (strncmp (buff + i - 3, "taa", 3) == 0
                   || strncmp (buff + i - 3, "tag", 3) == 0
                   || strncmp (buff + i - 3, "tga", 3) == 0))
          i --;  // Don't allow stop codons
     }

   for  (start = 0;  i <= len;  start ++, i ++)
     {
      prob = Full_Window_Prob_Dist (buff + start, MODEL [start % 3]);
      x = drand48 ();
      sum = prob [0];
      for  (j = 0;  j < alpha_len - 1 && x > sum;  j ++)
        sum += prob [j + 1];
      buff [i - 1] = alphabet [j];
      if  (i % 3 == 0
             && (strncmp (buff + i - 3, "taa", 3) == 0
                   || strncmp (buff + i - 3, "tag", 3) == 0
                   || strncmp (buff + i - 3, "tga", 3) == 0))
          i --;  // Don't allow stop codons
     }

   // Attach stop codon
   x = drand48 ();
   if  (x < 0.333333)
       strcpy (buff + len, "taa");
   else if  (x < 0.666666)
       strcpy (buff + len, "tag");
     else
       strcpy (buff + len, "tga");
   len += 3;
   buff [len] = '\0';

   reversed = (drand48 () < Reverse_Prob_g);
   if  (reversed)
       Reverse_Complement (buff - 1, len);
   
   return;
  }


void  Emit_IID_Sequence
    (char * alphabet, int alpha_len, double prob [], int lo, int hi,
     char * buff, int & len)

//  Put in  buff  a sequence of characters from  alphabet [0 .. (alpha_len - 1)]
//  where each character is chosen independently with cumulative
//  probability in  prob [] .  The length of the sequence is uniformly
//  distributed between  lo  and  hi  and put in  len .   buff  is
//  assumed to be large enough to hold the sequence.

  {
   int  i;

   len = lo + lrand48 () % (1 + hi - lo);

   for  (i = 0;  i < len;  i ++)
     {
      double  p;
      int  j;

      p = drand48 ();
      for  (j = 0;  j < alpha_len - 1 && p > prob [j];  j ++)
        ;
      buff [i] = alphabet [j];
     }
   buff [len] = '\0';

   return;
  }



void  Parse_Command_Line
    (int argc, char * argv [])

//  Process command-line options and parameters in argv [0 .. argc]

  {
   double  sum;
   int  ch;
   bool  errflg = false;
   char  * p;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "g:hi:l:p:r:")) != EOF))
     switch  (ch)
       {
        case  'g' :
          if  (sscanf (optarg, "%d-%d", & Lo_Gene_Len_g, & Hi_Gene_Len_g)
                       != 2
                 || Lo_Gene_Len_g <= 0 || Hi_Gene_Len_g < Lo_Gene_Len_g)
              {
               fprintf (stderr, "Bad gene length range \"%s\"\n",
                        optarg);
               errflg = true;
              }

          // round up and down, resp, to nearest multiples of 3
          Lo_Gene_Len_g = 3 * ((Lo_Gene_Len_g + 2) / 3);
          Hi_Gene_Len_g = 3 * (Hi_Gene_Len_g / 3);
          if  (Hi_Gene_Len_g < Lo_Gene_Len_g)
              {
               fprintf (stderr, "No multiple of 3 in range \"%s\"\n",
                        optarg);
               errflg = true;
              }
          break;
          
        case  'h' :
          errflg = true;
          break;

        case  'i' :
          if  (sscanf (optarg, "%d-%d", & Lo_Inter_Len_g, & Hi_Inter_Len_g)
                       != 2
                 || Lo_Inter_Len_g <= 0 || Hi_Inter_Len_g < Lo_Inter_Len_g)
              {
               fprintf (stderr, "Bad intergenic length range \"%s\"\n",
                        optarg);
               errflg = true;
              }
          break;
          
        case  'l' :
          Genome_Len_g = int (strtol (optarg, & p, 10));
          if  (p == optarg || Genome_Len_g <= 0)
              {
               fprintf (stderr, "Bad genome length value \"%s\"\n",
                        optarg);
               errflg = true;
              }
          break;
          
        case  'p' :
          if  (sscanf (optarg, "%lf,%lf,%lf,%lf",
                       Inter_Prob_g + 0, Inter_Prob_g + 1,
                       Inter_Prob_g + 2, Inter_Prob_g + 3)
                       != 4
                 || Inter_Prob_g [0] < 0.0
                 || Inter_Prob_g [1] < 0.0
                 || Inter_Prob_g [2] < 0.0
                 || Inter_Prob_g [3] < 0.0)
              {
               fprintf (stderr, "Bad intergenic probabilities \"%s\"\n",
                        optarg);
               errflg = true;
               break;
              }

          sum = Inter_Prob_g [0] + Inter_Prob_g [1] + Inter_Prob_g [2]
                  + Inter_Prob_g [3];
          if  (sum <= 0.0)
              {
               fprintf (stderr, "Error:  Intergenic probabilities sum to %.3f\n",
                        sum);
               errflg = true;
               break;
              }

          if  (sum < 99.0 || sum > 101.0)
              {
               fprintf (stderr, "Intergenic probabilities sum to %.3f  normalizing\n",
                        sum);
              }
          Inter_Prob_g [0] /= sum;
          Inter_Prob_g [1] /= sum;
          Inter_Prob_g [2] /= sum;
          Inter_Prob_g [3] /= sum;
          Inter_Probs_Set = true;

          break;
          
        case  'r' :
          Reverse_Prob_g = strtod (optarg, & p);
          if  (p == optarg)
              {
               fprintf (stderr, "Bad reverse prob value \"%s\"\n",
                        optarg);
               errflg = true;
              }
          break;
          
        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = true;
       }

   if  (errflg || optind != argc - 2)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Model_Path_g = argv [optind ++];
   Key_Path_g = argv [optind ++];
          

   return;
  }



void  Read_Probability_Model  (char * Param)

//  Read in the probability model indicated by  Param .

  {
   FILE  * fp;

   fp = File_Open (Param, "r");   // maybe rb ?

   Read_Scoring_Model (fp);

   fclose (fp);

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [options] <icm-model-file> <key-file>\n"
           "\n"
           "Create an artificial genome with gene sequences generated\n"
           "according to <icm-model_file>.  Output goes to stdout.\n"
           "<key-file> is the name of a file where the actual positions\n"
           "of the genes will be written\n"
           "\n"
           "Options:\n"
           " -g <lo>-<hi>  Set range of gene lengths\n"
           " -h            Print this message\n"
           " -i <lo>-<hi>  Set range of intergenic lengths\n"
           " -l <num>      Set approx genome length.  Will go beyond this\n"
           "                 until finish gene and intergenic region if\n"
           "                 necessary\n"
           " -p <n1>,<n2>,<n3>,<n4>  Set percent a,c,g,t in intergenic\n"
           "                 to <n1>,<n2>,<n3>,<n4> resp\n"
           " -r <num>      Set prob of reverse-compl gene to <num>\n"
           "\n",
           command);

   return;
  }



