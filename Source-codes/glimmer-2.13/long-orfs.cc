//
//  Programmer:  A. Delcher
//        File:  ~adelcher/Glimmer/long-orfs.c
//
//     Version:  1.01  29 Oct 97
//     Version:  1.02  25 Feb 98
//               Removed unused variables and parameters.
//               Report length of longest orf to stderr.
//               Allow 0 for min-overlap.
//     Version:  1.03   8 Feb 99
//               Make easier to modify start/stop codons
//     Version:  1.04  10 May 99
//               Add -l switch to turn off circularity of genome
//     Version:  1.1 April 2003 (S. Salzberg)
//               Compute the optimal length for minimum "long"
//               orfs, so that the program will return the largest
//               number of orfs possible.  The -g switch still works
//               if specified, but I don't know why anyone would want
//               to use that for a training set.
//               Also, change min overlap by default to be 0.
//    Copyright (c) 1997-2003 by Arthur Delcher, Steven Salzberg, Simon
//    Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  This program finds long open reading frames in the file named
//  on the command line,  eliminates those that overlap too much,
//  and prints out the resulting coordinates to the standard output.
//

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include  "delcher.h"
#include  "gene.h"
using namespace std;

const char ID_PREFIX [] = "T";
const double  DEFAULT_MIN_OLAP_FRACTION = 0.0;
const int  DEFAULT_MIN_OLAP = 0;
const int  DEFAULT_CHOOSE_FIRST_START_CODON = TRUE;
const int  ORF_SIZE_INCR = 1000;
const long int  PREV_INIT_VAL = LONG_MAX - 10;

const int  MAX_NAME_LEN = 256;


const unsigned int  OK = 0x0;
const unsigned int  SCORES_WORSE = 0x1;
const unsigned int  SHORTER = 0x2;
const unsigned int  SHADOWED_BY = 0x4;
const unsigned int  SHADOWS_ANOTHER = 0x8;
const unsigned int  REJECT_MASK = 0x3;
const char  SCORES_WORSE_CHAR = 'B';
const char  SHORTER_CHAR = 'S';
const char  SHADOWED_BY_CHAR = 'W';
const char  REJECT_CHAR = 'R';

struct  Gene_Ref
  {
   long int  Lo, Hi, Max_Hi, Min_Lo;
   int  Frame;
   unsigned int  Status;
   char  Reject_Code;
  };


void  Check_Previous  (Gene_Ref * &, long int, long int, long int, int);
void  Check_Wraparound  (Gene_Ref * &, long int);
void  Compare_Genes  (Gene_Ref * &, long int, long int, long int, long int,
                      long int);
void  Process_Options  (int, char * []);
void  Transfer  (char *, long int, int);
long  Compute_Optimal_Orf_Len (char *);

int  Choose_First_Start_Codon = DEFAULT_CHOOSE_FIRST_START_CODON;
char  * Data;
long int  Data_Len;
int  Genome_Is_Circular = TRUE;
long int  Longest_Orf = 0;
long int  Min_Gene_Len = -1;   // this either is set in Process_Options or in Compute_Optimal_Orf_Len
int  Min_Olap = DEFAULT_MIN_OLAP;
double  Min_Olap_Fraction = DEFAULT_MIN_OLAP_FRACTION;
char  * Orf_Buffer;
long int  Orf_Buffer_Len;


int  main  (int argc, char * argv [])

  {
   FILE  * fp;
   Gene_Ref  * Previous;
   char  Name [MAX_LINE];
   int  Frame;
   unsigned  Codon;
   long int  For_Prev [3] = {PREV_INIT_VAL, PREV_INIT_VAL, PREV_INIT_VAL};
   long int  Rev_Prev [3] = {PREV_INIT_VAL, PREV_INIT_VAL, PREV_INIT_VAL};
   long int  For_Start [3] = {0};
   long int  Rev_Start [3] = {0};
   long int  ID_Num = 0;
   long int  Virtual_Start, Virtual_End;
   long int  Lo, Hi;
   long int  i, k, Ct, Input_Size, Len, Gene_Len, Start;

   if  (argc < 2)
       {
        fprintf (stderr,
            "USAGE:  %s <genome-file> [options] \n",
            argv [0]);
        exit (-1);
       }

   Process_Options (argc, argv);

   fp = File_Open (argv [1], "r");

   Data = (char *) Safe_malloc (INIT_SIZE);
   Input_Size = INIT_SIZE;

   Read_String (fp, Data, Input_Size, Name, FALSE);
   fclose (fp);

   Data_Len = strlen (Data + 1);
   // convert input genome to lowercase and remove non-ACGT chars
   for  (i = 1;  i <= Data_Len;  i ++)
     Data [i] = Filter (tolower (Data [i]));


   Orf_Buffer_Len = ORF_SIZE_INCR;
   Orf_Buffer = (char *) Safe_malloc (Orf_Buffer_Len);
   Orf_Buffer [0] = ' ';

   Previous = (Gene_Ref *) Safe_malloc (sizeof (Gene_Ref));

   if (Min_Gene_Len < 0) {  // if it was not set on command line
     Min_Gene_Len = Compute_Optimal_Orf_Len (Data);  // compute optimal value to maximum number of ORFs
   }

   printf ("Minimum gene length = %ld\n", Min_Gene_Len);
   printf ("Minimum overlap length = %d\n", Min_Olap);
   printf ("Minimum overlap percent = %.1f%%\n", 100.0 * Min_Olap_Fraction);
   putchar ('\n');

   if  (Genome_Is_Circular)
     // Ch_Mask takes a character and coverts to a 4-bit representation, returns
     // a 32-bit unsigned int with last 4 bits holding our nucleotide
       Codon = Ch_Mask (Data [Data_Len - 1]) << 4
                            | Ch_Mask (Data [Data_Len]);
     else
       Codon = Ch_Mask ('g') << 4 | Ch_Mask ('g');
     // now Codon holds two bases in its last 8 bits.  This will be the first two
     // bases of our 1st codon with Data[1] being the 3rd codon position.

   Frame = 0;
   for  (i = 1;  Data [i] != '\0';  i ++)
     {
       // SHIFT_MASK erases 4 bits representing 1st base in codon
      Codon = (Codon & SHIFT_MASK) << 4;
      Codon |= Ch_Mask (Data [i]);  // store base in last 4 bits of Codon
      Frame = (Frame + 1) % 3;
      if  (Is_Forward_Stop (Codon))
          {
           if  (Genome_Is_Circular || For_Prev [Frame] != PREV_INIT_VAL)
               Len = i - For_Prev [Frame] - 3;
             else
               {
                Len = 3 * ((i / 3) - 1);
                For_Prev [Frame] = i % 3;  // seems unnecessary and wrong
               }
             
           if  (Len > Longest_Orf)
               Longest_Orf = Len;
           if  (Len >= Min_Gene_Len)
                if  (For_Start [Frame] != 0)
                    {
                     Gene_Len = 1 + i - 3 - For_Start [Frame];
                     Start = For_Start [Frame];
                     if  (Gene_Len >= Min_Gene_Len)
                         {
                          ID_Num ++;
                          Check_Previous (Previous, ID_Num, Start,
                                              i - 3, Frame + 1);
                         }
                    }
           For_Prev [Frame] = i;
           For_Start [Frame] = 0;
          }
      if  (Is_Forward_Start (Codon) && For_Start [Frame] == 0)
          For_Start [Frame] = i - 2;

      if  (Is_Reverse_Stop (Codon))
          {
           Len = i - Rev_Prev [Frame] - 3;
           if  (Len > Longest_Orf)
               Longest_Orf = Len;
           if  (Len >= Min_Gene_Len)
                if  (Rev_Start [Frame] != 0)
                    {
                     Gene_Len = Rev_Start [Frame] - Rev_Prev [Frame];
                     Start = Rev_Start [Frame];
                     Gene_Len = Start - Rev_Prev [Frame];  // seems unnecessary
                     if  (Gene_Len >= Min_Gene_Len)
                         {
                          ID_Num ++;
                          Check_Previous (Previous, ID_Num, Rev_Prev [Frame] + 1,
                                   Start, - Frame - 1);
                         }
                    }
           Rev_Prev [Frame] = i;
           Rev_Start [Frame] = 0;
          }
      if  (Is_Reverse_Start (Codon))
          Rev_Start [Frame] = i;
     }

   assert (Data_Len == i - 1);

   if  (! Genome_Is_Circular)
       for  (i = Data_Len;  Data_Len - i < 3;  i --)
         {
          Frame = i % 3;
          if  (Rev_Start [Frame] != 0
                 && (Gene_Len = Rev_Start [Frame] - Rev_Prev [Frame])
                        >= Min_Gene_Len)
              {
               Len = i - Rev_Prev [Frame];
               Start = Rev_Start [Frame];
               Gene_Len = Start - Rev_Prev [Frame];
               if  (Gene_Len >= Min_Gene_Len)
                   {
                    ID_Num ++;
                    Check_Previous (Previous, ID_Num, Rev_Prev [Frame] + 1,
                             Start, - Frame - 1);
                   }
              }
         }

   Previous [ID_Num] . Min_Lo = Previous [ID_Num] . Lo;
   for  (i = ID_Num - 1;  i > 0;  i --)
     Previous [i] . Min_Lo = Min (Previous [i] . Lo, Previous [i + 1] . Min_Lo);

   if  (Genome_Is_Circular)
       {
        Ct = 0;
        for  (i = 1;  Ct < 6 && Data [i] != '\0';  i ++)
          {
           Codon = (Codon & SHIFT_MASK) << 4;
           Codon |= Ch_Mask (Data [i]);
           Frame = (Frame + 1) % 3;
           if  (Is_Forward_Stop (Codon) && i < For_Prev [Frame])
               {
                Len = Data_Len + i - For_Prev [Frame] - 3;
                if  (For_Start [Frame] != 0)
                    {
                     if  (Len > Longest_Orf)
                         Longest_Orf = Len;
                     k = i - 3;
                     Virtual_End = k + Data_Len;
                     if  (i <= 3)
                         k += Data_Len;
                     Gene_Len = i - 2 - For_Start [Frame];
                     if  (For_Start [Frame] > i)
                         Gene_Len += Data_Len;
                     Start = For_Start [Frame];
                     Virtual_Start = Start;
                     Gene_Len = i - 2 - Start;
                     if  (Start > i)
                         Gene_Len += Data_Len;
                       else
                         Virtual_Start += Data_Len;
                     if  (Gene_Len >= Min_Gene_Len)
                         {
                          ID_Num ++;
                          Check_Previous (Previous, ID_Num, Virtual_Start,
                                             Virtual_End, Frame + 1);
                          Check_Wraparound (Previous, ID_Num);
                         }
                    }
                For_Prev [Frame] = i;
                For_Start [Frame] = 0;
                Ct ++;
               }
           if  (Is_Forward_Start (Codon) && For_Start [Frame] == 0)
               {
                if  (i > 2)
                    For_Start [Frame] = i - 2;
                  else
                    For_Start [Frame] = Data_Len + i - 2;
               }

           if  (Is_Reverse_Stop (Codon) && i < Rev_Prev [Frame])
               {
                Len = Data_Len + i - Rev_Prev [Frame] - 3;
                if  (Rev_Start [Frame] != 0)
                    {
                     if  (Len > Longest_Orf)
                         Longest_Orf = Len;
                     k = i - 3;
                     if  (i <= 3)
                         k += Data_Len;
                     Gene_Len = Rev_Start [Frame] - Rev_Prev [Frame];
                     if  (Rev_Start [Frame] < i)
                         Gene_Len += Data_Len;
                     Start = Rev_Start [Frame];
                     Virtual_Start = Start;
                     Gene_Len = Start - Rev_Prev [Frame];
                     if  (Start < i)
                         {
                          Gene_Len += Data_Len;
                          Virtual_Start += Data_Len;
                         }
                     if  (Gene_Len >= Min_Gene_Len)
                         {
                          ID_Num ++;
                          Check_Previous (Previous, ID_Num, Rev_Prev [Frame] + 1,
                                              Virtual_Start, - Frame - 1);
                          Check_Wraparound (Previous, ID_Num);
                         }
                    }
                Rev_Prev [Frame] = i;
                Rev_Start [Frame] = 0;
                Ct ++;
               }
           if  (Is_Reverse_Start (Codon))
               Rev_Start [Frame] = i;
          }
       }
     
   printf ("End = %ld\n", Data_Len);

   printf ("\n\nPutative Genes:\n");
   for  (i = 1;  i <= ID_Num;  i ++)
     {
      if  (Previous [i] . Reject_Code != 'R')
          {
           Lo = Previous [i] . Lo;
           if  (Lo > Data_Len)
               Lo -= Data_Len;
           Hi = Previous [i] . Hi;
           if  (Hi > Data_Len)
               Hi -= Data_Len;
           if  (Previous [i] . Status != 0)
               continue;
           sprintf (Name, "%s%ld", ID_PREFIX, i);
           if  (Previous [i] . Frame > 0)
               printf ("%6s %8ld %8ld", Name, Lo, Hi);
             else
               printf ("%6s %8ld %8ld", Name, Hi, Lo);
//           printf ("  <%04o>", Previous [i] . Status);
           putchar ('\n');
          }
     }

   fprintf (stderr, "Longest orf = %ld\n", Longest_Orf);

   return  0;
  }



void  Check_Previous  (Gene_Ref * & Previous, long int ID,
                       long int Lo, long int Hi, int Frame)

//  Add reference for new gene  ID  at positions  Lo .. Hi  in
//  frame  Frame  to  Previous .  Then check for overlaps with other
//  genes and score them and set their status.  Overlaps must be
//  at least  Min_Olap  long to be reported.   Lo  and  Hi  may be
//  greater than  Data_Len .

  {
   long int  i, Curr_Len;

   Previous = (Gene_Ref *) Safe_realloc (Previous,
                                    (1 + ID) * sizeof (Gene_Ref));

   Previous [ID] . Lo = Lo;
   Previous [ID] . Hi = Hi;
   Previous [ID] . Frame = Frame;
   Previous [ID] . Status = OK;
   Previous [ID] . Reject_Code = ' ';
   if  (ID == 1 || Previous [ID - 1] . Hi < Hi)
       Previous [ID] . Max_Hi = Hi;
     else
       Previous [ID] . Max_Hi = Previous [ID - 1] . Hi;

   Curr_Len = 1 + Hi - Lo;

   for  (i = ID - 1;  i > 0 && Min_Olap <= 1 + Previous [i] . Max_Hi - Lo;  i --)
     Compare_Genes (Previous, i, ID, Lo, Hi, Curr_Len);

   return;
  }



void  Check_Wraparound  (Gene_Ref * & Previous, long int ID)

//  Check for overlaps of the gene in  Previous [ID]  with
//  genes at the beginning of  Previous , i.e., at the beginning
//  of the genome.  Overlaps must be
//  at least  Min_Olap  long to be reported.   Lo  and  Hi  may
//  wraparound the end of the data string.

  {
   int  Frame;
   long int  i, Curr_Hi, Curr_Lo, Curr_Len;

   if  (Previous [ID] . Hi > Data_Len)
       Curr_Hi = Previous [ID] . Hi - Data_Len;
     else
       return;

   Curr_Lo = Previous [ID] . Lo;
   if  (Curr_Lo > Data_Len / 2)
       Curr_Lo -= Data_Len;
   Curr_Len = 1 + Curr_Hi - Curr_Lo;

   Frame = 1 + (Curr_Hi % 3);        // Frame may have changed from wraparound

   for  (i = 1;  i <= ID && 1 + Curr_Hi - Previous [i] . Min_Lo >= Min_Olap;  i ++)
     Compare_Genes (Previous, i, ID, Curr_Lo, Curr_Hi, Curr_Len);

   return;
  }



void  Compare_Genes  (Gene_Ref * & Previous, long int i, long int ID,
                      long int Lo, long int Hi, long int Curr_Len)

//  Check for overlaps between  Previous [i]  and  Previous [ID]
//  where the latter goes from positions  Lo .. Hi  with length
//   Curr_Len .  Overlaps must be at least  Min_Olap  long to be reported.
//   Lo  and  Hi  may be greater than  Data_Len .

  {
   long int  Olap, Prev_Len;
   int  Prev_On_Right;

   // If 100% overlaps are allowed, there is no shadowing
   if  (Min_Olap_Fraction >= 1.0)
       return;

   Prev_Len = 1 + Previous [i] . Hi - Previous [i] . Lo;
   if  (Lo <= Previous [i] . Lo && Hi >= Previous [i] . Hi)
       {
        Previous [i] . Status |= SHADOWED_BY;
        Previous [ID] . Status |= SHADOWS_ANOTHER;
       }
   else if  (Previous [i] . Lo < Lo && Hi < Previous [i] . Hi)
       {
        Previous [ID] . Status |= SHADOWED_BY;
        Previous [i] . Status |= SHADOWS_ANOTHER;
       }
     else
       {
        Prev_On_Right = Hi < Previous [i] . Hi;
        if  (Prev_On_Right)
            Olap = 1 + Hi - Previous [i] . Lo;
          else
            Olap = 1 + Previous [i] . Hi - Lo;
        if  (Olap < Min_Olap
               || Olap <= Min_Olap_Fraction * Curr_Len
                    && Olap <= Min_Olap_Fraction * Prev_Len)
            return;

        if  (Prev_Len <= Curr_Len)
            {
             Previous [i] . Status |= SHADOWED_BY;
             Previous [ID] . Status |= SHADOWS_ANOTHER;
            }
        else
            {
             Previous [ID] . Status |= SHADOWED_BY;
             Previous [i] . Status |= SHADOWS_ANOTHER;
            }
       }

   return;
  }


void  Process_Options  (int argc, char * argv [])

//  Process command-line options and set corresponding global switches
//  and parameters.
//
//    -g n   Set minimum gene length to n
//    -l     Assume genome is linear, i.e., not circular
//    -o n   Set minimum overlap length to n.  Overlaps shorter than this
//           are ignored.
//    -p n   Set minimum overlap percentage to n%.  Overlaps shorter than
//           or equal to this percentage of *both* strings are ignored.

  {
   char  * P;
   long int  W;
   double  X;
   int  i;

   for  (i = 2;  i < argc;  i ++)
     {
      switch  (argv [i] [0])
        {
         case  '-' :
           switch  (argv [i] [1])
             {
              case  'g' :       // minimum gene length
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum gene length %s\n", P);
                  else
                    Min_Gene_Len = W;
                assert (Min_Gene_Len > 0);
                break;
              case  'l' :       // non-circular genome
                Genome_Is_Circular = FALSE;
                break;
              case  'o' :       // minimum overlap length
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum overlap length %s\n", P);
                  else
                    Min_Olap = W;
                assert (Min_Olap >= 0);
                break;
              case  'p' :       // minimum overlap percent
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                X = strtod (P, NULL);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum overlap percent %s\n", P);
                  else
                    Min_Olap_Fraction = X / 100.0;
                assert (Min_Olap_Fraction >= 0.0 && Min_Olap_Fraction <= 1.0);
                break;
              default :
                fprintf (stderr, "Unrecognized option %s\n", argv [i]);
             }
           break;
         case  '+' :
           switch  (argv [i] [1])
             {
              default :
                fprintf (stderr, "Unrecognized option %s\n", argv [i]);
             }
           break;
         default :
           fprintf (stderr, "Unrecognized option %s\n", argv [i]);
        }
     }

   return;
  }



void  Transfer  (char * S, long int Start, int Len)

/* Transfer  |Len|  characters from  Data [Start ...]  to
*  S  and add null terminator.  If  Len > 0 go in forward direction;
*  otherwise, go in reverse direction and use complements.
*  Allow for wraparound at  Data_Len . */

  {
   long int  i, j;

   if  (Len > 0)
       {
        for  (i = 0;  i < Len;  i ++)
          {
           j = Start + i;
           if  (j > Data_Len)
               j -= Data_Len;
           else if  (j < 1)
               j += Data_Len;
           S [i] = Filter (tolower (Data [j]));
          }
        S [i] = '\0';
       }
     else
       {
        for  (i = 0;  i < - Len;  i ++)
          {
           j = Start - i;
           if  (j > Data_Len)
               j -= Data_Len;
           else if  (j < 1)
               j += Data_Len;
           S [i] = Filter (tolower (Complement (Data [j])));
          }
        S [i] = '\0';
       }

   return;
  }

// here's the function to sort our vector of pairs

bool compare (const pair<long, long>x, const pair<long, long>y)
{
  return x.first < y.first;
}

long Compute_Optimal_Orf_Len (char * Data)
/*  compute optimal value to maximum number of ORFs.  Resets the
    global variable Min_Gene_Len with whatever it finds.  Go through
    entire fasta file (Data), find all the orfs, and find lengths of
    the ORFs plus their overlapping ORFs.  An ORF is used only if it
    doesn't overlap another long orf.  So we need to calculate, for
    each ORF, the length of the longest ORF overlapping it.  Example:
    an orf of 500bp overlaps another orf of 400bp.  Then if
    Min_Gene_Len is anywhere from 401 to 500, this ORF will be
    counted.  For each ORF we need to compute an interval like this.
    Store the intervals in a vector, sort the vector and then find the
    length that maximizes the number of ORFs.  */
{  // begin function
   unsigned  Codon;
   int  Frame;
   long For_Prev [3] = {PREV_INIT_VAL, PREV_INIT_VAL, PREV_INIT_VAL};
   long Rev_Prev [3] = {PREV_INIT_VAL, PREV_INIT_VAL, PREV_INIT_VAL};
   long For_Start [3] = {0};
   long Rev_Start [3] = {0};
   long i, j, Len, Gene_Len, Start, leftmost, rightmost;
   long min_gene_len = 200;
   long max_len_overlapping_orf, myLen, hisLen, longest_gene;
   long deepest_coverage, deepest_coverage_length;
   pair<long int, long int>orf_pair;
   vector<pair<long, long> >orf_intervals;   // store locations of all ORFs
   vector<pair<long, long> >orf_lengths;     // store min, max lengths for each orf to 'count'
   vector<int>Coverage_Ct;    

   Codon = Ch_Mask ('g') << 4 | Ch_Mask ('g');
   // now Codon holds two bases in its last 8 bits.  This will be the first two
   // bases of our 1st codon with Data[1] being the 3rd codon position.
   Frame = 0;
   for  (i = 1;  Data [i] != '\0';  i ++)
     {
       Codon = (Codon & SHIFT_MASK) << 4;    // SHIFT_MASK erases 4 bits representing 1st base in codon
       Codon |= Ch_Mask (Data [i]);  // store base in last 4 bits of Codon
       Frame = (Frame + 1) % 3;
       if  (Is_Forward_Stop (Codon))
	 {
           if  (For_Prev [Frame] != PREV_INIT_VAL)  // if we've seen a stop before in this frame
	     Len = i - For_Prev [Frame] - 3;      // compute length from stop to stop
             else
               {
		 Len = 3 * ((i / 3) - 1);   // else length is just entire genome to this point
               }
	   if  (For_Start [Frame] != 0)  // if we have seen a start previously
	     {
	       Gene_Len = 1 + i - 3 - For_Start [Frame];  // get length from start to stop
	       if (Gene_Len > min_gene_len) {
		 orf_pair = make_pair(For_Start[Frame],i-3);  // make a pair out of interval
		 orf_intervals.push_back(orf_pair);   // put it on our vector
	       }
	       Start = For_Start [Frame];  // store position of start codon
	     }
           For_Prev [Frame] = i;   // record this stop codon as 'previous' stop
           For_Start [Frame] = 0;  // start over looking for a start codon in new orf
          }
      if  (Is_Forward_Start (Codon) && For_Start [Frame] == 0)
	For_Start [Frame] = i - 2;   // keep track of the last start codon

      if  (Is_Reverse_Start (Codon))
	Rev_Start [Frame] = i;     // keep track of start codons on other strand

      if  (Is_Reverse_Stop (Codon))
          {
           Len = i - Rev_Prev [Frame] - 3;
	   if  (Rev_Start [Frame] != 0)   // if we have seen a reverse start previously
	     {
	       Gene_Len = Rev_Start [Frame] - Rev_Prev [Frame];  // gene length, opposite strand
	       if (Gene_Len > min_gene_len) {
		 // make a pair from interval, always put lower value first
		 orf_pair = make_pair(Rev_Prev[Frame] + 1, Rev_Start[Frame]);  
		 orf_intervals.push_back(orf_pair);     // store the pair on our vector
	       }
	       Start = Rev_Start [Frame];
	     }
           Rev_Prev [Frame] = i;
           Rev_Start [Frame] = 0;
          }
     }
   sort(orf_intervals.begin(),orf_intervals.end(),compare);  // sort intervals, increasing by left end
   // now scan through all the intervals.  For each one, look backward
   // and forward in the sorted list to find all the overlapping intervals.
   // Only look as far as needed.  Then compute max length of those overlapping
   // intervals, and from this information build a new interval specifying the 
   // range in which this ORF will be 'accepted' by long-orfs.
   longest_gene = 200;
   for (i=0; i< (long) orf_intervals.size(); i++) {
     leftmost = i;
     rightmost = i;
     while (leftmost >= 0 && orf_intervals[i].first <= orf_intervals[leftmost].second) {
       leftmost--;
     }
     leftmost += 1;  // went one to far, now go forward
     // because the intervals are sorted by the 5' end but not the 3' end,
     // we have to go further back and check to see who else might overlap
     // Go back until we are 10,000bp to the left of the current ORF, this
     // will catch almost every gene.  If we find one way back to the left, 
     //it will be much longer than all the intervening ones and will contain them.  
     // So don't worry about excluding those when we're comparing lengths - they
     // don't really overlap but they are guaranteed to be shorter than one
     // that does.
     j = leftmost;
     while (j >= 0 && orf_intervals[j].first + 10000 > orf_intervals[i].first) {
       if (orf_intervals[i].first <= orf_intervals[j].second) {
     	   leftmost = j;
       }
       j--;
     }
     // going to the right is easy because orfs are sorted by 5' end
     while (orf_intervals[i].second >= orf_intervals[rightmost].first &&
	    rightmost <= (long) orf_intervals.size() ) {
       rightmost++;
     }
     rightmost -= 1;  // went one to far, now go back
     // now we've got our range of orfs to look, find min and max length
     myLen = orf_intervals[i].second - orf_intervals[i].first + 1;
     if (myLen > longest_gene) {
       longest_gene = myLen;
     }
     max_len_overlapping_orf = 200;  // shortest ORF we'll ever consider
     for (j = leftmost; j <= rightmost; j++) {
       if (j != i && j >= 0) { // ignore self
	 hisLen = orf_intervals[j].second - orf_intervals[j].first + 1;
	 if (hisLen > max_len_overlapping_orf) { 
	   max_len_overlapping_orf = hisLen;
	 }
       }
     }
     // count this ORF only if it's longer than the longest overlapping ORF
     if (max_len_overlapping_orf < myLen) { 
       orf_pair = make_pair(max_len_overlapping_orf + 1, myLen);
       orf_lengths.push_back(orf_pair);     
     }
   }
   for (i = 0; i < longest_gene + 1; i++) 
     Coverage_Ct.push_back(0);  // fill vector with 0's
   // now count how many intervals cover every orf length
   for (i = 0; i < (long) orf_lengths.size(); i++) {  
     for (j = orf_lengths[i].first; j <= orf_lengths[i].second; j++)
       Coverage_Ct[j]++;    // add 1 for each point that is covered by this ORF
   }
   deepest_coverage = 0;
   deepest_coverage_length = 0;
   for (i = 0; i < longest_gene + 1; i++) {
     if (Coverage_Ct[i] >= deepest_coverage) {  // use >= so the length will be as long as possible
       deepest_coverage = Coverage_Ct[i];
       deepest_coverage_length = i;
     }
   }
   return(deepest_coverage_length);
}

