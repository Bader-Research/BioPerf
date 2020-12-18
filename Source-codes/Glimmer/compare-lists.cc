//  A. L. Delcher
//
//        File:  Glimmer/compare-lists.c
//
//     Version:  2.03  9 Dec 2002
//
//  Compare lists of coordinates from files on command line
//  and report how well they matched.


#include  "delcher.h"


const int  MAX_LINE_LEN = 1024;
const int  SHORT_ANSWER_LEN = 90;

struct  Gene_Rec
  {
   char  * Tag, * Note;
   long int  Lo, Hi, Start;
   signed char  Matched, Start_Frame, Stop_Frame;
  };


int  Cmp  (const void * A, const void * B)
  {
   if  (((Gene_Rec *) A) -> Hi < ((Gene_Rec *) B) -> Hi)
       return  -1;
   else if  (((Gene_Rec *) A) -> Hi > ((Gene_Rec *) B) -> Hi)
       return  1;
     else
       return  0;
  }

void  Read_Genes  (FILE *, Gene_Rec * &, long int &);

long int  Genome_Len = 0;



int  main
    (int argc, char * argv [])

  {
   FILE  * Infile;
   Gene_Rec  * Predict, * Answer;
   long int  Predict_Ct, Answer_Ct, Match_Ct, Head_Len, Tail_Len;
   long int  Exact_Match_Ct = 0, Early_Start_Ct = 0, Early_Start_Sum = 0;
   long int  Late_Start_Ct = 0, Late_Start_Sum = 0;
   long int  Len, Long_Unmatched_Answer_Ct = 0;
   long int  Unmatched_Predict_Ct = 0, Unmatched_Answer_Ct = 0;
   long int  Extra_Match_Ct = 0, Matched_Predict_Ct = 0;
   long int  i, j, Stop;
   
   if  (argc < 4)
       {
        fprintf (stderr, "USAGE:  %s <predict-file> <answer-file> <genome-len>\n",
                    argv [0]);
        exit (-1);
       }

   for  (i = 4;  i < argc;  i ++)
     {
      if  (argv [i] [0] != '-')
          fprintf (stderr, "Unrecognized option %s\n", argv [i]);
        else
          {
           switch  (argv [i] [1])
             {
              default :
                fprintf (stderr, "Unrecognized option %s\n", argv [i]);
             }
          }
     }

   Genome_Len = strtol (argv [3], NULL, 10);
   fprintf (stderr, "Genome Length = %ld\n", Genome_Len);
   Infile = File_Open (argv [1], "r");
   Read_Genes (Infile, Predict, Predict_Ct);
   qsort (Predict, Predict_Ct, sizeof (Gene_Rec), Cmp);

#if  0
   printf ("\nPredicted Genes:\n");
   for  (i = 0;  i < Predict_Ct;  i ++)
     printf ("%8ld  %8ld  %8ld  %2d\n", Predict [i] . Lo, Predict [i] . Hi,
                 Predict [i] . Start, int (Predict [i] . Start_Frame));
#endif

   Infile = File_Open (argv [2], "r");
   Read_Genes (Infile, Answer, Answer_Ct);
   qsort (Answer, Answer_Ct, sizeof (Gene_Rec), Cmp);
   for  (i = 0;  i < Answer_Ct;  i ++)
     Answer [i] . Matched = 0;

#if  0
   printf ("\nAnswer Genes:\n");
   for  (i = 0;  i < Answer_Ct;  i ++)
     printf ("%8ld  %8ld  %8ld  %2d\n", Answer [i] . Lo, Answer [i] . Hi,
                 Answer [i] . Start, int (Answer [i] . Start_Frame));
#endif

   for  (i = 0;  i < Predict_Ct;  i ++)
     {
      Match_Ct = 0;
#if  0
      if  (strstr (Predict [i] . Note, "Shadowed") != NULL
//             || strstr (Predict [i] . Note, "Shorter") != NULL
//           || strstr (Predict [i] . Note, "Bad") != NULL
          )
          continue;
#endif

      if  (Predict [i] . Start_Frame != Predict [i] . Stop_Frame)
          fprintf (stderr,
                   "*** Predict [%ld] = %8ld %8ld %8ld %2d %2d has frame shift\n",
                   i, Predict [i] . Lo, Predict [i] . Hi,
                   Predict [i] . Start, int (Predict [i] . Start_Frame),
                   int (Predict [i] . Stop_Frame));

      for  (j = 0;  j < Answer_Ct;  j ++)
        {
         if  (Predict [i] . Lo > Answer [j] . Hi)
             /* Do nothing */ ;
         else if  (Predict [i] . Hi >= Answer [j] . Lo)
           {
            if  (Predict [i] . Start_Frame == Answer [j] . Start_Frame
                     || Predict [i] . Start_Frame == Answer [j] . Stop_Frame)
                {
                 Answer [j] . Matched ++;
                 Head_Len = Answer [j] . Lo - Predict [i] . Lo;
                 Tail_Len = Answer [j] . Hi - Predict [i] . Hi;
                 if  (Head_Len == 0)
                     {
                      if  (Tail_Len == 0)
                          {
                           Exact_Match_Ct ++;
                           Match_Ct ++;
                           continue;
                          }
                      else if  (Tail_Len < 0)
                          {
                           if  (Answer [j] . Start_Frame < 0)
                               {
                                Early_Start_Ct ++;
                                Early_Start_Sum += - Tail_Len;
                               }
                             else
                               printf ("*** %s *** contains stop of %s\n",
                                         Predict [i] . Tag, Answer [j] . Tag);
                           Match_Ct ++;
                           continue;
                          }
                        else
                          {
                           if  (Answer [j] . Start_Frame < 0)
                               {
                                Late_Start_Ct ++;
                                Late_Start_Sum += Tail_Len;
                               }
                             else
                               printf ("*** %s *** hits %s with stop\n",
                                         Predict [i] . Tag, Answer [j] . Tag);
                           Match_Ct ++;
                           continue;
                          }
                     }
                 else if  (Head_Len < 0)
                     {
                      if  (Tail_Len == 0)
                          {
                           if  (Answer [j] . Start_Frame < 0)
                               printf ("*** %s *** hits %s with stop\n",
                                         Predict [i] . Tag, Answer [j] . Tag);
                             else
                               {
                                Late_Start_Ct ++;
                                Late_Start_Sum -= Head_Len;
                               }
                           Match_Ct ++;
                           continue;
                          }
                      else if  (Tail_Len < 0)
                          {
                           if  (Answer [j] . Start_Frame < 0)
                               printf ("*** %s *** hits %s with stop\n",
                                         Predict [i] . Tag, Answer [j] . Tag);
                             else
                               printf ("*** %s *** contains stop of %s\n",
                                         Predict [i] . Tag, Answer [j] . Tag);
                           Match_Ct ++;
                           continue;
                          }
                        else
                          {
                           printf ("*** %s *** contained within %s\n",
                                     Predict [i] . Tag, Answer [j] . Tag);
                           Match_Ct ++;
                           continue;
                          }
                     }
                   else
                     {
                      if  (Tail_Len == 0)
                          {
                           if  (Answer [j] . Start_Frame < 0)
                               printf ("*** %s *** contains stop of %s\n",
                                         Predict [i] . Tag, Answer [j] . Tag);
                             else
                               {
                                Early_Start_Ct ++;
                                Early_Start_Sum += Head_Len;
                               }
                           Match_Ct ++;
                           continue;
                          }
                      else if  (Tail_Len < 0)
                          {
                           printf ("*** %s *** completely contains %s\n",
                                     Predict [i] . Tag, Answer [j] . Tag);
                           Match_Ct ++;
                           continue;
                          }
                        else
                          {
                           if  (Answer [j] . Start_Frame < 0)
                               printf ("*** %s *** contains stop of %s\n",
                                         Predict [i] . Tag, Answer [j] . Tag);
                             else
                               printf ("*** %s *** hits %s with stop\n",
                                         Predict [i] . Tag, Answer [j] . Tag);
                           Match_Ct ++;
                           continue;
                          }
                     }
                }
           }
        }
      if  (Match_Ct == 0)
          {
           Unmatched_Predict_Ct ++;
           printf ("%s unmatched  Tag = %s\n", Predict [i] . Tag,
                         Predict [i] . Note);
          }
        else
          Matched_Predict_Ct ++;
     }

   printf ("\nUnmatched Answers:\n");
   for  (i = 0;  i < Answer_Ct;  i ++)
     if  (Answer [i] . Matched == 0)
         {
          Unmatched_Answer_Ct ++;
          if  (Answer [i] . Start_Frame > 0)
              Stop = Answer [i] . Hi;
            else
              Stop = Answer [i] . Lo;
          Len = 1 + Answer [i] . Hi - Answer [i] . Lo;
          printf ("%-8s %8ld %8ld  Len = %4ld\n", Answer [i] . Tag,
                        Answer [i] . Start, Stop,
                        Len);
          if  (Len >= SHORT_ANSWER_LEN)
              Long_Unmatched_Answer_Ct ++;
         }
     else if  (Answer [i] . Matched > 1)
         {
          printf ("%s was matched %d times\n", Answer [i] . Tag,
                  int (Answer [i] . Matched));
          Extra_Match_Ct += Answer [i] . Matched - 1;
         }

   printf (" %5ld Exact Matches\n", Exact_Match_Ct);

   printf (" %5ld Early Starts", Early_Start_Ct);
   if  (Early_Start_Ct > 0)
       printf ("   Avg Length = %.1f\n", double (Early_Start_Sum) / Early_Start_Ct);
     else
       putchar ('\n');
                      
   printf (" %5ld Late Starts", Late_Start_Ct);
   if  (Late_Start_Ct > 0)
       printf ("   Avg Length = %.1f\n", double (Late_Start_Sum) / Late_Start_Ct);
     else
       putchar ('\n');

   printf (" %5ld Total predictions\n", Predict_Ct);
   printf (" %5ld Matched predictions\n", Matched_Predict_Ct);
   printf (" %5ld Unmatched predictions\n", Unmatched_Predict_Ct);
   printf (" %5ld Extra matches\n", Extra_Match_Ct);
   printf (" %5ld Unmatched answers\n", Unmatched_Answer_Ct);
   printf (" %5ld Unmatched answers at least %d long\n", Long_Unmatched_Answer_Ct,
                    SHORT_ANSWER_LEN);

   return  0;
  }



void  Read_Genes  (FILE * fp, Gene_Rec * & A, long int & Ct)

//  Read gene coordinates from  fp  into  A  and set  Ct  to how
//  many there are.

  {
   char  Line [MAX_LINE_LEN];
   char  * P, * Q;
   long int  Len, Start, Stop;

   A = (Gene_Rec *) Safe_malloc (sizeof (Gene_Rec));
   Ct = 0;

   while  (fgets (Line, MAX_LINE_LEN, fp) != NULL)
     {
      for  (P = Line;  isspace (* P);  P ++)   // Skip ID string
        ;
      Q = P;
      while  (isgraph (* P))
        P ++;
      A [Ct] . Tag = (char *) Safe_malloc (1 + P - Q);
      strncpy (A [Ct] . Tag, Q, P - Q);
      A [Ct] . Tag [P - Q] = '\0';

      Start = strtol (P, & Q, 10);
      if  (P == Q)
          {
           fprintf (stderr, "ERROR:  Bad number  Ct = %ld  P = %s\n", Ct, P);
           exit (-1);
          }

      P = Q;
      Stop = strtol (P, & Q, 10);
      if  (P == Q)
          {
           fprintf (stderr, "ERROR:  Bad number  Ct = %ld  P = %s\n", Ct, P);
           exit (-1);
          }

      while  (isspace (* Q))                   // Point to comment
        Q ++;
      A [Ct] . Note = strdup (Q);

      Len = abs (1 + Stop - Start);
      if  (Len > Genome_Len / 2)               // Wraparound
          {
           if  (Start < Stop)
               Start += Genome_Len;
             else
               Stop += Genome_Len;
          }

      A [Ct] . Start = Start;
      if  (Start < Stop)
          {
           A [Ct] . Lo = Start;
           A [Ct] . Hi = Stop;
           A [Ct] . Start_Frame = 1 + ((Start + 2) % 3);
           A [Ct] . Stop_Frame = 1 + (Stop % 3);
          }
        else
          {
           A [Ct] . Lo = Stop;
           A [Ct] . Hi = Start;
           A [Ct] . Start_Frame = - (1 + ((Start + 2) % 3));
           A [Ct] . Stop_Frame = - (1 + ((Stop + 1) % 3));
          }

      Ct ++;
      A = (Gene_Rec *) Safe_realloc (A, (Ct + 1) * sizeof (Gene_Rec));
     }
  }
