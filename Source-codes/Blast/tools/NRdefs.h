/* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: NRdefs.h

Authors: E. Michael Gertz, Alejandro Schaffer, Yi-Kuo Yu

Contents: Definitions for Newton's method used in compositional score
          matrix adjustment
******************************************************************************/
/*
 * $Log: NRdefs.h,v $
 * Revision 1.1  2005/05/16 16:11:41  papadopo
 * Initial revision
 *
 */

#ifndef NRDEFS
#define NRDEFS

#define Alphsize 20

/* Next five constants are the options for how relative entropy is
   specified */
enum{ SMITH_WATERMAN_ONLY   = (-1),
      KEEP_OLD_MATRIX       = 0,
      RE_NO_CONSTRAINT      = 1,
      RE_OLDMAT_NEWCONTEXT  = 2,
      RE_OLDMAT_OLDCONTEXT  = 3,
      RE_USER_SPECIFIED     = 4,
      NUM_RE_OPTIONS};

#define NR_ERROR_TOLERANCE     0.00000001   /* bound on error for Newton's
                                               method */
#define NR_ITERATION_LIMIT     2000
#define PROB_SUM_TOLERANCE     0.000000001  /* bound on error for sum of
                                               probabilities*/
#define SCORE_BOUND            0.0000000001

#define LAMBDA_STEP_FRACTION   0.5          /* default step fraction in
                                               Newton's method */
#define INITIAL_LAMBDA         1.0
#define LAMBDA_ITERATION_LIMIT 300
#define LAMBDA_ERROR_TOLERANCE 0.0000001    /* bound on error for estimating
                                               lambda */

typedef struct NRitems {
    Int4 flag;                  /* determines which of the optimization
                                   problems are solved */
    Nlm_FloatHiPtr PNTR mat_b;
    Nlm_FloatHiPtr PNTR score_old;
    Nlm_FloatHiPtr PNTR mat_final;
    Nlm_FloatHiPtr PNTR score_final;
    Nlm_FloatHi RE_final;       /* the relative entropy used, either
                                   re_o_implicit or re_o_newcontext */
    Nlm_FloatHi RE_o_implicit;  /* used for  RE_OLDMAT_OLDCONTEXT mode */
    Nlm_FloatHiPtr first_seq_freq;      /* freq vector of first seq */
    Nlm_FloatHiPtr second_seq_freq;     /* freq. vector for the second. */
    Nlm_FloatHiPtr first_standard_freq; /* background freq vector of first
                                           seq using matrix */
    Nlm_FloatHiPtr second_standard_freq;        /* background freq vector for
                                                   the second. */
    Nlm_FloatHiPtr first_seq_freq_wpseudo;      /* freq vector of first seq
                                                   w/pseudocounts */
    Nlm_FloatHiPtr second_seq_freq_wpseudo;     /* freq. vector for the
                                                   second seq w/pseudocounts */
} NRitems;
typedef NRitems PNTR NRitemsPtr;

NRitemsPtr allocate_NR_memory();
void free_NR_memory(NRitemsPtr NRrecord);

Nlm_FloatHi Get_RE(Nlm_FloatHi *A, Nlm_FloatHi *B);


Nlm_FloatHi 
RE_interface(char *matrixName,
             Int4 length1, Int4 length2,
             Nlm_FloatHi *probArray1, Nlm_FloatHi *probArray2,
             Int4 pseudocounts, Nlm_FloatHi specifiedRE, 
             NRitems *NRrecord);


Int4 initializeNRprobabilities(NRitems * NRrecord, char *matrixName);

Int4
score_matrix_direct_solve(NRitems * NRrecord, Nlm_FloatHi tol, Int4 maxits);

Nlm_FloatHi
compute_lambda(NRitemsPtr NRrecord,
               Int4 compute_re,
	       Nlm_FloatHiPtr lambdaToReturn);

#endif
