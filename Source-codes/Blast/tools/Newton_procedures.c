static char const rcsid[] = "$Id: Newton_procedures.c,v 1.1 2005/05/16 16:11:41 papadopo Exp $";

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

File name: Newton_procedures.c

Authors: Yi-Kuo Yu, Alejandro Schaffer, E. Michael Gertz

Contents: Highest level functions to solve the optimization problem 
          for compositional score matrix adjustment

******************************************************************************/
/*
 * $Log: Newton_procedures.c,v $
 * Revision 1.1  2005/05/16 16:11:41  papadopo
 * Initial revision
 *
 */

/* Functions that find the new joint probability matrix given 
   an original joint prob and new set(s) of background frequency(ies)
   The eta parameter is a lagrangian multiplier to
   fix the relative entropy between the new target and the new background

   RE_FLAG is used to indicate various choices as encoded in NRdefs.h
*/

#include <ncbi.h>
#include <NRdefs.h>
#include <re_newton.h>
#include <nlm_numerics.h>

/*allocate one record of type NRitems, allocate memory for its
  arrays, and return a pointer to the record*/
NRitemsPtr
allocate_NR_memory()
{
    NRitems PNTR NRrecord;      /* record to allocate and return */
    Int4 i;                     /* loop index */

    NRrecord = (NRitemsPtr) Nlm_MemNew(1 * sizeof(NRitems));
    NRrecord->first_standard_freq =
        (Nlm_FloatHiPtr) Nlm_MemNew(Alphsize * sizeof(Nlm_FloatHi));
    NRrecord->second_standard_freq =
        (Nlm_FloatHiPtr) Nlm_MemNew(Alphsize * sizeof(Nlm_FloatHi));
    NRrecord->first_seq_freq =
        (Nlm_FloatHiPtr) Nlm_MemNew(Alphsize * sizeof(Nlm_FloatHi));
    NRrecord->second_seq_freq =
        (Nlm_FloatHiPtr) Nlm_MemNew(Alphsize * sizeof(Nlm_FloatHi));
    NRrecord->first_seq_freq_wpseudo =
        (Nlm_FloatHiPtr) Nlm_MemNew(Alphsize * sizeof(Nlm_FloatHi));
    NRrecord->second_seq_freq_wpseudo =
        (Nlm_FloatHiPtr) Nlm_MemNew(Alphsize * sizeof(Nlm_FloatHi));

    NRrecord->score_old   = Nlm_DenseMatrixNew(Alphsize, Alphsize);
    NRrecord->score_final = Nlm_DenseMatrixNew(Alphsize, Alphsize);
    NRrecord->mat_final   = Nlm_DenseMatrixNew(Alphsize, Alphsize);
    NRrecord->mat_b       = Nlm_DenseMatrixNew(Alphsize, Alphsize);
    for(i = 0; i < Alphsize; i++) {
        NRrecord->first_standard_freq[i] = NRrecord->second_standard_freq[i] =
            0.0;
        NRrecord->first_seq_freq[i] = NRrecord->second_seq_freq[i] = 0.0;
        NRrecord->first_seq_freq_wpseudo[i] =
            NRrecord->second_seq_freq_wpseudo[i] = 0.0;
    }
    return (NRrecord);
}

/*free memory assoicated with a record of type NRitems*/
void
free_NR_memory(NRitemsPtr NRrecord)
{
    NRrecord->first_standard_freq =
        Nlm_MemFree(NRrecord->first_standard_freq);
    NRrecord->second_standard_freq =
        Nlm_MemFree(NRrecord->second_standard_freq);

    NRrecord->first_seq_freq  = Nlm_MemFree(NRrecord->first_seq_freq);
    NRrecord->second_seq_freq = Nlm_MemFree(NRrecord->second_seq_freq);

    NRrecord->first_seq_freq_wpseudo =
        Nlm_MemFree(NRrecord->first_seq_freq_wpseudo);
    NRrecord->second_seq_freq_wpseudo =
        Nlm_MemFree(NRrecord->second_seq_freq_wpseudo);

    NRrecord->score_old = 
        Nlm_DenseMatrixFree(NRrecord->score_old);
    NRrecord->score_final = 
        Nlm_DenseMatrixFree(NRrecord->score_final);
    NRrecord->mat_final = 
        Nlm_DenseMatrixFree(NRrecord->mat_final);
    NRrecord->mat_b = 
        Nlm_DenseMatrixFree(NRrecord->mat_b);

    NRrecord = Nlm_MemFree(NRrecord);
}

/*compute Lambda and if flag set according return re_o_newcontext,
  otherwise return 0.0, also test for the possibility of average
  score >= 0*/
Nlm_FloatHi
compute_lambda(NRitemsPtr NRrecord,
               Int4 compute_re,
	       Nlm_FloatHiPtr lambdaToReturn)
{
    Int4 iteration_count;       /* counter for number of iterations of
                                   Newton's method */
    Int4 i, j;                  /* loop indices */
    Nlm_FloatHi sum;            /* used to compute the sum for estimating
                                   lambda */
    Nlm_FloatHi lambda_error = 1.0;     /* error when estimating lambda */
    Nlm_FloatHi lambda;         /* scale parameter of the Extreme Value
                                   Distribution of scores */
    Nlm_FloatHi ave_score;      /* average score in new context */
    Nlm_FloatHi slope;          /* used to compute the derivative when
                                   estimating lambda */
    Nlm_FloatHi re_to_return = 0.0;     /* relative entropy if using old
                                           joint probabilities*/


    *lambdaToReturn = 1.0;

    if(RE_OLDMAT_NEWCONTEXT == NRrecord->flag) {
        ave_score = 0.0;
        for(i = 0; i < Alphsize; i++) {
            for(j = 0; j < Alphsize; j++) {
                ave_score +=
                    NRrecord->score_old[i][j] * NRrecord->first_seq_freq[i] *
                    NRrecord->second_seq_freq[j];
            }
        }
    }
    if((RE_OLDMAT_NEWCONTEXT == NRrecord->flag) &&
       (ave_score >= (-SCORE_BOUND))) {
        /* fall back to no constraint mode when ave score becomes global
           alignment-like */
        NRrecord->flag = RE_NO_CONSTRAINT;

        printf("scoring matrix has nonnegative average score %12.8f,"
               " reset to mode 0 \n", ave_score);
    }

    /* Need to find the relative entropy here. */
    if(compute_re) {
        slope = 0.0;
        lambda = INITIAL_LAMBDA;
        while(slope <= LAMBDA_ERROR_TOLERANCE) {
            /* making sure iteration starting point belongs to nontrivial
               fixed point */
            lambda = 2.0 * lambda;
            for(i = 0; i < Alphsize; i++) {
                for(j = 0; j < Alphsize; j++) {
                    if(RE_OLDMAT_NEWCONTEXT == NRrecord->flag)
                        slope +=
                            NRrecord->score_old[i][j] *
                            exp(NRrecord->score_old[i][j] * lambda) *
                            NRrecord->first_seq_freq[i] *
                            NRrecord->second_seq_freq[j];
                    else
                        slope +=
                            NRrecord->score_final[i][j] *
                            exp(NRrecord->score_final[i][j] * lambda) *
                            NRrecord->first_seq_freq[i] *
                            NRrecord->second_seq_freq[j];
                }
            }
        }

        iteration_count = 0;
        while((fabs(lambda_error) > LAMBDA_ERROR_TOLERANCE) &&
              (iteration_count < LAMBDA_ITERATION_LIMIT)) {
            sum = 0.0;
            slope = 0.0;
            for(i = 0; i < Alphsize; i++) {
                for(j = 0; j < Alphsize; j++) {
                    if(RE_OLDMAT_NEWCONTEXT == NRrecord->flag) {
                        sum +=
                            exp(NRrecord->score_old[i][j] * lambda) *
                            NRrecord->first_seq_freq[i] *
                            NRrecord->second_seq_freq[j];
                        slope +=
                            NRrecord->score_old[i][j] *
                            exp(NRrecord->score_old[i][j] * lambda) *
                            NRrecord->first_seq_freq[i] *
                            NRrecord->second_seq_freq[j];
                    } else {
                        if(RE_NO_CONSTRAINT == NRrecord->flag) {
                            sum +=
                                exp(NRrecord->score_final[i][j] * lambda) *
                                NRrecord->first_seq_freq[i] *
                                NRrecord->second_seq_freq[j];
                            slope +=
                                NRrecord->score_final[i][j] *
                                exp(NRrecord->score_final[i][j] * lambda) *
                                NRrecord->first_seq_freq[i] *
                                NRrecord->second_seq_freq[j];
                        }
                    }
                }
            }
            lambda_error = (1.0 - sum) / slope;
            lambda = lambda + LAMBDA_STEP_FRACTION * lambda_error;
            iteration_count++;
        }
	*lambdaToReturn = lambda;
        printf("Lambda iteration count %d\n", iteration_count );
        printf("the lambda value = %lf \t sum of jp = %12.10f \n", lambda,
               sum);
        re_to_return = 0.0;
        for(i = 0; i < Alphsize; i++) {
            for(j = 0; j < Alphsize; j++) {
                if(RE_OLDMAT_NEWCONTEXT == NRrecord->flag) {
                    re_to_return +=
                        lambda * NRrecord->score_old[i][j] *
                        exp(NRrecord->score_old[i][j] * lambda) *
                        NRrecord->first_seq_freq[i] *
                        NRrecord->second_seq_freq[j];
                } else {
                    if(RE_NO_CONSTRAINT == NRrecord->flag) {
                        re_to_return +=
                            lambda * NRrecord->score_final[i][j] *
                            exp(NRrecord->score_final[i][j] * lambda) *
                            NRrecord->first_seq_freq[i] *
                            NRrecord->second_seq_freq[j];
                    }
                }
            }
        }
    }
    return (re_to_return);
}

/** This procedures helps set up the direct solution of the new
 * score matrix by a Newtonian procedure and converts the results into a 
 * score matrix. This is the highest level procedure
 * shared by the code as it is used both inside and outside BLAST.
 * NRrecord keeps track of the variabales used in the Newtonian
 * optimization; tol is the tolerence to be used to declare
 * convergence to a local optimum; maxits is the maximum number of
 * iterations allowed*/

Int4
score_matrix_direct_solve(NRitems * NRrecord,
                              Nlm_FloatHi tol,
                              Int4 maxits)
{
    Int4 its; /*number of iterations used*/
    Nlm_FloatHi sum; /*sum of joint letter probabilities*/
    Int4 i, j; /*loop indices*/

    /*Is the relative entropy constrained? Behaves as boolean for now*/
    Int4 constrain_rel_entropy =
        RE_NO_CONSTRAINT != NRrecord->flag;

    its =
        OptimizeTargetFrequencies(&NRrecord->mat_final[0][0], Alphsize,
                                  &NRrecord->mat_b[0][0],
                                  NRrecord->first_seq_freq_wpseudo,
                                  NRrecord->second_seq_freq_wpseudo,
                                  constrain_rel_entropy, NRrecord->RE_final,
                                  tol, maxits);

    if(its <= maxits) {
        sum = 0.0;
        for(i = 0; i < Alphsize; i++) {
            for(j = 0; j < Alphsize; j++) {
                sum += NRrecord->mat_final[i][j];
            }
        }
        /*      printf("new joint probabilities sum to %16.12f ,"
           " but then normalized to 1\n",sum); */

        for(i = 0; i < Alphsize; i++) {
            for(j = 0; j < Alphsize; j++) {
                NRrecord->score_final[i][j] =
                    log(NRrecord->mat_final[i][j] / sum /
                        NRrecord->first_seq_freq_wpseudo[i] /
                        NRrecord->second_seq_freq_wpseudo[j]);
            }
        }
    }

    return its;
}

/*compute the symmetric form of the relative entropy of two
 *probability vectors 
 *In this software relative entropy is expressed in "nats", 
 *meaning that logarithms are base e. In some other scientific 
 *and engineering domains where entropy is used, the logarithms 
 *are taken base 2 and the entropy is expressed in bits.
*/
Nlm_FloatHi
Get_RE(Nlm_FloatHi * A,
       Nlm_FloatHi * B)
{
    Int4 i;                     /* loop index over letters */
    Nlm_FloatHi temp;           /* intermediate term */
    Nlm_FloatHi value = 0.0;    /* square of relative entropy */

    for(i = 0; i < Alphsize; i++) {
        temp = (A[i] + B[i]) / 2;
        if(temp > 0) {
            if(A[i] > 0) {
                value += A[i] * log(A[i] / temp) / 2;
            }
            if(B[i] > 0) {
                value += B[i] * log(B[i] / temp) / 2;
            }
        }
    }
    if(value < 0) {             /* must be numerical rounding error */
        value = 0;
    }

    return sqrt(value);
}
