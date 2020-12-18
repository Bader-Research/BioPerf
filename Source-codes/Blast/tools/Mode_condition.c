static char const rcsid[] = "$Id: Mode_condition.c,v 1.1 2005/05/16 16:11:41 papadopo Exp $";

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

File name: Mode_condition.c

Authors: Alejandro Schaffer, Yi-Kuo Yu

Contents: Functions to test whether conditional score matrix 
          adjustment should be applied for a pair of matching sequences.

******************************************************************************/
/*
 * $Log: Mode_condition.c,v $
 * Revision 1.1  2005/05/16 16:11:41  papadopo
 * Initial revision
 *
 */

#include <ncbi.h>
#include <NRdefs.h>
#include <Mode_condition.h>

double BLOSUM62_bg[Alphsize] =
    { 0.0742356686, 0.0515874541, 0.0446395713, 0.0536092024, 0.0246865086,
      0.0342500470, 0.0543174458, 0.0741431988, 0.0262119099, 0.0679331197,
      0.0989057232, 0.0581774322, 0.0249972837, 0.0473970070, 0.0385382904,
      0.0572279733, 0.0508996546, 0.0130298868, 0.0322925130, 0.0729201182
    };
  /* BLOSUM 62 is the correct bg */

#define HALF_CIRCLE_DEGREES 180
#define PI 3.1415926543
#define QUERY_MATCH_DISTANCE_THRESHOLD 0.16
#define LENGTH_RATIO_THRESHOLD 3.0
#define ANGLE_DEGREE_THRESHOLD 70

/* declaration of Htype function for future use
 *
 *     typedef int (*Condition) (int , int , int *, int *, char *);
 *
 * variable orders: Queryseq_length, Matchseq_length,
 *                  query_amino_count, match_amino_account, matrix_name
 */

Int4 TestToApplyREAdjustmentUnconditional(Int4,
                                          Int4,
                                          Nlm_FloatHi *,
                                          Nlm_FloatHi *,
                                          char *);


Int4 TestToApplyREAdjustmentConditional(Int4,
                                        Int4,
                                        Nlm_FloatHi *,
                                        Nlm_FloatHi *,
                                        char *);


/* If this function is used relative-entropy score adjustment is
 * always applied, with a fixed value as the target relative entropy*/
Int4
TestToApplyREAdjustmentUnconditional(Int4 Len_query,
                                     Int4 Len_match,
                                     Nlm_FloatHi * P_query,
                                     Nlm_FloatHi * P_match,
                                     char *matrix_name)
{
    return RE_USER_SPECIFIED;
}


/* Decide whether a relative-entropy score adjustment should be used
 * based on lengths and letter counts of the two matched sequences;
 * matrix_name is the underlying score matrix; for now only BLOSUM62
 * is supported */
Int4
TestToApplyREAdjustmentConditional(Int4 Len_query,
                                   Int4 Len_match,
                                   Nlm_FloatHi * P_query,
                                   Nlm_FloatHi * P_match,
                                   char *matrix_name)
{
    Int4 mode_value;            /* which relative entropy mode to return */
    Int4 i;                     /* loop indices */
    Nlm_FloatHi p_query[Alphsize], p_match[Alphsize];   /*letter probabilities
                                                         *for query and match*/
    Nlm_FloatHi *p_matrix;      /* letter probabilities used in constructing
                                 * matrix name*/
    Nlm_FloatHi D_m_mat, D_q_mat, D_m_q;        /* distances between
                                                 * match and original
                                                 * between query and
                                                 * original between
                                                 * match and query*/
    Nlm_FloatHi corr_factor = 0.0;      /* correlation between how
                                           p_query and p_match deviate
                                           from p_matrix */
    Nlm_FloatHi len_q, len_m;   /* lengths of query and matching
                                   sequence in floating point */
    Nlm_FloatHi len_large, len_small;   /* store the larger and smaller of
                                         * len_q and len_m */
    Nlm_FloatHi angle;          /* angle between query and match
                                   probabilities */

    p_matrix = Get_bg_freq(matrix_name);

    for(i = 0; i < Alphsize; i++) {
        p_query[i] = P_query[i];
        p_match[i] = P_match[i];
        corr_factor +=
            (p_query[i] - p_matrix[i]) * (p_match[i] - p_matrix[i]);
    }
    D_m_mat = Get_RE(p_match, p_matrix);
    D_q_mat = Get_RE(p_query, p_matrix);
    D_m_q = Get_RE(p_match, p_query);   /* distance between match and query */

    angle =
        acos((D_m_mat * D_m_mat + D_q_mat * D_q_mat -
              D_m_q * D_m_q) / 2.0 / D_m_mat / D_q_mat);
    /* convert from radians to degrees */
    angle = angle * HALF_CIRCLE_DEGREES / PI;

    len_q = 1.0 * Len_query;
    len_m = 1.0 * Len_match;
    if(len_q > len_m) {
        len_large = len_q;
        len_small = len_m;
    } else {
        len_large = len_m;
        len_small = len_q;
    }

    if((D_m_q > QUERY_MATCH_DISTANCE_THRESHOLD) &&
       (len_large / len_small > LENGTH_RATIO_THRESHOLD) &&
       (angle > ANGLE_DEGREE_THRESHOLD)) {
        mode_value = KEEP_OLD_MATRIX;
    } else {
        mode_value = RE_USER_SPECIFIED;
    }

    return mode_value;
}


/* Retrieve the background letter probabilities implicitly used in
 * constructing the score matrix matrix_name*/
Nlm_FloatHi *
Get_bg_freq(char *matrix_name)
{
    if(0 == strcmp(matrix_name, "BLOSUM62")) {
        return BLOSUM62_bg;
    } else {                    /* default */
        printf("matrix not supported, exit now! \n");
        exit(1);
    }
}


/* initialization of array of functions that can be used to decide
 * which optimization formulation should be used for score
 * adjustment */
Condition Cond_func[] ={ TestToApplyREAdjustmentConditional,
                         TestToApplyREAdjustmentUnconditional,
                         NULL };

/* Choose how the relative entropy should be constrained based on
 * properties of the two sequences to be aligned. length1 an length2
 * are the lengths of the two sequences; probArray1 and probArray2 are
 * arrays of probabilities of letters in each sequence, using the
 * 20-letter alphabet; matrixName is the name of the underlying 20x20
 * score matrix; testFunctionIndex allows different rules to be tested
 * for the relative entropy decision. */
Int4
chooseMode(Int4 length1,
           Int4 length2,
           Nlm_FloatHi * probArray1,
           Nlm_FloatHi * probArray2,
           char *matrixName,
           Int4 testFunctionIndex)
{
    return
        Cond_func[testFunctionIndex] (length1,    length2,
                                      probArray1, probArray2, matrixName);
}
