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

File name: Mode_condition.h

Authors: Alejandro Schaffer, Yi-Kuo Yu

Contents: Definitions used only in Mode_condition.c

******************************************************************************/
/*
 * $Log: Mode_condition.h,v $
 * Revision 1.1  2005/05/16 16:11:41  papadopo
 * Initial revision
 *
 */
#ifndef MODE_CONDITION
#define MODE_CONDITION

#define Mode_1_per  0.3
#define Mode_unchange_per 0.6
#define RE_mode_1_limit 0.18

double *Get_bg_freq(char *matrix_name);

/* declaration of function type for future use 
 *
 * variable orders: Queryseq_length, Matchseq_length, query_amino_count,
 * match_amino_account, matrix_name
 *
 * return values for both Test_0 and Test_1
 *      -1: no adjustment; 0: mode 0 (unconstrained);
 *       1: mode 1 (with RE in new context)
 */
typedef Int4 (*Condition) (Int4 , Int4 ,
                           Nlm_FloatHi *, Nlm_FloatHi *, char *);

Int4
TestToApplyREAdjustmentUnconditional(Int4 Len_query,
                                     Int4 Len_match,
                                     Nlm_FloatHi * P_query,
                                     Nlm_FloatHi * P_match,
                                     char *matrix_name);
Int4
TestToApplyREAdjustmentConditional(Int4 Len_query,
                                   Int4 Len_match,
                                   Nlm_FloatHi * P_query,
                                   Nlm_FloatHi * P_match,
                                   char *matrix_name);

Int4
chooseMode(Int4 length1, Int4 length2,
           Nlm_FloatHi * probArray1, Nlm_FloatHi * probArray2,
           char *matrixName, Int4 testFunctionIndex);

#endif
