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

File name: nlm_numerics.h

Author: E. Michael Gertz

Contents: Definitions used in nlm_numerics.c

******************************************************************************/
/*
 * $Log: nlm_numerics.h,v $
 * Revision 1.1  2005/05/16 16:11:41  papadopo
 * Initial revision
 *
 */
#ifndef NLMNUMERICS
#define NLMNUMERICS

#include <ncbistd.h>

Nlm_FloatHiPtr PNTR Nlm_DenseMatrixNew(Int4 nrows, Int4 ncols);
Nlm_FloatHiPtr PNTR Nlm_LtriangMatrixNew(Int4 n);
Nlm_FloatHiPtr PNTR Nlm_DenseMatrixFree(Nlm_FloatHiPtr PNTR mat);

void Nlm_FactorLtriangPosDef(Nlm_FloatHiPtr PNTR A, Int4 n);
void Nlm_SolveLtriangPosDef(Nlm_FloatHiPtr x, Int4 n,
                            Nlm_FloatHiPtr PNTR L );

Nlm_FloatHi Nlm_EuclideanNorm(const Nlm_FloatHi PNTR v, Int4 n);

void Nlm_AddVectors(Nlm_FloatHiPtr y, Int4 n, Nlm_FloatHi alpha,
                    const Nlm_FloatHi PNTR x);

Nlm_FloatHi Nlm_StepBound(const Nlm_FloatHi PNTR x, Int4 n,
                          const Nlm_FloatHi PNTR step_x, Nlm_FloatHi max);

#endif
