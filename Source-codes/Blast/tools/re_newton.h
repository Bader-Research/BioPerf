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

File name: re_newton.h

Author: E. Michael Gertz

Contents: Exports for re_newton.c
          Mid-level functions that directly solve the optimization 
          problem for compositional score matrix adjustment.  
          Used in conjunction with Newton_procedures.c and nlm_numerics

******************************************************************************/
/*
 * $Log: re_newton.h,v $
 * Revision 1.1  2005/05/16 16:11:41  papadopo
 * Initial revision
 *
 */

#ifndef RE_NEWTON
#define RE_NEWTON

Int4
OptimizeTargetFrequencies(Nlm_FloatHiPtr x,
                          Int4 alphsize,
                          const Nlm_FloatHi PNTR q,
                          const Nlm_FloatHi PNTR row_sums,
                          const Nlm_FloatHi PNTR col_sums,
                          Int4 constrain_rel_entropy,
                          Nlm_FloatHi relative_entropy,
                          Nlm_FloatHi tol,
                          Int4 maxits);
    
#endif
