static char const rcsid[] = "$Id: nlm_numerics.c,v 1.1 2005/05/16 16:11:41 papadopo Exp $";

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

File name: nlm_numerics.c

Author: E. Michael Gertz

Contents: Basic matrix and vector operations for use in conjunction
          with higher-level procedures in re_newton.c

******************************************************************************/
/*
 * $Log: nlm_numerics.c,v $
 * Revision 1.1  2005/05/16 16:11:41  papadopo
 * Initial revision
 *
 */

#include <ncbi.h>
#include <nlm_numerics.h>

/**
 * Create and return a new, dense matrix.  Elements of the matrix A
 * may be accessed as A[i][j]
 *
 * @param nrows     the number of rows for the new matrix.
 * @param ncols     the number of columns for the new matrix.
 */
Nlm_FloatHiPtr PNTR
Nlm_DenseMatrixNew(Int4 nrows,
                   Int4 ncols)
{
    Int4 i;                    /* iteration index */
    Nlm_FloatHiPtr PNTR mat;   /* the new matrix */

    mat = (Nlm_FloatHiPtr PNTR) Nlm_Calloc(nrows, sizeof(Nlm_FloatHiPtr));

    mat[0] =
        (Nlm_FloatHiPtr) Nlm_MemNew((size_t) nrows *
                                    (size_t) ncols * sizeof(Nlm_FloatHi));
    for(i = 1; i < nrows; i++) {
        mat[i] = &mat[0][i * ncols];
    }

    return mat;
}


/**
 * Create and return a new, dense, lower-triangular matrix.  Elements
 * of the matrix A may be accessed as A[i][j] for i <= j.
 *
 * @param n         the dimension of the matrix.
 */
Nlm_FloatHiPtr PNTR
Nlm_LtriangMatrixNew(Int4 n)
{
    Int4 i;                     /* iteration index */
    Nlm_FloatHiPtr PNTR L;      /* the new, lower triangular matrix */
    size_t nelts;               /* the number of elements in
                                   the matrix */

    nelts = ((size_t) n * (n + 1))/2;

    L    = (Nlm_FloatHiPtr PNTR) Nlm_Calloc(n, sizeof(Nlm_FloatHi *));
    L[0] = (Nlm_FloatHiPtr) Nlm_MemNew(nelts * sizeof(Nlm_FloatHi) );

    for( i = 1; i < n; i++ ) {
        L[i] = L[i - 1] + i;
    }

    return L;
}


/**
 * Free a matrix created by Nlm_DenseMatrixNew or
 * Nlm_LtriangMatrixNew.
 *
 * @param mat       the matrix to be freed
 * @return          always NULL
 */
Nlm_FloatHiPtr PNTR
Nlm_DenseMatrixFree(Nlm_FloatHiPtr PNTR mat)
{
    mat[0] = (Nlm_FloatHiPtr) Nlm_MemFree(mat[0]);
    mat    = (Nlm_FloatHiPtr PNTR) Nlm_MemFree(mat);

    return NULL;
}


/**
 * Accessing only the lower triangular elements of the symmetric,
 * positive definite matrix A, compute a lower triangular matrix L
 * such that A = L L^T (Cholesky factorization.)  Overwrite the lower
 * triangle of A with L.
 *
 * This routine may be used with the Nlm_SolveLtriangPosDef routine to
 * solve systems of equations.
 *
 * @param A         the lower triangle of a symmetric, positive-definite
 *                  matrix
 * @param n         the size of A
 */
void
Nlm_FactorLtriangPosDef(Nlm_FloatHiPtr PNTR A, Int4 n)
{
    Int4 i, j, k;               /* iteration indices */
    Nlm_FloatHi temp;           /* temporary variable for intermediate
                                   values in a computation */

    for( i = 0; i < n; i++ ) {
        for( j = 0; j < i; j++ ) {
            temp = A[i][j];
            for( k = 0; k < j; k++ ) {
                temp -= A[i][k] * A[j][k];
            }
            A[i][j] = temp/A[j][j];
        }
        temp = A[i][i];
        for(k = 0; k < i; k++ ) {
            temp -= A[i][k] * A[i][k];
        }
        A[i][i] = sqrt(temp);
    }
}


/**
 * Solve the linear system L L\T y = b, where L is a non-singular
 * lower triangular matrix, usually computed using
 * the Nlm_FactorLtriangPosDef routine.
 *
 * @param x         on entry, the right hand size of the linear system
 *                  L L\T y = b; on exit the solution
 * @param n         the size of x
 * @param L         a non-singular lower triangular matrix
 */
void Nlm_SolveLtriangPosDef(Nlm_FloatHiPtr x, Int4 n,
                            Nlm_FloatHiPtr PNTR L )
{
    Int4 i, j;                  /* iteration indices */
    Nlm_FloatHi temp;           /* temporary variable for intermediate
                                   values in a computation */

    /* At point x = b in the equation L L\T y = b */

    /* Forward solve; L z = b */
    for( i = 0; i < n; i++ ) {
        temp = x[i];
        for( j = 0; j < i; j++ ) {
            temp -= L[i][j] * x[j];
        }
        x[i] = temp/L[i][i];
    }
    /* Now x = z */

    /* Back solve; L\T y = z */
    for( j = n - 1; j >= 0; j-- ) {
        x[j] /= L[j][j];
        for( i = 0; i < j; i++ ) {
            x[i] -= L[j][i] * x[j];
        }
    }
    /* Now x = y, the solution to  L L\T y = b */
}


/**
 * Compute the Euclidean norm (2-norm) of a vector.
 *
 * This routine is based on the (freely available) BLAS routine dnrm2,
 * which handles the scale of the elements of v in a stable fashion.
 *
 * @param v      a vector
 * @param n      the length of v
 */
Nlm_FloatHi
Nlm_EuclideanNorm(const Nlm_FloatHi PNTR v, Int4 n)
{
    Nlm_FloatHi sum   = 1.0;   /* sum of squares of elements in v */
    Nlm_FloatHi scale = 0.0;   /* a scale factor for the elements in v */
    Int4 i;                    /* iteration index */

    for( i = 0; i < n; i++ ) {
        if( v[i] != 0.0 ) {
            Nlm_FloatHi absvi = ABS(v[i]);
            if( scale < absvi ) {
                sum = 1.0 + sum * (scale/absvi) * (scale/absvi);
                scale = absvi;
            } else {
                sum += (absvi/scale) * (absvi/scale);
            }
        }
    }

    return scale * sqrt(sum);
}


/**
 * Let y = y + alpha * x
 *
 * @param y         a vector
 * @param x         another vector
 * @param n         the length of x and y
 */
void Nlm_AddVectors(Nlm_FloatHiPtr y, Int4 n, Nlm_FloatHi alpha,
                    const Nlm_FloatHi PNTR x )
{
    Int4 i;                     /* iteration index */

    for( i = 0; i < n; i++ ) y[i] += alpha * x[i];
}


/**
 * Given a nonnegative vector x and a nonnegative scalar max, returns
 * the largest value in [0, max] for which x + alpha * step_x >= 0.
 *
 * @param x         a vector with nonnegative elements
 * @param step_x    another vector
 * @param n         the size of x and step_x
 * @param max       a nonnegative scalar
 */
Nlm_FloatHi
Nlm_StepBound(const Nlm_FloatHi PNTR x, Int4 n,
              const Nlm_FloatHi PNTR step_x, Nlm_FloatHi max )
{
    Int4 i;                     /* iteration index */
    Nlm_FloatHi alpha;          /* current largest permitted step */

    alpha = max;

    for( i = 0; i < n; i++ ) {
        Nlm_FloatHi alpha_i;      /* a step to the boundary for the
                                   current i */

        alpha_i = -x[i] / step_x[i];
        if( alpha_i >= 0 && alpha_i < alpha ) {
            alpha = alpha_i;
        }
    }

    return alpha;
}
