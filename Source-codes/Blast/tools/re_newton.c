static char const rcsid[] = "$Id: re_newton.c,v 1.1 2005/05/16 16:11:41 papadopo Exp $";

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

File name: re_newton.c

Authors: E. Michael Gertz, Alejandro Schaffer

Contents: Mid-level functions that directly solve the optimization 
          problem for compositional score matrix adjustment.  
          Used in conjunction with Newton_procedures.c and nlm_numerics

******************************************************************************/
/*
 * $Log: re_newton.c,v $
 * Revision 1.1  2005/05/16 16:11:41  papadopo
 * Initial revision
 *
 */
#include <ncbi.h>
#include <nlm_numerics.h>
#include <re_newton.h>

/**
 * @file re_newton.c
 *
 * Author E. Michael Gertz
 * 
 * Routines for finding an optimal set of target frequencies for the
 * purpose of generating a compositionally adjusted score matrix.  The
 * function for performing this optimization is named
 * OptimizeTargetFrequencies (see below).
 *
 * The optimal target frequencies x minimize the Kullback-Liebler
 * "distance"
 *
 *       sum_k x[k] * ln(x[k]/q[k])
 *
 * from the set q of target frequencies from a standard matrix.  They
 * also satisfy the constraints
 *
 * \sum_{i = 0...alphsize - 1} x[i * alphsize + j] = col_sums[j]
 *      for j = 0...alphsize - 1
 *
 * \sum_{j = 0...alphsize - 1} x[i * alphsize + j] = row_sums[j]
 *      for i = 1...alphsize - 1
 *
 * where col_sums and row_sums are sets of background frequencies.
 * Note that the second set of equations, the first index for i is one
 * (the i = 0 case holds automatically if the other constraints are
 * satisfied).
 *
 * The coefficient matrix A of these linear equations is used
 * implicitly in this file.  The routines MultiplyByA,
 * MultiplyByAtranspose and ScaledSymmetricProductA all perform
 * operations with A, without explicitly forming the matrix A.  The
 * indices i and j are related to the column indices k of A by the
 * formula
 *
 *    k = i * alphsize + j.
 *
 * Thus A[j][k] = 1 and A[i + alphsize - 1][k] = 1 if i > 0 and
 * A[ell][k] = 0 otherwise.
 *
 * The target frequencies are also usually subject to a constraint on
 * the relative entropy.
 *
 *  relative_entropy =
 *      sum_{ij} x[i * alphsize + j] *
 *                    log(x[i * alphsize + j]/(row_sums[i] * col_sums[j]))
 *
 * but this constraint is optional.
 *
 * REFERENCES
 *
 * Yi-Kuo Yu, John C Wootton, Stephen F. Altschul (2003) The
 * compositional adjustment of amino-acid substitution matrices. Proc
 * Natl Acad Sci USA. 100, 15688-93.
 *
 * Stephen F. Altschul, John C. Wootton, E. Michael Gertz, Richa
 * Agarwala, Aleksandr Morgulis, Alejandro Schaffer and Yi-Kuo Yu
 * (2005) Protein Database Searches Using Compositionally Adjusted
 * Substitution Matrices.  Submitted to FEBS Journal.
 */

/**
 * Compute the symmetric product A D A^T, where A is the matrix of
 * linear constraints from the problem of generating optimal target
 * frequencies. The matrix A is constant, and is used implicitly in
 * this routine.
 *
 * The symmetric product
 *
 *     A D A\T = \sum_{k = 0}^{n - 1} d_k * A[:][k] * A[:][k]\T,
 *
 * where n = alphsize * alphsize and A[:][k] is column k of A.
 *
 * @param alphsize     the size of the alphabet for this minimization
 *                     problem
 * @param W            the product, a matrix of size 2 * alphsize - 1
 * @param diagonal     a vector that represents the diagonal of D, of
 *                     length alphsize * alphsize
 */
static void
ScaledSymmetricProductA(Nlm_FloatHiPtr PNTR W,
                        Nlm_FloatHiPtr diagonal,
                        Int4 alphsize)
{
    Int4 rowW, colW;   /* iteration indices over the rows and columns of W */
    Int4 i, j;         /* iteration indices over characters in the alphabet */
    Int4 m;            /* The number of rows in A; also the size of W */
    
    m = 2 * alphsize - 1;

    for(rowW = 0; rowW < m; rowW++) {
        for(colW = 0; colW <= rowW; colW++) {
            W[rowW][colW] = 0.0;
        }
    }
    for(i = 0; i < alphsize; i++) {
        for(j = 0; j < alphsize; j++) {
            Nlm_FloatHi dd;     /* an individual diagonal element */

            dd = diagonal[i * alphsize + j];

            W[j][j] += dd;
            if(i > 0) {
                W[i + alphsize - 1][j] += dd;
                W[i + alphsize - 1][i + alphsize - 1] += dd;
            }
        }
    }
}


/**
 * Compute the product y = beta * y + alpha * A * x, where A is the
 * constraint matrix for the problem of generating optimal target
 * frequencies.  If beta == 0.0, then y need not be initialized before
 * calling routine.
 *
 * @param beta         a scalar
 * @param alphsize     the alphabet size for this minimization problem
 * @param y            a vector of size 2 * alphsize - 1
 * @param alpha        a scalar
 * @param x            a vector of size alphsize * alphsize
 */
static void
MultiplyByA(Nlm_FloatHi beta,
            Nlm_FloatHiPtr y,
            Int4 alphsize,
            Nlm_FloatHi alpha,
            const Nlm_FloatHi PNTR x)
{
    Int4 i, j;     /* iteration indices over characters in the alphabet */
    if(beta == 0.0) {
        /* Initialize y to zero, without reading any elements from y */
        for(i = 0; i < 2 * alphsize - 1; i++) {
            y[i] = 0.0;
        }
    } else if(beta != 1.0) {
        /* rescale y */
        for(i = 0; i < 2 * alphsize - 1; i++) {
            y[i] *= beta;
        }
    }
    for(i = 0; i < alphsize; i++) {
        for(j = 0; j < alphsize; j++) {
            y[j] += alpha * x[i * alphsize + j];
        }
    }
    for(i = 1; i < alphsize; i++) {
        for(j = 0; j < alphsize; j++) {
            y[i + alphsize - 1] += alpha * x[i * alphsize + j];
        }
    }
}


/**
 * Compute the product y = beta * y + alpha * A^T * x, where A^T is
 * the transpose of the constraint matrix for the problem of
 * generating optimal target frequencies.  If beta == 0.0, then y need
 * not be initialized before calling routine.
 *
 * @param beta         a scalar
 * @param alphsize     the alphabet size for this minimization problem
 * @param y            a vector of size alphsize * alphsize
 * @param alpha        a scalar
 * @param x            a vector of size 2 * alphsize - 1
 */
static void
MultiplyByAtranspose(Nlm_FloatHi beta,
                     Nlm_FloatHiPtr y,
                     Int4 alphsize,
                     Nlm_FloatHi alpha,
                     const Nlm_FloatHi PNTR x)
{
    Int4 i, j;     /* iteration indices over characters in the alphabet */
    Int4 k;        /* index of a row of A transpose (a column of A); also
                      an index into y */

    if(beta == 0.0) {
        /* Initialize y to zero, without reading any elements from y */
        for(k = 0; k < alphsize * alphsize; k++) {
            y[k] = 0.0;
        }
    } else if(beta != 1.0) {
        /* rescale y */
        for(k = 0; k < alphsize * alphsize; k++) {
            y[k] *= beta;
        }
    }
    for(i = 0; i < alphsize; i++) {
        for(j = 0; j < alphsize; j++) {
            k = i * alphsize + j;

            y[k] += alpha * x[j];
            if(i > 0) {
                y[k] += alpha * x[i + alphsize - 1];
            }
        }
    }
}


/**
 * Calculate the residuals of the linear constraints of the score
 * matrix optimization problem.
 *
 * @param rA           the residual vector of the linear constraints
 * @param alphsize     the alphabet size for this minimization problem
 * @param x            the substitution probabilities
 * @param row_sums     row sums of the substitution probabilities
 * @param col_sums     column sums of the substitution probabilities
 */
static void
ResidualsLinearConstraints(Nlm_FloatHiPtr rA,
                           Int4 alphsize,
                           const Nlm_FloatHi PNTR x,
                           const Nlm_FloatHi PNTR row_sums,
                           const Nlm_FloatHi PNTR col_sums)
{
    Int4 i;             /* iteration index */

    for(i = 0; i < alphsize; i++) {
        rA[i] = col_sums[i];
    }
    for(i = 1; i < alphsize; i++) {
        rA[i + alphsize - 1] = row_sums[i];
    }
    MultiplyByA(1.0, rA, alphsize, -1.0, x);
}


/**
 * Compute the dual residual vector of the optimization problem.  The
 * dual residual vector is also known as the gradient of the
 * Lagrangian function.
 *
 * @param resids_x     the dual residual vector
 * @param alphsize     the alphabet size for this optimization problem
 * @param grads        the gradient of the objective function, an
 *                     possibly the nonlinear relative entropy constraint.
 * @param constrain_rel_entropy    if true, then the relative entropy
 *                                 constraint is used in this optimization
 *                                 problem.
 */
static void
DualResiduals(Nlm_FloatHiPtr resids_x,
              Int4 alphsize,
              Nlm_FloatHiPtr PNTR grads,
              const Nlm_FloatHi PNTR z,
              Int4 constrain_rel_entropy)
{
    Int4 i;                        /* iteration index */
    Int4 n = alphsize * alphsize;  /* size of resids_x */

    if(constrain_rel_entropy) {
        Nlm_FloatHi eta;           /* dual variable for the relative
                                      entropy constraint */
        eta = z[2 * alphsize - 1];
        for(i = 0; i < n; i++) {
            resids_x[i] = -grads[0][i] + eta * grads[1][i];
        }
    } else {
        for(i = 0; i < n; i++) {
            resids_x[i] = -grads[0][i];
        }
    }
    MultiplyByAtranspose(1.0, resids_x, alphsize, 1.0, z);
}


/**
 * Calculate the primal and dual residuals the optimization problem,
 * and their combined Euclidean norm.
 *
 * @param rnorm        the Euclidean norm of the residuals
 * @param resids_x     the dual residual vector (gradient of the
 *                     Lagrangian)
 * @param alphsize     size of the alphabet for this minimization
 *                     problem
 * @param resids_z     the primal residual vector (residuals of the
 *                     original constraints)
 * @param values       a vector containing the values of the nonlinear
 *                     functions
 * @param grads        a matrix whose rows are the gradients of the
 *                     nonlinear functions
 * @param row_sums     row sums of the substitution probabilities
 * @param col_sums     column sums of the substitution probabilities
 * @param relative_entropy          target relative entropy
 * @param x            the substitution probabilities
 * @param z            the dual variables (Lagrange multipliers)
 * @param constrain_rel_entropy    if true, the relative entropy
 *                                 constraint is used
 *
 */
static void
CalculateResiduals(Nlm_FloatHiPtr rnorm,
                   Nlm_FloatHiPtr resids_x,
                   Int4 alphsize,
                   Nlm_FloatHiPtr resids_z,
                   const Nlm_FloatHi PNTR values,
                   Nlm_FloatHiPtr PNTR grads,
                   const Nlm_FloatHi PNTR row_sums,
                   const Nlm_FloatHi PNTR col_sums,
                   const Nlm_FloatHi PNTR x,
                   const Nlm_FloatHi PNTR z,
                   Int4 constrain_rel_entropy,
                   Nlm_FloatHi relative_entropy)
{
    /* Euclidean norms of the primal and dual residuals */
    Nlm_FloatHi norm_resids_z, norm_resids_x;

    DualResiduals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    norm_resids_x = Nlm_EuclideanNorm(resids_x, alphsize * alphsize);

    ResidualsLinearConstraints(resids_z, alphsize, x, row_sums, col_sums);

    if(constrain_rel_entropy) {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];

        norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize);
    } else {
        norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize - 1);
    }
    *rnorm =
        sqrt(norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z);
}


/**
 * Represents the inverse of the matrix of the linear system that must
 * be solved at every iteration of Newton's method.  If J is the
 * Jacobian of all constraints, then the system has the form
 *
 *     (D     J^T)
 *     (J     0  ),
 *
 * where D is a diagonal matrix.  This matrix may be block reduced to
 * the form.
 *
 *     (D    J^T          )
 *     (0    -J D^{-1} J^T)
 *
 * Neither the matrix nor its inverse are stored explicitly. Rather, a
 * factorization of -J D^{-1} J^T, and the sufficient information to
 * backsolve using this factorization are stored.
 */
struct ReNewtonSystem {
    Int4 alphsize;              /*< the size of the alphabet */
    Int4 constrain_rel_entropy; /*< if true, use the relative entropy
                                    constraint for this optimization
                                    problem */
    Nlm_FloatHiPtr PNTR W;      /*< A lower-triangular matrix
                                    representing a factorization of
                                    the (2,2) block, -J D^{-1} J^T, of
                                    the condensed linear system */
    Nlm_FloatHiPtr Dinv;        /*< The diagonal elements of the
                                    inverse of the necessarily
                                    diagonal (1,1) block of the linear
                                    system */
    Nlm_FloatHiPtr grad_re;     /*< the gradient of the
                                    relative-entropy constraint, if
                                    this constraint is used. */
};
typedef struct ReNewtonSystem ReNewtonSystem;
typedef ReNewtonSystem PNTR ReNewtonSystemPtr;


/**
 * Create a new uninitialized ReNewtonSystem; the fields are
 * initialized by the FactorReNewtonSystem procedure.
 * ReNewtonSystemNew and FactorReNewtonSystem are called from only the
 * newt procedure.
 *
 * @param alphsize    the size of the alphabet for this optimization
 *                     problem.
 */
static ReNewtonSystemPtr
ReNewtonSystemNew(Int4 alphsize)
{
    ReNewtonSystemPtr newton_system;  /* the new ReNewtonSystem */

    newton_system = (ReNewtonSystem *) Nlm_MemNew(sizeof(ReNewtonSystem));

    newton_system->alphsize              = alphsize;
    newton_system->constrain_rel_entropy = 1;
    newton_system->W                     = Nlm_LtriangMatrixNew(2 * alphsize);

    newton_system->Dinv =
        (Nlm_FloatHiPtr) Nlm_MemNew(alphsize * alphsize * sizeof(Nlm_FloatHi));
    newton_system->grad_re =
        (Nlm_FloatHiPtr) Nlm_MemNew(alphsize * alphsize * sizeof(Nlm_FloatHi));

    return newton_system;
}


/**
 * Free the memory associated with a ReNewtonSystem.
 *
 * @param newton_system      on entry *newton_system points to the
 *                           system to be freed.  On exit, *newton_system
 *                           is set to NULL.
 */
static void
ReNewtonSystemFree(ReNewtonSystemPtr PNTR newton_system)
{
    (*newton_system)->W       = Nlm_DenseMatrixFree((*newton_system)->W);
    (*newton_system)->Dinv    =
        (Nlm_FloatHiPtr) Nlm_MemFree((*newton_system)->Dinv);
    (*newton_system)->grad_re = 
        (Nlm_FloatHiPtr) Nlm_MemFree((*newton_system)->grad_re);

    *newton_system = (ReNewtonSystemPtr) Nlm_MemFree(*newton_system);
}


/**
 * Factor the linear system to be solved in this iteration of Newton's
 * method.
 *
 * @param newton_system      holds the factorization
 * @param x                  the primal variables, representing the
 *                           target frequencies.
 * @param z                  the dual variables (Lagrange multipliers)
 * @param grads              gradient of the objective function and
 *                           (if used) the relative entropy constraint
 * @param constrain_rel_entropy    if true, then the relative entropy
 *                                 constraint is used for this optimization
 *                                 problem.
 */
static void
FactorReNewtonSystem(ReNewtonSystemPtr newton_system,
                     const Nlm_FloatHi PNTR x,
                     const Nlm_FloatHi PNTR z,
                     Nlm_FloatHiPtr PNTR grads,
                     Int4 constrain_rel_entropy)
{
    Int4 i;                     /* iteration index */
    Int4 n;                     /* the length of x */
    Int4 m;                     /* the length of z */

    /* Pointers to fields in newton_systems; the names of the local
     * variables match the names of the fields. */
    Nlm_FloatHiPtr PNTR W   = newton_system->W;
    Int4 alphsize           = newton_system->alphsize;
    Nlm_FloatHiPtr Dinv     = newton_system->Dinv;
    Nlm_FloatHiPtr grad_re  = newton_system->grad_re;

    n = alphsize * alphsize;
    m = constrain_rel_entropy ? 2 * alphsize : 2 * alphsize - 1;

    newton_system->constrain_rel_entropy = constrain_rel_entropy;

    /* The original system has the form
     *     
     *     (D     J^T)
     *     (J     0  ).
     *
     * We block reduce the system to 
     *
     *     (D    J^T          )
     *     (0    -J D^{-1} J^T).
     *
     * First we find the inverse of the diagonal matrix D. */
    
     if(constrain_rel_entropy) {
        Nlm_FloatHi eta;        /* dual variable for the relative
                                   entropy constraint */
        eta = z[m - 1];
        for(i = 0; i < n; i++) {
            Dinv[i] = x[i] / (1 - eta);
        }
    } else {
        Nlm_MemCpy(Dinv, x, n * sizeof(Nlm_FloatHi));
    }

    /* Then we compute J D^{-1} J^T; First fill in the part that corresponds
     * to the linear constraints */
    ScaledSymmetricProductA(W, Dinv, alphsize);

    if(constrain_rel_entropy) {
        Nlm_FloatHiPtr work;      /* a vector for intermediate computations */

        /* Save the gradient of the relative entropy constraint. */
        Nlm_MemCpy(grad_re, grads[1], n * sizeof(Nlm_FloatHi));

        /* Fill in the part of J D^{-1} J^T that corresponds to the relative
         * entropy constraint. */
        work = (Nlm_FloatHiPtr) Nlm_MemNew(n * sizeof(Nlm_FloatHi));

        W[m - 1][m - 1] = 0.0;
        for(i = 0; i < n; i++) {
            work[i] = Dinv[i] * grad_re[i];

            W[m - 1][m - 1] += grad_re[i] * work[i];
        }
        MultiplyByA(0.0, &W[m - 1][0], alphsize, 1.0, work);

        work = (Nlm_FloatHiPtr) Nlm_MemFree(work);
    }
    /* Factor J D^{-1} J^T and save the result in W. */
    Nlm_FactorLtriangPosDef(W, m);
}


/**
 * Solve the linear system for this iteration of Newton's method, using
 * the matrix factored by the FactorReNewtonSystem routine.
 *
 * @param x               on entry, the dual residuals; on exit, the
 *                        step in the primal variables.
 * @param z               on entry, the primal residuals; on exit, the
 *                        step in the dual variables.
 * @param newton_system   the factored matrix for the Newton system.
 */
static void
SolveReNewtonSystem(Nlm_FloatHiPtr x,
                    Nlm_FloatHiPtr z,
                    const ReNewtonSystem PNTR newton_system)
{
    Int4 i;                    /* iteration index */
    Int4 n;                    /* the size of x */
    Int4 mA;                   /* the number of linear constraints */
    Int4 m;                    /* the size of z */
    Nlm_FloatHiPtr work;       /* vector for intermediate calculations */

    /* Local variables that represent fields of newton_system */
    Nlm_FloatHiPtr PNTR W       = newton_system->W;
    Nlm_FloatHiPtr Dinv         = newton_system->Dinv;
    Nlm_FloatHiPtr grad_re      = newton_system->grad_re;
    Int4 alphsize               = newton_system->alphsize;
    Int4 constrain_rel_entropy  = newton_system->constrain_rel_entropy;

    n  = alphsize * alphsize;
    mA = 2 * alphsize - 1;
    m  = constrain_rel_entropy ? mA + 1 : mA;

    work = (Nlm_FloatHiPtr) Nlm_MemNew(n * sizeof(Nlm_FloatHi));

    /* Apply the same block reduction to the right-hand side as was
     * applied to the matrix:
     *
     *     rzhat = rz - J D^{-1} rx
     */
    for(i = 0; i < n; i++) {
        work[i] = x[i] * Dinv[i];
    }
    MultiplyByA(1.0, z, alphsize, -1.0, work);

    if(constrain_rel_entropy) {
        for(i = 0; i < n; i++) {
            z[m - 1] -= grad_re[i] * work[i];
        }
    }

    /* Solve for step in z, using the inverse of J D^{-1} J^T */
    Nlm_SolveLtriangPosDef(z, m, W);

    /* Backsolve for the step in x, using the newly-computed step in z.
     *
     *     x = D^{-1) (rx + J\T z)
     */
    if(constrain_rel_entropy) {
        for(i = 0; i < n; i++) {
            x[i] += grad_re[i] * z[m - 1];
        }
    }
    MultiplyByAtranspose(1.0, x, alphsize, 1.0, z);

    for(i = 0; i < n; i++) {
        x[i] *= Dinv[i];
    }
    work = (Nlm_FloatHiPtr) Nlm_MemFree(work);
}


/**
 * Evaluate the nonlinear functions and derivatives associated with
 * the problem of generating optimal target frequencies.
 *
 * @param values        values[0] is the value of the objective function
 *                      and values[1] is the relative entropy
 * @param grads         grads[0] is the gradient of the objective function
 *                      and grad[1] is the gradient of the relative entropy
 *                      constraint
 * @param x             the primal variables
 * @param q             target frequencies of the standard matrix
 * @param scores        scores as computed using the target frequencies
 *                      of the standard matrix, but the composition
 *                      of the sequences being compared (scaled to
 *                      lambda = 1).
 * @param constrain_rel_entropy     if true, the relative entropy constraint
 *                                  is used in this optimization problem
 */
static void
EvaluateReFunctions(Nlm_FloatHiPtr values,
                    Nlm_FloatHiPtr PNTR grads,
                    Int4 alphsize,
                    const Nlm_FloatHi PNTR x,
                    const Nlm_FloatHi PNTR q,
                    const Nlm_FloatHi PNTR scores,
                    Int4 constrain_rel_entropy)
{
    Int4 k;        /* iteration index over elements of x, q and scores */
    Nlm_FloatHi temp;   /* holds intermediate values in a computation */

    values[0] = 0.0; values[1] = 0.0;
    for(k = 0; k < alphsize * alphsize; k++) {
        temp = log(x[k] / q[k]);

        values[0]   += x[k] * temp;
        grads[0][k]  = temp + 1;

        if(constrain_rel_entropy) {
            temp += scores[k];

            values[1]   += x[k] * temp;
            grads[1][k]  = temp + 1;
        }
    }
}


/**
 * Compute a set of scores (scaled to lambda = 1) from a set of target
 * frequencies and a set of background frequencies.  The target
 * frequencies need not be consistent with the background
 * frequencies.
 *
 * @param scores        the resulting vector of scores, interpreted
 *                      as a matrix stored in row major order
 * @param alphsize      the size of the alphabet
 * @param target_freqs  target frequencies, interpreted as a matrix
 *                      stored in row-major order
 * @param row_freqs     background frequencies of one sequence
 * @param col_freqs     background frequencies of the other sequence
 */
static void
ComputeScoresFromProbs(Nlm_FloatHiPtr scores,
                       Int4 alphsize,
                       const Nlm_FloatHi PNTR target_freqs,
                       const Nlm_FloatHi PNTR row_freqs,
                       const Nlm_FloatHi PNTR col_freqs)
{
    Int4 i, j;     /* iteration indices over characters in the alphabet */
    Int4 k;        /* index into scores and target_freqs */

    for(i = 0; i < alphsize; i++) {
        for(j = 0; j < alphsize; j++) {
            k = i * alphsize + j;

            scores[k] = log(target_freqs[k] / (row_freqs[i] * col_freqs[j]));
        }
    }
}


/**
 * Find an optimal set of target frequencies for the purpose of
 * generating a compositionally adjusted score matrix.
 *
 * @param x           On exit, the optimal set of target frequencies,
 *                    interpreted as a two dimensional array in
 *                    row-major order.  x need not be initialized on
 *                    entry; any initial value will be ignored.
 * @param alphsize    the size of the alphabet for this optimization
 *                    problem.
 * @param q           a set of target frequencies from a standard
 *                    matrix
 * @param row_sums    the required row sums for the target frequencies;
 *                    the composition of one of the sequences being compared.
 * @param col_sums    the required column sums for the target frequencies;
 *                    the composition of the other sequence being compared.
 * @param constrain_rel_entropy   if true, constrain the relative
 *                                entropy of the optimal target
 *                                frequencies to equal
 *                                relative_entropy
 * @param relative_entropy  if constrain_rel_entropy is true, then this
 *                          is the required relative entropy for the
 *                          optimal target frequencies.  Otherwise,
 *                          this argument is ignored.
 * @param maxits    the maximum number of iterations permitted for the
 *                  optimization algorithm; a good value is 2000.
 * @param tol       the solution tolerance; the residuals of the optimization 
 *                  program must have Euclidean norm <= tol for the
 *                  algorithm to terminate.
 *
 * @returns         if an optimal set of target frequencies is
 *                  found, then the number of iterations used by the
 *                  optimization algorithm; otherwise maxits + 1.
 */
Int4
OptimizeTargetFrequencies(Nlm_FloatHiPtr x,
                          Int4 alphsize,
                          const Nlm_FloatHi PNTR q,
                          const Nlm_FloatHi PNTR row_sums,
                          const Nlm_FloatHi PNTR col_sums,
                          Int4 constrain_rel_entropy,
                          Nlm_FloatHi relative_entropy,
                          Nlm_FloatHi tol,
                          Int4 maxits)
{
    Int4 its;       /* number of iterations that have been performed */
    Int4 n;         /* number of target frequencies; the size of x */
    Int4 mA;        /* number of linear constraints */
    Int4 m;         /* total number of constraints */

    Nlm_FloatHi values[2];      /* values of the nonlinear functions
                                   at this iterate */
    Nlm_FloatHiPtr PNTR grads;  /* gradients of the nonlinear
                                   functions at this iterate */

    ReNewtonSystemPtr newton_system;    /* factored matrix of the
                                           linear system to be solved
                                           at this iteration */
    Nlm_FloatHiPtr z;           /* dual variables (Lagrange multipliers) */
    Nlm_FloatHiPtr resids_x;    /* dual residuals (gradient of Lagrangian) */
    Nlm_FloatHiPtr resids_z;    /* primal (constraint) residuals */
    Nlm_FloatHiPtr old_scores;  /* a scoring matrix, with lambda = 1,
                                   generated from q, row_sums and
                                   col_sums */
    Int4 converged;             /* true if Newton's method converged
                                   to a *minimizer* (strong
                                   second-order point) */
    n  = alphsize * alphsize;
    mA = 2 * alphsize - 1;
    m  = constrain_rel_entropy ? mA + 1 : mA;

    newton_system = ReNewtonSystemNew(alphsize);

    resids_x = (Nlm_FloatHiPtr) Nlm_MemNew(n * sizeof(Nlm_FloatHi));
    resids_z = (Nlm_FloatHiPtr) Nlm_MemNew((mA + 1) * sizeof(Nlm_FloatHi));
    /* z must be initialized to zero */
    z        = (Nlm_FloatHiPtr) Nlm_Calloc( mA + 1,   sizeof(Nlm_FloatHi));

    old_scores = (Nlm_FloatHiPtr) Nlm_MemNew(n * sizeof(Nlm_FloatHi));
    ComputeScoresFromProbs(old_scores, alphsize, q, row_sums, col_sums);

    grads = Nlm_DenseMatrixNew(2, n);

    /* Use q as the initial value for x */
    Nlm_MemCpy(x, q, n * sizeof(Nlm_FloatHi));
    its = 0;        /* Initialize the iteration count. Note that we may
                       converge in zero iterations if the initial x is 
                       optimal. */
    while(its <= maxits) {
        Nlm_FloatHi rnorm;      /* norm of the residuals for this iterate */

        /* Compute the residuals */
        EvaluateReFunctions(values, grads, alphsize, x, q, old_scores,
                            constrain_rel_entropy);
        CalculateResiduals(&rnorm, resids_x, alphsize, resids_z, values,
                           grads, row_sums, col_sums, x, z,
                           constrain_rel_entropy, relative_entropy);

        /* and check convergence */
        if(rnorm < tol) {
            /* We converged at the current iterate */
            break;
        } else {
            /* we did not converge, so increment the iteration counter
               and start a new iteration */
            if(++its <= maxits) {
                /* We have not exceeded the maximum number of iterations;
                   take a Newton step. */
                Nlm_FloatHi alpha;  /* a positive number used to scale the
                                       Newton step. */

                FactorReNewtonSystem(newton_system, x, z, grads,
                                     constrain_rel_entropy);
                SolveReNewtonSystem(resids_x, resids_z, newton_system);

                /* Calculate a value of alpha that ensure that x is
                   positive */
                alpha = Nlm_StepBound(x, n, resids_x, 1.0 / .95);
                alpha *= 0.95;

                Nlm_AddVectors(x, n, alpha, resids_x);
                Nlm_AddVectors(z, m, alpha, resids_z);
            }
        }
    }

    converged = 0;
    if( its <= maxits ) {
        /* Newton's iteration converged */
        if( !constrain_rel_entropy || z[m - 1] < 1 ) {
            /* and the final iterate is a minimizer */
            converged = 1;
        }
    }
            
    grads      = Nlm_DenseMatrixFree(grads);
    old_scores = (Nlm_FloatHiPtr) Nlm_MemFree(old_scores);
    z          = (Nlm_FloatHiPtr) Nlm_MemFree(z);
    resids_z   = (Nlm_FloatHiPtr) Nlm_MemFree(resids_z);
    resids_x   = (Nlm_FloatHiPtr) Nlm_MemFree(resids_x);

    ReNewtonSystemFree(&newton_system);

    return converged ? its : maxits + 1;
}
