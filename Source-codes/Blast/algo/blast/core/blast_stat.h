/*  $Id: blast_stat.h,v 1.67 2005/04/27 19:49:17 dondosha Exp $
 * ===========================================================================
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
 * ===========================================================================
 *
 * Author:  Tom Madden
 *
 */

/** @file blast_stat.h
 * Definitions and prototypes used by blast_stat.c to calculate BLAST
 * statistics. @todo FIXME: needs doxygen comments
 */

#ifndef __BLAST_STAT__
#define __BLAST_STAT__

#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_message.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Defines for the matrix 'preferences' (as specified by S. Altschul).
 * Used in arrays such as blosum45_prefs in blast_stat.c
 */
#define BLAST_MATRIX_NOMINAL 0   /**< acceptable values, not recommended. */
#define BLAST_MATRIX_PREFERRED 1 /**< These value are preferred over others. */
#define BLAST_MATRIX_BEST 2      /**< This is the best value, only one per matrix. */


#define BLASTMAT_DIR "/usr/ncbi/blast/matrix" /**< Default location for blast databases. */

/**
  Structure to hold the Karlin-Altschul parameters.
*/
typedef struct Blast_KarlinBlk {
      double   Lambda; /**< Lambda value used in statistics */
      double   K; /**< K value used in statistics */
      double   logK; /**< natural log of K value used in statistics */
      double   H; /**< H value used in statistics */
      double   paramC;  /**< for use in seed. */
   } Blast_KarlinBlk;




/********************************************************************

   Structures relating to scoring or the BlastScoreBlk

********************************************************************/

#define BLAST_SCORE_MIN INT2_MIN   /**< minimum allowed score (for one letter comparison). */
#define BLAST_SCORE_MAX INT2_MAX   /**< maximum allowed score (for one letter comparison). */


/** Holds score frequencies used in calculation
of Karlin-Altschul parameters for an ungapped search.
*/
typedef struct Blast_ScoreFreq {
    Int4         score_min; /**< lowest allowed scores */
    Int4         score_max; /**< highest allowed scores */
    Int4         obs_min;   /**< lowest observed (actual) scores */
    Int4         obs_max;   /**< highest observed (actual) scores */
    double       score_avg; /**< average score, must be negative for local alignment. */
    double*      sprob0;    /**< arrays for frequency of given score */
    double*      sprob;     /**< arrays for frequency of given score, shifted down by score_min. */
} Blast_ScoreFreq;

/** Scoring matrix used in BLAST */
typedef struct SBlastScoreMatrix {
    int** data;         /**< actual scoring matrix data, stored in row-major
                          form */
    size_t ncols;       /**< number of columns */
    size_t nrows;       /**< number of rows */
} SBlastScoreMatrix;

/** Scoring matrix data used in PSI-BLAST */
typedef struct SPsiBlastScoreMatrix {
    SBlastScoreMatrix* pssm;    /**< position-specific score matrix */
    double** freq_ratios;       /**< PSSM's frequency ratios, dimensions are
                                  specified in pssm data above */
    Blast_KarlinBlk* kbp;       /**< Karlin-Altschul block associated with this
                                  PSSM */
} SPsiBlastScoreMatrix;

/** Allocates a new SPsiBlastScoreMatrix structure of dimensions ncols by
 * BLASTAA_SIZE.
 * @param ncols number of columns (i.e.: query length) [in]
 * @return NULL in case of memory allocation failure, else new
 * SPsiBlastScoreMatrix structure
 */
SPsiBlastScoreMatrix*
SPsiBlastScoreMatrixNew(size_t ncols);

/** Deallocates a SPsiBlastScoreMatrix structure.
 * @param matrix structure to deallocate [in]
 * @return NULL
 */
SPsiBlastScoreMatrix*
SPsiBlastScoreMatrixFree(SPsiBlastScoreMatrix* matrix);

/** Structure used for scoring calculations.
*/
typedef struct BlastScoreBlk {
   Boolean     protein_alphabet; /**< TRUE if alphabet_code is for a 
protein alphabet (e.g., ncbistdaa etc.), FALSE for nt. alphabets. */
   Uint1    alphabet_code; /**< NCBI alphabet code. */
   Int2     alphabet_size;  /**< size of alphabet. */
   Int2     alphabet_start;  /**< numerical value of 1st letter. */
   char*    name;           /**< name of scoring matrix. */
   ListNode*   comments;    /**< Comments about scoring matrix. */
   SBlastScoreMatrix* matrix;   /**< scoring matrix data */
   SPsiBlastScoreMatrix* psi_matrix;    /**< PSSM and associated data. If this
                                         is not NULL, then the BLAST search is
                                         position specific (i.e.: PSI-BLAST) */
   Int4  loscore;   /**< Min.  substitution scores */
   Int4  hiscore;   /**< Max. substitution scores */
   Int4  penalty;   /**< penalty for mismatch in blastn. */
   Int4  reward;    /**< reward for match in blastn. */
        double  scale_factor; /**< multiplier for all cutoff and dropoff scores */
   Boolean     read_in_matrix; /**< If TRUE, matrix is read in, otherwise
               produce one from penalty and reward above. @todo should this be
                an allowed way of specifying the matrix to use? */
   Blast_ScoreFreq** sfp;  /**< score frequencies for scoring matrix. */
   /* kbp & kbp_gap are ptrs that should be set to kbp_std, kbp_psi, etc. */
   Blast_KarlinBlk** kbp;  /**< Karlin-Altschul parameters. Actually just a placeholder. */
   Blast_KarlinBlk** kbp_gap; /**< K-A parameters for gapped alignments.  Actually just a placeholder. */
   /* Below are the Karlin-Altschul parameters for non-position based ('std')
   and position based ('psi') searches. */
   Blast_KarlinBlk **kbp_std,  /**< K-A parameters for ungapped alignments */
                    **kbp_psi,       /**< K-A parameters for position-based alignments. */
                    **kbp_gap_std,  /**< K-A parameters for std (not position-based) alignments */
                    **kbp_gap_psi;  /**< K-A parameters for psi alignments. */
   Blast_KarlinBlk*  kbp_ideal;  /**< Ideal values (for query with average database composition). */
   Int4 number_of_contexts;   /**< Used by sfp and kbp, how large are these*/
   Uint1*   ambiguous_res; /**< Array of ambiguous res. (e.g, 'X', 'N')*/
   Int2     ambig_size, /**< size of array above. FIXME: not needed here? */
         ambig_occupy;  /**< How many occupied? */
} BlastScoreBlk;

/** 
Stores the letter frequency of a sequence or database.
*/
typedef struct Blast_ResFreq {
    Uint1   alphabet_code;    /**< indicates alphabet. */
    double* prob;       /**< letter probs, (possible) non-zero offset. */
    double* prob0;            /**< probs, zero offset. */
} Blast_ResFreq;

/**
 *   Allocates and initializes BlastScoreBlk
 * @param alphabet either BLASTAA_SEQ_CODE or BLASTNA_SEQ_CODE [in]
 * @param number_of_contexts how many strands or sequences [in]
 * @return BlastScoreBlk*
*/
BlastScoreBlk* BlastScoreBlkNew (Uint1 alphabet, Int4 number_of_contexts);

/** Deallocates BlastScoreBlk as well as all associated structures.
 * @param sbp BlastScoreBlk to be deallocated [in]
 * @return NULL pointer.
 */
BlastScoreBlk* BlastScoreBlkFree (BlastScoreBlk* sbp);

/** Set the ambiguous residue (e.g, 'N', 'X') in the BlastScoreBlk*.
 *  Convert from ncbieaa to sbp->alphabet_code (i.e., ncbistdaa) first.
 *
 * @param sbp the object to be modified [in|out]
 * @param ambiguous_res the residue to be set on the BlastScoreBlk
 * @return zero on success, others on error
 */
Int2 BLAST_ScoreSetAmbigRes (BlastScoreBlk* sbp, char ambiguous_res);


/** Calculate and fill the ungapped Karlin-Altschul parameters in the
 * BlastScoreBlk structure.
 * @param program BLAST program type, needed to decide whether to substitute
 *                ideal values. [in]
 * @param sbp Scoring block to work with [in] [out]
 * @param query Buffer containing (concatenated) query sequence [in]
 * @param query_info Information about offsets of concatenated queries [in]
 * @return 0 if ungapped Karlin-Altschul parameters could be 
 *               calculated for at least one context; 1 otherwise.
 */
Int2
Blast_ScoreBlkKbpUngappedCalc(EBlastProgramType program, 
                              BlastScoreBlk* sbp, Uint1* query, 
                              BlastQueryInfo* query_info);

/** This function fills in the BlastScoreBlk structure.  
 * Tasks are:
 * -read in the matrix
 * -set maxscore
 * @param sbp Scoring block [in] [out]
 * @param matrix Full path to the matrix in the directory structure [in]
*/
Int2 Blast_ScoreBlkMatrixFill (BlastScoreBlk* sbp, char* matrix);
 
/** Callocs a Blast_KarlinBlk
 * @return pointer to the Blast_KarlinBlk
*/
Blast_KarlinBlk* Blast_KarlinBlkNew (void);

/** Copies contents of one Karlin block to another. Both must be allocated
 * before this function is called.
 * @param kbp_to Karlin block to copy values to [in] [out]
 * @param kbp_from Karlin block to copy values from [in]
 * @return 0 on success; -1 if either argument is NULL on input.
 */ 
Int2 Blast_KarlinBlkCopy(Blast_KarlinBlk* kbp_to, Blast_KarlinBlk* kbp_from);


/** Deallocates the KarlinBlk
 * @param kbp KarlinBlk to be deallocated [in]
 * @return NULL
*/
Blast_KarlinBlk* Blast_KarlinBlkFree(Blast_KarlinBlk* kbp);

/** Fills in lambda, H, and K values, as calculated by Stephen Altschul 
 *  in Methods in Enzy. (vol 266, page 474).
 * @param kbp object to be filled in [in|out]
 * @param gap_open cost of gap existence [in]
 * @param gap_extend cost to extend a gap one letter [in]
 * @param decline_align cost of declining to align a letter [in]
 * @param matrix_name name of the matrix to be used [in]
 * @param error_return filled in with error message if needed [out]
 * @return zero on success
 */
Int2 Blast_KarlinBlkGappedCalc (Blast_KarlinBlk* kbp, Int4 gap_open, 
        Int4 gap_extend, Int4 decline_align, const char* matrix_name, 
        Blast_Message** error_return);


/** Calculates the Karlin-Altschul parameters assuming standard residue
 * compositions for the query and subject sequences. It populates the kbp_ideal
 * field of its sbp argument. This is used if the query is translated and the
 * calculated (real) Karlin parameters are bad, as they're calculated for 
 * non-coding regions.
 * @param sbp ScoreBlk used to calculate "ideal" values. [in|out]
 * @return 0 on success, 1 on failure
*/
Int2 Blast_ScoreBlkKbpIdealCalc(BlastScoreBlk* sbp);

/** Attempts to fill KarlinBlk for given gap opening, extensions etc.
 *
 * @param kbp object to be filled in [in|out]
 * @param gap_open gap existence cost [in]
 * @param gap_extend gap extension cost [in]
 * @param decline_align cost to not align part of a sequence [in]
 * @param matrix_name name of the matrix used [in]
 * @return  -1 if matrix_name is NULL;
 *          1 if matrix not found
 *           2 if matrix found, but open, extend etc. values not supported.
*/
Int2 Blast_KarlinBlkGappedLoadFromTables(Blast_KarlinBlk* kbp, Int4 gap_open, Int4 gap_extend, Int4 decline_align, const char* matrix_name);

/** Prints a messages about the allowed matrices, BlastKarlinBlkGappedFill should return 1 before this is called. 
 * @param matrix the matrix to print a message about [in]
 * @return the message
 */
char* BLAST_PrintMatrixMessage(const char *matrix);

/** Prints a messages about the allowed open etc values for the given matrix, 
 * BlastKarlinBlkGappedFill should return 2 before this is called. 
 * @param matrix name of the matrix [in]
 * @param gap_open gap existence cost [in]
 * @param gap_extend cost to extend a gap by one [in]
 * @param decline_align cost of declining to align [in]
 * @return message
 */
char* BLAST_PrintAllowedValues(const char *matrix, Int4 gap_open, Int4 gap_extend, Int4 decline_align);

/** Calculates the parameter Lambda given an initial guess for its value */
double
Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess);

/** Calculates the Expect value based upon the search space and some Karlin-Altschul 
 * parameters.  It is "simple" as it does not use sum-statistics.
 * @param S the score of the alignment. [in]
 * @param kbp the Karlin-Altschul parameters. [in]
 * @param searchsp total search space to be used [in]
 * @return the expect value
 */
double BLAST_KarlinStoE_simple (Int4 S, Blast_KarlinBlk* kbp, Int8  searchsp);

/** Compute a divisor used to weight the evalue of a collection of
 * "nsegs" distinct alignments.  These divisors are used to compensate
 * for the effect of choosing the best among multiple collections of
 * alignments.  See
 *
 * Stephen F. Altschul. Evaluating the statitical significance of
 * multiple distinct local alignments. In Suhai, editior, Theoretical
 * and Computational Methods in Genome Research, pages 1-14. Plenum
 * Press, New York, 1997.
 *
 * The "decayrate" parameter of this routine is a value in the
 * interval (0,1). Typical values of decayrate are .1 and .5.
 * @param decayrate adjusts for (multiple) tests of number of HSP sum groups [in]
 * @param nsegs the number of HSPs in the sum group [in]
 * @return divisor used to compensate for multiple tests
 */
double BLAST_GapDecayDivisor(double decayrate, unsigned nsegs );

/** Calculate the cutoff score from the expected number of HSPs or vice versa.
 * @param S The calculated score [in] [out]
 * @param E The calculated e-value [in] [out]
 * @param kbp The Karlin-Altschul statistical parameters [in]
 * @param searchsp The effective search space [in]
 * @param dodecay Use gap decay feature? [in]
 * @param gap_decay_rate Gap decay rate to use, if dodecay is set [in]
 */
Int2 BLAST_Cutoffs (Int4 *S, double* E, Blast_KarlinBlk* kbp, 
                    Int8 searchsp, Boolean dodecay, double gap_decay_rate);

/** Calculates the e-value for alignments with "small" gaps (typically
 *  under fifty residues/basepairs) following ideas of Stephen Altschul's.
 * @param start_points the number of starting points permitted between
 *    adjacent alignments; max_overlap + max_gap + 1 [in]
 * @param num the number of distinct alignments in this collection [in]
 * @param xsum the sum of the scores of these alignments each individually
 *    normalized using an appropriate value of Lambda and logK [in]
 * @param query_length effective len of the query seq [in]
 * @param subject_length effective len of the subject seq [in]
 * @param searchsp_eff effective size of the search space [in]
 * @param weight_divisor a divisor used to weight the e-value
 *    when multiple collections of alignments are being considered by 
 *    the calling routine [in]
 * @return the expect value 
 */
double BLAST_SmallGapSumE (Int4 start_points, Int2 num,  double xsum,
                           Int4 query_length, Int4 subject_length,
                           Int8 searchsp_eff, double weight_divisor);

/** Calculates the e-value of a collection multiple distinct
 *   alignments with asymmetric gaps between the alignments. The gaps
 *   in one (protein) sequence are typically small (like in
 *   BLAST_SmallGapSumE) gap an the gaps in the other (translated DNA)
 *   sequence are possibly large (up to 4000 bp.)  This routine is used
 *   for linking HSPs representing exons in the DNA sequence that are
 *   separated by introns.
 * @param query_start_points the number of starting points in
 *     the query sequence permitted between adjacent alignments;
 *     max_overlap + max_gap + 1. [in]
 * @param subject_start_points  the number of starting points in
 *     the subject sequence permitted between adjacent alignments [in]
 * @param num number of distinct alignments in one collection [in]
 * @param xsum the sum of the scores of these alignments each individually
 *    normalized using an appropriate value of Lambda and logK [in]
 * @param query_length effective length of query [in]
 * @param subject_length effective length of subject [in]
 * @param searchsp_eff effective size of the search space [in]
 * @param weight_divisor  a divisor used to weight the e-value
 *    when multiple collections of alignments are being considered by
 *    the calling routine [in]
 * @return sum expect value.
 */
double BLAST_UnevenGapSumE (Int4 query_start_points, Int4 subject_start_points,
                            Int2 num, double xsum,
                            Int4 query_length, Int4 subject_length,
                            Int8 searchsp_eff, double weight_divisor);

/** Calculates the e-value if a collection of distinct alignments with
 *   arbitrarily large gaps between the alignments
 * @param num number of distinct alignments in the collection [in]
 * @param xsum the sum of the scores of these alignments each individually
 *     normalized using an appropriate value of Lambda and logK [in]
 * @param query_length effective length of query sequence [in]
 * @param subject_length effective length of subject sequence [in]
 * @param searchsp_eff effective size of the search space [in]
 * @param weight_divisor  a divisor used to weight the e-value
 *    when multiple collections of alignments are being considered by the
 *    calling routine [in]
 * @return sum expect value.
 */
double BLAST_LargeGapSumE (Int2 num,  double xsum,
                           Int4 query_length, Int4 subject_length,
                           Int8 searchsp_eff, double weight_divisor );

/** Extract the alpha and beta settings for this matrixName, and these
 *  gap open and gap extension costs
 * @param matrixName name of the matrix used [in]
 * @param alpha Karlin-Altschul parameter to be set [out]
 * @param beta Karlin-Altschul parameter to be set [out]
 * @param gapped TRUE if a gapped search [in]
 * @param gap_open existence cost of a gap [in]
 * @param gap_extend extension cost of a gap [in]
*/
void BLAST_GetAlphaBeta (const char* matrixName, double *alpha,
                    double *beta, Boolean gapped, Int4 gap_open, Int4 gap_extend);


/** Rescale the PSSM, using composition-based statistics, for use
 *  with RPS BLAST. This function produces a PSSM for a single RPS DB
 *  sequence (of size db_seq_length) and incorporates information from 
 *  the RPS blast query. Each individual database sequence must call this
 *  function to retrieve its own PSSM. The matrix is returned (and must
 *  be freed elsewhere). posMatrix is the portion of the complete 
 *  concatenated PSSM that is specific to this DB sequence 
 * @todo revise to use existing code
 * @param scalingFactor used to rescale Lambda [in]
 * @param rps_query_length length of query sequence [in]
 * @param rps_query_seq the query sequence [in]
 * @param db_seq_length Length of the database sequence [in]
 * @param posMatrix matrix (actual) values to be used [in]
 * @param matrix_name Name of the score matrix underlying the RPS search [in]
 * @return rescaled pssm 
 */
Int4 ** RPSRescalePssm(double scalingFactor, Int4 rps_query_length, 
                   const Uint1 * rps_query_seq, Int4 db_seq_length, 
                   Int4 **posMatrix, const char *matrix_name);


/** 
 * Computes the adjustment to the lengths of the query and database sequences
 * that is used to compensate for edge effects when computing evalues. 
 *
 * The length adjustment is an integer-valued approximation to the fixed
 * point of the function
 *
 *    f(ell) = beta + 
 *               (alpha/lambda) * (log K + log((m - ell)*(n - N ell)))
 *
 * where m is the query length n is the length of the database and N is the
 * number of sequences in the database. The values beta, alpha, lambda and
 * K are statistical, Karlin-Altschul parameters.
 * 
 * The value of the length adjustment computed by this routine, A, 
 * will always be an integer smaller than the fixed point of
 * f(ell). Usually, it will be the largest such integer.  However, the
 * computed length adjustment, A, will also be so small that 
 *
 *    K * (m - A) * (n - N * A) > MAX(m,n).
 *
 * Moreover, an iterative method is used to compute A, and under
 * unusual circumstances the iterative method may not converge. 
 *
 * @param K      the statistical parameter K [in]
 * @param logK   the natural logarithm of K [in]
 * @param alpha_d_lambda    the ratio of the statistical parameters 
 *                          alpha and lambda (for ungapped alignments, the
 *                          value 1/H should be used) [in]
 * @param beta              the statistical parameter beta (for ungapped
 *                          alignments, beta == 0) [in]
 * @param query_length      the length of the query sequence [in]
 * @param db_length         the length of the database [in]
 * @param db_num_seqs       the number of sequences in the database [in]
 * @param length_adjustment the computed value of the length adjustment [out]
 *
 * @return   0 if length_adjustment is known to be the largest integer less
 *           than the fixed point of f(ell); 1 otherwise.
 */
Int4
BLAST_ComputeLengthAdjustment(double K,
                              double logK,
                              double alpha_d_lambda,
                              double beta,
                              Int4 query_length,
                              Int8 db_length,
                              Int4 db_num_seqs,
                              Int4 * length_adjustment);


/** Allocates a new Blast_ResFreq structure and fills in the prob element
    based upon the contents of sbp.
 * @param sbp The BlastScoreBlk* used to init prob [in]
*/
Blast_ResFreq* Blast_ResFreqNew(const BlastScoreBlk* sbp);

/** Deallocates Blast_ResFreq and prob0 element.
 * @param rfp the Blast_ResFreq to be deallocated.
*/
Blast_ResFreq* Blast_ResFreqFree(Blast_ResFreq* rfp);


/** Calculates residues frequencies given a standard distribution.
 * @param sbp the BlastScoreBlk provides information on alphabet.
 * @param rfp the prob element on this Blast_ResFreq is used.
 * @return zero on success
*/
Int2 Blast_ResFreqStdComp(const BlastScoreBlk* sbp, Blast_ResFreq* rfp);

/** Creates a new structure to keep track of score frequencies for a scoring
 * system.
 * @param score_min Minimum score [in]
 * @param score_max Maximum score [in]
 * @return allocated and initialized pointer to Blast_ScoreFreq
 */
Blast_ScoreFreq*
Blast_ScoreFreqNew(Int4 score_min, Int4 score_max);

/** Deallocates the score frequencies structure 
 * @param sfp the structure to deallocate [in]
 * @return NULL
 */
Blast_ScoreFreq*
Blast_ScoreFreqFree(Blast_ScoreFreq* sfp);

/** Fills a buffer with the 'standard' alphabet 
 * (given by STD_AMINO_ACID_FREQS[index].ch).
 *
 * @param alphabet_code specifies alphabet [in]
 * @param residues buffer to be filled in [in|out]
 * @param residue_size size of "residues" buffer [in]
 * @return Number of residues in alphabet or negative returns upon error.
 */
Int2
Blast_GetStdAlphabet(Uint1 alphabet_code, Uint1* residues, 
                     Uint4 residue_size);

/** Computes the parameters lambda, H K for use in calculating the
 * statistical significance of high-scoring segments or subalignments (see
 * comment on blast_stat.c for more details).
 * @param kbp object containing Lambda, H, and K as well as scoring information [in|out]
 * @param sfp array of probabilities for all scores [in]
 * @return zero on success, 1 on error.
 */
Int2
Blast_KarlinBlkUngappedCalc(Blast_KarlinBlk* kbp, Blast_ScoreFreq* sfp);

/**  Given a sequence of 'length' amino acid residues, compute the
 *   probability of each residue and put that in the array resProb
 *   Excludes ambiguity characters.
 *
 * @param sequence the sequence to be computed upon [in]
 * @param length the length of the sequence [in]
 * @param resProb the object to be filled in [in|out]
 */
void
Blast_FillResidueProbability(const Uint1* sequence, Int4 length, double * resProb);

/** Fill in the matrix for blastn using the penaly and rewards
 * The query sequence alphabet is blastna, the subject sequence
 * is ncbi2na.  The alphabet blastna is defined in blast_stat.h
 * and the first four elements of blastna are identical to ncbi2na.
 * if sbp->matrix==NULL, it is allocated.
 * @param sbp the BlastScoreBlk on which reward, penalty, and matrix will be set
 [in|out]
 * @return zero on success.
*/
Int2 BlastScoreBlkNuclMatrixCreate(BlastScoreBlk* sbp);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_STAT__ */
