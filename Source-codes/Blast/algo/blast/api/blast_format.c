/* $Id: blast_format.c,v 1.88 2005/06/02 20:41:52 dondosha Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
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
* Author: Ilya Dondoshansky
*
*/

/** @file blast_format.c
 * Formatting of BLAST results (SeqAlign)
 */

static char const rcsid[] = "$Id: blast_format.c,v 1.88 2005/06/02 20:41:52 dondosha Exp $";

#include <algo/blast/api/blast_format.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/api/blast_returns.h>
#include <algo/blast/api/blast_seq.h>
#include <readdb.h>
#include <txalign.h>
#include <blfmtutl.h>
#include <xmlblast.h>

extern void PrintTabularOutputHeader PROTO((char* blast_database, 
               Bioseq* query_bsp, SeqLoc* query_slp, char* blast_program, 
               Int4 iteration, Boolean believe_query, FILE *outfp));
extern void BlastPrintTabulatedResults PROTO((SeqAlign* seqalign, 
               Bioseq* query_bsp, SeqLoc* query_slp, Int4 num_alignments, 
               char* blast_program, Boolean is_ungapped, 
               Boolean believe_query, Int4 q_shift, Int4 s_shift, FILE *fp, 
               Boolean print_query_info));

extern void MBXmlClose(MBXml* xmlp, ValNode* other_returns, 
                       Boolean ungapped);

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

Int2 BlastFormattingOptionsNew(EBlastProgramType program_number, 
                               EAlignView align_view, 
                               BlastFormattingOptions** format_options_ptr)
{
   BlastFormattingOptions* format_options; 

   /* For translated RPS blast, the query will be a
      nucleotide sequence and subject sequences will
      be protein; thus the formatting must behave as
      if a blastx search was performed */

   if (program_number == eBlastTypeRpsTblastn)
      program_number = eBlastTypeBlastx;

   if ((format_options = (BlastFormattingOptions*)
        calloc(1, sizeof(BlastFormattingOptions))) == NULL)
      return -1;
   format_options->believe_query = FALSE;

   if ((program_number == eBlastTypeBlastx || 
        program_number == eBlastTypeTblastx)
        && align_view > 0 && align_view < eAlignViewXml)
   {
          ErrPostEx(SEV_FATAL, 1, 0, "This alignment view option is not "
                    "supported for blastx or tblastx");
          sfree(format_options);
          return -1;
   }


   format_options->print_options = 0;
   format_options->align_options = 0;
   format_options->align_options += TXALIGN_COMPRESS;
   format_options->align_options += TXALIGN_END_NUM;
   
   if (align_view == eAlignViewPairwise)
   {
      format_options->align_options += TXALIGN_MATRIX_VAL;
      format_options->align_options += TXALIGN_SHOW_QS;
   }
   else
   {
       format_options->align_options += TXALIGN_MASTER;
       if (align_view == eAlignViewQueryAnchoredIdent || 
           align_view == eAlignViewFlatQueryAnchoredIdent)
            format_options->align_options += TXALIGN_MISMATCH;
       if (align_view == eAlignViewFlatQueryAnchoredIdent || 
           align_view == eAlignViewFlatQueryAnchoredNoIdent ||
           align_view == eAlignViewFlatQueryAnchoredBluntEnds)
            format_options->align_options += TXALIGN_FLAT_INS;
        if (align_view == eAlignViewQueryAnchoredBluntEnds || 
            align_view == eAlignViewFlatQueryAnchoredBluntEnds)
            format_options->align_options += TXALIGN_BLUNT_END;
   }

   if (program_number == eBlastTypeBlastx)
      format_options->align_options += TXALIGN_BLASTX_SPECIAL;
   
   format_options->align_view = align_view;

   *format_options_ptr = format_options;
   return 0;
}

Int2 
BlastPrintLogReport(const BlastFormattingInfo* format_info)
{
     if (!format_info || !format_info->outfp)
         return -1;

     /* Index of the next query to be processed is equal to the number of
        queries processed so far. */
     fprintf(format_info->outfp, 
             "Mega BLAST run finished, processed %ld queries\n", 
             (long) format_info->num_formatted);

     return 0;
}

#define BXML_INCLUDE_QUERY 0x1

/** Creates the header part of an XML report for a BLAST search.
 * @param program Program name [in]
 * @param database Database name [in]
 * @param query_loc Query Seq-loc [in]
 * @param flags Flag to indicate whether query sequence should be included in
 *              the output. [in]
 * @param search_params Search parameters [in]
 */
static BlastOutput* 
s_CreateBlastOutputHead(const char* program, const char* database, 
                        SeqLoc* query_loc, Int4 flags, 
                        const Blast_SearchParams* search_param)
{
    BlastOutput* boutp;
    Char buffer[1024];
    SeqPort* spp;
    Boolean is_aa = FALSE;
    Int4 i;
    char* program_to_use = NULL;
    
    if((boutp = BlastOutputNew()) == NULL)
        return FALSE;
    
    if (strcmp(program, "rpsblast") == 0)
        program_to_use = strdup("blastp");
    else if (strcmp(program, "rpstblastn") == 0)
        program_to_use = strdup("blastx");
    else
        program_to_use = strdup(program);

    /* For optimization BLOSUM62 may be loaded ones */
    if (query_loc) {
       SeqId* sip = SeqLocId(query_loc);
       Bioseq* bsp;
       SeqIdWrite(sip, buffer, PRINTID_FASTA_LONG, sizeof(buffer));
       boutp->query_ID = strdup(buffer);

       bsp = BioseqLockById(sip);
       if(!bsp ||
          (boutp->query_def = strdup(BioseqGetTitle(bsp))) == NULL) {
          boutp->query_def = strdup("No definition line found");
       }
       BioseqUnlock(bsp);

       boutp->query_len = SeqLocLen(query_loc);

       if(flags & BXML_INCLUDE_QUERY) {
          boutp->query_seq = MemNew(boutp->query_len+1);
          is_aa = (strcmp(program_to_use, "blastp") || 
                   strcmp(program_to_use, "tblastn"));
          spp = SeqPortNewByLoc(query_loc, 
                   (Uint1)((is_aa) ? Seq_code_ncbieaa : Seq_code_iupacna));
          
          for (i = 0; i < boutp->query_len; i++) {
             boutp->query_seq[i] = SeqPortGetResidue(spp);
          }
          spp = SeqPortFree(spp);
       } else {
          boutp->query_seq = NULL;    /* Do we need sequence here??? */
       }
    }
    /* Program name. Use the local version of the program. No need to copy it
       since it was locally allocated. */
    boutp->program = program_to_use;

    /* Database name */
    if (database)
        boutp->db = strdup(database);

    /* Version text */
    sprintf(buffer, "%s %s [%s]", program_to_use, BlastGetVersionNumber(), 
            BlastGetReleaseDate());
    boutp->version = strdup(buffer);

    /* Reference */
    boutp->reference = BlastGetReference(FALSE);

    /* Filling parameters */
    boutp->param = ParametersNew();    
    boutp->param->expect = search_param->expect;
    boutp->param->gap_open = search_param->gap_open;
    boutp->param->gap_extend = search_param->gap_extension;
    if (search_param->matrix)
        boutp->param->matrix = strdup(search_param->matrix);
    boutp->param->sc_match = search_param->match;
    boutp->param->sc_mismatch = search_param->mismatch;
    boutp->param->include = search_param->ethresh;
    if (search_param->filter_string)
        boutp->param->filter = strdup(search_param->filter_string);
    
    return boutp;
}

/** Finds ASN.1 type for printing appropriate parts of the XML output. */
#define MACRO_atp_find(atp,name)\
        if((atp = AsnTypeFind(amp, #name))==NULL){\
                ErrPostEx(SEV_ERROR,0,0,\
                        "Could not find type <%s>", #name);\
                return NULL; \
        }


/** Initializes the XML output structure and prints the XML output header. 
 * @param aip ASN.1 output stream [in]
 * @param program BLAST program name [in]
 * @param database BLAST database name [in]
 * @param query_loc Query Seq-loc [in]
 * @param flags Flag for query sequence inclusion [in]
 * @param search_params BLAST search parameters [in]
 * @return Initialized structure for printing XML output.
 */
static MBXml* 
s_MBXmlInit(AsnIo* aip, char* program, char* database, 
          SeqLoc* query_loc, Int4 flags, Blast_SearchParams* search_params)
{
    MBXml* xmlp;
    AsnModule* amp;
    DataVal av;
    AsnType* atp;
    Boolean retval = FALSE;

    AsnType*       BLASTOUTPUT;
    AsnType*       BLASTOUTPUT_program;
    AsnType*       BLASTOUTPUT_version;
    AsnType*       BLASTOUTPUT_reference;
    AsnType*       BLASTOUTPUT_db;
    AsnType*       BLASTOUTPUT_query_ID;
    AsnType*       BLASTOUTPUT_query_def;
    AsnType*       BLASTOUTPUT_query_len;
    AsnType*       BLASTOUTPUT_query_seq;
    AsnType*       BLASTOUTPUT_param;
    AsnType*       BLASTOUTPUT_iterations;
    AsnType*       BLASTOUTPUT_iterations_E;
    AsnType*       BLASTOUTPUT_mbstat;

    if (strcmp(program, "rpsblast") == 0)
       program = "blastp";
    else if (strcmp(program, "rpstblastn") == 0)
       program = "blastx";

    xmlp = (MBXml*) MemNew(sizeof(MBXml));
    
    xmlp->aip = aip;
    
    if (! bxmlobjAsnLoad()) {
        return NULL;
    }
    
    amp = AsnAllModPtr();

    MACRO_atp_find(BLASTOUTPUT,BlastOutput);
    MACRO_atp_find(BLASTOUTPUT_program,BlastOutput.program);
    MACRO_atp_find(BLASTOUTPUT_version,BlastOutput.version);
    MACRO_atp_find(BLASTOUTPUT_reference,BlastOutput.reference);
    MACRO_atp_find(BLASTOUTPUT_db,BlastOutput.db);
    MACRO_atp_find(BLASTOUTPUT_query_ID,BlastOutput.query-ID);
    MACRO_atp_find(BLASTOUTPUT_query_def,BlastOutput.query-def);
    MACRO_atp_find(BLASTOUTPUT_query_len,BlastOutput.query-len);
    MACRO_atp_find(BLASTOUTPUT_query_seq,BlastOutput.query-seq);
    MACRO_atp_find(BLASTOUTPUT_param,BlastOutput.param);
    MACRO_atp_find(BLASTOUTPUT_iterations,BlastOutput.iterations);
    MACRO_atp_find(BLASTOUTPUT_iterations_E,BlastOutput.iterations.E);
    MACRO_atp_find(BLASTOUTPUT_mbstat,BlastOutput.mbstat);

    /* Start of iterations structure */
    xmlp->atp = BLASTOUTPUT_iterations_E;

    /* Head of all BlastOutput structure */
    xmlp->BlastOutput = BLASTOUTPUT;
    
    /* Head of iterations strucure */
    xmlp->BlastOutput_iterations = BLASTOUTPUT_iterations;

    /* Head of the final statistics for Mega BLAST */
    xmlp->BlastOutput_mbstat = BLASTOUTPUT_mbstat;
    
    xmlp->boutp = s_CreateBlastOutputHead(program, database, query_loc, flags, search_params);
    
    atp = AsnLinkType(NULL, BLASTOUTPUT);   /* link local tree */
    
    if (atp == NULL) {
        return NULL;
    }
    
    if (! AsnOpenStruct(xmlp->aip, atp, (Pointer) xmlp->boutp)) {
        return NULL;
    }

    if (xmlp->boutp->program != NULL) {
        av.ptrvalue = xmlp->boutp -> program;
        retval = AsnWrite(xmlp->aip, BLASTOUTPUT_program,  &av);
    }
    
    if (xmlp->boutp->version != NULL) {
        av.ptrvalue = xmlp->boutp->version;
        retval = AsnWrite(xmlp->aip, BLASTOUTPUT_version,  &av);
    }
    
    if (xmlp->boutp->reference != NULL) {
        av.ptrvalue = xmlp->boutp->reference;
        retval = AsnWrite(xmlp->aip, BLASTOUTPUT_reference,  &av);
    }

    if (xmlp->boutp -> db != NULL) {
        av.ptrvalue = xmlp->boutp->db;
        retval = AsnWrite(xmlp->aip, BLASTOUTPUT_db,  &av);
    }

    if (xmlp->boutp -> query_ID != NULL) {
        av.ptrvalue = xmlp->boutp->query_ID;
        retval = AsnWrite(xmlp->aip, BLASTOUTPUT_query_ID,  &av);
    }

    if (xmlp->boutp->query_def != NULL) {
        av.ptrvalue = xmlp->boutp->query_def;
        retval = AsnWrite(xmlp->aip, BLASTOUTPUT_query_def,  &av);
    }

    av.intvalue = xmlp->boutp->query_len;
    retval = AsnWrite(xmlp->aip, BLASTOUTPUT_query_len,  &av);

    if (xmlp->boutp->query_seq != NULL) {
        av.ptrvalue = xmlp->boutp->query_seq;
        retval = AsnWrite(xmlp->aip, BLASTOUTPUT_query_seq,  &av);
    }

    if (xmlp->boutp->param != NULL) {
        if (!ParametersAsnWrite(xmlp->boutp->param, 
                                xmlp->aip, BLASTOUTPUT_param)) {
            return NULL;
        }
    }

    if(!AsnOpenStruct(xmlp->aip, BLASTOUTPUT_iterations, NULL))
        return NULL;
    
    AsnIoFlush(xmlp->aip);
    
    return xmlp;
}

/** Fills the Statistics ASN.1 generated structure for the XML report.
 * @param sum_returns Search summary return data [in]
 * @param ungapped Is this an ungapped search? [in]
 * @return Populated structure.
 */
static Statistics*
s_XMLBuildStatistics(Blast_SummaryReturn* sum_returns, Boolean ungapped)
{
   Statistics* stat;
   Blast_DatabaseStats* db_stats = sum_returns->db_stats;

   stat = StatisticsNew();

   stat->eff_space = db_stats->eff_searchsp;
   stat->hsp_len = db_stats->hsp_length;
   stat->db_num= db_stats->dbnum;
   stat->db_len = db_stats->dblength;

   if(ungapped) {
      if(sum_returns->ka_params != NULL) {
         stat->lambda = sum_returns->ka_params->Lambda;
         stat->kappa = sum_returns->ka_params->K;
         stat->entropy = sum_returns->ka_params->H;
      }
   } else {
      if(sum_returns->ka_params_gap != NULL) {
         stat->lambda = sum_returns->ka_params_gap->Lambda;
         stat->kappa = sum_returns->ka_params_gap->K;
         stat->entropy = sum_returns->ka_params_gap->H;
      }
   }

   return stat;
}

/** Fills the Iteration ASN.1 structure, for part of the BLAST XML report
 * corresponding to one query.
 * @param seqalign Seq-align list with results [in]
 * @param sum_returns Search summary data [in]
 * @param is_ooframe Was out-of-frame gapping used in this search? [in]
 * @param ungapped Was this an ungapped search? [in]
 * @param iter_num Index of this "iteration" (query). [in]
 * @param message Error or warning message [in]
 * @param query Query Bioseq [in]
 * @param mask_loc List of masking locations [in]
 * @return Populated structure.
 */
static Iteration* 
s_XMLBuildOneQueryIteration(SeqAlign* seqalign,
                            Blast_SummaryReturn* sum_returns,
                            Boolean is_ooframe, Boolean ungapped,
                            Int4 iter_num, char* message,
                            Bioseq* query, ValNode* mask_loc)
{
    Iteration* iterp = IterationNew();
    iterp->iter_num = iter_num;

    if (query) {
       char buffer[1024];
       SeqIdWrite(query->id, buffer, PRINTID_FASTA_LONG, sizeof(buffer));
       iterp->query_ID = strdup(buffer);

       if((iterp->query_def = strdup(BioseqGetTitle(query))) == NULL) {
          iterp->query_def = strdup("No definition line found");
       }

       iterp->query_len = query->length;
    }

    if(seqalign != NULL) {
       iterp->hits =
          BXMLSeqAlignToHits(seqalign, ungapped, is_ooframe, mask_loc);
    }

    iterp->stat = s_XMLBuildStatistics(sum_returns, ungapped);

    if (message)
        iterp->message = strdup(message);

    return iterp;
}

Int2 BlastFormattingInfoNew(EAlignView align_view, 
                            const SBlastOptions* search_options,
                            const char* program_name, const char* db_name, 
                            const char* outfile_name, 
                            BlastFormattingInfo* *info_ptr)
{
    BlastFormattingInfo* info;
    
    if (!outfile_name || !info_ptr || !search_options || 
        align_view >= eAlignViewMax)
        return -1;

    info = *info_ptr = 
        (BlastFormattingInfo*) calloc(1, sizeof(BlastFormattingInfo));

    BlastFormattingOptionsNew(search_options->program, align_view, 
                              &info->format_options);
    info->search_options = search_options;
    info->program_name = strdup(program_name);
    if (db_name)
        info->db_name = strdup(db_name);

   if (align_view != eAlignViewXml && align_view < eAlignViewAsnText) {
      FILE* outfp;
      if ((outfp = FileOpen(outfile_name, "w")) == NULL)
         return -1;
      info->outfp = outfp;
   } else {
      char write_mode[3];
      /* AsnIoOpen requires a non-const argument, so use a local variable. */
      char* filename_copy = strdup(outfile_name);
      if (align_view == eAlignViewXml)
          strcpy(write_mode, "wx");
      else if (align_view == eAlignViewAsnText)
          strcpy(write_mode, "w"); 
      else
          strcpy(write_mode, "wb"); 
      
      if ((info->aip = AsnIoOpen(filename_copy, write_mode)) == NULL)
         return -1;
      sfree(filename_copy);
   }
   return 0;
}

Int2 BlastFormattingInfoNewBasic(EAlignView align_view, 
                                 const BLAST_SummaryOptions* basic_options,
                                 SeqLoc* query_seqloc,
                                 const char* outfile_name, 
                                 SBlastOptions* *advanced_options,
                                 BlastFormattingInfo* *info_ptr)
{
    Blast_SummaryReturn* extra_returns = Blast_SummaryReturnNew();
    char* program_name = NULL;

    Blast_SearchOptionsFromSummaryOptions(basic_options, query_seqloc,
                                          extra_returns, advanced_options,
                                          &program_name);

    Blast_SummaryReturnFree(extra_returns);

    return BlastFormattingInfoNew(align_view, *advanced_options, program_name, 
                                  NULL, outfile_name, info_ptr);
    sfree(program_name);
}

void 
BlastFormattingInfoSetUpOptions(BlastFormattingInfo* info, Int4 num_descriptions,
                                Int4 num_alignments, Boolean html, 
                                Boolean is_megablast, Boolean show_gi,
                                Boolean believe_defline)
{
    ASSERT(info && info->format_options && info->search_options);

    info->format_options->number_of_descriptions = num_descriptions;
    info->format_options->number_of_alignments = num_alignments;
    info->format_options->html = html;
    if (html) {
        info->format_options->align_options += TXALIGN_HTML;
        info->format_options->print_options += TXALIGN_HTML;
    }
    info->format_options->is_megablast = is_megablast;
    if (show_gi) {
        info->format_options->align_options += TXALIGN_SHOW_GI;
        info->format_options->print_options += TXALIGN_SHOW_GI;
    }

    if (!info->search_options->score_options->gapped_calculation)
        info->format_options->print_options += TXALIGN_SHOW_NO_OF_SEGS;

    info->format_options->believe_query = believe_defline;
}

BlastFormattingInfo* BlastFormattingInfoFree(BlastFormattingInfo* info)
{
    if (info->outfp)
        FileClose(info->outfp);
    /* Free the XML or ASN.1 i/o structure. Note that info->aip is set for
       XML too, but it should not be freed twice. */
    if (info->xmlp) {
        Boolean ungapped = 
            !info->search_options->score_options->gapped_calculation;
        MBXmlClose(info->xmlp, NULL, ungapped);
    } else if (info->aip) {
        AsnIoClose(info->aip);
    }
    sfree(info->program_name);
    sfree(info->db_name);
    sfree(info->format_options);
    /* The remaining fields are not owned by the structure. */
    sfree(info);

    return NULL;
}

/** Returns the integer alignment type value used in the C toolkit 
 * formatting code, given the EBlastProgramType enumeration value.
 * @param prog The BLAST program enumeration value [in]
 * @param db_is_na Is database nucleotide or protein for this program? [out]
 * @return Integer alignment type for use in old formatting code.
 */
static Uint1 
s_GetOldAlignType(EBlastProgramType prog, Boolean* db_is_na)
{
    Uint1 align_type = 0;

    switch (prog) {
    case eBlastTypeBlastp: case eBlastTypeRpsBlast: case eBlastTypePhiBlastp:
    case eBlastTypePsiBlast:
        align_type = 2;
        *db_is_na = FALSE;
        break;
    case eBlastTypeBlastn: case eBlastTypePhiBlastn: 
        align_type = 1;
        *db_is_na = TRUE;
        break;
    case eBlastTypeBlastx: case eBlastTypeRpsTblastn:
        align_type = 3;
        *db_is_na = FALSE;
        break;
    case eBlastTypeTblastn:
        align_type = 4;
        *db_is_na = TRUE;
        break;
    case eBlastTypeTblastx:
        align_type = 5;
        *db_is_na = TRUE;
        break;        
    default:
        break;
    }
    return align_type;
}

Int2 BLAST_FormatResults(SeqAlign* head, Int4 num_queries, 
        SeqLoc* query_slp, SeqLoc* mask_loc_head, 
        BlastFormattingInfo* format_info,
        Blast_SummaryReturn* sum_returns)
{  
   SeqAlign* seqalign = head;
   SeqAlign* sap;
   SeqAlign* next_seqalign;
   SeqLoc* mask_loc;
   SeqLoc* next_mask_loc = NULL;
   SeqLoc* tmp_loc = NULL;
   Uint1 align_type;
   Boolean db_is_na;
   Int4 query_index;
   SeqLoc* slp;
   SeqLoc* mask_slp;
   SeqId* query_id;
   Bioseq* bsp = NULL;
   AsnIo* aip = NULL;
   BlastPruneSapStruct* prune;
   SeqAnnot* seqannot;
   MBXml* xmlp = NULL;
   FILE *outfp = NULL;
   BlastFormattingOptions* format_options;
   EAlignView align_view;

   ASSERT(format_info && format_info->format_options && 
          format_info->search_options && query_slp);


   format_options = format_info->format_options;
   align_view = format_options->align_view;

   if (align_view == eAlignViewXml) {
       xmlp = format_info->xmlp;
       if (!xmlp) {
           xmlp = format_info->xmlp = 
               s_MBXmlInit(format_info->aip, format_info->program_name, 
                           format_info->db_name, query_slp, 0, 
                           sum_returns->search_params);
       }
   } else if (align_view == eAlignViewAsnText || 
              align_view == eAlignViewAsnBinary)
       aip = format_info->aip; 
   else 
       outfp = format_info->outfp;

   align_type = 
       s_GetOldAlignType(format_info->search_options->program, &db_is_na);

   if (format_info->db_name)
      ReadDBBioseqFetchEnable ("blast", format_info->db_name, db_is_na, TRUE);

   query_index = 0;
   slp = query_slp;
   mask_loc = mask_loc_head;

   if (!seqalign && align_view < eAlignViewXml) {
      /* Acknowledge all queries and report the "no hits found" message. */
      for ( ; slp; slp = slp->next) {
         bsp = BioseqLockById(SeqLocId(slp));
         init_buff_ex(70);
         AcknowledgeBlastQuery(bsp, 70, outfp, 
            format_options->believe_query, format_options->html);
         free_buff();
         BioseqUnlock(bsp);
         fprintf(outfp, " ***** No hits found ******\n\n\n");
      }
   }
      
   while (seqalign) {
      /* Find which query the current SeqAlign is for */
      query_id = TxGetQueryIdFromSeqAlign(seqalign);
      for ( ; slp && query_index < num_queries; 
            ++query_index, slp = slp->next) {
         if (SeqIdComp(query_id, SeqLocId(slp)) == SIC_YES) {
            break;
         } else if (align_view < eAlignViewXml) {
            /* This query has no results; acknowledge the query and report a "no
               results" message. */
            bsp = BioseqLockById(SeqLocId(slp));
            init_buff_ex(70);
            AcknowledgeBlastQuery(bsp, 70, outfp, 
               format_options->believe_query, format_options->html);
            free_buff();
            BioseqUnlock(bsp);
            fprintf(outfp, " ***** No hits found ******\n\n\n");
         }
      }
      /* The following should never happen! We mustn't have a SeqAlign left
         without a corresponding query */
      if (!slp || query_index >= num_queries)
         break;

      slp = slp->next;

      /* Find the masking location for this query. Initialize next_mask_loc
	 to the current start of the chain, in case nothing for this query 
	 will be found. */
      next_mask_loc = mask_loc;
      for ( ; mask_loc; mask_loc = mask_loc->next) {
         mask_slp = (SeqLoc*) mask_loc->data.ptrvalue;
         if (SeqIdComp(query_id, SeqLocId(mask_slp)) == SIC_YES) {
            break;
         }
      }
      /* Unlink the masking location for this query and save the next one */
      if (mask_loc) {
         for (next_mask_loc = mask_loc; next_mask_loc->next; 
              next_mask_loc = next_mask_loc->next) {
            mask_slp = (SeqLoc*) next_mask_loc->next->data.ptrvalue;
            if (SeqIdComp(query_id, SeqLocId(mask_slp))
                != SIC_YES) {
               break;
            }
         }
         tmp_loc = next_mask_loc;
         next_mask_loc = next_mask_loc->next;
         tmp_loc->next = NULL;
      }

      /* On the next iteration we can start from the next query */
      ++query_index;

      /* Find where this query's SeqAlign chain ends */
      for (sap = seqalign; sap->next; sap = sap->next) {
         if (SeqIdComp(query_id, TxGetQueryIdFromSeqAlign(sap->next)) != 
             SIC_YES)
            break;
      }
      /* Unlink this query's SeqAlign chain and save the start of the 
         next chain*/
      next_seqalign = sap->next;
      sap->next = NULL;

      /* Retrieve this query's Bioseq */
      bsp = BioseqLockById(query_id);

      if (align_view < eAlignViewXml) {
         init_buff_ex(70);
         AcknowledgeBlastQuery(bsp, 70, outfp, 
            format_options->believe_query, format_options->html);
         free_buff();
      }
      if (align_view == eAlignViewTabular || 
          align_view == eAlignViewTabularWithComments) {
         if (align_view == eAlignViewTabularWithComments)
            PrintTabularOutputHeader(format_info->db_name, bsp, NULL, 
                                     format_info->program_name,
                                     0, format_options->believe_query, outfp);
         
         BlastPrintTabulatedResults(seqalign, bsp, NULL, 
            format_options->number_of_alignments, format_info->program_name, 
            format_options->ungapped, format_options->believe_query, 0, 0, 
            outfp, (Boolean)(align_view == eAlignViewTabularWithComments));
      } else if(align_view == eAlignViewXml) {
         Iteration* iterp;
         
         ASSERT(xmlp && xmlp->aip);
         /* The index of this "query iteration" is the query_index in the 
            current formatting round, plus the number of previously formatted
            queries. */
         iterp = 
             s_XMLBuildOneQueryIteration(seqalign, sum_returns, FALSE, 
                                         format_options->ungapped, 
                                         query_index+format_info->num_formatted,
                                         NULL, bsp, mask_loc);
         IterationAsnWrite(iterp, xmlp->aip, xmlp->atp);
         AsnIoFlush(xmlp->aip);
         IterationFree(iterp);
      } else {
         seqannot = SeqAnnotNew();
         seqannot->type = 2;
         AddAlignInfoToSeqAnnot(seqannot, align_type);
         seqannot->data = seqalign;
         if (aip) {
            SeqAnnotAsnWrite((SeqAnnot*) seqannot, aip, NULL);
            AsnIoReset(aip);
         } 
         if (outfp) {
             Int4** matrix = NULL;
            ObjMgrSetHold();
            init_buff_ex(85);
            PrintDefLinesFromSeqAlignEx2(seqalign, 80, outfp, 
               format_options->print_options, FIRST_PASS, NULL,
               format_options->number_of_descriptions, NULL, NULL);
            free_buff();
            
            /** @todo FIXME: note that by calling BlastPruneHitsFromSeqAlign
             * we're making a COPY of the seqalign to print it out! Clearly
             * this could use a better design */
            prune = BlastPruneHitsFromSeqAlign(seqalign, 
                       format_options->number_of_alignments, NULL);
            seqannot->data = prune->sap;

             /** @todo FIXME: matrix must be filled here if it is not default 
                (BLOSUM62). Matrix name can be retrieved from 
                sum_returns->search_params->matrix. */
            if(format_info->search_options->score_options->is_ooframe) {
               OOFShowBlastAlignment(prune->sap, mask_loc, outfp, 
                                     format_options->align_options, NULL);
            } else if (align_view != eAlignViewPairwise) {
               ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, 
                  format_options->align_options, matrix, mask_loc, NULL);
            } else {
               ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, 
                  format_options->align_options, matrix, mask_loc, 
                  FormatScoreFunc);
            }
            seqannot->data = seqalign;
            prune = BlastPruneSapStructDestruct(prune);
            ObjMgrClearHold();
         }
         /* Set data to NULL, because we do not free Seq-align here. */
         seqannot->data = NULL;
         seqannot = SeqAnnotFree(seqannot);
      }
      BioseqUnlock(bsp);
      bsp = NULL;
      seqalign = next_seqalign;
      /* Re-link the chain of Seq-aligns, to prevent memory leak. */
      sap->next = next_seqalign;
      /* Relink the mask locations so chain can be freed in the end.
       The 'tmp_loc' variable points to the location that was unlinked. */
      if (tmp_loc)
          tmp_loc->next = next_mask_loc;
      
      mask_loc = next_mask_loc;
      ObjMgrFreeCache(0);

   } /* End loop on seqaligns for different queries */

   if (format_info->db_name)
      ReadDBBioseqFetchDisable();

   /* Update the count of the formatted queries. */
   format_info->num_formatted += num_queries;

   return 0;
}

/** Creates a list of SeqLoc structures with data about PHI BLAST pattern 
 * occurrences, to be used as features on Query Seq-locs.
 * @param pattern_info Pattern information structure. [in]
 * @param query_seqloc Query SeqLoc, needed to retrieve Seq-id. [in]
 * @param seed_seqloc_ptr List of SeqLoc's with pattern data. [out]
 */
static Int2
s_PHIBlastCreateSeedSeqLoc(const SPHIQueryInfo* pattern_info, 
                           SeqLoc* query_seqloc, 
                           SeqLoc** seed_seqloc_ptr)
{
    Int4 index;
    for (index = 0; index < pattern_info->num_patterns; ++index) {
        const SPHIPatternInfo* this_occurrence = 
            &pattern_info->occurrences[index];
        SeqInt* si = SeqIntNew();
        si->id = SeqIdDup(SeqLocId(query_seqloc));
        si->from = this_occurrence->offset;
        si->to = this_occurrence->offset + this_occurrence->length - 1;
        ValNodeAddPointer(seed_seqloc_ptr, SEQLOC_INT, si);
    }
    return 0;
}

/** Produces part of the text report describing the PHI BLAST pattern 
 * occurrences.
 * @param sum_returns Search summary return data [in]
 * @param outfp Output file stream [in]
 */ 
static Int2
s_PHIBlastFormatPatternInfo(Blast_SummaryReturn* sum_returns, FILE* outfp)
{
    Int4 index;
    SPHIQueryInfo* pattern_info = sum_returns->pattern_info;
    double lenXprob = 
        sum_returns->db_stats->eff_dblength * pattern_info->probability;
    Int8 db_patterns;

    /* Get the pattern count in database from the diagnostics structure - it is
       equal to the number of "lookup" hits. */
    ASSERT(sum_returns->diagnostics && sum_returns->diagnostics->ungapped_stat);
    db_patterns = sum_returns->diagnostics->ungapped_stat->lookup_hits;

    fprintf(outfp, "\n%d occurrence(s) of pattern in query\n", 
            pattern_info->num_patterns);
    for (index = 0; index < pattern_info->num_patterns; ++index) {
        fprintf(outfp, "\n pattern %s\n at position %d of query sequence\n",
                sum_returns->search_params->pattern, 
                pattern_info->occurrences[index].offset + 1);
        fprintf(outfp, "effective database length=%.1e\n", 
                (double)sum_returns->db_stats->eff_dblength);
        fprintf(outfp, " pattern probability=%.1e\nlengthXprobability=%.1e\n", 
                pattern_info->probability, lenXprob);
        fprintf(outfp, 
                "\nNumber of occurrences of pattern in the database is %lld\n",
                db_patterns);
        fprintf(outfp, "WARNING: There may be more matching sequences with "
                "e-values below the threshold of %f\n", 
                sum_returns->search_params->expect);

    }
    return 0;
}

Int2 PHIBlastFormatResults(ValNode* phivnps, SeqLoc* query_slp,
                           const BlastFormattingInfo* format_info, 
                           Blast_SummaryReturn* sum_returns)
{  
   Boolean db_is_na;
   Bioseq* query_bsp = NULL;
   FILE *outfp = NULL;
   ValNode* pruneSeed = NULL;
   Uint1 featureOrder[FEATDEF_ANY];
   Uint1 groupOrder[FEATDEF_ANY];
   SeqLoc* seed_seqloc = NULL; /* SeqLoc containing pattern locations. */ 
   EBlastProgramType program;
   BlastFormattingOptions* format_options;

   if (!format_info || !format_info->outfp || !query_slp)
      return -1;

   format_options = format_info->format_options;
   program = format_info->search_options->program;

   ASSERT(Blast_ProgramIsPhiBlast(program));

   outfp = format_info->outfp;

   s_PHIBlastFormatPatternInfo(sum_returns, outfp);

   /* Old toolkit might have different values for program numbers, so 
      use old toolkit function to determine alignment type. */
   if (program == eBlastTypePhiBlastn)
       db_is_na = TRUE;
   else
       db_is_na = FALSE;

   if (format_info->db_name)
      ReadDBBioseqFetchEnable ("blast", format_info->db_name, db_is_na, TRUE);

   pruneSeed = 
       SeedPruneHitsFromSeedReturn(phivnps,
                                   format_options->number_of_descriptions);

   s_PHIBlastCreateSeedSeqLoc(sum_returns->pattern_info, query_slp, 
                              &seed_seqloc);

   PrintDefLinesExtra(pruneSeed, 80, outfp, format_options->print_options,
                      FIRST_PASS, NULL, seed_seqloc);

   if (format_options->number_of_alignments < 
       format_options->number_of_descriptions) {
       pruneSeed = 
           SeedPruneHitsFromSeedReturn(phivnps,
                                       format_options->number_of_alignments);
   }

   query_bsp = BioseqLockById(SeqLocId(query_slp));
   memset(featureOrder, 0, sizeof(featureOrder));
   memset(groupOrder, 0, sizeof(groupOrder));
   featureOrder[FEATDEF_REGION] = 1;
   groupOrder[FEATDEF_REGION] = 1;

   if (format_options->align_view != eAlignViewPairwise) {
       ShowTextAlignFromAnnotExtra(query_bsp, pruneSeed, seed_seqloc, 60, outfp,
                                   featureOrder, groupOrder,
                                   format_options->align_options, NULL, NULL, 
                                   NULL);
   } else {
       ShowTextAlignFromAnnotExtra(query_bsp, pruneSeed, seed_seqloc, 60, outfp,
                                   featureOrder, groupOrder,
                                   format_options->align_options, NULL, NULL, 
                                   FormatScoreFunc);
   }

   SeqLocSetFree(seed_seqloc);

   if (format_info->db_name)
      ReadDBBioseqFetchDisable();

   return 0;
}

ValNode* PHIBlastResultsFree(ValNode* phivnps)
{
    ValNode* vnp;

    for (vnp = phivnps; vnp; vnp = vnp->next)
        SeqAlignSetFree((SeqAlign*) vnp->data.ptrvalue);
    ValNodeFree(phivnps);
    return NULL;
}

Int2 Blast_PrintOutputFooter(const BlastFormattingInfo* format_info, 
                             const Blast_SummaryReturn* sum_returns) 
{
   FILE *outfp;
   char* params_buffer=NULL;
   BlastFormattingOptions* format_options;
   EBlastProgramType program;

   if (!format_info || !format_info->outfp)
      return -1;

   outfp = format_info->outfp;

   format_options = format_info->format_options;
   program = format_info->search_options->program;

   if (format_options->align_view >= eAlignViewXml)
      return 0;

   if(format_options->html) 
      fprintf(outfp, "<PRE>\n");
   init_buff_ex(85);

   if (format_info->db_name) {
       const Boolean kDbIsProt = 
           (program == eBlastTypeBlastp    ||
            program == eBlastTypeBlastx    ||
            program == eBlastTypePhiBlastp ||
            program == eBlastTypePsiBlast  ||
            program == eBlastTypeRpsBlast  ||
            program == eBlastTypeRpsTblastn);
       ReadDBFILE* rdfp = readdb_new(format_info->db_name, kDbIsProt);
       TxDfDbInfo* dbinfo_head,* dbinfo;
       dbinfo_head = Blast_GetDbInfo(rdfp);
       rdfp = readdb_destruct(rdfp);

       for (dbinfo = dbinfo_head; dbinfo; dbinfo = dbinfo->next) {
           PrintDbReport((TxDfDbInfo*) dbinfo, 70, outfp);
       }
       dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
   }

      
   if (sum_returns->ka_params)
   {
      BLAST_KAParameters* ka_params=sum_returns->ka_params;
      PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, 
                           outfp, FALSE);
   }

 
   if (sum_returns->ka_params_gap && sum_returns->search_params->gapped_search)
   {
      BLAST_KAParameters* ka_params=sum_returns->ka_params_gap;
      if (Blast_ProgramIsPhiBlast(program)) {
          PrintKAParametersExtra(ka_params->Lambda, ka_params->K, ka_params->H, 
                                 ka_params->C, 70, outfp, TRUE);
      } else {
          PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, 
                            outfp, TRUE);
      }
   }
   params_buffer = 
      Blast_GetParametersBuffer(program, sum_returns);
   PrintTildeSepLines(params_buffer, 70, outfp);
   sfree(params_buffer);

   if(format_options->html)
      fprintf(outfp, "</PRE>\n</BODY>\n</HTML>\n");

   free_buff();

   return 0;
}

Int2 
BLAST_PrintOutputHeader(const BlastFormattingInfo* format_info)
{
   Int2 status = 0;
   BlastFormattingOptions* format_options;

   if (!format_info)
       return -1;
   
   format_options = format_info->format_options;

   if (format_options->align_view < eAlignViewXml) {
      if (format_options->html) {
         fprintf(format_info->outfp, 
                 "<HTML>\n<TITLE>BLAST Search Results</TITLE>\n");
         fprintf(format_info->outfp, 
                 "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
                 "VLINK=\"#660099\" ALINK=\"#660099\">\n");
         fprintf(format_info->outfp, "<PRE>\n");
      }
      init_buff_ex(90);
      if (format_options->is_megablast) {
         BlastPrintVersionInfo("MEGABLAST", format_options->html, format_info->outfp);
         fprintf(format_info->outfp, "\n");
         MegaBlastPrintReference(format_options->html, 90, format_info->outfp);
      }
      else if (strcmp(format_info->program_name, "rpsblast") == 0 ||
               strcmp(format_info->program_name, "rpstblastn") == 0) {
         BlastPrintVersionInfo("RPS-BLAST", format_options->html, format_info->outfp);
      }
      else {
         BlastPrintVersionInfo(format_info->program_name, format_options->html, 
                               format_info->outfp);
         fprintf(format_info->outfp, "\n");
         BlastPrintReference(format_options->html, 90, format_info->outfp);
      }

      if(format_info->db_name) { 
          EBlastProgramType prog = format_info->search_options->program;
          const Boolean kDbIsProt = 
              (prog == eBlastTypeBlastp || prog == eBlastTypeBlastx ||
               prog == eBlastTypeRpsBlast || prog == eBlastTypeRpsTblastn ||
               prog == eBlastTypePhiBlastp);
         status = 
            PrintDbInformation(format_info->db_name, kDbIsProt, 70, 
                               format_info->outfp, format_options->html);
      }
      free_buff();
	
#ifdef OS_UNIX
      fprintf(format_info->outfp, "%s", "Searching\n\n");
#endif
   }

   return status;
}

#ifdef BUFFER_LENGTH
#undef BUFFER_LENGTH
#endif

/** Buffer length for printing sequence ids. */
#define BUFFER_LENGTH 256

void 
Blast_SeqIdGetDefLine(SeqId* sip, char** buffer_ptr, Boolean ncbi_gi, 
                      Boolean accession_only)
{
   char* seqid_buffer = NULL;
   Int4 gi;
   Boolean numeric_id_type = FALSE;

   *buffer_ptr = NULL;

   if (sip == NULL)
	return;

   /* Check for ad hoc ID's generated by formatdb if the user does not provide 
      any. */
   if (sip->choice != SEQID_GENERAL ||
       StringCmp(((Dbtag*)sip->data.ptrvalue)->db, "BL_ORD_ID"))  
   {
      if ((!accession_only && !ncbi_gi) || sip->choice == SEQID_LOCAL) {
         seqid_buffer = (char*) malloc(BUFFER_LENGTH + 1);
         SeqIdWrite(sip, seqid_buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
      } else if (accession_only) {
         seqid_buffer = (char*) malloc(BUFFER_LENGTH + 1);
         SeqIdWrite(SeqIdFindBestAccession(sip), seqid_buffer, 
                    PRINTID_TEXTID_ACC_VER, BUFFER_LENGTH);
      } else if (ncbi_gi) {
         numeric_id_type = 
            GetAccessionFromSeqId(SeqIdFindBest(sip, SEQID_GI), 
                                  &gi, &seqid_buffer);
      } else {
         numeric_id_type = 
	    GetAccessionFromSeqId(SeqIdFindBestAccession(sip), 
                                  &gi, &seqid_buffer);
      }
   }

   if (numeric_id_type && gi > 0) {
      seqid_buffer = (char*) malloc(16);
      sprintf(seqid_buffer, "%ld", (long) gi);
   }   
   if (!seqid_buffer) {
      /* If it's still NULL make a last ditch effort to get info. */
      char* title=NULL;
      Bioseq* bsp = BioseqLockById(sip);
      if (bsp)
         title = strdup(BioseqGetTitle(bsp));
      BioseqUnlock(bsp);
      
      if (title) /* Use first token as id. */
         seqid_buffer = StringTokMT(title, " \t\n\r", &title);  
   }
   *buffer_ptr = seqid_buffer;

}

Int2 
Blast_SummaryReturnsPostError(Blast_SummaryReturn* sum_return)
{
    Blast_Message* error = NULL;

    if (!sum_return)
        return -1;

    /* If there was no error, there is nothing to post. */
    if (!sum_return->error)
        return 0;

    error = sum_return->error;

    ErrPostEx(error->severity, error->code, error->subcode, error->message);

    return 0;
}
/* @} */

