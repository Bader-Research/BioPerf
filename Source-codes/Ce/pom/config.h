/* frequently used variables */
#ifndef H_CONFIG
#define H_CONFIG

#define DBDIR	"/projects/ard/db1/shindyal/rcsb/data/pdb/db/pom"

#ifndef DBPATH
#define DBPATH	"/projects/ard/db1/shindyal/rcsb/data/pdb/db/pom/pdbplus"
#endif

// the next section is somewhat obsolete and should be removed
#define MIRRORLOG	"/users/science/pdb/mirror/sdsc.distr/logs/recent.log"
#define OBSDIR		"/users/science/pdb/obsolete"
#define CURRDIR		"/users/science/pdb/data"
#define OBSDIRREL	"/obsolete"
#define CURRDIRREL	"/data"
#define WWWROOT		""
#define URL	"http://www.rcsb.org/"

// PROPERTY GROUPS

#define ASSIGN_ENV		       0
#define ASSIGN_KS		       0
#define MOOSE_PATTERNS_GROUP           0
#define MOOSE_REPORT_GROUP             1
#define FDS_GROUP                      0
#define OBS_COL_GROUP                  1


// PROPERTY OBJECT SELECTION TABLE (1-selected, 0-not selected)

#define P_CODE3_MON                    1
#define P_CODE1_MON                    1
#define P_N_BOND_MON                   1
#define P_BOND_MON                     1
#define P_N_ATOM_MON                   1
#define P_ATOM_MON                     1
#define P_PREV_MON                     1
#define P_NEXT_MON                     1
#define P_TYPE_MON                     1

#define P_ID_COM                       1
#define P_FILE_COM                     1
#define P_STATUS_COM                   1
#define P_TITLE_COM                    1
#define P_COMPND_COM                   1
#define P_SOURCE_COM                   1
#define P_DATE_TEX_COM                 1
#define P_DATE_INT_COM                 1
#define P_HEADER_COM                   1
#define P_AUTH_COM                     1
#define P_JRNL_COM                     1
#define P_EXPDTA_COM                   1
#define P_EXPDTA_TXT_COM               1
#define P_EC_COM              MOOSE_REPORT_GROUP 
#define P_RES_COM                      1
#define P_N_ENP_COM                    1
#define P_I_ENP_COM                    1
#define P_N_ENC_COM                    1
#define P_I_ENC_COM                    1
#define P_OBS_COM                      1
#define P_OBS_IDS_COM                  1
#define P_OBS_DAT_COM                  1
#define P_RELDAT_COM                   1
#define P_SPR_COM                      1
#define P_SPR_IDS_COM                  1
#define P_SPR_DAT_COM                  1
#define P_CURRENT_COM                  1
#define P_SPCGRP_COM                   1
#define P_UNITCELL_COM                 1
#define P_ZVAL_COM                     1
#define P_RVAL_COM                     1
#define P_SSBOND_COM                   1
#define P_SITE_COM                     1
#define P_NDBMAP_COM                   1

#define P_CHI1_COM                     0
#define P_FDS_COM                 FDS_GROUP 

#define P_NAME_ENC                     1
#define P_I_COM_ENC                    1
#define P_I_ENP_ENC                    1
#define P_N_SE_ENC                     1
#define P_SE_ENC                       1
#define P_SEN_PDB_ENC                  1
#define P_N_XYZ_ENC                    1
#define P_XYZ_ENC                      1
#define P_BFAC_ENC                     1
#define P_SE_XYZ_ENC                   1
#define P_XYZ_SE_ENC                   1

#define P_NAME_ENP                     1
#define P_I_COM_ENP                    1
#define P_I_ENC_ENP                    1
#define P_TYPE_ENP               MOOSE_REPORT_GROUP
#define P_MW_ENP                 MOOSE_REPORT_GROUP 
#define P_N_SE_ENP                     1
#define P_ALPHA_C_ENP            MOOSE_REPORT_GROUP
#define P_BETA_C_ENP             MOOSE_REPORT_GROUP
#define P_ALPHA_N_ENP            MOOSE_REPORT_GROUP 
#define P_BETA_N_ENP             MOOSE_REPORT_GROUP
#define P_SS_SEG_ENP             MOOSE_REPORT_GROUP 
#define P_SEQ_ENP                      1
#define P_SEQ_FLT_ENP                  1
#define P_ISEQ_ENP                     1
#define P_ISEQ_FLT_ENP                 1
#define P_C_A_ENP                      1
#define P_K_S_ENP                      ASSIGN_KS
#define P_EXP_ENP                      ASSIGN_ENV
#define P_EXP_FLT_ENP                  ASSIGN_ENV
#define P_POL_ENP                      ASSIGN_ENV
#define P_POL_FLT_ENP                  ASSIGN_ENV
#define P_BFAC_FLT_ENP                 ASSIGN_ENV
#define P_SE_TYPE_ENP                  1
#define P_FDS_ENP                 FDS_GROUP
#define P_CHI1_ENP                     0

#define P_EXP_1B_ENP              MOOSE_PATTERNS_GROUP        
#define P_POL_1B_ENP              MOOSE_PATTERNS_GROUP        
#define P_BFAC_1B_ENP             MOOSE_PATTERNS_GROUP     
#define P_SEXP_1B_ENP             MOOSE_PATTERNS_GROUP     
#define P_SPOL_1B_ENP             MOOSE_PATTERNS_GROUP     
#define P_SHYD_1B_ENP             MOOSE_PATTERNS_GROUP     
#define P_SVOL_1B_ENP             MOOSE_PATTERNS_GROUP     
#define P_SISO_1B_ENP             MOOSE_PATTERNS_GROUP

#define P_EXP_IMM                 MOOSE_PATTERNS_GROUP
#define P_POL_IMM                 MOOSE_PATTERNS_GROUP
#define P_BFAC_IMM                MOOSE_PATTERNS_GROUP
#define P_SEXP_IMM                MOOSE_PATTERNS_GROUP
#define P_SPOL_IMM                MOOSE_PATTERNS_GROUP
#define P_SHYD_IMM                MOOSE_PATTERNS_GROUP
#define P_SVOL_IMM                MOOSE_PATTERNS_GROUP
#define P_SISO_IMM                MOOSE_PATTERNS_GROUP

#define P_EXP_IP6                 MOOSE_PATTERNS_GROUP
#define P_POL_IP6                 MOOSE_PATTERNS_GROUP
#define P_BFAC_IP6                MOOSE_PATTERNS_GROUP
#define P_SEXP_IP6                MOOSE_PATTERNS_GROUP
#define P_SPOL_IP6                MOOSE_PATTERNS_GROUP
#define P_SHYD_IP6                MOOSE_PATTERNS_GROUP
#define P_SVOL_IP6                MOOSE_PATTERNS_GROUP
#define P_SISO_IP6                MOOSE_PATTERNS_GROUP

#define P_DEPYEAR_COL             OBS_COL_GROUP
#define P_RELYEAR_COL             OBS_COL_GROUP
#define P_OBSYEAR_COL             OBS_COL_GROUP
#define P_SPRYEAR_COL             OBS_COL_GROUP
#define P_CHAINS_COL              OBS_COL_GROUP
#define P_RES_COL                 OBS_COL_GROUP
#define P_FDS_COL                 FDS_GROUP    
#define P_MISC_COL                OBS_COL_GROUP

#define P_EXP_PP6                 MOOSE_PATTERNS_GROUP
#define P_POL_PP6                 MOOSE_PATTERNS_GROUP
#define P_BFAC_PP6                MOOSE_PATTERNS_GROUP
#define P_EXP_PP8                 MOOSE_PATTERNS_GROUP
#define P_POL_PP8                 MOOSE_PATTERNS_GROUP
#define P_BFAC_PP8                MOOSE_PATTERNS_GROUP
#define P_KS_PP3                       1
#define P_PROP8_AVE               MOOSE_PATTERNS_GROUP
#define P_PROP8_STS               MOOSE_PATTERNS_GROUP

#define P_LOCK_PID_MISC                1

// END OF -- PROPERTY OBJECT SELECTION TABLE

#endif
