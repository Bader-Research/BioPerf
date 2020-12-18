/*
  copyright (c) 1999 William R. Pearson
*/

/* $Name: fa34t25d2 $ - $Id: mm_file.h,v 1.20 2004/11/10 21:06:04 wrp Exp $ */

/*
  mm_file.h - defines m_file_str for mmap()ed files 
*/

#ifndef USE_FSEEKO
#define FSEEK fseek
#define FTELL ftell
typedef long fseek_t;
#else
#define FSEEK fseeko
#define FTELL ftello
typedef off_t fseek_t;
#endif
#define FSEEK_T_DEF

#ifdef HAS_INTTYPES
#include <inttypes.h>
#else
typedef long int64_t;
typedef unsigned long uint64_t;
#endif
#ifdef BIG_LIB64
typedef int64_t MM_OFF;
#else
typedef long MM_OFF;
#endif

#ifdef MYSQL_DB
#include <mysql.h>
#endif
#ifdef PGSQL_DB
#include <libpq-fe.h>
#endif

struct lmf_str {
  FILE *libf;		/* sequence file being read */
  FILE *hfile;		/* BLAST2.0 description file */
  char lb_name[120];	/* file name */
  int lb_type;		/* library type */
  int *sascii;		/* ascii -> sq mapping */

  /* used by flat files */
  char *lline;		/* last line read */
  unsigned char *cpsave;	/* position in line for lgetlib() */
  fseek_t lpos;			/* position in file */

  /* Genbank Flat files */
  int lfflag;		/* flag for CRLF in EMBL CDROM files */

  /* stuff for GCG format files (5,6) */
  int gcg_binary;	/* flag for binary gcg format */
  long gcg_len;		/* length of GCG sequence */

  /* for ncbl2 */
  int bl_lib_pos;

  /* text after filename */
  char opt_text[MAX_FN];

  /* blast formatdb version */
  int bl_format_ver;

  /* used when memory mapping */
  int mm_flg;		/* mmap worked */
  int mmap_fd;		/* mmap_fd */
  char *mmap_base;	/* base */
  char *mmap_addr;	/* current pos */
  long st_size;		/* file size */

  MM_OFF *d_pos_arr;	/* pointer to desc. offsets */
  MM_OFF *s_pos_arr;	/* pointer to seq. offsets */
  MM_OFF *a_pos_arr;	/* pointer to aux offsets */

  /* currently available only for memory mapped files */
  int max_cnt;		/* # database entries */
  int64_t tot_len;		/* total residue length */
  long max_len;		/* maximum sequence lengh */
  int lib_aa;		/* 0 = DNA, 1 = prot */

  char *sql_db, *sql_query, *sql_getdesc, *sql_getseq;
  int sql_reopen;
  char **sql_uid_arr;	/* indexed by lpos */
  /* used to get sequence data */
  char *sql_seqp;

#ifdef MYSQL_DB
  /* used to open the database */
  MYSQL *mysql_conn;
  MYSQL_RES *mysql_res;
  MYSQL_ROW mysql_row;
#endif

#ifdef PGSQL_DB
  /* used to open the database */
  PGconn *pgsql_conn;
  PGresult *pgsql_res;
#endif

  int (*getlib)();
  void (*ranlib)();
};

