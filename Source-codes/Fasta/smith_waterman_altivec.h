
int
smith_waterman_altivec_word(unsigned char *     query_sequence,
                            unsigned short *    query_profile_word,
                            int                 query_length,
                            unsigned char *     db_sequence,
                            int                 db_length,
                            unsigned short      bias,
                            unsigned short      gap_open,
                            unsigned short      gap_extend,
                            struct f_struct *   f_str);


int
smith_waterman_altivec_byte(unsigned char *     query_sequence,
                            unsigned char *     query_profile_byte,
                            int                 query_length,
                            unsigned char *     db_sequence,
                            int                 db_length,
                            unsigned char       bias,
                            unsigned char       gap_open,
                            unsigned char       gap_extend,
                            struct f_struct *   f_str);

