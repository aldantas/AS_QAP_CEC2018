#ifndef _QAP_LS_H
#define _QAP_LS_H
#include "ls.h"
enum ls_types { LS_NONE = 0,
                LS_first_2_opt,
                LS_best_2_opt,
                LS_tabu_search_short,
                LS_tabu_search_long,
                LS_UNKNOWN };
/* Must be synced with the ls_types */
extern const char * const strings_ls_types  [];

void first_2_opt_symmetric ( long int *q );
void first_2_opt_asymmetric ( long int * q );
void best_2_opt_symmetric ( long int * q );
void best_2_opt_asymmetric ( long int * q );
void tabu_search( long int *s, int tabu_search_length);
#endif
