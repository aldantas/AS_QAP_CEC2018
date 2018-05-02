#ifndef _TSP_LS_H
#define _TSP_LS_H

enum ls_types { LS_NONE = 0,
                LS_TWO_OPT_FIRST,     /* 2-opt local search */
                LS_TWO_H_OPT_FIRST,   /* 2.5-opt local search */
                LS_THREE_OPT_FIRST,   /* 3-opt local search */
                LS_UNKNOWN };
/* Must be synced with the ls_types */
extern const char * const strings_ls_types  [];

#include "ls.h"
void two_opt_first( long int *tour );
void two_h_opt_first( long int *tour );
void three_opt_first( long int *tour );

#endif

