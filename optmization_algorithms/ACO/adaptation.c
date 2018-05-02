#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "adaptation.h"
#include "parse.h"
#include "ants.h"
#include "utilities.h"


extern long int end_n_ants;
extern double end_beta;
extern double end_rho;
extern double end_q0;

#define CLAMP(VALUE, MIN_VALUE, MAX_VALUE)                   \
    do {                                                     \
    if (MIN_VALUE < MAX_VALUE) {                             \
        if (VALUE < MIN_VALUE) VALUE = MIN_VALUE;            \
        else if (VALUE > MAX_VALUE) VALUE = MAX_VALUE;       \
    } else {                                                 \
        if (VALUE > MIN_VALUE) VALUE = MIN_VALUE;            \
        else if (VALUE < MAX_VALUE) VALUE = MAX_VALUE;       \
    }} while (0)

static double cum_n_ants;

void
adapt_parameters_init(void)
{
#define ADAPT_PARAMETER_INIT(PARAM)                      \
    do {                                                 \
        switch(opt_var_##PARAM) {                        \
        case VAR_NONE:                                   \
            PARAM = opt_##PARAM;                         \
            break;                                       \
                                                         \
        case VAR_DELTA:                                                        \
            PARAM = opt_##PARAM;                                               \
            break;                                                             \
                                                                               \
        case VAR_SWITCH:                                                       \
            PARAM = opt_##PARAM;                                               \
            break;                                                             \
                                                                               \
        default: abort();                                                      \
        }                                                                      \
        DEBUG(fprintf (stderr, #PARAM " = %g\n", (double)PARAM));             \
    } while(0)

    cum_n_ants = 0;
    ADAPT_PARAMETER_INIT(n_ants);
    ADAPT_PARAMETER_INIT(beta);
    ADAPT_PARAMETER_INIT(rho);
    ADAPT_PARAMETER_INIT(q0);
}

void
adapt_parameters_next_iteration(void)
{
#define ADAPT_PARAMETER_NEXT_REAL(PARAM) \
    do {                                 \
        switch(opt_var_##PARAM) {        \
                                         \
        case VAR_NONE: break;                        \
                                                     \
        case VAR_DELTA:                              \
            PARAM += (opt_##PARAM > end_##PARAM) \
                ? -delta_##PARAM : delta_##PARAM;                              \
            CLAMP (PARAM, opt_##PARAM, end_##PARAM);                           \
            break;                                                             \
                                                                               \
        case VAR_SWITCH:                                                       \
            if (iteration == switch_##PARAM)                                   \
                PARAM = end_##PARAM;                                           \
            break;                                                             \
                                                                               \
        default: abort();                                                      \
        }                                                                      \
        DEBUG(fprintf (stderr, #PARAM " = %g\n", PARAM));                     \
    }  while (0)

    
    switch(opt_var_n_ants) {
        
    case VAR_NONE: break;

    case VAR_DELTA:
    {
        /* delta_n_ants may not be an integer. In such a case, we
           accumulate in cum_n_ants until we reach something >= 1, and
           then we change n_ants.  */
        cum_n_ants += (opt_n_ants > end_n_ants) ? -delta_n_ants : delta_n_ants;
        int trunc_n_ants = trunc(cum_n_ants);
        if (abs(trunc_n_ants) >= 1) {
            n_ants += trunc_n_ants;
            cum_n_ants -= trunc_n_ants;
            CLAMP (n_ants, opt_n_ants, end_n_ants);
        }
        break;
    }
    case VAR_SWITCH: 
        if (iteration == switch_n_ants) 
            n_ants = end_n_ants;
        break;

    default: abort();
    }
    DEBUG(fprintf (stderr, "n_ants = %ld\n", n_ants));

    ADAPT_PARAMETER_NEXT_REAL(beta);
    ADAPT_PARAMETER_NEXT_REAL(rho);
    ADAPT_PARAMETER_NEXT_REAL(q0);
}
