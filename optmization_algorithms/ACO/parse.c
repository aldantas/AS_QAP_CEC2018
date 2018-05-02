/* This file was originally generated with opag 0.6.4, but it has been manually
   modified afterwards.  */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdlib.h>


#include "InOut.h"
#include "utilities.h"
#include "ants.h"
#include "ls.h"
#include "parse.h"

struct options {

/* Set to 1 if option --tsplibfile (-i) has been specified.  */
unsigned int opt_tsplibfile : 1;

/* Set to 1 if option --optimum (-o) has been specified.  */
unsigned int opt_optimum : 1;

/* Set to 1 if option --ants (-m) has been specified.  */
unsigned int opt_ants : 1;

/* Set to 1 if option --elitistants (-c) has been specified.  */
unsigned int opt_elitistants : 1;

/* Set to 1 if option --localsearch (-l) has been specified.  */
unsigned int opt_localsearch : 1;

/* Set to 1 if option --dlb (-d) has been specified.  */
unsigned int opt_dlb : 1;

/* Set to 1 if option --as has been specified.  */
unsigned int opt_as : 1;

/* Set to 1 if option --eas has been specified.  */
unsigned int opt_eas : 1;

/* Set to 1 if option --ras has been specified.  */
unsigned int opt_ras : 1;

/* Set to 1 if option --mmas has been specified.  */
unsigned int opt_mmas : 1;

/* Set to 1 if option --bwas has been specified.  */
unsigned int opt_bwas : 1;

/* Set to 1 if option --acs has been specified.  */
unsigned int opt_acs : 1;

/* Argument to option --tsplibfile (-i).  */
const char *arg_tsplibfile;

/* Argument to option --optimum (-o).  */
const char *arg_optimum;

/* Argument to option --ants (-m).  */
const char *arg_ants;

/* Argument to option --elitistants (-c).  */
const char *arg_elitistants;

/* Argument to option --localsearch (-l).  */
const char *arg_localsearch;

/* Argument to option --dlb (-d).  */
const char *arg_dlb;
#define DEFINE_OPTIONS_PARAMETER(A, B, C, D, TYPE) DEFINE_PARAMETER(A,B,C,D,TYPE, 0, 0)
#define DEFINE_PARAMETER(VAR, B, C, D, TYPE, F, G)                    \
    /* Set to 1 if option VAR has been specified.  */                 \
    unsigned int opt_##VAR : 1;                                       \
    /* Argument to option VAR.  */                                    \
    const char *arg_##VAR;                                            \

#include "aco-parameters.def"
};

#define DEFINE_BOOLEAN_PARAMETER(PARAM)                                        \
    bool PARAM;
#define DEFINE_DOUBLE_PARAMETER(PARAM) \
    double PARAM;
#define DEFINE_INTEGER_PARAMETER(PARAM) \
    long int PARAM;
#define DEFINE_OPTIONS_PARAMETER(PARAM, B, C, D, TYPE) \
    enum  enum_##TYPE PARAM;
#define DEFINE_OPTION(E,S) S
#define DEFINE_OPTIONS(NAME, ...)                                              \
    const char * const strings_##NAME [] = { __VA_ARGS__, NULL };
#define DEFINE_PARAMETER(PARAM,B,C,D, TYPE, X, Y) TYPE(PARAM)
#include "aco-parameters.def"
long int opt_n_ants;


#ifndef STR_ERR_UNKNOWN_LONG_OPT
# define STR_ERR_UNKNOWN_LONG_OPT   "%s: unrecognized option `--%s'\n"
#endif

#ifndef STR_ERR_LONG_OPT_AMBIGUOUS
# define STR_ERR_LONG_OPT_AMBIGUOUS "%s: option `--%s' is ambiguous\n"
#endif

#ifndef STR_ERR_MISSING_ARG_LONG
# define STR_ERR_MISSING_ARG_LONG   "%s: option `--%s' requires an argument\n"
#endif

#ifndef STR_ERR_UNEXPEC_ARG_LONG
# define STR_ERR_UNEXPEC_ARG_LONG   "%s: option `--%s' doesn't allow an argument\n"
#endif

#ifndef STR_ERR_UNKNOWN_SHORT_OPT
# define STR_ERR_UNKNOWN_SHORT_OPT  "%s: unrecognized option `-%c'\n"
#endif

#ifndef STR_ERR_MISSING_ARG_SHORT
# define STR_ERR_MISSING_ARG_SHORT  "%s: option `-%c' requires an argument\n"
#endif


#define STR_HELP_TSPLIBFILE \
  "  -i, --tsplibfile     f inputfile (TSPLIB format necessary)\n"

#define STR_HELP_OPTIMUM \
  "  -o, --optimum        # stop if tour better or equal optimum is found\n"

#define STR_HELP_ANTS \
  "  -m, --ants           # number of ants\n"

#define STR_HELP_ELITISTANTS \
  "  -c, --elitistants    # number of elitist ants\n"

#define STR_HELP_LOCALSEARCH \
  "  -l, --localsearch    0: no local search   1: 2-opt   2: 2.5-opt   3: 3-opt\n"

#define STR_HELP_DLB \
  "  -d, --dlb            1 use don't look bits in local search\n"

#define STR_HELP_AS \
  "      --as              apply basic Ant System\n"

#define STR_HELP_EAS \
  "      --eas             apply elitist Ant System\n"

#define STR_HELP_RAS \
  "      --ras             apply rank-based version of Ant System\n"

#define STR_HELP_MMAS \
  "      --mmas            apply MAX-MIN ant system\n"

#define STR_HELP_BWAS \
  "      --bwas            apply best-worst ant system\n"

#define STR_HELP_ACS \
  "      --acs             apply ant colony system\n"

static const char *const STR_HELP[] = {
    STR_HELP_TSPLIBFILE ,
    STR_HELP_OPTIMUM ,
    STR_HELP_ANTS ,
    STR_HELP_ELITISTANTS ,
    STR_HELP_LOCALSEARCH ,
    STR_HELP_DLB ,
    STR_HELP_AS ,
    STR_HELP_EAS ,
    STR_HELP_RAS ,
    STR_HELP_MMAS ,
    STR_HELP_BWAS ,
    STR_HELP_ACS ,
    NULL
};

static const char *const STR_HELP_NEW[][3] = {
#define DEFINE_OPTIONS_PARAMETER(A, B, C, D, TYPE) DEFINE_PARAMETER(A,B,C,D,TYPE,0,0)
#define DEFINE_PARAMETER(PARAM, SHORT, LONG, HELP, TYPE, X, Y)                 \
    {SHORT, LONG, HELP},
#include "aco-parameters.def"
    { NULL, NULL, NULL}
};

static const char *const optstr__tsplibfile = "tsplibfile";
static const char *const optstr__optimum = "optimum";
static const char *const optstr__ants = "ants";
static const char *const optstr__elitistants = "elitistants";
static const char *const optstr__localsearch = "localsearch";
static const char *const optstr__dlb = "dlb";
static const char *const optstr__as = "as";
static const char *const optstr__eas = "eas";
static const char *const optstr__ras = "ras";
static const char *const optstr__mmas = "mmas";
static const char *const optstr__bwas = "bwas";
static const char *const optstr__acs = "acs";
#define DEFINE_OPTIONS_PARAMETER(A, B, C, D, TYPE) DEFINE_PARAMETER(A,B,C,D,TYPE,0,0)
#define DEFINE_PARAMETER(VAR, SHORT, LONG, D, TYPE, X, Y)                      \
    static const char *const optchar__##VAR = SHORT;                           \
    static const char *const optstr__##VAR = LONG;
#include "aco-parameters.def"

#define HANDLE_SHORT_OPTION(A)                                                 \
    if (option [1] != '\0') {                                                  \
        options->arg_##A = option + 1;                                         \
    } else if (++i < argc) {                                                   \
        options->arg_##A = argv [i];                                           \
    } else {                                                                   \
        goto error_missing_arg_short;                                          \
    }                                                                          \
    option = "\0";                                                             \
    options->opt_##A = 1;                                                      \
    break;


#define HANDLE_SHORT_BOOLEAN_OPTION(A)                                         \
    options->opt_##A = true;                                                   \
    return i + 1;


#define HANDLE_LONG_OPTION(A)                                                  \
    if (strncmp (option + 1, optstr__##A + 1, option_len - 1) == 0)            \
    {                                                                          \
        if (option_len <= 1)                                                   \
            goto error_long_opt_ambiguous;                                     \
        if (argument != 0)                                                     \
            options->arg_##A = argument;                                       \
        else if (++i < argc)                                                   \
            options->arg_##A = argv [i];                                       \
        else                                                                   \
        {                                                                      \
            option = optstr__##A;                                              \
            goto error_missing_arg_long;                                       \
        }                                                                      \
        options->opt_##A = 1;                                                  \
        break;                                                                 \
    } else

/* FIXME: How to automatically choose between the boolean and the standard
   HANDLE_LONG_OPTION */
#define HANDLE_LONG_BOOLEAN_OPTION(A)                                          \
    if (strncmp (option + 1, optstr__##A + 1, option_len - 1) == 0)            \
    {                                                                          \
        if (argument != 0)                                                     \
        {                                                                      \
            option = optstr__##A;                                              \
            goto error_unexpec_arg_long;                                       \
        }                                                                      \
        options->opt_##A = 1;                                                  \
        break;                                                             \
    } else

/* Parse command line options.  Return index of first non-option argument,
   or -1 if an error is encountered.  */
static int
parse_options (struct options *const options, const char *const program_name,
               const int argc, char **const argv)
{
  int i = 0;
#define DEFINE_OPTIONS_PARAMETER(A, B, C, D, TYPE) DEFINE_PARAMETER(A,B,C,D,TYPE,0,0)
#define DEFINE_PARAMETER(VAR, B, C, D, TYPE, X, Y)                             \
  options->opt_##VAR = 0;                                                      \
  options->arg_##VAR = 0;
#include "aco-parameters.def"

  options->opt_tsplibfile = 0;
  options->opt_optimum = 0;
  options->opt_ants = 0;
  options->opt_elitistants = 0;
  options->opt_localsearch = 0;
  options->opt_dlb = 0;
  options->opt_as = 0;
  options->opt_eas = 0;
  options->opt_ras = 0;
  options->opt_mmas = 0;
  options->opt_bwas = 0;
  options->opt_acs = 0;
  options->arg_tsplibfile = 0;
  options->arg_optimum = 0;
  options->arg_ants = 0;
  options->arg_elitistants = 0;
  options->arg_localsearch = 0;
  options->arg_dlb = 0;

  while (++i < argc)
  {
    const char *option = argv [i];
    if (*option != '-')
      return i;
    else if (*++option == '\0')
      return i;
    else if (*option == '-')
    {
      const char *argument;
      size_t option_len;
      ++option;
      if ((argument = strchr (option, '=')) == option)
        goto error_unknown_long_opt;
      else if (argument == 0)
        option_len = strlen (option);
      else
        option_len = argument++ - option;
      switch (*option)
      {
       case '\0':
        return i + 1;
       case 'a':
        if (strncmp (option + 1, optstr__acs + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
          {
            option = optstr__acs;
            goto error_unexpec_arg_long;
          }
          options->opt_acs = 1;
          break;
        }
        else if (strncmp (option + 1, optstr__ants + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
            options->arg_ants = argument;
          else if (++i < argc)
            options->arg_ants = argv [i];
          else
          {
            option = optstr__ants;
            goto error_missing_arg_long;
          }
          options->opt_ants = 1;
          break;
        }
        else if (strncmp (option + 1, optstr__as + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
          {
            option = optstr__as;
            goto error_unexpec_arg_long;
          }
          options->opt_as = 1;
          break;
        } else HANDLE_LONG_OPTION(alpha)
        goto error_unknown_long_opt;
       case 'b':
        if (strncmp (option + 1, optstr__bwas + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
          {
            option = optstr__bwas;
            goto error_unexpec_arg_long;
          }
          options->opt_bwas = 1;
          break;
        } else HANDLE_LONG_OPTION(opt_beta)
        goto error_unknown_long_opt;
       case 'd':
        if (strncmp (option + 1, optstr__dlb + 1, option_len - 1) == 0)
        {
          if (argument != 0)
            options->arg_dlb = argument;
          else if (++i < argc)
            options->arg_dlb = argv [i];
          else
          {
            option = optstr__dlb;
            goto error_missing_arg_long;
          }
          options->opt_dlb = 1;
          break;
        }
        HANDLE_LONG_OPTION(delta_n_ants)
        HANDLE_LONG_OPTION(delta_beta)
        HANDLE_LONG_OPTION(delta_rho)
        HANDLE_LONG_OPTION(delta_q0)
        goto error_unknown_long_opt;
       case 'e':
        if (strncmp (option + 1, optstr__eas + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
          {
            option = optstr__eas;
            goto error_unexpec_arg_long;
          }
          options->opt_eas = 1;
          break;
        }
        else if (strncmp (option + 1, optstr__elitistants + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
            options->arg_elitistants = argument;
          else if (++i < argc)
            options->arg_elitistants = argv [i];
          else
          {
            option = optstr__elitistants;
            goto error_missing_arg_long;
          }
          options->opt_elitistants = 1;
          break;
        }
        HANDLE_LONG_OPTION(end_n_ants)
        HANDLE_LONG_OPTION(end_beta)
        HANDLE_LONG_OPTION(end_rho)
        HANDLE_LONG_OPTION(end_q0)
        goto error_unknown_long_opt;
       case 'h':
        HANDLE_LONG_BOOLEAN_OPTION(flag_help)
        goto error_unknown_long_opt;
       case 'l':
        if (strncmp (option + 1, optstr__localsearch + 1, option_len - 1) == 0)
        {
          if (argument != 0)
            options->arg_localsearch = argument;
          else if (++i < argc)
            options->arg_localsearch = argv [i];
          else
          {
            option = optstr__localsearch;
            goto error_missing_arg_long;
          }
          options->opt_localsearch = 1;
          break;
        }
        goto error_unknown_long_opt;
       case 'm':
        if (strncmp (option + 1, optstr__mmas + 1, option_len - 1) == 0)
        {
          if (argument != 0)
          {
            option = optstr__mmas;
            goto error_unexpec_arg_long;
          }
          options->opt_mmas = 1;
          break;
        }
        goto error_unknown_long_opt;
       case 'n':
           HANDLE_LONG_OPTION(nn_ants)
           HANDLE_LONG_OPTION(nn_ls)
        goto error_unknown_long_opt;
       case 'o':
        if (strncmp (option + 1, optstr__optimum + 1, option_len - 1) == 0)
        {
          if (argument != 0)
            options->arg_optimum = argument;
          else if (++i < argc)
            options->arg_optimum = argv [i];
          else
          {
            option = optstr__optimum;
            goto error_missing_arg_long;
          }
          options->opt_optimum = 1;
          break;
        }
        goto error_unknown_long_opt;
       case 'p':
           HANDLE_LONG_OPTION(flag_ph_limits)
           HANDLE_LONG_OPTION(p_dec)
           goto error_unknown_long_opt;

       case 'q':
            HANDLE_LONG_BOOLEAN_OPTION(quiet_flag)
            HANDLE_LONG_OPTION(opt_q0)
            goto error_unknown_long_opt;
       case 'r':
        if (strncmp (option + 1, optstr__ras + 1, option_len - 1) == 0)
        {
          if (option_len < 3)
            goto error_long_opt_ambiguous;
          if (argument != 0)
          {
            option = optstr__ras;
            goto error_unexpec_arg_long;
          }
          options->opt_ras = 1;
          break;
        }
        else HANDLE_LONG_OPTION(ras_ranks)
        HANDLE_LONG_OPTION(opt_rho)
        HANDLE_LONG_OPTION(flag_restart)
        HANDLE_LONG_OPTION(restart_branch_factor)
        HANDLE_LONG_OPTION(restart_avg_distance)
        HANDLE_LONG_OPTION(min_iters_after_restart_best)
        goto error_unknown_long_opt;
      case 's':
          HANDLE_LONG_OPTION(seed)
          HANDLE_LONG_OPTION(schedule_length)
          HANDLE_LONG_OPTION(switch_n_ants)
          HANDLE_LONG_OPTION(switch_beta)
          HANDLE_LONG_OPTION(switch_rho)
          HANDLE_LONG_OPTION(switch_q0)
          goto error_unknown_long_opt;
       case 't':
           HANDLE_LONG_OPTION(max_time)
           HANDLE_LONG_OPTION(max_tours)
           HANDLE_LONG_OPTION(max_tries)

        if (strncmp (option + 1, optstr__tsplibfile + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
            options->arg_tsplibfile = argument;
          else if (++i < argc)
            options->arg_tsplibfile = argv [i];
          else
          {
            option = optstr__tsplibfile;
            goto error_missing_arg_long;
          }
          options->opt_tsplibfile = 1;
          break;
        }
      case 'v':
          HANDLE_LONG_OPTION(opt_var_n_ants)
          HANDLE_LONG_OPTION(opt_var_beta)
          HANDLE_LONG_OPTION(opt_var_rho)
          HANDLE_LONG_OPTION(opt_var_q0)
          goto error_unknown_long_opt;

      case 'x':
          HANDLE_LONG_OPTION(xi)
          goto error_unknown_long_opt;

      default:
       error_unknown_long_opt:
        fprintf (stderr, STR_ERR_UNKNOWN_LONG_OPT, program_name, option);
        return -1;
       error_long_opt_ambiguous:
        fprintf (stderr, STR_ERR_LONG_OPT_AMBIGUOUS, program_name, option);
        return -1;
       error_missing_arg_long:
        fprintf (stderr, STR_ERR_MISSING_ARG_LONG, program_name, option);
        return -1;
       error_unexpec_arg_long:
        fprintf (stderr, STR_ERR_UNEXPEC_ARG_LONG, program_name, option);
        return -1;
      }
    }
    else
      do
      {
        switch (*option)
        {
         case 'a':
             HANDLE_SHORT_OPTION(alpha)
         case 'b':
             HANDLE_SHORT_OPTION(opt_beta)
         case 'c':
          if (option [1] != '\0')
            options->arg_elitistants = option + 1;
          else if (++i < argc)
            options->arg_elitistants = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_elitistants = 1;
          break;
         case 'd':
          if (option [1] != '\0')
            options->arg_dlb = option + 1;
          else if (++i < argc)
            options->arg_dlb = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_dlb = 1;
          break;
         case 'e':
             HANDLE_SHORT_OPTION(opt_rho)
         case 'f':
             HANDLE_SHORT_OPTION(ras_ranks)
         case 'g':
             HANDLE_SHORT_OPTION(nn_ants)
         case 'h':
             HANDLE_SHORT_BOOLEAN_OPTION(flag_help)
         case 'i':
          if (option [1] != '\0')
            options->arg_tsplibfile = option + 1;
          else if (++i < argc)
            options->arg_tsplibfile = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_tsplibfile = 1;
          break;
         case 'k':
             HANDLE_SHORT_OPTION(nn_ls)
         case 'l':
          if (option [1] != '\0')
            options->arg_localsearch = option + 1;
          else if (++i < argc)
            options->arg_localsearch = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_localsearch = 1;
          break;
         case 'm':
          if (option [1] != '\0')
            options->arg_ants = option + 1;
          else if (++i < argc)
            options->arg_ants = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_ants = 1;
          break;
         case 'o':
          if (option [1] != '\0')
            options->arg_optimum = option + 1;
          else if (++i < argc)
            options->arg_optimum = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_optimum = 1;
          break;
         case 'q':
             HANDLE_SHORT_OPTION(opt_q0)
         case 'r':
             HANDLE_SHORT_OPTION(max_tries)
         case 's':
             HANDLE_SHORT_OPTION(max_tours)
         case 't':
             HANDLE_SHORT_OPTION(max_time)
         default:
          fprintf (stderr, STR_ERR_UNKNOWN_SHORT_OPT, program_name, *option);
          return -1;
         error_missing_arg_short:
          fprintf (stderr, STR_ERR_MISSING_ARG_SHORT, program_name, *option);
          return -1;
        }
      } while (*++option != '\0');
  }
  return i;
}

static void
check_out_of_range (double value, double minval, double maxval,
                    const char *optionName)
/*
      FUNCTION: check whether parameter values are within allowed range
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    if (value < minval || value > maxval) {
        fprintf (stderr,"Error: Option `--%s' out of range [%g, %g]\n",
                 optionName, minval, maxval);
        exit(1);
    }
}

int parse_commandline (int argc, char *argv [])
{
    int i;
    const char *progname;
    struct options options;

    progname = argv [0] != NULL && *(argv [0]) != '\0'
	? argv [0]
	: "acotsp";

    i = parse_options (&options, progname, argc, argv);

    if (i < 2)
    {
	fprintf (stderr, "No options are specified\n");
	fprintf (stderr, "Try `%s --help' for more information.\n",
		 progname);
	exit(1);
    }

    if (options.opt_flag_help)
    {
        int k;
	printf ("Usage: %s [OPTION]... [ARGUMENT]...\n"
		"Options:\n", progname);
        for (k = 0; STR_HELP[k]; k++)
            printf ("%s", STR_HELP[k]);
        for (k = 0; STR_HELP_NEW[k][0]; k++) {
            if (strlen(STR_HELP_NEW[k][0]) == 0) {
                printf ("     --%-22s    %s\n", STR_HELP_NEW[k][1], STR_HELP_NEW[k][2]);
            } else {
                printf (" -%s, --%-22s    %s\n", STR_HELP_NEW[k][0], STR_HELP_NEW[k][1], STR_HELP_NEW[k][2]);
            }
        }
	exit(0);
    }

    /* puts ("\t OPTIONS:"); */



    if ( options.opt_optimum )
    {
	optimal = atol(options.arg_optimum);
 	/* fputs ("  -o  --optimum ", stdout); */
	/* if (options.arg_optimum != NULL) */
	    /* printf ("with argument \"%ld\"\n", optimal); */
    } else {
        /* fprintf(stdout,"\tNote: optimal solution value is set to default %ld\n", */
        /*         optimal); */
    }



    if ( options.opt_tsplibfile )
    {
        if (strlen(options.arg_tsplibfile) >= LINE_BUF_LEN) {
            fprintf (stderr, "error: too long input filename '%s', maximum length is %d",
                     options.arg_tsplibfile, LINE_BUF_LEN);
            exit (1);
        }
	strcpy (name_buf, options.arg_tsplibfile);
 	/* fputs ("  -i  --tsplibfile ", stdout); */
	/* if (options.arg_tsplibfile != NULL) */
	    /* printf ("with argument \"%s\"\n", name_buf ); */
    }


    if (options.opt_as + options.opt_eas + options.opt_ras + options.opt_mmas
        + options.opt_bwas + options.opt_acs > 1) {
        fprintf (stderr, "error: more than one ACO algorithm enabled in the command line");
        exit (1);
    } else if (options.opt_as + options.opt_eas + options.opt_ras + options.opt_mmas
               + options.opt_bwas + options.opt_acs == 1)  {
        as_flag = eas_flag = ras_flag = mmas_flag = bwas_flag = acs_flag = FALSE;
    }

    if (options.opt_as || as_flag) {
        as_flag = TRUE;
        set_default_as_parameters();
        /* fprintf(stdout,"as_flag is set to 1, run Ant System\n"); */
    }

    if (options.opt_eas || eas_flag) {
        eas_flag = TRUE;
        set_default_eas_parameters();
        fprintf(stdout,"eas_flag is set to 1, run elitist Ant System\n");
    }

    if (options.opt_ras || ras_flag) {
        ras_flag = TRUE;
        set_default_ras_parameters();
        fprintf(stdout,"ras_flag is set to 1, run rank-based Ant System\n");
    }

    if (options.opt_mmas || mmas_flag) {
        mmas_flag = TRUE;
        set_default_mmas_parameters();
        /* fprintf(stdout,"mmas_flag is set to 1, run MAX-MIN Ant System\n"); */
    }

    if (options.opt_bwas || bwas_flag) {
        bwas_flag = TRUE;
        set_default_bwas_parameters();
        fprintf(stdout,"bwas_flag is set to 1, run Best-Worst Ant System\n");
    }

    if ( options.opt_acs || acs_flag ) {
        acs_flag = TRUE;
        set_default_acs_parameters();
        fprintf(stdout,"acs_flag is set to 1, run Ant Colony System\n");
    }

    /* Problem-specific parameters settings may override the algorithm-specific
       settings. */
    problem_set_default_parameters ();

    if ( options.opt_localsearch ) {
        ls_flag = atol(options.arg_localsearch);
        /* fputs ("  -l  --localsearch ", stdout); */
        /* if (options.arg_localsearch != NULL) */
        /*     printf ("with argument \"%u\"\n", (unsigned int) ls_flag); */
        check_out_of_range(ls_flag, 0, LS_MAX, "ls_flag");
    } else {
        /* fprintf(stdout,"\tNote: local search flag is set to default %u (%s)\n", */
        /*         (unsigned int) ls_flag, ls_type_to_string(ls_flag)); */
    }

    if (ls_flag) {
        set_default_ls_parameters();
    }

#define ASSIGN_ENUM_PARAMETER(PARAM, TYPE)                                     \
    do {                                                                       \
        if ( options.opt_##PARAM )                                             \
        {                                                                      \
            int opt = -1;                                                      \
            int k;                                                             \
            for (k = 0; strings_##TYPE[k]; k++) {                              \
                if (strncmp (options.arg_##PARAM, strings_##TYPE[k],           \
                             strlen (strings_##TYPE[k])) == 0) {               \
                    opt = k; break;                                            \
                }                                                              \
            }                                                                  \
            if (opt < 0) {                                                     \
                printf ("Unknown value \"%s\" for option --%s\n",              \
                        options.arg_##PARAM, optstr__##PARAM);                 \
                exit(1);                                                       \
            }                                                                  \
            PARAM = opt;                                                       \
        }                                                                      \
    } while(0)                                                                 \


#define ASSIGN_DOUBLE_PARAMETER(PARAM, MIN_RANGE, MAX_RANGE)                   \
    do {                                                                       \
        if ( options.opt_##PARAM ) {                                           \
            PARAM = atof(options.arg_##PARAM);                                 \
        }                                                                      \
    } while(0)


#define ASSIGN_INTEGER_PARAMETER(PARAM, MIN_RANGE, MAX_RANGE)                  \
    do {                                                                       \
        if ( options.opt_##PARAM ) {                                           \
            PARAM = atol(options.arg_##PARAM);                                 \
        }                                                                      \
    } while(0)

#define ASSIGN_BOOLEAN_PARAMETER(PARAM)                                        \
    do {                                                                       \
        if (options.opt_##PARAM) {                                             \
            PARAM = TRUE;                                                      \
        }                                                                      \
    } while(0)

#define DEFINE_BOOLEAN_PARAMETER(PARAM, MIN_RANGE, MAX_RANGE) ASSIGN_BOOLEAN_PARAMETER(PARAM);
#define DEFINE_DOUBLE_PARAMETER(PARAM, MIN_RANGE, MAX_RANGE) ASSIGN_DOUBLE_PARAMETER(PARAM, MIN_RANGE, MAX_RANGE);
#define DEFINE_INTEGER_PARAMETER(PARAM, MIN_RANGE, MAX_RANGE) ASSIGN_INTEGER_PARAMETER(PARAM, MIN_RANGE, MAX_RANGE);
#define DEFINE_OPTIONS_PARAMETER(PARAM, B, C, D, TYPE) ASSIGN_ENUM_PARAMETER(PARAM,TYPE);
#define DEFINE_PARAMETER(PARAM,B,C,D, TYPE, MIN_RANGE, MAX_RANGE) TYPE(PARAM, MIN_RANGE, MAX_RANGE)
#include "aco-parameters.def"

    if (options.opt_schedule_length && !mmas_flag) {
        fprintf (stderr, "error: --schedule-length only has effect with --mmas\n");
        exit(1);
    }

    if (flag_ph_limits && acs_flag) {
        fprintf (stderr, "error: --ph-limits 1 has not effect with --acs\n");
        exit(1);
    }

    if ( options.opt_ants ) {
	opt_n_ants = atol(options.arg_ants);
	/* fputs ("  -m  --ants ", stdout); */
	/* if (options.arg_ants != NULL) */
	/*     printf ("with argument \"%ld\"\n", opt_n_ants); */
	/* check_out_of_range( opt_n_ants, -1, MAX_ANTS-1, "opt_ants"); */
    } else {
        /* if (opt_n_ants < 0) */
        /*     fprintf(stdout,"\tNote: number of ants is set to default n\n"); */
        /* else */
        /*     fprintf(stdout,"\tNote: number of ants is set to default %ld\n", */
        /*             opt_n_ants); */
    }


    if ( options.opt_elitistants ) {
	elitist_ants = atol(options.arg_elitistants);
	/* fputs ("  -c  --elitistants ", stdout); */
	/* if (options.arg_elitistants != NULL) */
	/*     printf ("with argument \"%ld\"\n", elitist_ants); */
        /* check_out_of_range(elitist_ants, 0, LONG_MAX, "elitistants"); */
    } else {
        /* if (elitist_ants <= 0) { */
        /*     fprintf(stdout,"\tNote: number of elitist ants is set to default n\n"); */
        /* } else { */
        /*     fprintf(stdout,"\tNote: number of elitist ants is set to default %ld\n", elitist_ants); */
        /* } */
    }

    if ( options.opt_dlb ) {
	dlb_flag = atol(options.arg_dlb);
	/* fputs ("  -d  --dlb ", stdout); */
	/* if (options.arg_dlb != NULL) */
	/*     printf ("with argument \"%ld\"\n", dlb_flag); */
	/* check_out_of_range( dlb_flag, 0, 1, "dlb_flag"); */
    /* } else { */
    /*     fprintf(stdout,"\tNote: dlb flag is set to default %d (%s don't look bits)\n", */
    /*             dlb_flag ? 1 : 0, dlb_flag ? "use" : "not use"); */
    }

    /* fputs ("Non-option arguments:", stderr); */
    while (i < argc) {
	fprintf (stderr,"  \"%s\"\n", argv [i++]);
	fprintf (stderr,"\nThere were non-option arguments\n");
	fprintf (stderr,"I suspect there is something wrong, maybe wrong option name; exit\n");
	exit(1);
    }

    return 0;
}
