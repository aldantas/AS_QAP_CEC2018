/*

       AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
      AA  AA  CC     OO  OO    TT    SS      PP  PP
      AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
      AA  AA  CC     OO  OO    TT        SS  PP
      AA  AA   CCCC   OOOO     TT    SSSSS   PP

######################################################
##########    ACO algorithms for the TSP    ##########
######################################################

      Version: 1.0
      File:    InOut.c
      Author:  Thomas Stuetzle
      Purpose: mainly input / output / statistic routines
      Check:   README and gpl.txt
      Copyright (C) 2002  Thomas Stuetzle
*/

/***************************************************************************

    Program's name: acotsp

    Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the
    symmetric TSP

    Copyright (C) 2004  Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    email: stuetzle no@spam ulb.ac.be
    mail address: Universite libre de Bruxelles
                  IRIDIA, CP 194/6
                  Av. F. Roosevelt 50
                  B-1050 Brussels
		  Belgium

***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "InOut.h"
#include "timer.h"
#include "utilities.h"
#include "ants.h"
#include "problem.h"
#include "ls.h"
#include "parse.h"

#ifndef VERSION
#define VERSION "0.99"
#endif

long int *best_in_try;
long int *best_found_at;
double   *time_best_found;
double   *time_total_run;

long int n_try;
long int n_tours;
long int iteration;         /* iteration counter */
long int restart_iteration; /* remember iteration when restart was done if any */
double   restart_time;      /* remember time when restart was done if any */
long int max_tries;         /* maximum number of independent tries */
long int max_tours;         /* maximum number of tour constructions in one try */
long int seed;

const double   lambda = 0.05; /* Parameter to determine branching factor */
int     restart_freq;

double   max_time;          /* maximal allowed run time of a try  */
double   time_used;         /* time used until some given event */
double   time_passed;       /* time passed until some moment*/
long int optimal;           /* optimal solution or bound to find */


double   mean_ants;         /* average tour length */
double   stddev_ants;       /* stddev of tour lengths */
double   branching_factor;  /* average node branching factor when searching */
double   found_branching;   /* branching factor when best solution is found */

long int found_best;        /* iteration in which best solution is found */
long int restart_found_best;/* iteration in which restart-best solution is found */

/* ------------------------------------------------------------------------ */

FILE *report, *comp_report, *stat_report;

char name_buf[LINE_BUF_LEN];



void write_report( void )
/*
      FUNCTION: output some info about trial (best-so-far solution quality, time)
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    /* printf("best %ld, iteration: %ld, time %.2f\n", */
    /*        best_so_far_ant->tour_length, iteration, time_used); */
    /* if (comp_report) */
    /*     fprintf(comp_report, */
    /*             "best %ld\t iteration %ld\t tours %ld\t time %.3f\n", */
    /*             best_so_far_ant->tour_length, iteration, n_tours, time_used); */
}


static void
fprintf_parameters (FILE *stream)
{
    fprintf(stream,"max_tours\t\t %ld\n", max_tours);
    fprintf(stream,"max_time\t\t %.2f\n", max_time);
    fprintf(stream,"seed\t\t %ld\n", seed);
    fprintf(stream,"optimum\t\t\t %ld\n", optimal);
    fprintf(stream,"n_ants\t\t\t %ld\n", n_ants);
    fprintf(stream,"nn_ants\t\t\t %ld\n", nn_ants);
    fprintf(stream,"alpha\t\t\t %.2f\n", alpha);
    fprintf(stream,"beta\t\t\t %.2f\n", beta);
    fprintf(stream,"rho\t\t\t %.2f\n", rho);
    fprintf(stream,"q_0\t\t\t %.2f\n", q0);
    fprintf(stream,"elitist_ants\t\t %ld\n", elitist_ants);
    fprintf(stream,"ras_ranks\t\t %ld\n", ras_ranks);
    fprintf(stream,"ls_flag\t\t\t %u\n", (unsigned int) ls_flag);
    fprintf(stream,"nn_ls\t\t\t %ld\n", nn_ls);
    fprintf(stream,"dlb_flag\t\t %ld\n", dlb_flag);
    fprintf(stream,"as_flag\t\t\t %ld\n", as_flag);
    fprintf(stream,"eas_flag\t\t %ld\n", eas_flag);
    fprintf(stream,"ras_flag\t\t %ld\n", ras_flag);
    fprintf(stream,"mmas_flag\t\t %ld\n", mmas_flag);
    fprintf(stream,"bwas_flag\t\t %ld\n", bwas_flag);
    fprintf(stream,"acs_flag\t\t %ld\n", acs_flag);

#define DEFINE_BOOLEAN_PARAMETER(PARAM)                                        \
    fprintf(stream, "%-20s %10ld\n", #PARAM, (long int) PARAM);
#define DEFINE_DOUBLE_PARAMETER(PARAM) \
    fprintf(stream, "%-20s %13.2f\n", #PARAM, PARAM);
#define DEFINE_INTEGER_PARAMETER(PARAM) \
    fprintf(stream, "%-20s %10ld\n", #PARAM, (long int) PARAM);
#define DEFINE_OPTIONS_PARAMETER(PARAM, B, C, D, TYPE) \
    fprintf(stream, "%-20s %10s\n", #PARAM, strings_##TYPE[PARAM]);
#define DEFINE_OPTION(E,S) S
#define DEFINE_OPTIONS(NAME, ...)                                              \
    const char * const strings_##NAME [] = { __VA_ARGS__, NULL };
#define DEFINE_PARAMETER(PARAM,B,C,D, TYPE, X, Y) TYPE(PARAM)
#include "aco-parameters.def"
}

void print_default_parameters()
/*
      FUNCTION: output default parameter settings
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    fprintf(stderr,"\nDefault parameter settings are:\n\n");
    fprintf_parameters (stderr);
}


void
set_default_as_parameters (void)
{
    assert (as_flag);
    opt_n_ants     = -1;    /* number of ants (-1 means instance size) */
    alpha          = 1.0;
    opt_beta       = 2.0;
    opt_rho        = 0.5;
    opt_q0         = 0.0;
    ras_ranks      = 0;
    elitist_ants   = 0;
}

void
set_default_eas_parameters (void)
{
    assert (eas_flag);
    opt_n_ants     = -1;    /* number of ants (-1 means instance size) */
    alpha          = 1.0;
    opt_beta       = 2.0;
    opt_rho        = 0.5;
    opt_q0         = 0.0;
    ras_ranks      = 0;
    elitist_ants   = opt_n_ants;
}

void
set_default_ras_parameters (void)
{
    assert (ras_flag);
    opt_n_ants     = -1;    /* number of ants (-1 means instance size) */
    alpha          = 1.0;
    opt_beta       = 2.0;
    opt_rho        = 0.1;
    opt_q0         = 0.0;
    ras_ranks      = 6;
    elitist_ants   = 0;
}

void
set_default_bwas_parameters (void)
{
    assert (bwas_flag);
    opt_n_ants     = -1;    /* number of ants (-1 means instance size) */
    alpha          = 1.0;
    opt_beta       = 2.0;
    opt_rho        = 0.1;
    opt_q0         = 0.0;
    ras_ranks      = 0;
    elitist_ants   = 0;
    flag_restart   = RESTART_AVG_DISTANCE;
}

void
set_default_mmas_parameters (void)
{
    assert (mmas_flag);
    opt_n_ants     = 5;    /* number of ants (-1 means instance size) */
    alpha          = 1.0;
    opt_beta       = 1.0;
    opt_rho        = 0.2;
    rho        = 0.2;
    opt_q0         = 0.0;
    ras_ranks      = 0;
    elitist_ants   = 0;
    flag_restart   = RESTART_BRANCH_FACTOR;
    flag_ph_limits = true;
    quiet_flag     = true;
}

void
set_default_acs_parameters (void)
{
    assert (acs_flag);
    opt_n_ants     = 10;    /* number of ants (-1 means instance size) */
    alpha          = 1.0;
    opt_beta       = 2.0;
    opt_rho        = 0.1;
    opt_q0         = 0.9;
    ras_ranks      = 0;
    elitist_ants   = 0;
    xi             = 0.1;   /* as suggested by Gambardella and Dorigo for the TSP. */
}

void set_default_ls_parameters(void)
{
    assert (ls_flag);
    alpha          = 1.0;
    opt_q0         = 0.0;

    if (acs_flag) {
        opt_q0 = 0.98;
    }
    problem_set_default_ls_parameters ();
}

void set_default_parameters(void)
/*
      FUNCTION: set default parameter settings
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    n_ants     = 5;     /* number of ants */
    opt_n_ants     = 5;     /* number of ants */
    end_n_ants     = opt_n_ants;
    alpha          = 1.0;
    beta           = 1.0;
    rho            = 2.0;
    opt_beta       = 1.0;
    end_beta       = opt_beta;
    opt_rho        = 0.2;
    end_rho        = opt_rho;
    xi             = 0.0;   /* only used for ACS */
    opt_q0         = 0.0;
    end_q0         = opt_q0;
    max_tries      = 1;
    max_tours      = 0;
    seed           = (long int) time(NULL);
    max_time       = 600;
    optimal        = 1;
    as_flag        = false;
    eas_flag       = false;
    ras_flag       = false;
    mmas_flag      = true;
    bwas_flag      = false;
    acs_flag       = false;
    ras_ranks      = 0;
    elitist_ants   = 0;
    flag_ph_limits = false;
    flag_restart   = RESTART_NEVER;
    restart_avg_distance = 0.05;
    quiet_flag     = true;
}



void population_statistics ( void )
/*
      FUNCTION:       compute some population statistics like average tour length,
                      standard deviations, average distance, branching-factor and
		      output to a file gathering statistics
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  none
*/
{
    long int j, k;
    long int *l;
    double   pop_mean, pop_stddev, avg_distance = 0.0;

    l = malloc(n_ants * sizeof(long int));
    for( k = 0 ; k < n_ants ; k++ ) {
	l[k] = ant[k].tour_length;
    }

    pop_mean = mean( l, n_ants);
    pop_stddev = std_deviation( l, n_ants, pop_mean );
    branching_factor = node_branching(lambda);

    for ( k = 0 ; k < n_ants-1 ; k++ )
	for ( j = k+1 ; j < n_ants ; j++) {
	    avg_distance += (double)distance_between_ants( &ant[k], &ant[j]);
	}
    avg_distance /= ((double)n_ants * (double)(n_ants-1) / 2.);

    /* if (stat_report) { */
    /*     fprintf(stat_report,"%ld\t%.1f\t%.5f\t%.7f\t%.5f\t%.1f\t%.1f\t%.5f\n", */
    /*             iteration, pop_mean, pop_stddev, pop_stddev / pop_mean, */
    /*             branching_factor, (branching_factor - 1.) * (double)n, */
    /*             avg_distance, avg_distance / (double)n); */
    /* } */
    free (l);
}



void output_solution( void )
/*
      FUNCTION:       output a solution together with node coordinates
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       not used in the default implementation but may be useful anyway
*/
{
    long int i;
    if (stat_report) {
        for ( i = 0 ; i < n ; i++ ) {
            fprintf(stat_report," %ld",
                    best_so_far_ant->tour[i]);
/*
            fprintf(stat_report," %ld %f %f\n",
                    best_so_far_ant->tour[i],
                    instance.nodeptr[best_so_far_ant->tour[i]].x,
                    instance.nodeptr[best_so_far_ant->tour[i]].y);
*/
        }
        fprintf(stat_report, "\n");
    }
}



void exit_try( long int ntry )
/*
      FUNCTION:       save some statistical information on a trial once it finishes
      INPUT:          trial number
      OUTPUT:         none
      COMMENTS:
*/
{
  checkTour( best_so_far_ant->tour );
/*    printTourFile( best_so_far_ant->tour ); */

  /* printf("\n Best Solution in try %ld is %ld\n",ntry, best_so_far_ant->tour_length); */
  /* if (report) */
  /*     fprintf(report, "Best: %ld\t Iterations: %6ld\t B-Fac %.5f\t Time %.2f\t Tot.time %.2f\n", */
  /*             best_so_far_ant->tour_length, found_best, found_branching, */
  /*             time_used, elapsed_time( VIRTUAL )); */
  /* printf(" Best Solution was found after %ld iterations\n", found_best); */

  best_in_try[ntry] = best_so_far_ant->tour_length;
  best_found_at[ntry] = found_best;
  time_best_found[ntry] = time_used;
  time_total_run[ntry] = elapsed_time( VIRTUAL );
  printf("Best %ld, found at iteration %ld, found at time %f\n",
		  best_in_try[ntry], best_found_at[ntry],
		  time_best_found[ntry]);

  if (comp_report) fprintf(comp_report,"end try %ld\n\n",ntry);
  if (stat_report) fprintf(stat_report,"end try %ld\n\n",ntry);
#if TRACE
  output_solution();
#endif
  if (report) fflush(report);
  if (comp_report) fflush(comp_report);
  if (stat_report) fflush(stat_report);
}



void exit_program( void )
/*
      FUNCTION:       save some final statistical information on a trial once it finishes
      INPUT:          none
      OUTPUT:         none
      COMMENTS:
*/
{
  long int best_tour_length, worst_tour_length;
  double   t_avgbest, t_stdbest, t_avgtotal, t_stdtotal;
  double   avg_sol_quality = 0., avg_cyc_to_bst = 0., stddev_best, stddev_iterations;

  best_tour_length = best_of_vector( best_in_try ,max_tries );
  worst_tour_length = worst_of_vector( best_in_try , max_tries );

  avg_cyc_to_bst = mean( best_found_at , max_tries );
  stddev_iterations = std_deviation( best_found_at, max_tries, avg_cyc_to_bst );

  avg_sol_quality = mean( best_in_try , max_tries );
  stddev_best = std_deviation( best_in_try, max_tries, avg_sol_quality);

  t_avgbest = meanr( time_best_found, max_tries );
  printf(" t_avgbest = %.4f\n", t_avgbest );
  t_stdbest = std_deviationr( time_best_found, max_tries, t_avgbest);

  t_avgtotal = meanr( time_total_run, max_tries );
  printf(" t_avgtotal = %.4f\n", t_avgtotal );
  t_stdtotal = std_deviationr( time_total_run, max_tries, t_avgtotal);

  if (report) {
      fprintf(report,"\nAverage-Best: %.2f\t Average-Iterations: %.2f", avg_sol_quality, avg_cyc_to_bst);
      fprintf(report,"\nStddev-Best: %.2f \t Stddev Iterations: %.2f", stddev_best, stddev_iterations);
      fprintf(report,"\nBest try: %ld\t\t Worst try: %ld\n", best_tour_length, worst_tour_length);
      fprintf(report,"\nAvg.time-best: %.2f stddev.time-best: %.2f\n", t_avgbest, t_stdbest);
      fprintf(report,"\nAvg.time-total: %.2f stddev.time-total: %.2f\n", t_avgtotal, t_stdtotal);

      if ( optimal > 0 ) {
          fprintf(report,
                  " excess best = %f, excess average = %f, excess worst = %f\n",
                  (double)(best_tour_length - optimal) / (double)optimal,
                  (double)(avg_sol_quality - optimal) / (double)optimal,
                  (double)(worst_tour_length - optimal) / (double)optimal);
      }
  }

  if (comp_report)
      fprintf(comp_report,"end problem %s\n", get_instance_name(&instance));
}



void init_program( long int argc, char *argv[] )
/*
      FUNCTION:       initialize the program,
      INPUT:          program arguments, needed for parsing commandline
      OUTPUT:         none
      COMMENTS:
*/
{

  char temp_buffer[LINE_BUF_LEN];

  /* printf("\n%s, version %s\n", PROG_ID_STR, VERSION); */
  set_default_parameters();
  setbuf(stdout,NULL);
  parse_commandline(argc, argv);
  /* print_default_parameters(); */

  assert (max_tries <= MAXIMUM_NO_TRIES);

  best_in_try = calloc(max_tries, sizeof(long int));
  best_found_at = calloc(max_tries, sizeof(long int));
  time_best_found = calloc(max_tries, sizeof(double));
  time_total_run = calloc(max_tries, sizeof(double));

  /* trace_print("read problem data  `%s` ...\n", name_buf); */
  read_instance(name_buf, &instance);
  /* trace_print("... done\n\n"); */

  if (opt_n_ants < 0) opt_n_ants = n;
  /* default setting for elitist_ants is 0; if EAS is applied and
     option elitist_ants is not used, we set the default to
     elitist_ants = n, that is, the instance size.  */
  /* FIXME: the parameter should be a factor (0, large] of the number of ants.  */
  if (eas_flag && elitist_ants <= 0) elitist_ants = n;

  nn_ls = MIN(n - 1, nn_ls);
  nn_ants = (nn_ants <= 0) ? n : MIN(n - 1, nn_ants);

  /* Set maximum value from option setting. */
  max_n_ants = MAX(opt_n_ants, end_n_ants);

  /* Set runtime value from option setting. */
  n_ants = opt_n_ants;
  beta = opt_beta;
  rho = opt_rho;
  q0 = opt_q0;

  assert (max_n_ants < MAX_ANTS-1);
  assert (nn_ants < MAX_NEIGHBOURS);
  assert (nn_ants > 0);
  assert (nn_ls > 0);

  assert (n_ants < MAX_ANTS-1);
  assert (nn_ants < MAX_NEIGHBOURS);
  assert (nn_ants > 0);
  assert (nn_ls > 0);

  const char * instance_name = get_instance_name(&instance);
  if (!quiet_flag) {
      sprintf(temp_buffer,"best.%s", instance_name);
      trace_print("%s\n",temp_buffer);
      report = fopen(temp_buffer, "w");
      sprintf(temp_buffer,"cmp.%s", instance_name);
      trace_print("%s\n",temp_buffer);
      comp_report = fopen(temp_buffer, "w");
      sprintf(temp_buffer,"stat.%s", instance_name);
      trace_print("%s\n",temp_buffer);
      stat_report = fopen(temp_buffer, "w");
      printf("oi\n");
  } else {
      report = NULL;
      comp_report = NULL;
      stat_report = NULL;
  }

  /* write_params(); */
  /* if (comp_report) */
  /*     fprintf(comp_report,"begin problem %s\n",name_buf); */
  /* trace_print("allocate ants' memory ...\n"); */
  allocate_ants();
  /* trace_print("... done\n"); */
}






void printTrail(void)
/*
      FUNCTION:       print pheromone trail values
      INPUT:          none
      OUTPUT:         none
*/
{
  long int i,j;

  printf("pheromone Trail matrix, iteration: %ld\n\n",iteration);
  for ( i = 0 ; i < n ; i++) {
    printf("From %ld:  ",i);
    for ( j = 0 ; j < n ; j++ ) {
      printf(" %.10f ", pheromone[i][j]);
      if ( pheromone[i][j] > 1.0 )
	printf("XXXXX\n");
    }
    printf("\n");
  }
  printf("\n");
}



void printTotal(void)
/*
      FUNCTION:       print values of pheromone times heuristic information
      INPUT:          none
      OUTPUT:         none
*/
{
  long int i, j;

  printf("combined pheromone and heuristic info\n\n");
  for (i=0; i < n; i++) {
    for (j = 0; j < n - 1 ; j++) {
      printf(" %.15f &", total[i][j]);
      if ( total[i][j] > 1.0 )
	printf("XXXXX\n");
    }
    printf(" %.15f\n", total[i][n-1]);
    if ( total[i][n-1] > 1.0 )
      printf("XXXXX\n");
  }
  printf("\n");
}



void printProbabilities(void)
/*
      FUNCTION:       prints the selection probabilities as encountered by an ant
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       this computation assumes that no choice has been made yet.
*/
{
  long int i, j;
  double   *p;
  double   sum_prob;

  printf("Selection Probabilities, iteration: %ld\n",iteration);
  p = calloc( n, sizeof(double) );

  for (i=0; i < n; i++) {
    printf("From %ld:  ",i);
    sum_prob = 0.;
    for ( j = 0 ; j < n ; j++) {
      if ( i == j )
	p[j] = 0.;
      else
	p[j] = total[i][j];
      sum_prob += p[j];
    }
    for ( j = 0 ; j < n ; j++) {
      p[j] = p[j] / sum_prob;
    }
    for ( j = 0 ; j < n-1 ; j++) {
      printf(" %.5f ", p[j]);
    }
    printf(" %.5f\n", p[n-1]);
    if (!(j % 26)) {
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
  free ( p );
}



void fprintTour(FILE * stream, const long int *t )
/*
      FUNCTION:       print the tour *t
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    int   i;
    fprintf(stream, "\n");
    for( i = 0 ; i <= n ; i++ ) {
	if (!i%25)
	    fprintf(stream, "\n");
	fprintf(stream, "%ld ", t[i]);
    }
    fprintf(stream, "\n");
    fprintf(stream, "Tour length = %ld\n\n", compute_tour_length( t ));
}

void printTour(const long int *t )
/*
      FUNCTION:       print the tour *t
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    fprintTour(stdout, t);
}



void printTourFile( long int *t )
/*
      FUNCTION:       print the tour *t to cmp.tsplibfile
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    long int   i;
    if (!comp_report) return;

    fprintf(comp_report,"begin solution\n");
    for( i = 0 ; i < n ; i++ ) {
	fprintf(comp_report,"%ld ", t[i]);
    }
    fprintf(comp_report,"\n");
    fprintf(comp_report,"Tour length %ld\n",compute_tour_length( t ));
    fprintf(comp_report,"end solution\n");
}



void checkTour( long int *t )
/*
      FUNCTION:       make a simple check whether tour *t can be feasible
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    long int   i, sum=0;

    for( i = 0 ; i < n ; i++ ) {
	sum += t[i];
    }
    if ( sum != (n-1) * n / 2 ) {
	fprintf(stderr,"Next solution must be flawed !!\n");
        check_solution (t);
	fprintTour(stderr, t);
	exit(1);
    }
}



void write_params( void )
/*
      FUNCTION:       writes chosen parameter settings in standard output and in
                      report files
      INPUT:          none
      OUTPUT:         none
*/
{
    fprintf(stdout, "\nParameter settings: \n\n");
    fprintf_parameters (stdout);
    fprintf(stdout, "\n");

  if (report) {
      fprintf(report,"\nParameter settings: \n\n");
      fprintf_parameters (report);
      fprintf(report,"\n");
  }

  if (comp_report) {
      fprintf(comp_report,"\n%s, version %s\n", PROG_ID_STR, VERSION);
      fprintf(comp_report,"\nParameter settings: \n\n");
      fprintf_parameters (comp_report);
      fprintf(comp_report,"\n");
  }
}
