#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "problem.h"
#include "ants.h"
#include "utilities.h"
#include "InOut.h"
#include "timer.h"
#include "aco.h"
#include "ls.h"
#include "adaptation.h"

/* This is a short-hand for instance->nn_list, allowing to have different
   implementations of "instance" per problem.  */
static long int **nn_list;

static long int termination_condition( void )
/*
      FUNCTION:       checks whether termination condition is met
      INPUT:          none
      OUTPUT:         0 if condition is not met, number neq 0 otherwise
      (SIDE)EFFECTS:  none
*/
{
  return ( ((n_tours >= max_tours) || (elapsed_time( VIRTUAL ) >= max_time)));
}

void apply_local_search( void )
/*
      FUNCTION:       manage the local search phase; apply local search to ALL ants; in
                      dependence of ls_flag one of 2-opt, 2.5-opt, and 3-opt local search
		      is chosen.
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants of the colony have locally optimal tours
      COMMENTS:       typically, best performance is obtained by applying local search
                      to all ants. It is known that some improvements (e.g. convergence
		      speed towards high quality solutions) may be obtained for some
		      ACO algorithms by applying local search to only some of the ants.
		      Overall best performance is typically obtained by using 3-opt.
*/
{
    long int k;

    /* trace_print("apply local search to all ants\n"); */

    for ( k = 0 ; k < n_ants ; k++ ) {
        local_search (ls_flag, ant[k].tour);
	ant[k].tour_length = compute_tour_length( ant[k].tour );
        if (termination_condition()) return;
    }
}


static void init_try( long int ntry )
/*
      FUNCTION: initialize variables appropriately when starting a trial
      INPUT:    trial number
      OUTPUT:   none
      COMMENTS: none
*/
{
    /* trace_print("INITIALIZE TRIAL\n"); */

    start_timers();
    time_used = elapsed_time( VIRTUAL );
    time_passed = time_used;

    if (comp_report) {
        /* fprintf(comp_report,"seed %ld\n",seed); */
        fflush(comp_report);
    }

    adapt_parameters_init();

    /* Initialize variables concerning statistics etc. */
    n_tours      = 0;
    iteration    = 1;
    restart_iteration = 1;
    best_so_far_ant->tour_length = INFTY;
    found_best   = 0;
    u_gb         = INFTY;

    /* FIXME: Why save_seed? */
    long int save_seed = seed;
    /* Initialize the Pheromone trails */
    double nn_tour_length = nn_tour(best_so_far_ant);
    trail_max = 1. / (rho * nn_tour_length);
    trail_min = trail_max / ( 2. * n );
    seed = save_seed;

    if ( acs_flag || bwas_flag ) {
	trail_0 = 1. / ( (double) n * nn_tour_length);
    } else if ( mmas_flag ) {
	trail_0 = trail_max;
    } else {
	trail_0 = trail_max;
	/* in the original papers on Ant System, Elitist Ant System, and
	   Rank-based Ant System it is not exactly defined what the
	   initial value of the pheromones is. Here we set it to some
	   small constant, analogously as done in MAX-MIN Ant System.
	*/
    }

    init_pheromone_trails( trail_0 );

    /* Calculate combined information pheromone times heuristic information */
    compute_total_information();

    /* printf("Init trial end \n"); */
    if (comp_report) fprintf(comp_report,"begin try %li \n",ntry);
    if (stat_report) fprintf(stat_report,"begin try %li \n",ntry);
}


static void pheromone_reinit (double trail_value)
{
    /* printf("INIT TRAILS!!!\n"); */
    init_pheromone_trails( trail_value );
    restart_best_ant->tour_length = INFTY;
    restart_iteration = iteration;
    restart_time = elapsed_time( VIRTUAL );
}


static bool restart_condition(void)
{
    switch (flag_restart) {
      case RESTART_NEVER:
          return false;
      case RESTART_ALWAYS:
          return true;
      case RESTART_BRANCH_FACTOR:
          return branching_factor < restart_branch_factor;
      case RESTART_AVG_DISTANCE:
      {
          int iteration_worst_ant = find_worst();
          long int distance_best_worst = distance_between_ants(best_so_far_ant, &ant[iteration_worst_ant]);
          vector_long_fprint (stderr, best_so_far_ant->tour, n);
          fprintf (stderr, "\n");
          vector_long_fprint (stderr, ant[iteration_worst_ant].tour, n);
          fprintf (stderr, "\n");
          /* trace_print("distance_best_worst %ld (tour length worst %ld) < %ld ?\n", */
          /*             distance_best_worst, ant[iteration_worst_ant].tour_length, */
          /*             (long int) (restart_avg_distance * (double) n)); */
          return distance_best_worst < (long int) (restart_avg_distance * (double) n);
      }
      default:
          printf("restart_condition: invalid flag_restart value!\n");
          abort();
    }
}

static void search_control_and_statistics( void )
/*
      FUNCTION:       occasionally compute some statistics and check whether algorithm
                      has converged
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min
                      and trail_max used by MMAS may be updated
*/
{
    /* trace_print("SEARCH CONTROL AND STATISTICS\n"); */

    bool reinit_done = false;

    if (!(iteration % restart_freq)) {
	population_statistics();
	branching_factor = node_branching(lambda);
	/* trace_print("best so far %ld, iteration: %ld, time %.2f, b_fac %.5f\n", */
                    /* best_so_far_ant->tour_length, iteration, elapsed_time(VIRTUAL), */
                    /* branching_factor); */

        /* MAX-MIN Ant System was the first ACO algorithm to use pheromone
           trail re-initialisation as implemented here.  Other ACO algorithms
           may also profit from this mechanism.
        */
	if (iteration - restart_found_best > min_iters_after_restart_best
            && restart_condition())
        {
            /* trace_print("pheromone reinit: " */
            /*             "iteration (%ld) - restart_found_best (%ld), " */
            /*             "branching factor = %g\n", */
            /*             iteration, restart_found_best, */
            /*             branching_factor); */
            pheromone_reinit (mmas_flag ? trail_max : trail_0);
            reinit_done = true;
        }

	/* printf("try %li, iteration %li, b-fac %f \n\n", */
               /* n_try, iteration, branching_factor); */
    }
    if (bwas_flag && !reinit_done)
        bwas_pheromone_mutation();
}

static double update_trail_max (void)
{
    return 1. / ( rho * best_so_far_ant->tour_length );
}

static void update_statistics( void )
/*
      FUNCTION:       manage some statistical information about the trial, especially
                      if a new best solution (best-so-far or restart-best) is found and
                      adjust some parameters if a new best solution is found
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min
                      and trail_max used by MMAS may be updated
*/
{
    long int iteration_best_ant; /* MANUEL: make it global. */

    iteration_best_ant = find_best();

    if ( ant[iteration_best_ant].tour_length < best_so_far_ant->tour_length ) {

	time_used = elapsed_time( VIRTUAL ); /* best sol found after time_used */
	copy_from_to( &ant[iteration_best_ant], best_so_far_ant );
	copy_from_to( &ant[iteration_best_ant], restart_best_ant );

	found_best = iteration;
	restart_found_best = iteration;
	found_branching = node_branching(lambda);
	branching_factor = found_branching;

        trail_max = update_trail_max ();
        trail_min = update_trail_min (trail_max);
        if ( mmas_flag ) {
            trail_0 = trail_max;
        }
        /* trace_print("trail_min %.15lf trail_max %.15lf\n", trail_min, trail_max); */

        /* if (time_used < 1.) time_used = 1.0; */
        write_report();
    }
    if ( ant[iteration_best_ant].tour_length < restart_best_ant->tour_length ) {
	copy_from_to( &ant[iteration_best_ant], restart_best_ant );
	restart_found_best = iteration;
	/* printf("restart best: %ld, restart_found_best %ld, time %.2f\n", */
               /* restart_best_ant->tour_length, restart_found_best, elapsed_time ( VIRTUAL )); */
    }
}

static void as_update( void )
/*
      FUNCTION:       manage global pheromone deposit for Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    trace_print("Ant System pheromone deposit\n");

    for ( k = 0 ; k < n_ants ; k++ )
	global_update_pheromone( &ant[k] );
}



static void eas_update( void )
/*
      FUNCTION:       manage global pheromone deposit for Elitist Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants plus elitist ant deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    trace_print("Elitist Ant System pheromone deposit\n");

    for ( k = 0 ; k < n_ants ; k++ )
	global_update_pheromone( &ant[k] );
    global_update_pheromone_weighted( best_so_far_ant, elitist_ants );
}



static void ras_update( void )
/*
      FUNCTION:       manage global pheromone deposit for Rank-based Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the ras_ranks-1 best ants plus the best-so-far ant deposit pheromone
                      on matrix "pheromone"
      COMMENTS:       this procedure could be implemented slightly faster, but it is
                      anyway not critical w.r.t. CPU time given that ras_ranks is
		      typically very small.
*/
{
    long int i, k, b, target;
    long int *help_b;

    trace_print("Rank-based Ant System pheromone deposit\n");

    help_b = malloc( n_ants  * sizeof(long int) );
    for ( k = 0 ; k < n_ants ; k++ )
	help_b[k] = ant[k].tour_length;

    for ( i = 0 ; i < ras_ranks-1 ; i++ ) {
	b = LONG_MAX; target = -1;
	for ( k = 0 ; k < n_ants ; k++ ) {
	    if ( help_b[k] < b ) {
		b = help_b[k]; target = k;
	    }
	}
	if (target == -1) {
            break;
	}
        help_b[target] = LONG_MAX;
        global_update_pheromone_weighted( &ant[target], ras_ranks-i-1 );
    }
    /* FIXME: use update_ant = update_schedule(); */
    global_update_pheromone_weighted( best_so_far_ant, ras_ranks );
    free ( help_b );
}

static void bwas_update( void )
/*
      FUNCTION:       manage global pheromone deposit for Best-Worst Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  either the iteration-best or the best-so-far ant deposit pheromone
                      on matrix "pheromone"
*/
{
    int   iteration_worst_ant;

    trace_print("Best-worst Ant System pheromone deposit\n");
    /* FIXME: Use     update_ant = update_schedule(); */
    global_update_pheromone( best_so_far_ant );

    iteration_worst_ant = find_worst();
    bwas_worst_ant_update( &ant[iteration_worst_ant], best_so_far_ant );
    /* bwas_pheromone_mutation () is called later when deciding whether to do
       pheromone_reinit. */
}



static void acs_global_update( void )
/*
      FUNCTION:       manage global pheromone deposit for Ant Colony System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the best-so-far ant deposits pheromone on matrix "pheromone"
      COMMENTS:       global pheromone deposit in ACS is done per default using
                      the best-so-far ant; Gambardella & Dorigo examined also iteration-best
		      update (see their IEEE Trans. on Evolutionary Computation article),
		      but did not use it for the published computational results.
*/
{
    trace_print("Ant colony System global pheromone deposit\n");
    /* FIXME: Use  update_ant = update_schedule(); */
    global_acs_pheromone_update( best_so_far_ant );
}



static void pheromone_trail_update( void )
/*
      FUNCTION:       manage global pheromone trail update for the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited
                      according to the rules defined by the various ACO algorithms.
*/
{
    /* Simulate the pheromone evaporation of all pheromones; this is not necessary
       for ACS (see also ACO Book) */
    if ( as_flag || eas_flag || ras_flag || bwas_flag || mmas_flag ) {
	if ( nn_list && ls_flag) {
	    /* evaporate only pheromones on arcs of candidate list to make the
	       pheromone evaporation faster for being able to tackle large TSP
	       instances. For MMAS additionally check lower pheromone trail limits.
	    */
	    if (flag_ph_limits )
		mmas_evaporation_nn_list();
	    else
		evaporation_nn_list();
	} else {
	    /* if no local search is used, evaporate all pheromone trails */
	    evaporation();
	}
    }

    /* Next, apply the pheromone deposit for the various ACO algorithms */
    if ( as_flag )
	as_update();
    else if ( eas_flag )
	eas_update();
    else if ( ras_flag )
	ras_update();
    else if ( mmas_flag )
	mmas_update();
    else if ( bwas_flag )
	bwas_update();
    else if ( acs_flag )
	acs_global_update();

  /* check pheromone trail limits (default for MMAS); not necessary if nn_list
     and local search are used, because in that case lower pheromone trail
     limits are checked in procedure mmas_evaporation_nn_list above */
    if ( as_flag || eas_flag || ras_flag || mmas_flag || bwas_flag ) {
        if (flag_ph_limits && !(nn_list && ls_flag))
            check_pheromone_trail_limits();

  /* Compute combined information pheromone times heuristic info after
     the pheromone update for all ACO algorithms except ACS; in the ACS case
     this is already done in the pheromone update procedures of ACS */
	if ( nn_list && ls_flag ) {
	    compute_nn_list_total_information();
	} else {
	    compute_total_information();
	}
    }
}




/* --- main program ------------------------------------------------------ */

void run_aco(int argc, char *argv[], long *cost, double *best_time, double
		*total_time)
/*
      FUNCTION:       main control for running the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  none
      COMMENTS:       this function controls the run of "max_tries" independent trials

*/
{
	start_timers();

	init_program(argc, argv);

	nn_list = compute_nn_lists(&instance);
	pheromone = generate_double_matrix( n, n );
	total = generate_double_matrix( n, n );

	time_used = elapsed_time( VIRTUAL );
	/* printf("Initialization took %.10f seconds\n", time_used); */
	n_try = 0;


	init_try(n_try);
	/* write_report(); /1* we print the initial heuristic solution.  *1/ */

	while (!termination_condition()) {

		construct_solutions();

		if (ls_flag)
			apply_local_search();

		update_statistics();

		pheromone_trail_update();

		iteration++;

		search_control_and_statistics();

		adapt_parameters_next_iteration();
	}
	exit_try(n_try);
	exit_program();
	*cost = best_in_try[n_try];
	*best_time = time_best_found[n_try];
	*total_time = time_total_run[n_try];

	free_instance( &instance );
	free( pheromone );
	free( total );
	free( best_in_try );
	free( best_found_at );
	free( time_best_found );
	free( time_total_run );
	free_ants();
}
