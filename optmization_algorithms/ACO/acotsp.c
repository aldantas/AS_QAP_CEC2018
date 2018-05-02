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
      File:    main.c
      Author:  Thomas Stuetzle
      Purpose: main routines and control for the ACO algorithms
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


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ants.h"
#include "utilities.h"
#include "InOut.h"
#include "TSP.h"
#include "timer.h"
#include "tsp-ls.h"
#include "adaptation.h"
#include "aco.h"

const char const * PROG_ID_STR = "ACO algorithms for the TSP";

void problem_set_default_ls_parameters(void)
{
    dlb_flag       = TRUE;  /* apply don't look bits in local search */
    nn_ls          = 20;    /* use fixed radius search in the 20 nearest neighbours */
    nn_ants        = 20;    /* number of nearest neighbours in tour construction */

    opt_n_ants     = 25;    /* number of ants */
    opt_beta       = 2.0;
    opt_rho        = 0.5;

    if (mmas_flag) {
        opt_n_ants = 25;
        opt_rho = 0.2;
    } else if (acs_flag) {
        opt_n_ants = 10;
        opt_rho = 0.1;
    } else if (eas_flag) {
        elitist_ants = opt_n_ants;
    }
}

/* These override any algorithm-specific settings. */
void problem_set_default_parameters(void)
{
    ls_flag        = LS_THREE_OPT_FIRST;     /* per default run 3-opt*/
    nn_ls          = 20;    /* use fixed radius search in the 20 nearest neighbours */
    nn_ants        = 20;    /* number of nearest neighbours in tour construction */

    p_dec          = 0.05;
    schedule_length    = 250;
    min_iters_after_restart_best = 250;
    restart_freq        = 100;
    restart_branch_factor = 1.00001;
}


void construct_solutions( void )
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    trace_print("construct solutions for all ants\n");

    /* Mark all cities as unvisited */
    for ( k = 0 ; k < n_ants ; k++) {
	ant_empty_memory( &ant[k] );
    }
    
    step = 0; 
    /* Place the ants on same initial city */
    for ( k = 0 ; k < n_ants ; k++ )
	place_ant( &ant[k], step);

    while ( step < n-1 ) {
	step++;
	for ( k = 0 ; k < n_ants ; k++ ) {
	    choose_and_move_to_next( &ant[k], step);
	    if ( acs_flag )
		local_acs_pheromone_update( &ant[k], step );
	}
    }

    step = n;
    for ( k = 0 ; k < n_ants ; k++ ) {
	ant[k].tour[n] = ant[k].tour[0];
	ant[k].tour_length = compute_tour_length( ant[k].tour );
	if ( acs_flag )
	    local_acs_pheromone_update( &ant[k], step );
    }
    n_tours += n_ants;
}


double node_branching(double l) 
/*    
      FUNCTION:       compute the average node lambda-branching factor 
      INPUT:          lambda value 
      OUTPUT:         average node branching factor 
      (SIDE)EFFECTS:  none
      COMMENTS:       see the ACO book for a definition of the average node 
                      lambda-branching factor 
*/
{
  long int  i, m;
  double    min, max, cutoff;
  double    avg;
  double    *num_branches;

  num_branches = calloc(n, sizeof(double));

  for ( m = 0 ; m < n ; m++ ) {
    /* determine max, min to calculate the cutoff value */
    min = pheromone[m][instance.nn_list[m][0]];
    max = pheromone[m][instance.nn_list[m][0]];
    for ( i = 1 ; i < nn_ants ; i++ ) {
      if ( pheromone[m][instance.nn_list[m][i]] > max )
	max = pheromone[m][instance.nn_list[m][i]];
      if ( pheromone[m][instance.nn_list[m][i]] < min )
	min = pheromone[m][instance.nn_list[m][i]];
    }
    cutoff = min + l * (max - min);
    
    for ( i = 0 ; i < nn_ants ; i++ ) {    
      if ( pheromone[m][instance.nn_list[m][i]] > cutoff )
	num_branches[m] += 1.;
    }
  }
  avg = 0.;
  for ( m = 0 ; m < n ; m++ ) {
    avg += num_branches[m];
  }
  free ( num_branches );
  /* Norm branching factor to minimal value 1 */
  return ( avg / (double) (n * 2)  );
}

/* The convergence factor gives an indication about how far the algorithm is
   from convergence.

   C. Blum and M. Dorigo. 2004. The Hyper-Cube Framework for Ant Colony
   Optimization. IEEE Transactions on Systems, Man and Cybernetics, Part B
   (Cybernetics) 34 (2): 1161–72.
*/
double compute_convergence_factor()
{
    double cf = 0.0;
    long int i,m;
    for (i = 0; i < n; i++) {
        for (m = 0; m < nn_ants; m++) {
            long int j = instance.nn_list[i][m];
            cf += MAX (trail_max - pheromone[i][j],
                       pheromone[i][j] - trail_min);
        }
    }
    cf /= (n * nn_ants * (trail_max - trail_min));
    cf = 2.0 * (cf - 0.5);
    return cf;
}



void mmas_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for MAX-MIN Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  either the iteration-best or the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
*/
{
    /* we use default upper pheromone trail limit for MMAS and hence we
       do not have to worry regarding keeping the upper limit */

    trace_print("MAX-MIN Ant System pheromone deposit\n");

    if ( iteration % u_gb ) {
        long int iteration_best_ant = find_best();
	global_update_pheromone( &ant[iteration_best_ant] );
    }
    /* This means that for !ls_flag, best_so_far_ant is never used. */ 
    else if ( u_gb == 1 && iteration - restart_found_best > 50)
        global_update_pheromone( best_so_far_ant );
    else 
        global_update_pheromone( restart_best_ant );

    if ( ls_flag ) {
	/* implement the schedule for u_gb as defined in the 
	   Future Generation Computer Systems article or in Stuetzle's PhD thesis.
	   This schedule is only applied if local search is used.
	*/
        long int iterations_since_restart = iteration - restart_iteration;

	if ( iterations_since_restart < schedule_length / 8 )
	    u_gb = 25;
	else if ( iterations_since_restart < schedule_length / 4 )
	    u_gb = 5;
	else if ( iterations_since_restart < schedule_length / 2 )
	    u_gb = 3;
	else if ( iterations_since_restart < schedule_length )
	    u_gb = 2;
	else 
	    u_gb = 1;
    } else
	u_gb = 25;
}

