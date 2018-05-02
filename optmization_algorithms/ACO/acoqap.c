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
#include "QAP.h"
#include "timer.h"
#include "aco.h"
#include "qap-ls.h"
#include "adaptation.h"

const char const * PROG_ID_STR = "ACO algorithms for the QAP";

void problem_set_default_ls_parameters(void)
{
    dlb_flag       = FALSE;  /* don't apply don't look bits in local search */
    nn_ls          = 1;    /* unused in the QAP */
    nn_ants        = 0;      /* no candidate list */
}

/* These override any algorithm-specific settings. */
void problem_set_default_parameters(void)
{
    ls_flag        = LS_best_2_opt;
    nn_ls          = 1;    /* unused in the QAP */
    nn_ants        = 0;    /* no candidate list */

    opt_n_ants     = 5;
    opt_rho        = 0.2; /* Default in the published paper on MMAS for the
                             QAP */
    opt_beta       = 1.0;

    p_dec              = 0.005;
    schedule_length    = 20;
    min_iters_after_restart_best = 5;
    restart_freq        = 1;
    restart_branch_factor = 1.1;
}


void construct_solutions( void )
/*
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution
*/
{
    long int k;        /* counter variable for ants */
    long int j;        /* counter variable */
    long int *objects;

    /* trace_print("construct solutions for all ants\n"); */

    /* Different from the TSP, here ants construct solutions sequentially.  */
    for ( k = 0 ; k < n_ants ; k++ ) {
      /* trace_print("construct solution for ant %ld\n",k); */
      /* free ants for next construction */
      ant_empty_memory( &ant[k] );
      /* construct for one ant after another (sequential) */
      /* random order of objects to be assigned to locations */
      objects = generate_random_permutation( n );
      for ( j = 0 ; j < n ; j++ ) {
          choose_and_move_to_next( &ant[k], objects[j] );
          if ( acs_flag )
	      local_acs_pheromone_update( &ant[k], j );
      }
      ant[k].tour_length = compute_tour_length( ant[k].tour );
      checkTour( ant[k].tour );
      free(objects);
      n_tours++;
    }
    /* trace_print("end construction\n"); */
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
    min = pheromone[m][0];
    max = min;
    for ( i = 1 ; i < n ; i++ ) {
      if ( pheromone[m][i] > max )
          max = pheromone[m][i];
      else if ( pheromone[m][i] < min )
          min = pheromone[m][i];
    }
    cutoff = min + l * (max - min);

    for ( i = 0 ; i < n ; i++ ) {
      if ( pheromone[m][i] > cutoff )
          num_branches[m]++;
    }
  }
  avg = 0.;
  for ( m = 0 ; m < n ; m++ ) {
    avg += num_branches[m];
  }
  free ( num_branches );
  /* Norm branching factor to minimal value 1 */
  return ( avg / (double) n );
}

/* The convergence factor gives an indication about how far the algorithm is
   from convergence.

   C. Blum and M. Dorigo. 2004. The Hyper-Cube Framework for Ant Colony
   Optimization. IEEE Transactions on Systems, Man and Cybernetics, Part B
   (Cybernetics) 34 (2): 1161-72.
*/
double compute_convergence_factor()
{
    double cf = 0.0;
    long int i,j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            cf += MAX (trail_max - pheromone[i][j],
                       pheromone[i][j] - trail_min);
        }
    }

    cf /= (n * n * (trail_max - trail_min));
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

    long int iterations_since_restart = iteration - restart_iteration;
    long int iteration_best_ant;

    /* trace_print("MAX-MIN Ant System pheromone deposit\n"); */
    /* This is what the old MMASQAP code did : */

    if ( iterations_since_restart < 5 ) {
	iteration_best_ant = find_best();
	global_update_pheromone( &ant[iteration_best_ant] );
        /* trace_print("pheromone update: iteration-best: iteration_since_restart = %ld\n", */
        /*             iterations_since_restart); */
    }
    /* ??? Thus, it doesn't make sense that schedule_length is < 5  */
    else if ( iterations_since_restart < schedule_length ) {
        if (iterations_since_restart % u_gb) {
            iteration_best_ant = find_best();
            global_update_pheromone( &ant[iteration_best_ant] );
            /* trace_print("pheromone update: iteration-best: " */
            /*             "iteration_since_restart = %ld, u_gb = %ld\n", */
            /*             iterations_since_restart, u_gb); */
        } else {
            global_update_pheromone( restart_best_ant );
            /* trace_print("pheromone update: restart-best: " */
            /*             "iteration_since_restart = %ld, u_gb = %ld\n", */
            /*             iterations_since_restart, u_gb); */
        }
    } else {
        global_update_pheromone( best_so_far_ant );
        /* trace_print("pheromone update: global-best: " */
        /*             "iteration_since_restart = %ld\n", */
        /*             iterations_since_restart); */
    }

    if (ls_flag == LS_tabu_search_short || ls_flag == LS_tabu_search_long) {
	if ( iterations_since_restart < schedule_length )
	    u_gb = 2;
	else
	    u_gb = 1;
    }
    else if ( iterations_since_restart > schedule_length )
        u_gb = 1;
    else if ( iterations_since_restart > schedule_length / 2)
        u_gb = 2;
    else
        u_gb = 3;
}

