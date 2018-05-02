/*

       AAAA    CCCC   OOOO    OOOO     AAAA   PPPPP
      AA  AA  CC     OO  OO  OO  OO   AA  AA  PP  PP
      AAAAAA  CC     OO  OO  OO  OO   AAAAAA  PPPPP
      AA  AA  CC     OO  OO  OO  OO   AA  AA  PP
      AA  AA   CCCC   OOOO    OOOO O  AA  AA  PP

######################################################
##########    ACO algorithms for the QAP    ##########
######################################################

      Version: 1.0
      File:    ants.c
      Author:  Thomas Stuetzle & Manuel Lopez-Ibanez
      Purpose: implementation of procedures for ants' behaviour
      Check:   README.txt and legal.txt
      Copyright (C) 2014  Thomas Stuetzle & Manuel Lopez-Ibanez
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
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>

#include "InOut.h"
#include "QAP.h"
#include "ants.h"
#include "ls.h"
#include "utilities.h"
#include "timer.h"


/************************************************************
 ************************************************************
Procedures for pheromone manipulation
 ************************************************************
 ************************************************************/

void evaporation( void )
/*
      FUNCTION:      implements the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
*/
{
    long int    i, j;

    /* trace_print("pheromone evaporation\n"); */

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < n ; j++ ) {
	    pheromone[i][j] = (1. - rho) * pheromone[i][j];
            /* We need this because we update with extremely small values, thus
               pheromone may become extremely small and underflow. If that
               happens, we may get invalid solutions.  */
            if (pheromone[i][j] < trail_absolute_min)
                pheromone[i][j] = trail_absolute_min;
	}
    }
}


void evaporation_nn_list( void )
/*
      FUNCTION:      simulation of the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
      REMARKS:       if local search is used, this evaporation procedure
                     only considers links between a city and those cities
		     of its candidate list
*/
{
    printf("pheromone evaporation nn_list NOT implemented for QAP\n");
    abort();
}


void update_pheromone (const long int *s, double d_tau)
{
    int i;
    for ( i = 0 ; i < n ; i++ ) {
	pheromone[s[i]][i] += d_tau;
    }
}



/* THOMAS: CHANGED; COULD BE REMOVED BUT THEN MORE EFFORT TO REMOVE total AND RENAME IN CODE */
void compute_total_information( void )
/*
      FUNCTION: calculates heuristic info times pheromone for each arc
      INPUT:    none
      OUTPUT:   none
*/
{
    long int     i, j;

    /* trace_print("compute total information\n"); */

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < n ; j++ ) {
            total[i][j] = powx(pheromone[i][j], alpha);  /* * pow(HEURISTIC(i,j),beta); THOMAS: NO HEURISTIC IN CODE */
	}
    }
}

void compute_nn_list_total_information( void )
/*
      FUNCTION: calculates heuristic info times pheromone for arcs in nn_list
      INPUT:    none
      OUTPUT:   none
*/
{
    printf("compute total information nn_list NOT implemented for QAP\n");
    abort();
}



/****************************************************************
 ****************************************************************
Procedures implementing solution construction and related things
 ****************************************************************
 ****************************************************************/



void choose_and_move_to_next(  ant_struct *a, long int item /* long int individual, long int item */)
/*
      FUNCTION:      choose a location where to put the item
      INPUT:         number indicating the ant and the item to be placed
      OUTPUT:        none
      (SIDE)EFFECTS: places the item on a location loc
*/
{
    long int  i;
    double     rnd, partial_sum = 0., sum_prob = 0.;
    long int  loc;
    double  * prob = prob_of_selection;

    if ( (q0 > 0.0) && (ran01( &seed ) < q0)  ) {
	/* with a probability q0 make the best possible choice
	   according to pheromone trails and heuristic information */
	/* we first check whether q0 > 0.0, to avoid the very common case
	   of q0 = 0.0 to have to compute a random number, which is
	   expensive computationally */
	choose_best_next(a, item);
	return;
    }

    DEBUG( assert ( item >= 0 && item < n ) );
    for ( i = 0 ; i < n ; i++) {
	if ( a->visited[i] )
	    prob[i] = 0.0;
	else {
	    prob[i] = total[item][i];
	    sum_prob += prob[i];
	}
    }

    /* now select to which location assign next an object */
    rnd = ran01(&seed) * sum_prob;
    i = 0;
    partial_sum = prob[i];
    while ((partial_sum <= rnd) && ( i < n-1 )) {
	i++;
	partial_sum += prob[i];
    }
    loc = i;
    DEBUG( assert ( loc >= 0 && loc < n ) );
    a->visited[loc] = TRUE; /* colony[individual].occupied[loc] = TRUE; */
    a->tour[loc] = item; /* colony[individual].s[loc] = item; */
}



void choose_best_next( ant_struct *a, long int item )
/*
      FUNCTION:      chooses for an ant as the next city the one with
                     maximal value of heuristic information times pheromone
      INPUT:         pointer to ant and the construction step
      OUTPUT:        none
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{
    long int   i, loc;
    double   value_best;

    DEBUG( assert ( item >= 0 && item < n ));
    loc = -1;
    value_best = -1.;             /* values in total matrix are always >= 0.0 */
    for ( i = 0 ; i < n ; i++ ) {
	if ( a->visited[i] )
	    ; /* city already visited, do nothing */
	else {
	    if ( total[item][i] > value_best ) {
		loc = i;
		value_best = total[item][i];
	    }
	}
    }
    DEBUG( assert ( 0 <= loc && loc < n);
           assert ( value_best > 0.0 );
           assert ( a->visited[loc] == FALSE ));
    a->tour[loc] = item;
    a->visited[loc] = TRUE;
}




/****************************************************************
 ****************************************************************
Procedures specific to MAX-MIN Ant System
 ****************************************************************
****************************************************************/



void mmas_evaporation_nn_list( void )
/*
      FUNCTION:      simulation of the pheromone trail evaporation for MMAS
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
      REMARKS:       if local search is used, this evaporation procedure
                     only considers links between a city and those cities
		     of its candidate list
*/
{
    printf("mmas specific evaporation on nn_lists NOT implemented for QAP\n");
    abort();
}

double update_trail_min (double tau_max)
{
    double tau_min;
/*
    if (ls_flag ) {
        tau_min = tau_max / ( 2. * n );
    } else
*/ {
        /* More or less equal to pow(p_dec, 1.0 / n) */
        double p_x = exp(log(p_dec)/n);
        tau_min = (1. - p_x) / (p_x * (double)((nn_ants + 1) / 2));
        tau_min = tau_max * tau_min;
    }
    return tau_min;
}

void check_pheromone_trail_limits(void )
/*
      FUNCTION:      only for MMAS without local search:
                     keeps pheromone trails inside trail limits
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are forced to interval [trail_min,trail_max]
*/
{
    int i, j;
    const double tau_min = trail_min;
    const double tau_max = trail_max;

    /* trace_print("check pheromone trail limits\n"); */

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < n ; j++ ) {
	    if ( pheromone[i][j] < tau_min )
		pheromone[i][j] = tau_min;
        }
    }
    if (iteration == 1) {
        for ( i = 0 ; i < n ; i++ ) {
            for ( j = 0 ; j < n ; j++ ) {
                if ( pheromone[i][j] > tau_max ) {
                    pheromone[i][j] = tau_max;
                }
            }
        }
    }
}



/****************************************************************
 ****************************************************************
Procedures specific to Ant Colony System
 ****************************************************************
****************************************************************/

void global_acs_pheromone_update( ant_struct *a )
/*
      FUNCTION:      reinforces the edges used in ant's solution as in ACS
      INPUT:         pointer to ant that updates the pheromone trail
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{
    long int i;
    double   d_tau;

    trace_print("acs specific: global pheromone update\n");

    d_tau = 1.0 / (double) a->tour_length;

    for ( i = 0 ; i < n ; i++ ) {
        int j = a->tour[i];
	pheromone[j][i] = (1. - rho) * pheromone[j][i] + rho * d_tau;
	total[a->tour[i]][i] = powx(pheromone[a->tour[i]][i], alpha);
    }
}



void local_acs_pheromone_update( ant_struct *a, long int j )
/*
      FUNCTION:      removes some pheromone on edge just passed by the ant
      INPUT:         pointer to ant and number of constr. phase
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are decreased
*/
{
    int h = a->tour[j];
    pheromone[h][j] = (1. - xi) * pheromone[h][j] + xi * trail_0;
    total[h][j] = powx(pheromone[h][j], alpha);
}



/****************************************************************
 ****************************************************************
Procedures specific to Best-Worst Ant System
 ****************************************************************
****************************************************************/

void bwas_worst_ant_update( ant_struct *worst, ant_struct *best)
/*
      FUNCTION:      uses additional evaporation on the arcs of iteration worst
                     ant that are not shared with the global best ant
      INPUT:         pointer to the worst (a1) and the best (a2) ant
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones on some arcs undergo additional evaporation
*/
{
    long int    i;
    trace_print("bwas specific: best-worst pheromone update\n");

    for ( i = 0 ; i < n ; i++ ) {
        /* edge (i,h) does not occur in ant best */
	if (worst->tour[i] != best->tour[i])
            pheromone[worst->tour[i]][i] *= (1 - rho);
    }
}


void bwas_pheromone_mutation( void )
/*
      FUNCTION: implements the pheromone mutation in Best-Worst Ant System
      INPUT:    none
      OUTPUT:   none
*/
{
    long int     i, j, k;
    long int     num_mutations;
    double       avg_trail = 0.0, mutation_strength = 0.0, mutation_rate = 0.3;

    trace_print("bwas specific: pheromone mutation\n");

    /* compute average pheromone trail on edges of global best solution */
    for ( i = 0 ; i < n ; i++ ) {
        /* MANUEL: This could be
           sol_comp c = solution_component(best_so_far_ant, i);
           avg_trail += get_pheromone(c);
         */
        int j = best_so_far_ant->tour[i];
        int k = i;
	avg_trail += pheromone[j][k];
    }
    avg_trail /= (double) n;

    /* determine mutation strength of pheromone matrix */
    /* FIXME: we add a small value to the denominator to avoid any
       potential division by zero. This may not be fully correct
       according to the original BWAS. */
    /* FIXME: Remove this error. */
    if ( max_time > 0.1 )
	mutation_strength = 4. * avg_trail * (elapsed_time(VIRTUAL) - restart_time)
            / (max_time - restart_time + 0.0001);
    else if ( max_tours > 100 )
	mutation_strength = 4. * avg_trail * (iteration - restart_iteration)
            / (max_tours - restart_iteration + 1);
    else
	printf("bwas_pheromone_mutation: apparently no termination condition applied!!\n");

    /* finally use fast version of matrix mutation */
    num_mutations = mutation_rate * (double) n;

    if ( restart_iteration < 2 )
	num_mutations = 0;

    for ( i = 0 ; i < num_mutations ; i++ ) {
	j =   (long int) (ran01( &seed ) * (double) n);
	k =   (long int) (ran01( &seed ) * (double) n);
	if ( ran01( &seed ) < 0.5 ) {
	    pheromone[j][k] += mutation_strength;
	}
	else {
	    pheromone[j][k] -= mutation_strength;
	    if ( pheromone[j][k] < trail_absolute_min ) {
		pheromone[j][k] = trail_absolute_min;
	    }
	}
    }
}



/**************************************************************************
 **************************************************************************
Procedures specific to the ant's tour manipulation other than construction
***************************************************************************
 **************************************************************************/



void copy_from_to(const ant_struct *a1, ant_struct *a2)
{
/*
      FUNCTION:       copy solution from ant a1 into ant a2
      INPUT:          pointers to the two ants a1 and a2
      OUTPUT:         none
      (SIDE)EFFECTS:  a2 is copy of a1
*/
    int   i;

    a2->tour_length = a1->tour_length;
    for ( i = 0 ; i < n ; i++ ) {
	a2->tour[i] = a1->tour[i];
    }
    a2->tour[n] = a2->tour[0];
}



long int nn_tour(ant_struct * nn_ant)
/*
      FUNCTION:       generate random solution and compute tour length
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  needs ant colony and one statistic ants
*/

{
    long int i;
    long int *objects;

    ant_empty_memory( nn_ant );
    objects = generate_random_permutation( n );

    for ( i = 0 ; i < n ; i++ ) {
        nn_ant->tour[i] = objects[i];
    }
    n_tours++;
    nn_ant->tour_length = compute_tour_length( nn_ant->tour );

    free (objects);
    return nn_ant->tour_length;
}


long int distance_between_ants( ant_struct *a1, ant_struct *a2)
/*
      FUNCTION: compute the distance between the tours of ant a1 and a2
      INPUT:    pointers to the two ants a1 and a2
      OUTPUT:   distance between ant a1 and a2
*/
{
    return num_different_positions (a1->tour, a2->tour, n);
}


