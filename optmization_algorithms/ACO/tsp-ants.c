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
      File:    ants.c
      Author:  Thomas Stuetzle
      Purpose: implementation of procedures for ants' behaviour
      Check:   README.txt and legal.txt
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
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>

#include "InOut.h"
#include "TSP.h"
#include "ants.h"
#include "tsp-ls.h"
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

    trace_print("pheromone evaporation\n");

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j <= i ; j++ ) {
	    pheromone[i][j] = (1 - rho) * pheromone[i][j];
	    pheromone[j][i] = pheromone[i][j];
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
    long int    i, j, help_city;

    trace_print("pheromone evaporation nn_list\n");

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < nn_ants ; j++ ) {
	    help_city = instance.nn_list[i][j];
	    pheromone[i][help_city] = (1 - rho) * pheromone[i][help_city];
	}
    }
}


void update_pheromone (const long int *s, double d_tau)
{
    int i;
    for ( i = 0 ; i < n ; i++ ) {
	int j = s[i];
	int h = s[i+1];
	pheromone[j][h] += d_tau;
	pheromone[h][j] = pheromone[j][h];
    }
}

void compute_total_information( void )
/*    
      FUNCTION: calculates heuristic info times pheromone for each arc
      INPUT:    none  
      OUTPUT:   none
*/
{
    long int     i, j;

    trace_print("compute total information\n");

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < i ; j++ ) {
	    total[i][j] = powx(pheromone[i][j],alpha) * powx(HEURISTIC(i,j),beta);
	    total[j][i] = total[i][j];
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
    long int    i, j, h;

    trace_print("compute total information nn_list\n");

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < nn_ants ; j++ ) {
	    h = instance.nn_list[i][j];
	    if ( pheromone[i][h] < pheromone[h][i] )
		/* force pheromone trails to be symmetric as much as possible */
		pheromone[h][i] = pheromone[i][h];
	    total[i][h] = powx(pheromone[i][h], alpha) * powx(HEURISTIC(i,h),beta);
	    total[h][i] = total[i][h];
	}
    }
}



/****************************************************************
 ****************************************************************
Procedures implementing solution construction and related things
 ****************************************************************
 ****************************************************************/



void place_ant( ant_struct *a , long int step )
/*    
      FUNCTION:      place an ant on a randomly chosen initial city
      INPUT:         pointer to ant and the number of construction steps 
      OUTPUT:        none
      (SIDE)EFFECT:  ant is put on the chosen city
*/
{
    long int     rnd;

    rnd = (long int) (ran01( &seed ) * (double) n); /* random number between 0 .. n-1 */
    a->tour[step] = rnd; 
    a->visited[rnd] = TRUE;
}



void choose_best_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next city the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int city, current_city, next_city;
    double   value_best;

    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ) );
    current_city = a->tour[phase-1];
    value_best = -1.;             /* values in total matrix are always >= 0.0 */    
    for ( city = 0 ; city < n ; city++ ) {
	if ( a->visited[city] ) 
	    ; /* city already visited, do nothing */
	else {
	    if ( total[current_city][city] > value_best ) {
		next_city = city;
		value_best = total[current_city][city];
	    }
	} 
    }
    DEBUG( assert ( 0 <= next_city && next_city < n);
           assert ( value_best > 0.0 );
           assert ( a->visited[next_city] == FALSE ) );
    a->tour[phase] = next_city;
    a->visited[next_city] = TRUE;
}



void neighbour_choose_best_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next city the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int i, current_city, next_city, help_city;
    double   value_best, help;
  
    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ) );
    current_city = a->tour[phase-1];
    DEBUG ( assert ( 0 <= current_city && current_city < n ) );
    value_best = -1.;             /* values in total matix are always >= 0.0 */    
    for ( i = 0 ; i < nn_ants ; i++ ) {
	help_city = instance.nn_list[current_city][i];
	if ( a->visited[help_city] ) 
	    ;   /* city already visited, do nothing */
	else {
	    help = total[current_city][help_city];
	    if ( help > value_best ) {
		value_best = help;
		next_city = help_city;
	    }
	}
    }
    if ( next_city == n )
	/* all cities in nearest neighbor list were already visited */
	choose_best_next( a, phase );
    else {
	DEBUG( assert ( 0 <= next_city && next_city < n);
               assert ( value_best > 0.0 );
               assert ( a->visited[next_city] == FALSE ));
	a->tour[phase] = next_city;
	a->visited[next_city] = TRUE;
    }
}



static void choose_closest_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      Chooses for an ant the closest city as the next one 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int city, current_city, next_city, min_distance;
  
    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ) );
    current_city = a->tour[phase-1];
    min_distance = INFTY;             /* Search shortest edge */    
    for ( city = 0 ; city < n ; city++ ) {
	if ( a->visited[city] ) 
	    ; /* city already visited */
	else {
	    if ( instance.distance[current_city][city] < min_distance) {
		next_city = city;
		min_distance = instance.distance[current_city][city];
	    }
	} 
    }
    DEBUG( assert ( 0 <= next_city && next_city < n) );
    a->tour[phase] = next_city;
    a->visited[next_city] = TRUE;
}



void choose_and_move_to_next( ant_struct *a , long int phase )
/*    
      FUNCTION:      Choose for an ant probabilistically a next city among all 
                     unvisited cities in the current city's candidate list. 
		     If this is not possible, choose the closest next
      INPUT:         pointer to ant the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int i, help; 
    long int current_city;
    double   rnd, partial_sum = 0., sum_prob = 0.0;
    /* stores the selection probabilities of the nearest neighbor cities */
    double  * prob = prob_of_selection;

    if ( (q0 > 0.0) && (ran01( &seed ) < q0)  ) {
	/* with a probability q0 make the best possible choice
	   according to pheromone trails and heuristic information */
	/* we first check whether q0 > 0.0, to avoid the very common case
	   of q0 = 0.0 to have to compute a random number, which is
	   expensive computationally */
	neighbour_choose_best_next(a, phase);
	return;
    }

    current_city = a->tour[phase-1]; /* current_city city of ant k */
    DEBUG( assert ( current_city >= 0 && current_city < n ) );
    for ( i = 0 ; i < nn_ants ; i++ ) {
        help = instance.nn_list[current_city][i];
	if ( a->visited[help] ) 
	    prob[i] = 0.0;   /* city already visited */
	else {
	    DEBUG( assert ( help >= 0 && help < n ));
	    prob[i] = total[current_city][help];
	    sum_prob += prob[i];
	} 
    }

    if (sum_prob <= 0.0) {
	/* All cities from the candidate set are tabu */
	choose_best_next( a, phase );
    }     
    else {  
	/* at least one neighbor is eligible, chose one according to the
	   selection probabilities */
	rnd = ran01( &seed ) * sum_prob;
	i = 0;
	partial_sum = prob[i];
        /* This loop always stops because prob[nn_ants] == HUGE_VAL  */
        while (partial_sum <= rnd) {
            i++;
            DEBUG(assert (i <= nn_ants));
            partial_sum += prob[i];
        }
        /* This may very rarely happen because of rounding if rnd is
           close to 1.  */
        if (i == nn_ants) {
            neighbour_choose_best_next(a, phase);
            return;
        }
	DEBUG( assert ( 0 <= i && i < nn_ants) );
	DEBUG( assert ( prob[i] >= 0.0) );
	help = instance.nn_list[current_city][i];
	DEBUG( assert ( help >= 0 && help < n ) );
	DEBUG( assert ( a->visited[help] == FALSE ) );
	a->tour[phase] = help; /* instance.nn_list[current_city][i]; */
	a->visited[help] = TRUE;
    }
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
    long int    i, j, help_city;

    trace_print("mmas specific evaporation on nn_lists\n");

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < nn_ants ; j++ ) {
	    help_city = instance.nn_list[i][j];
	    pheromone[i][help_city] = (1 - rho) * pheromone[i][help_city];
	    if ( pheromone[i][help_city] < trail_min )
		pheromone[i][help_city] = trail_min;
	}
    }
}


double update_trail_min (double tau_max)
{
    double tau_min;

    if (ls_flag ) {
        tau_min = tau_max / ( 2. * n );
    } else {
        /* More or less equal to pow(p_dec, 1.0 / n) */
        double p_x = exp(log(p_dec)/n);
        tau_min = (1. - p_x) / (p_x * (double)((nn_ants + 1) / 2));
        tau_min = tau_max * tau_min; 
    }
    return tau_min;
}

void check_pheromone_trail_limits( void )
/*    
      FUNCTION:      only for MMAS without local search: 
                     keeps pheromone trails inside trail limits
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are forced to interval [trail_min,trail_max]
*/
{ 
    long int    i, j;

    trace_print("mmas specific: check pheromone trail limits\n");

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < i ; j++ ) {
	    if ( pheromone[i][j] < trail_min ) {
		pheromone[i][j] = trail_min;
		pheromone[j][i] = trail_min;
	    } else if ( pheromone[i][j] > trail_max ) {
		pheromone[i][j] = trail_max;
		pheromone[j][i] = trail_max;
	    }
	}
    }
}



void check_nn_list_pheromone_trail_limits( void )
/*    
      FUNCTION:      only for MMAS with local search: keeps pheromone trails 
                     inside trail limits
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are forced to interval [trail_min,trail_max]
      COMMENTS:      currently not used since check for trail_min is integrated
                     mmas_evaporation_nn_list and typically check for trail_max 
		     is not done (see FGCS paper or ACO book for explanation)
*/
{ 
    long int    i, j, help_city;

    trace_print("mmas specific: check pheromone trail limits nn_list\n");

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < nn_ants ; j++ ) {
	    help_city = instance.nn_list[i][j];
	    if ( pheromone[i][help_city] < trail_min )
		pheromone[i][help_city] = trail_min;
	    if ( pheromone[i][help_city] > trail_max )
		pheromone[i][help_city] = trail_max;
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
    long int i, j, h;
    double   d_tau;

    trace_print("acs specific: global pheromone update\n");

    d_tau = 1.0 / (double) a->tour_length;

    for ( i = 0 ; i < n ; i++ ) {
	j = a->tour[i];
	h = a->tour[i+1];

	pheromone[j][h] = (1. - rho) * pheromone[j][h] + rho * d_tau;
	pheromone[h][j] = pheromone[j][h];

	total[h][j] = powx(pheromone[h][j], alpha) * powx(HEURISTIC(h,j),beta);
	total[j][h] = total[h][j];
    }
}



void local_acs_pheromone_update( ant_struct *a, long int phase )
/*    
      FUNCTION:      removes some pheromone on edge just passed by the ant
      INPUT:         pointer to ant and number of constr. phase
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    long int  h, j;
    
    DEBUG ( assert ( phase > 0 && phase <= n ) );
    j = a->tour[phase];
    h = a->tour[phase-1];
    DEBUG ( assert ( 0 <= j && j < n ) );
    DEBUG ( assert ( 0 <= h && h < n ) );
    pheromone[h][j] = (1. - xi) * pheromone[h][j] + xi * trail_0;
    pheromone[j][h] = pheromone[h][j];
    total[h][j] = powx(pheromone[h][j], alpha) * powx(HEURISTIC(h,j),beta);
    total[j][h] = total[h][j];
}



/****************************************************************
 ****************************************************************
Procedures specific to Best-Worst Ant System
 ****************************************************************
****************************************************************/



void bwas_worst_ant_update( ant_struct *a1, ant_struct *a2)
/*    
      FUNCTION:      uses additional evaporation on the arcs of iteration worst
                     ant that are not shared with the global best ant
      INPUT:         pointer to the worst (a1) and the best (a2) ant
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones on some arcs undergo additional evaporation
*/
{  
    long int    i, j, h, pos, pred;
    long int    *pos2;        /* positions of cities in tour of ant a2 */ 

    trace_print("bwas specific: best-worst pheromone update\n");

    pos2 = malloc(n * sizeof(long int));
    for ( i = 0 ; i < n ; i++ ) {
	pos2[a2->tour[i]] = i;
    }
 
    for ( i = 0 ; i < n ; i++ ) {
	j = a1->tour[i];
	h = a1->tour[i+1];
	pos = pos2[j];
	if (pos - 1 < 0)
	    pred = n - 1;
	else 
	    pred = pos - 1;
	if (a2->tour[pos+1] == h || a2->tour[pred] == h)
	    ; /* do nothing, edge is common with a2 (best solution found so far) */
	else {   /* edge (j,h) does not occur in ant a2 */       
	    pheromone[j][h] *= (1. - rho);
	    pheromone[h][j] = pheromone[j][h];
	}
    }
    free ( pos2 );
}



void bwas_pheromone_mutation( void )
/*    
      FUNCTION: implements the pheromone mutation in Best-Worst Ant System
      INPUT:    none  
      OUTPUT:   none
*/
{
    int     i, j, k;
    int     num_mutations;
    double  avg_trail = 0.0, mutation_strength = 0.0, mutation_rate = 0.3;

    trace_print("bwas specific: pheromone mutation\n");

    /* compute average pheromone trail on edges of global best solution */
    for ( i = 0 ; i < n ; i++ ) {
        int j = best_so_far_ant->tour[i];
        int k = best_so_far_ant->tour[i+1];
	avg_trail +=  pheromone[j][k];
    }
    avg_trail /= (double) n;
  
    /* determine mutation strength of pheromone matrix */ 
    /* FIXME: we add a small value to the denominator to avoid any
       potential division by zero. This may not be fully correct
       according to the original BWAS. */
    if ( max_time > 0.1 )
	mutation_strength = 4. * avg_trail * (elapsed_time(VIRTUAL) - restart_time) / (max_time - restart_time + 0.0001);
    else if ( max_tours > 100 )
	mutation_strength = 4. * avg_trail * (iteration - restart_iteration) 
            / (max_tours - restart_iteration + 1);
    else
	printf("bwas_pheromone_mutation: apparently no termination condition applied!!\n");

    /* finally use fast version of matrix mutation */
    num_mutations = mutation_rate * (double) n;
    /* / 2 because of adjustment for symmetry of pheromone trails */
    num_mutations /= 2;   
 
    if ( restart_iteration < 2 )
	num_mutations = 0; 

    for ( i = 0 ; i < num_mutations ; i++ ) {
	j =   (long int) (ran01( &seed ) * (double) n);
	k =   (long int) (ran01( &seed ) * (double) n);
	if ( ran01( &seed ) < 0.5 ) {
	    pheromone[j][k] += mutation_strength;
	    pheromone[k][j] = pheromone[j][k];
	}
	else {
	    pheromone[j][k] -= mutation_strength;
	    if ( pheromone[j][k] < trail_absolute_min ) {
		pheromone[j][k] = trail_absolute_min;
	    }
	    pheromone[k][j] = pheromone[j][k]; 
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
      FUNCTION:       generate some nearest neighbor tour and compute tour length
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  needs ant colony and one statistic ants
*/
{
    long int phase;

    ant_empty_memory(nn_ant);

    phase = 0; /* counter of the construction steps */
    place_ant( nn_ant, phase);

    while ( phase < n-1 ) {
	phase++;
	choose_closest_next( nn_ant,phase);
    }
    phase = n;
    nn_ant->tour[n] = nn_ant->tour[0];
    if ( ls_flag ) {
	two_opt_first( nn_ant->tour );
    }
    n_tours++;
    nn_ant->tour_length = compute_tour_length( nn_ant->tour );

    return nn_ant->tour_length;
}


long int distance_between_ants( ant_struct *a1, ant_struct *a2)
/*    
      FUNCTION: compute the distance between the tours of ant a1 and a2
      INPUT:    pointers to the two ants a1 and a2
      OUTPUT:   distance between ant a1 and a2
*/
{  
    return num_different_edges (a1->tour, a2->tour, n);
}
