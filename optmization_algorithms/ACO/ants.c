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

ant_struct *ant;
ant_struct *best_so_far_ant;
ant_struct *restart_best_ant;

double   **pheromone;
double   **total;

double   *prob_of_selection;

long int n_ants;      /* number of ants */
long int max_n_ants;
long int nn_ants;     /* length of nearest neighbor lists for the ants'
			 solution construction */

double rho;           /* parameter for evaporation */
double xi;           /* parameter for ACS local pheromone update */
double alpha;         /* importance of trail */
double beta;          /* importance of heuristic evaluate */
double q0;            /* probability of best choice in tour construction */


long int as_flag;     /* ant system */
long int eas_flag;    /* elitist ant system */
long int ras_flag;    /* rank-based version of ant system */
long int mmas_flag;   /* MAX-MIN ant system */
long int bwas_flag;   /* best-worst ant system */
long int acs_flag;    /* ant colony system */

long int elitist_ants;    /* additional parameter for elitist
			     ant system, no. elitist ants */

long int ras_ranks;       /* additional parameter for rank-based version
			     of ant system */

double   trail_max;       /* maximum pheromone trail in MMAS */
double   trail_min;       /* minimum pheromone trail in MMAS */
long int u_gb;            /* every u_gb iterations update with best-so-far ant */
double   trail_0;         /* initial pheromone level in ACS and BWAS */


static void calloc_ant (ant_struct * a)
{
    a->tour        = calloc(n+1, sizeof(long int));
    a->visited     = calloc(n, sizeof(bool));
}

static void free_ant (ant_struct * a)
{
    free(a->tour);
    free(a->visited);
}

void allocate_ants ( void )
/*
      FUNCTION:       allocate the memory for the ant colony, the best-so-far and
                      the iteration best ant
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  allocation of memory for the ant colony and two ants that
                      store intermediate tours

*/
{
    long int i;

    if ((ant = malloc(sizeof(ant_struct) * max_n_ants)) == NULL) {
	printf("Out of memory, exit.");
	exit(1);
    }
    for ( i = 0 ; i < n_ants ; i++ ) {
        calloc_ant (ant + i);
    }

    if((best_so_far_ant = malloc(sizeof( ant_struct ) )) == NULL){
	printf("Out of memory, exit.");
	exit(1);
    }
    calloc_ant (best_so_far_ant);

    if((restart_best_ant = malloc(sizeof( ant_struct ) )) == NULL){
	printf("Out of memory, exit.");
	exit(1);
    }
    calloc_ant (restart_best_ant);

    if ((prob_of_selection = malloc(sizeof(double) * (nn_ants + 1))) == NULL) {
	printf("Out of memory, exit.");
	exit(1);
    }
    /* Ensures that we do not run over the last element in the random wheel.  */
    prob_of_selection[nn_ants] = HUGE_VAL;
}

void free_ants (void)
/*
      FUNCTION:       free the memory for the ant colony, the best-so-far and
                      the iteration best ant. This should do the same as
                      allocate_ants in reverse order.
      INPUT:          none
      OUTPUT:         none
*/
{
    int i;

    free( prob_of_selection );
    free_ant(restart_best_ant);
    free( restart_best_ant);
    free_ant( best_so_far_ant);
    free( best_so_far_ant);
    for ( i = 0 ; i < n_ants ; i++ ) {
	free_ant( ant + i);
    }
    free( ant );
}

long int find_best( void )
/*
      FUNCTION:       find the best ant of the current iteration
      INPUT:          none
      OUTPUT:         index of struct containing the iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    long int   min;
    long int   k, k_min;

    min = ant[0].tour_length;
    k_min = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
	if( ant[k].tour_length < min ) {
	    min = ant[k].tour_length;
	    k_min = k;
	}
    }
    return k_min;
}



long int find_worst( void )
/*
      FUNCTION:       find the worst ant of the current iteration
      INPUT:          none
      OUTPUT:         pointer to struct containing iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    long int   max;
    long int   k, k_max;

    max = ant[0].tour_length;
    k_max = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
	if( ant[k].tour_length > max ) {
	    max = ant[k].tour_length;
	    k_max = k;
	}
    }
    return k_max;
}



/************************************************************
 ************************************************************
Procedures for pheromone manipulation
 ************************************************************
 ************************************************************/



void init_pheromone_trails( double initial_trail )
/*
      FUNCTION:      initialize pheromone trails
      INPUT:         initial value of pheromone trails "initial_trail"
      OUTPUT:        none
      (SIDE)EFFECTS: pheromone matrix is reinitialized
*/
{
    long int i, j;

    trace_print(" init trails with %.15f\n",initial_trail);

    /* Initialize pheromone trails */
    for ( i = 0 ; i < n ; i++ ) {
	for ( j =0 ; j <= i ; j++ ) {
	    pheromone[i][j] = initial_trail;
	    pheromone[j][i] = initial_trail;
	    total[i][j] = initial_trail;
	    total[j][i] = initial_trail;
	}
    }
}

void global_update_pheromone( ant_struct *a )
/*
      FUNCTION:      reinforces edges used in ant k's solution
      INPUT:         pointer to ant that updates the pheromone trail
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{
    double   d_tau;

    /* trace_print("global pheromone update\n"); */

    d_tau = 1.0 / (double) a->tour_length;
    /* trace_print("deposit pheromone: quality: %ld amount pheromone %.15lf\n", a->tour_length, d_tau); */
    update_pheromone (a->tour, d_tau);
}


void global_update_pheromone_weighted( ant_struct *a, long int weight )
/*
      FUNCTION:      reinforces edges of the ant's tour with weight "weight"
      INPUT:         pointer to ant that updates pheromones and its weight
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in the ant's tour are increased
*/
{
    double        d_tau;

    /* trace_print("global pheromone update weighted\n"); */

    d_tau = (double) weight / (double) a->tour_length;
    /* trace_print("deposit pheromone: quality: %ld amount pheromone %.15lf\n", a->tour_length,d_tau); */
    update_pheromone (a->tour, d_tau);
}

/****************************************************************
 ****************************************************************
Procedures implementing solution construction and related things
 ****************************************************************
 ****************************************************************/



void ant_empty_memory( ant_struct *a )
/*
      FUNCTION:       empty the ants's memory regarding visited cities
      INPUT:          ant identifier
      OUTPUT:         none
      (SIDE)EFFECTS:  vector of visited cities is reinitialized to FALSE
*/
{
    long int   i;

    for( i = 0 ; i < n ; i++ ) {
	a->visited[i]=FALSE;
    }
}



/****************************************************************
 ****************************************************************
Procedures specific to MAX-MIN Ant System
 ****************************************************************
****************************************************************/




/****************************************************************
 ****************************************************************
Procedures specific to Ant Colony System
 ****************************************************************
****************************************************************/





/****************************************************************
 ****************************************************************
Procedures specific to Best-Worst Ant System
 ****************************************************************
****************************************************************/









/**************************************************************************
 **************************************************************************
Procedures specific to the ant's tour manipulation other than construction
***************************************************************************
 **************************************************************************/







