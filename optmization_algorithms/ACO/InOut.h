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
      File:    InOut.h
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


#include <stdbool.h>

#define LINE_BUF_LEN     255

extern long int *best_in_try;
extern long int *best_found_at;
extern double   *time_best_found;
extern double   *time_total_run;

extern long int n_try;        /* number of runs counter */
extern long int n_tours;      /* constructed tours counter */
extern long int restart_iteration;    /* iteration counter */
extern double   restart_time; /* remember when restart was done if any */

extern long int max_tries;    /* maximum number of independent runs */
extern long int max_tours;    /* maximum number of tour constructions in one run */

extern const double lambda;       /* Parameter to determine branching factor */
extern int restart_freq; /* population statistics and pheromone reinitialization
                           conditions are only tested every reinit_freq iterations */

extern double   max_time;     /* maximal allowed run time of a run  */
extern double   time_used;    /* time used until some given event */
extern double   time_passed;  /* time passed until some moment*/
extern long int optimal;      /* optimal solution value or bound to find */

extern double mean_ants;      /* average tour length */
extern double stddev_ants;    /* stddev of tour lengths */
extern double branching_factor; /* average node branching factor when searching */
extern double found_branching;  /* branching factor when best solution is found */

extern long int found_best;         /* iteration in which best solution is found */
extern long int restart_found_best; /* iteration in which restart-best solution is found */

extern FILE *report, *comp_report, *stat_report;

extern char name_buf[LINE_BUF_LEN];
extern bool quiet_flag;

void write_report( void );

void print_default_parameters();

void set_default_as_parameters(void);
void set_default_eas_parameters(void);
void set_default_ras_parameters(void);
void set_default_mmas_parameters(void);
void set_default_bwas_parameters(void);
void set_default_acs_parameters(void);
void set_default_ls_parameters(void);
/* These are problem-specific.  */
extern void problem_set_default_ls_parameters(void);
extern void problem_set_default_parameters(void);
extern double node_branching(double l);

void set_default_parameters();

void output_solution( void );

void exit_try( long int ntry );

void exit_program( void );

void init_program( long int argc, char *argv[] );

void printDist(void);

void printHeur(void);

void printTrail(void);

void printTotal(void);

void printProbabilities(void);

void printTour(const long int *t );

void printTourFile( long int *t );

void checkTour( long int *t );

void population_statistics ( void );

void write_params ( void );
