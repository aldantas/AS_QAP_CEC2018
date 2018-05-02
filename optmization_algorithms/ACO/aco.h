#ifndef _ACO_H_
#define _ACO_H_

/* Problem specific */
void construct_solutions( void );
void mmas_update( void );
long int compute_tour_length(const long int *t );
void run_aco(int argc, char *argv[], long *cost, double *best_time, double
		*total_time);
#endif
