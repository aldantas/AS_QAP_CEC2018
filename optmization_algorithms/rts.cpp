/*****************************************************************
 Implementation of the robust taboo search of: E. Taillard
 "Robust taboo search for the quadratic assignment problem",
 Parallel Computing 17, 1991, 443-455.

 Data file format:
  n, optimum solution value
 (nxn) flow matrix,
 (nxn) distance matrix

 Copyright : E. Taillard, 1990-2004
 Standard C version with slight improvement regarding to
 1991 version: non uniform tabu duration
E. Taillard, 14.03.2006

 Compatibility: Unix and windows gcc, g++, bcc32.
 This code can be freely used for non-commercial purpose.
 Any use of this implementation or a modification of the code
 must acknowledge the work of E. Taillard

****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include "common.h"

using namespace std;

const int infinite = 999999999;
const int FALSE = 0;
const int TRUE = 1;

int opt;
double  somme_sol = 0.0;
long evaluations=0, linear_evals=0, const_evals=0;
long best_evals, total_evals;

/*********** return an integer between low and high *************/
int unif(int low, int high)
 {return low + rand() % (high + 1 - low) ;}

void transpose(int *a, int *b) {int temp = *a; *a = *b; *b = temp;}

int min(int a, int b) {if (a < b) return(a); else return(b);}

double cube(double x) {return x*x*x;}

/*--------------------------------------------------------------*/
/*       compute the cost difference if elements i and j        */
/*         are transposed in permutation (solution) p           */
/*--------------------------------------------------------------*/
int compute_delta(int n, type_matrix a, type_matrix b,
                   type_vector p, int i, int j)
 {int d; int k;
  d = (a[i][i]-a[j][j])*(b[p[j]][p[j]]-b[p[i]][p[i]]) +
      (a[i][j]-a[j][i])*(b[p[j]][p[i]]-b[p[i]][p[j]]);
  for (k = 0; k < n; k = k + 1) if (k!=i && k!=j)
    d = d + (a[k][i]-a[k][j])*(b[p[k]][p[j]]-b[p[k]][p[i]]) +
            (a[i][k]-a[j][k])*(b[p[j]][p[k]]-b[p[i]][p[k]]);
  linear_evals++;
  return(d);
 }

/*--------------------------------------------------------------*/
/*      Idem, but the value of delta[i][j] is supposed to       */
/*    be known before the transposition of elements r and s     */
/*--------------------------------------------------------------*/
int compute_delta_part(type_matrix a, type_matrix b,
                        type_vector p, type_matrix delta,
                        int i, int j, int r, int s)
{
	const_evals++;
	return(delta[i][j]+(a[r][i]-a[r][j]+a[s][j]-a[s][i])*
			(b[p[s]][p[i]]-b[p[s]][p[j]]+b[p[r]][p[j]]-b[p[r]][p[i]])+
			(a[i][r]-a[j][r]+a[j][s]-a[i][s])*
			(b[p[i]][p[s]]-b[p[j]][p[s]]+b[p[j]][p[r]]-b[p[i]][p[r]]) );
}

void tabu_search(int n,                  /* problem size */
                 type_matrix a,          /* flows matrix */
                 type_matrix b,          /* distance matrix */
                 type_vector best_sol,   /* best solution found */
                 int *best_cost,         /* cost of best solution */
                 int tabu_duration,      /* parameter 1 (< n^2/2) */
                 int aspiration,         /* parameter 2 (> n^2/2)*/
                 int max_iterations,
		 clock_t &best_time,
		 clock_t &total_time)      /* number of iterations */
{
	type_vector p;                        /* current solution */
	type_matrix delta;                    /* store move costs */
	type_matrix tabu_list;                /* tabu status */
	int current_iteration;                /* current iteration */
	int current_cost;                     /* current sol. value */
	int i, j, k, i_retained, j_retained;  /* indices */
	int min_delta;                        /* retained move cost */
	int autorized;                        /* move not tabu? */
	int aspired;                          /* move forced? */
	int already_aspired;                  /* in case many moves forced */
	int iter_without_improvement = 0;

        clock_t start = clock();

	/***************** dynamic memory allocation *******************/
	p = new int[n];

	delta = new long* [n];
	tabu_list = new long* [n];
	for (int i = 0; i < n; i++) {
		delta[i] = new long[n];
		tabu_list[i] = new long[n];
	}

	/************** current solution initialization ****************/
	for (i = 0; i < n; i = i + 1) p[i] = best_sol[i];

	/********** initialization of current solution value ***********/
	/**************** and matrix of cost of moves  *****************/
	current_cost = 0;
	for (i = 0; i < n; i = i + 1) {
		for (j = 0; j < n; j = j + 1)
		{
			current_cost = current_cost + a[i][j] * b[p[i]][p[j]];
			if (i < j) {delta[i][j] = compute_delta(n, a, b, p, i, j);};
		}
	}
	evaluations += 1;
	best_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n);
	*best_cost = current_cost;

	/****************** tabu list initialization *******************/
	for (i = 0; i < n; i = i + 1) for (j = 0; j < n; j = j+1)
		tabu_list[i][j] = -(n*i + j);

	/******************** main tabu search loop ********************/
	if (!has_zero_matrix)  {
		for (current_iteration = 1; current_iteration <= max_iterations && *best_cost > opt;
				current_iteration = current_iteration + 1)
		{/** find best move (i_retained, j_retained) **/
			i_retained = infinite;       /* in case all moves are tabu */
			j_retained = infinite;
			min_delta = infinite;
			already_aspired = FALSE;

			for (i = 0; i < n-1; i = i + 1)
				for (j = i+1; j < n; j = j+1)
				{autorized = (tabu_list[i][p[j]] < current_iteration) ||
					(tabu_list[j][p[i]] < current_iteration);

					aspired =
						(tabu_list[i][p[j]] < current_iteration-aspiration)||
						(tabu_list[j][p[i]] < current_iteration-aspiration)||
						(current_cost + delta[i][j] < *best_cost);

					if ((aspired && !already_aspired) || /* first move aspired*/
							(aspired && already_aspired &&    /* many move aspired*/
							 (delta[i][j] < min_delta)   ) || /* => take best one*/
							(!aspired && !already_aspired &&  /* no move aspired yet*/
							 (delta[i][j] < min_delta) && autorized))
					{i_retained = i; j_retained = j;
						min_delta = delta[i][j];
						if (aspired) {already_aspired = TRUE;};
					};
				};

			if (i_retained == infinite) {
				/* printf("All moves are tabu! \n"); */
			}
			else
			{/** transpose elements in pos. i_retained and j_retained **/
				transpose(&p[i_retained], &p[j_retained]);
				/* update solution value*/
				current_cost = current_cost + delta[i_retained][j_retained];
				/* forbid reverse move for a random number of iterations*/
				tabu_list[i_retained][p[j_retained]] =
					current_iteration + (int)(rand() % tabu_duration);
				tabu_list[j_retained][p[i_retained]] =
					current_iteration + (int)(rand() % tabu_duration);

				/* best solution improved ?*/
				if (current_cost < *best_cost) {
					*best_cost = current_cost;
					for (k = 0; k < n; k = k+1) best_sol[k] = p[k];
					best_time = clock();
					best_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n);
					iter_without_improvement = 0;
				} else {
					iter_without_improvement++;
				}

				/* update matrix of the move costs*/
				for (i = 0; i < n-1; i = i+1) for (j = i+1; j < n; j = j+1)
					if (i != i_retained && i != j_retained &&
							j != i_retained && j != j_retained)
					{delta[i][j] =
						compute_delta_part(a, b, p, delta,
								i, j, i_retained, j_retained);}
					else
					{delta[i][j] = compute_delta(n, a, b, p, i, j);};
			}
			if(ceil((clock() - start)/static_cast<double>(CLOCKS_PER_SEC))>=600) break;
		}
		best_time -= start;
		total_time = clock() - start;
	} else {
		best_time = clock() - start;
		total_time = best_time;
	}

	/* free memory*/
	free(p);
	for (i=0; i < n; i = i+1) free(delta[i]); free(delta);
	for (i=0; i < n; i = i+1) free(tabu_list[i]);
	free(tabu_list);
} /* tabu*/

void generate_random_solution(int n, type_vector  p)
{
	int i;
	for (i = 0; i < n;   i++) p[i] = i;
	for (i = 0; i < n-1; i++) transpose(&p[i], &p[unif(i, n-1)]);
}

int main(int argc, char *argv[])
{
	int n;                    /* problem size */
	type_matrix a, b;         /* flows and distances matrices*/
	type_vector solution;     /* solution (permutation) */
	int cost;                 /* solution cost */
	int i, j;

	readParameter(argc, argv);

	load_problem(n, a, b);

	solution = new int[n];

	clock_t best_time;
        clock_t total_time;

	srand(seed);

	if (!is_iteration_set) {
		max_iterations = 2000 * n;
	}

	generate_random_solution(n, solution);

	tabu_search(n, a, b,                     /* problem data */
			solution, &cost,              /* tabu search results */
			8*n, n*n*5,                   /* parameters */
			max_iterations,
			best_time,
			total_time);               /* number of iterations */

	total_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n);
	cout<<"Solution cost "<<cost<<" "<<(best_time)/static_cast<double>(CLOCKS_PER_SEC)<<" "<<(total_time)/static_cast<double>(CLOCKS_PER_SEC)<<endl;
	cout << best_evals << " " << total_evals << endl;

	write_results(cost, best_time, total_time, best_evals, total_evals);

	free(solution);
	for (i = 0; i < n; i = i+1) free(b[i]);
	free(b);
	for (i = 0; i < n; i = i+1) free(a[i]);
	free(a);
	fflush(stdin);
	return EXIT_SUCCESS;
}
