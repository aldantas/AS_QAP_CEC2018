/*****************************************************************/
// Adapted from the BLS source code.
//
// Copyright : Una Benlic, 2012
// This code can be freely used for non-commercial purpose.
// Any use of this implementation or a modification of the code
// must acknowledge the work of U. Benlic and J.K. Hao
/****************************************************************/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <map>
#include <set>
#include <vector>
#include <ctime>
#include <cmath>

#include "common.h"

/* #define mem_size 100000 */

using namespace std;

const long infinite = 999999999;

long iteration = 0;

void transpose(int & a, int & b) {long temp = a; a = b; b = temp;}


/*--------------------------------------------------------------*/
/*       compute the cost difference if elements i and j        */
/*         are transposed in permutation (solution) p           */
/*--------------------------------------------------------------*/
long compute_delta(int n, type_matrix & a, type_matrix & b,
                   type_vector & p, int i, int j)
 {long d; int k;
  d = (a[i][i]-a[j][j])*(b[p[j]][p[j]]-b[p[i]][p[i]]) +
      (a[i][j]-a[j][i])*(b[p[j]][p[i]]-b[p[i]][p[j]]);
  for (k = 0; k < n; k = k + 1) if (k!=i && k!=j)
    d = d + (a[k][i]-a[k][j])*(b[p[k]][p[j]]-b[p[k]][p[i]]) +
            (a[i][k]-a[j][k])*(b[p[j]][p[k]]-b[p[i]][p[k]]);
  return(d);
 }

/*--------------------------------------------------------------*/
/*      Idem, but the value of delta[i][j] is supposed to       */
/*    be known before the transposition of elements r and s     */
/*--------------------------------------------------------------*/
long compute_delta_part(type_matrix & a, type_matrix & b,
                        type_vector & p, type_matrix & delta,
                        int i, int j, int r, int s)
{
	return(delta[i][j]+(a[r][i]-a[r][j]+a[s][j]-a[s][i])*
			(b[p[s]][p[i]]-b[p[s]][p[j]]+b[p[r]][p[j]]-b[p[r]][p[i]])+
			(a[i][r]-a[j][r]+a[j][s]-a[i][s])*
			(b[p[i]][p[s]]-b[p[j]][p[s]]+b[p[j]][p[r]]-b[p[i]][p[r]]) );
}
void update_matrix_of_move_cost(int i_retained, int j_retained,long n, type_matrix & delta, type_vector & p, type_matrix & a, type_matrix & b)
{
     int i, j;
     for (i = 0; i < n-1; i = i+1) for (j = i+1; j < n; j = j+1)
        if (i != i_retained && i != j_retained &&
            j != i_retained && j != j_retained)
         {delta[i][j] =
            compute_delta_part(a, b, p, delta,
                               i, j, i_retained, j_retained);}
        else
         {delta[i][j] = compute_delta(n, a, b, p, i, j);};
}

void apply_move(type_vector & p,long n, type_matrix & delta, long &
		current_cost, type_matrix a, type_matrix b,int i_retained, int
		j_retained)
{
	if(i_retained!=-1 && j_retained!=-1) // apply the selected perturbation move
	{
		swap(p[i_retained], p[j_retained]);
		current_cost = current_cost + delta[i_retained][j_retained];
		update_matrix_of_move_cost(i_retained, j_retained, n, delta, p, a, b);
	}
}

void generate_random_solution(long n, type_vector  & p)
{
	int i;
	for (i = 0; i < n; i = i+1) p[i] = i;
	for (i = 1; i < n; i = i+1) transpose(p[i], p[i + rand()%(n-i)]);
}

void allocate_memory_and_initialize(type_vector & p, long n, type_matrix &
		delta, long & current_cost, type_matrix a, type_matrix b)
{
	int i, j;
	/***************** dynamic memory allocation *******************/
	delta = new long* [n+1];
	for (i = 0; i < n; i = i+1)
	{
		delta[i] = new long[n+1];
	}
	/********** initialization of current solution value ***********/
	/**************** and matrix of cost of moves  *****************/
	current_cost = 0;
	for (i = 0; i < n; i = i + 1) {
		for (j = 0; j < n; j = j + 1) {
			current_cost = current_cost + a[i][j] * b[p[i]][p[j]];
			if (i < j) {
				delta[i][j] = compute_delta(n, a, b, p, i, j);
			}
		}
	}
}

bool best_improvement_move(type_vector &p, long n, type_matrix & delta, long &
		current_cost, type_matrix a, type_matrix b)
{
	int i, j, i_retained, j_retained;
	long min_delta = infinite;   // retained move cost
	//select the best swap move for the descent local search phase
	for (i = 0; i < n-1; i = i + 1) {
		for (j = i+1; j < n; j = j+1) {
			if(delta[i][j] < min_delta) {
				i_retained = i; j_retained = j;
				min_delta = delta[i][j];
			}
		}
	}
	// apply the above selected move if it leads to a solution better than the current one
	if((current_cost + delta[i_retained][j_retained])<current_cost) {
		apply_move(p, n, delta, current_cost, a, b, i_retained, j_retained);
		return 1;
	} else {
		return 0;
	}
}

void two_opt_best(type_vector & solution, long & cost, long n, type_matrix a,
		type_matrix b, clock_t & best_time, clock_t & total_time)
{
	type_matrix delta;
	clock_t start = clock();

	allocate_memory_and_initialize(solution, n, delta, cost, a, b);
	bool descent = best_improvement_move(solution, n, delta, cost, a, b);
	while(descent) // while an improving solution in the neighborhood exists
	{
		descent = best_improvement_move(solution, n, delta, cost, a, b);
	}
	best_time = clock() - start;
	total_time = best_time;

	// free memory
	for (int i=0; i < n; i = i+1)
	{
		delete[] delta[i];
	}
	delete[] delta;
}

int n;                    // problem size
type_matrix a, b;         // flows and distances matrices
type_vector solution;     // solution (permutation)
long cost;                // solution cost
int i, j; //variables for iterations
int main(int argc, char *argv[])
{
	readParameter(argc, argv);

        load_problem(n, a, b);

	srand(seed);
	solution = new int[n];

	clock_t best_time;
        clock_t total_time;

	generate_random_solution(n, solution);

	two_opt_best(solution, cost, n, a, b, best_time, total_time);

	for (int i =0; i < n; i++) {
		cout << solution[i] << " ";
	}
	cout << endl;
	cout << cost;

	delete[]solution;

	return 0;
}
