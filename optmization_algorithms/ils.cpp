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

using namespace std;

const long infinite = 999999999;

long iteration = 0;
long evaluations=0, linear_evals=0, const_evals=0;
long best_evals, total_evals;


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
  linear_evals++;
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
	const_evals++;
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

void perturbe(type_vector & p, long n, type_matrix & delta, long & current_cost,
		type_matrix & a, type_matrix & b, long k)
{
	set <int> idx_history_set;
	vector <int> idx_history;
	int i = rand()%n;
	idx_history_set.insert(i);
	idx_history.push_back(i);
	int u, x, y;
	for (u = 1; u < k; u++) {
		int j;
		do {
			j = rand()%n;
		} while (idx_history_set.find(j) != idx_history_set.end());
		int idx = rand()%idx_history.size();
		i = idx_history[idx];
		idx_history_set.insert(j);
		idx_history.push_back(j);
		if(i > j) swap(i, j);
		apply_move(p, n, delta, current_cost, a, b, i, j);
	}
}

void allocate_memory_and_initialize(type_vector & p, long n, type_matrix &
		delta, type_matrix & old_delta, long & current_cost, type_matrix a,
		type_matrix b, type_vector & best_sol, long & best_cost)
{
	int i, j;
	/***************** dynamic memory allocation *******************/
	p = new int[n+1];
	delta = new long* [n+1];
	old_delta = new long* [n+1];
	for (i = 0; i < n; i = i+1)
	{
		delta[i] = new long[n+1];
		old_delta[i] = new long[n+1];
	}
	/************** current solution initialization ****************/
	for (i = 0; i < n; i = i + 1)
		p[i] = best_sol[i];
	/********** initialization of current solution value ***********/
	/**************** and matrix of cost of moves  *****************/
	current_cost = 0;
	for (i = 0; i < n; i = i + 1) {
		for (j = 0; j < n; j = j + 1) {
			current_cost = current_cost + a[i][j] * b[p[i]][p[j]];
			if (i < j) {
				delta[i][j] = compute_delta(n, a, b, p, i, j);
				old_delta[i][j] = delta[i][j];
			}
		}
	}
	evaluations++;
	best_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n);
	best_cost = current_cost;
}

void calculate_new_delta(type_vector & p, long n, type_matrix &
		delta, long & current_cost, type_matrix a,
		type_matrix b) {
	current_cost = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			current_cost = current_cost + a[i][j] * b[p[i]][p[j]];
			if (i < j) {
				delta[i][j] = compute_delta(n, a, b, p, i, j);
			}
		}
	}
}

bool best_improvement_move(type_vector & p, long n, type_matrix & delta, long &
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

bool best_improvement_local_search(type_vector & p, long n, type_matrix & delta,
		long & current_cost, type_matrix a, type_matrix b)
{
	bool descent = best_improvement_move(p, n, delta, current_cost, a, b);
	while(descent) // while an improving solution in the neighborhood exists
	{
		descent = best_improvement_move(p, n, delta, current_cost, a, b);
	}
}

void ILS(long n,                  // problem size
		type_matrix  a,         // flows matrix
		type_matrix  b,         // distance matrix
		type_vector & best_sol,  // best solution found
		long & best_cost,      // cost of best solution
                long num_iterations, long best_objective, clock_t & best_time,
		clock_t & total_time)

{
	type_vector p;                        // current solution
	type_matrix delta, old_delta;

	long current_cost;
	int iter_without_improvement = 0;
	long kmin = 3;
	long kmax = static_cast<long>(0.9 * n);
	long k = kmin;
	long restart_it = static_cast<long>(2.5 * kmax);

	clock_t start = clock();
	allocate_memory_and_initialize(p, n, delta, old_delta, current_cost, a, b, best_sol, best_cost);

	if (!has_zero_matrix ) {
		/******************** main ILS loop ********************/
		for (iteration = 1; iteration <= num_iterations; iteration++)
		{
			//best improvement descent procedure
			best_improvement_local_search(p, n, delta, current_cost, a, b);

			if(current_cost < best_cost) {
				best_cost = current_cost;
				copy_vector(best_sol, p, n);
				copy_matrix(old_delta, delta, n);
				best_time = clock();
				iter_without_improvement = 0;
				k = kmin;
				best_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n);
			} else {
				k++;
				if (k >= kmax) k = kmin;
				copy_vector(p, best_sol, n);
				copy_matrix(delta, old_delta, n);
				for (int m = 0; m < n; m++) p[m] = best_sol[m];
				current_cost = best_cost;
				iter_without_improvement++;
			}
			if((clock() - start)/static_cast<double>(CLOCKS_PER_SEC)>=600)
				break;
			if (iter_without_improvement >= restart_it) {
				generate_random_solution(n, p);
				calculate_new_delta(p, n, delta, current_cost, a, b);
				iter_without_improvement = 0;
			} else {
				perturbe(p, n, delta, current_cost, a, b, k);
			}
		}
		best_time -= start;
		total_time = clock() - start;
	} else {
		best_time = clock() - start;
		total_time = best_time;
	}

	// free memory
	delete[] p;
	for (int i=0; i < n; i = i+1)
	{
		delete[] delta[i];
	}
	delete[] delta;
}

int n;                    // problem size
type_matrix a, b;         // flows and distances matrices
type_vector solution;     // solution (permutation)
long cost, best_objective;                // solution cost
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

	if (!is_iteration_set) {
		max_iterations = 500;
	}

	ILS(n, a, b, solution, cost, max_iterations, best_objective, best_time, total_time);

	/* total_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n); */
	/* cout<<"Solution cost "<<cost<<" "<<(best_time)/static_cast<double>(CLOCKS_PER_SEC)<<" "<<(total_time)/static_cast<double>(CLOCKS_PER_SEC)<<endl; */
	/* cout << best_evals << " " << total_evals << endl; */
	for (int i =0; i < n; i++) {
		cout << solution[i] << " ";
	}
	cout << endl;
	cout << cost;

	/* write_results(cost, best_time, total_time, best_evals, total_evals); */

	delete[]solution;

	return 0;
}
