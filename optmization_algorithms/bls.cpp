/*****************************************************************/
// Implementation of the breakout local search of:
// "Breakout local search for the quadratic assignment problem",
// Applied Mathematics and Computation 219(9), 1991, 4800-4815.
//
// Data file format:
//  n,
// (nxn) flow matrix,
// (nxn) distance matrix
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
#include <ctime>
#include <cmath>

#include "common.h"

/* #define mem_size 100000 */

using namespace std;

const long infinite = 999999999;

long iteration = 0;

//BLS parameters
double r1 = 0.7,r2 = 0.2; // Parameters used for directed perturbation. Number of moves elapsed before perturb.
double init_pert_str = 0.15; //initial perturbation strength
double perturb_str; // actual perturbation strength
int T = 2500;
double P0=0.75, Q=0.3;
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

void apply_move(type_vector & p,long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol,
		int i_retained, int j_retained, clock_t & best_time)
{
	if(i_retained!=-1 && j_retained!=-1) // apply the selected perturbation move
	{
		last_swaped[i_retained][j_retained] = iteration;
                last_swaped[j_retained][i_retained] = iteration;
		swap(p[i_retained], p[j_retained]);
		current_cost = current_cost + delta[i_retained][j_retained];
		update_matrix_of_move_cost(i_retained, j_retained, n, delta, p, a, b);
		if(current_cost<best_cost)
		{
			best_cost = current_cost;
                        best_time = clock();
			iter_without_improvement = 0;
			for (int m = 0; m < n; m = m+1) best_sol[m] = p[m];
			best_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n);
		}
	}
	iteration++;
}

void tabu_search_perturb(type_vector & p, long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, long init_cost,
                clock_t & best_time)
{
	int i_retained=-1, j_retained=-1,  i, j, min_delta;
	min_delta = infinite;
	min_delta = 999999999;
        for (i = 0; i < n-1; i = i + 1)
        {
           for (j = i+1; j < n; j = j+1)
           {
             if((current_cost + delta[i][j])!=init_cost && delta[i][j] < min_delta && ((last_swaped[i][j] +n*r1+ rand()%(static_cast<int>(n*r2)))<iteration or (current_cost + delta[i][j])<best_cost))
             {
                i_retained = i; j_retained = j;
                min_delta = delta[i][j];
             };
           };
        };
	apply_move(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, i_retained, j_retained, best_time);
}
void recency_based_perturb(type_vector & p, long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, clock_t & best_time)
{
	int i, j, i_retained, j_retained, min = infinite;

        for (i = 0; i < n-1; i = i + 1)
        {
           for (j = i+1; j < n; j = j+1)
           {
              if(last_swaped[i][j]<min)
              {
                 i_retained = i; j_retained = j;
                 min = last_swaped[i][j];
              };
           };
        };
	apply_move(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, i_retained, j_retained, best_time);
}
void random_perturb(type_vector & p, long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, long init_cost, clock_t & best_time)
{
	int i_retained = rand()%n;
        int j_retained = rand()%n;
        if(i_retained>j_retained)
             swap(i_retained, j_retained);
        while(i_retained==j_retained || (current_cost + delta[i_retained][j_retained])==init_cost)
        {
             j_retained = rand()%n;
             if(i_retained>j_retained)
               swap(i_retained, j_retained);
        }
	apply_move(p, n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, i_retained, j_retained, best_time);
}
void perturbe(type_vector & p,long n, type_matrix & delta, long & current_cost, type_matrix & a, type_matrix & b,
             type_matrix & last_swaped,  int & iter_without_improvement, type_vector & best_sol, long & best_cost, clock_t & best_time)
{
   int i_retained=-1, j_retained=-1, k, i, j, min_delta, min;
   long cost = current_cost;

   double d = static_cast<double>(iter_without_improvement)/T;
   double e = pow(2.718, -d);
   int perturbation_bit = 0;
   if(e<P0)
      e=P0;
   if(e>(static_cast<double>(rand()%101)/100.0))
        perturbation_bit=1;

   for(k =0; k<(perturb_str); k++)
   {
      i_retained=-1; j_retained=-1;
      if(perturbation_bit==1)
      {
           tabu_search_perturb(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, cost, best_time);
      }
       else
       {
         if(Q>(static_cast<double>(rand()%101)/100.0))
           recency_based_perturb(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, best_time);
         else
            random_perturb(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, cost, best_time);
       }

   }
}

void allocate_memory_and_initialize(type_vector & p, long n, type_matrix & delta, type_matrix & last_swaped, long & current_cost,
		type_matrix a, type_matrix b, type_vector & best_sol, long & best_cost)
{
	int i, j;
	/***************** dynamic memory allocation *******************/
	p = new int[n+1];
	delta = new long* [n+1];
	last_swaped = new long* [n+1];
	for (i = 0; i < n; i = i+1)
	{
		delta[i] = new long[n+1];
		last_swaped[i] = new long [n+1];
	}
	/************** current solution initialization ****************/
	for (i = 0; i < n; i = i + 1)
		p[i] = best_sol[i];
	/************** tabu list initialization ****************/
	for (i = 0; i < n; i = i+1)
	{
		for(j=0; j<n; j++)
		{
			last_swaped[i][j] = 0;
		}
	}
	/********** initialization of current solution value ***********/
	/**************** and matrix of cost of moves  *****************/
	current_cost = 0;
	for (i = 0; i < n; i = i + 1)
		for (j = 0; j < n; j = j + 1)
		{
			current_cost = current_cost + a[i][j] * b[p[i]][p[j]];
			if (i < j)
			{
				delta[i][j] = compute_delta(n, a, b, p, i, j);
			};
		};
	evaluations++;
	best_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n);
	best_cost = current_cost;
}

bool best_improvement_move(type_vector & p,long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, clock_t & best_time)
{
    int i, j, i_retained, j_retained;
    long min_delta = infinite;   // retained move cost
    //select the best swap move for the descent local search phase
    for (i = 0; i < n-1; i = i + 1)
    {
		for (j = i+1; j < n; j = j+1)
		{
			if(delta[i][j] < min_delta)
			{
				i_retained = i; j_retained = j;
				min_delta = delta[i][j];
			};
		};
    };
    // apply the above selected move if it leads to a solution better than the current one
    if((current_cost + delta[i_retained][j_retained])<current_cost)
    {
	apply_move(p, n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, i_retained, j_retained, best_time);
		return 1;
    }else
		return 0;
}
void determine_jump_magnitude(int & iter_without_improvement, int descent_num, long previous_cost,
		long current_cost, long n)
{
	//the following lines determine the number of perturbation moves (jump magnitude)
	if(iter_without_improvement>T) // if the best found solution is not improved during the last T descent phases
	{
		// strong pertubation required
		iter_without_improvement = 0;
		perturb_str = n*(0.4 + (static_cast<double>(rand()%20)/100.0));
	}
	else if((descent_num!=0 && previous_cost != current_cost) || (descent_num!=0  && previous_cost != current_cost))
		// Escaped from the previous local optimum.
	{
		iter_without_improvement++;
		perturb_str=ceil(n*init_pert_str); if(perturb_str<5) perturb_str = 5;
	}
	else if(previous_cost == current_cost) // search returned to the previous local optimum
	{
		perturb_str+= 1;
	}
}
void BLS(long n,                  // problem size
		type_matrix  a,         // flows matrix
		type_matrix  b,         // distance matrix
		type_vector & best_sol,  // best solution found
		long & best_cost,      // cost of best solution
		long num_iterations, long best_objective, clock_t & best_time,
		clock_t & total_time)

{
	type_vector p;                        // current solution
	type_matrix delta;
	type_matrix last_swaped; // keeps track of the iteration number when a move was last performed

	long current_iteration;
	long current_cost, previous_cost;    // current sol. value and previous sol. value
	int descent_num;
	int iter_without_improvement = 0; // counter of the number of consecutive descent phases with no improvement of the best solution
	perturb_str = ceil(init_pert_str*n); // initialize the number of perturbation moves
	bool descent;

	clock_t start = clock();
	allocate_memory_and_initialize(p, n, delta, last_swaped, current_cost, a, b, best_sol, best_cost);
	previous_cost = current_cost;

	/******************** main BLS loop ********************/
	if (!has_zero_matrix) {
		for (current_iteration = 0; current_iteration < num_iterations; current_iteration++)
		{
			//best improvement descent procedure
			descent = best_improvement_move(p, n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, best_time);
			descent_num=0;
			while(descent) // while an improving solution in the neighborhood exists
			{
				descent = best_improvement_move(p, n, delta, current_cost, a, b, last_swaped, iter_without_improvement,
						best_cost, best_sol, best_time);
				descent_num++;
			}

			// a this point, a local optimum is attained and BLS enters a perturbation phase
			determine_jump_magnitude(iter_without_improvement, descent_num, previous_cost, current_cost, n);

			//update the objective value of the previous local optimum
			previous_cost = current_cost;

			perturbe(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_sol, best_cost, best_time);
			if(ceil((clock() - start)/static_cast<double>(CLOCKS_PER_SEC))>=600)
				break;
		}
		best_time -= start;
		total_time = clock() - start;
	} else {
		best_time = clock() - start;
		total_time = best_time;
	}
	printf("iters: %ld\n", current_iteration);

	// free memory
	delete[] p;
	for (int i=0; i < n; i = i+1)
	{
		delete[] delta[i];
		delete [] last_swaped[i];
	}
	delete[] delta;  delete [] last_swaped;
} // BLS

void generate_random_solution(long n, type_vector  & p)
{
	int i;
	for (i = 0; i < n; i = i+1) p[i] = i;
	for (i = 1; i < n; i = i+1) transpose(p[i], p[i + rand()%(n-i)]);
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
		max_iterations = 100 * n;
	}

	BLS(n, a, b, solution, cost, max_iterations, best_objective, best_time, total_time);

	total_evals = evaluations + static_cast<long>((linear_evals + static_cast<long>(const_evals/n))/n);
	cout<<"Solution cost "<<cost<<" "<<(best_time)/static_cast<double>(CLOCKS_PER_SEC)<<" "<<(total_time)/static_cast<double>(CLOCKS_PER_SEC)<<endl;
	cout << best_evals << " " << total_evals << endl;

	write_results(cost, best_time, total_time, best_evals, total_evals);

	delete[]solution;

	return 0;
}
