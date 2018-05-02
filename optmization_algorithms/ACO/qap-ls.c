#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>

#include "qap-ls.h"
#include "InOut.h"
#include "QAP.h"
#include "ants.h"
#include "utilities.h"

const char * const strings_ls_types [] = {
    "disable",
    "first 2-opt",
    "best 2-opt",
    "short tabu search",
    "long tabu search",
    NULL };

#define is_symmetric() ((d_symmetric_flag && f_symmetric_flag && null_diagonal_flag) \
                        || make_symmetric_flag)


const unsigned int LS_MAX = LS_UNKNOWN - 1;

unsigned int ls_flag;
long int nn_ls;  /* Unused here */
long int dlb_flag; /* set to one if don't look bits should be used */

long int null_diagonal_flag;
long int d_symmetric_flag;
long int f_symmetric_flag;
long int make_symmetric_flag;

static long int   **move_values;       /* move values in best improvement local search */
static long int   **tabu_values;       /* entries of matrix give iteration up to which an attribute is forbidden */
static long int  tabu_list_length;

const char * ls_type_to_string(unsigned int ls_type)
{
    if (ls_type > LS_MAX)
        abort();
    return strings_ls_types[ls_type];
}

void local_search(unsigned int ls_type, long int *tour)
{
    switch (ls_type) {
      case LS_first_2_opt:
          if (is_symmetric())
              first_2_opt_symmetric( tour );
          else
              first_2_opt_asymmetric( tour );
          break;
      case LS_best_2_opt:
          if (is_symmetric())
              best_2_opt_symmetric( tour );
          else
              best_2_opt_asymmetric( tour );
          break;
      case LS_tabu_search_short:
          tabu_search( tour, /*tabu_search_length =*/ 4);
          break;
      case LS_tabu_search_long:
          tabu_search( tour, /*tabu_search_length =*/ 10);
          break;
      default:
          fprintf(stderr,"type of local search procedure not correctly specified\n");
          exit(1);
    }
}

void first_2_opt_symmetric ( long int *q )
{
/*
      FUNCTION:      first improvement 2-opt local search for symmetric instances
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
*/

    long int   improvement = TRUE;
    long int   improve_item = FALSE;
    long int   u, v, k;
    long int   tmp;
    long int   original_symmetric_factor; /* = 2: original symmetric instance
					     = 1: original asymmetric instance
					  */
    long int  *dlb;                 /* don't look bits, just one vector needed */
    long int **d = instance.distance;
    long int **f = instance.flow;

    trace_print("first imp, symmetric case\n");
    dlb = calloc(n, sizeof(long int));

    if ( make_symmetric_flag )
	original_symmetric_factor = 1; /* compensation because of not dividing matrix by 2 */
    else
	original_symmetric_factor = 2;
    improvement   = TRUE;
    while ( improvement ) {
	improvement = FALSE;
	for ( u = 0 ; u < n ; u++ ) {
	    if ( dlb_flag && dlb[u] )
		continue;
	    improve_item = FALSE;
	    for ( v = 0 ; v < n ; v++ ) {
		if (u == v)
		    continue;
		tmp = 0;
		for ( k = 0 ; k < n ; k++ ) {
		    if ( (k != u) && (k != v) ) {
			tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
		    }
		}
		tmp *= original_symmetric_factor;
		if (tmp < 0) {
		    improvement = TRUE;
		    improve_item = TRUE;
		    dlb[u] = FALSE;
		    dlb[v] = FALSE;
		    swap(q, u, v);
		    trace_print("improvement %ld\n",tmp);
	        }
	    }
	    if ( improve_item )
		dlb[u] = FALSE;
	    else
		dlb[u] = TRUE;
	}
    }
    free (dlb);
}

void first_2_opt_asymmetric ( long int * q )
{
/*
      FUNCTION:      first improvement 2-opt local search for asymmetric instances
      INPUT:         pointer to initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
*/

    long int   improvement = TRUE;
    long int   improve_item = FALSE;
    long int   u, v, k;
    long int   tmp;
    long int  *dlb;                 /* don't look bits, just one vector needed */
    long int ** d = instance.distance;
    long int ** f = instance.flow;

    trace_print("first imp, asymmetric case\n");
    dlb = malloc( n * sizeof(long int) );
    for ( k = 0 ; k < n ; k++ ) {
        dlb[k] = FALSE;
    }

    while ( improvement ) {
        improvement = FALSE;
        for ( u = 0 ; u < n ; u++) {
            if ( dlb_flag && dlb[u] )
                continue;
            improve_item = FALSE;
            for ( v = 0 ; v < n ; v++) {
                if (u == v)
                    continue;
                tmp = 0;
                for ( k = 0 ; k < n ; k++ ) {
                    if ( (k != u) && (k != v) ) {
                        tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) +
                            d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) +
                            d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) +
                            d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
                    }
                }
                tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
                    d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
                    d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+
                    d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
                if (tmp < 0) {
                    improvement = TRUE;
                    improve_item = TRUE;
                    dlb[u] = FALSE;
                    dlb[v] = FALSE;
                    swap(q, u, v);
                }
            }
            if ( !improve_item )
                dlb[u] = TRUE;
        }
    }
    free (dlb);
}

void best_2_opt_symmetric ( long int * q )
{
/*
      FUNCTION:      best improvement 2-opt local search for symmetric instances
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENTS:      faster evaluation of moves with additional matrix move_values
*/

  long int   improvement = TRUE;
  long int   u, v, k;
  long int   tmp;
  long int   original_symmetric_factor;
  long int   **move_values;
  long int   first_it_flag = TRUE;
  long int   max_decrease;
  long int   rchosen = n, schosen = n;
  long int   r, s;
  long int ** d = instance.distance;
  long int ** f = instance.flow;

  trace_print("best imp, symmetric case\n");
  if ( make_symmetric_flag )
    original_symmetric_factor = 1;
  else
    original_symmetric_factor = 2;

  /* allocate and prepare matrix with move_values */
  move_values = generate_int_matrix (n, n);
  r = rchosen;
  s = schosen;

  while ( improvement ) {
    improvement = FALSE;
    max_decrease = LONG_MAX;
    /* in the first iteration the full neighborhood has to be evaluated */
    if (first_it_flag) {
      first_it_flag = FALSE;
      for ( u = 0 ; u < n-1 ; u++) {
	for ( v = u+1 ; v < n ; v++) {
	  tmp = 0;
	  for ( k = 0 ; k < n ; k++ ) {
	    if ( (k != u) && (k != v) ) {
	      tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	    }
	  }
	  tmp *= original_symmetric_factor;
	  move_values[u][v] = tmp;
	  if (tmp < max_decrease) {
	    max_decrease = tmp;
	    rchosen = u;
	    schosen = v;
/*    	    printf(" max-decr = %ld\n",tmp); */
	  }
	}
      }
    } else {
      for ( u = 0 ; u < n-1 ; u++) {
	for ( v = u+1 ; v < n ; v++) {
	  if (u == r || v == s || u == s || v == r) {
	    tmp = 0.;
	    for (k = 0 ; k < n ; k++) {
	      if ( (k != u) && (k != v) ) {
		tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	      }
	    }
	    tmp *= original_symmetric_factor;
	    move_values[u][v] = tmp;
	  } else { /* change derived from prev iteration, u and v differ from rchosen or schosen */
	    tmp = original_symmetric_factor * ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
	      ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] );
	    tmp += move_values[u][v];
	    move_values[u][v] = tmp;
	  }
	  if (tmp < max_decrease) {
	    max_decrease = tmp;
/*    	    printf(" max-decr = %ld\n",tmp); */
	    rchosen = u; /* memorize move */
	    schosen = v; /* memorize move */
	  }
	}
      }
    }
    if ( max_decrease < 0 ) {      /* Obj. function value can be improved */
      assert (rchosen < schosen);
      improvement = TRUE;
      swap(q, rchosen, schosen);
/*        printf(" max-decr = %ld\n",max_decrease); */
      r = rchosen; /* memorize previously done move */
      s = schosen; /* memorize previously done move */

      trace_print("improvement %ld, exchange %ld and %ld\n",max_decrease, rchosen, schosen);
    }
  }
  trace_print("best 2-opt, value after %ld\n",compute_evaluation_function(q));

  free ( move_values );
}


void best_2_opt_asymmetric ( long int * q )
/*
      FUNCTION:      best improvement 2-opt local search for asymmetric instances
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENTS:      fast evaluation of moves with additional matrix move_values
                     the first local search iteration is slow
                     local search looks for best neighboring solution in each iteration
*/
{
  long int   improvement = TRUE;
  long int   u, v, k;
  long int   tmp;
  long int   **move_values;             /* matrix of move values in previous iteration
					   allows for fast evaluation of neighbourhood
					*/
  long int   first_it_flag = TRUE;      /* first iteration of local search: TRUE */
  long int   max_decrease;              /* largest decrease found so far in neighbourhood scan */
  long int   rchosen = n, schosen = n;  /* memorize which is best move in current iteration */
  long int   r, s;                      /* memorize which is best move in previous iteration */
  long int ** d = instance.distance;
  long int ** f = instance.flow;

  trace_print("best imp, asymmetric case\n");

  move_values = generate_int_matrix (n, n);
  r = rchosen;
  s = schosen;

  while ( improvement ) {
    improvement = FALSE;
    max_decrease = LONG_MAX;
    /* in the first local search iteration the full neighborhood has to be evaluated */
    if (first_it_flag) {
      first_it_flag = FALSE;
      for ( u = 0 ; u < n-1 ; u++) {
	for ( v = u+1 ; v < n ; v++) {
	  tmp = 0;
	  for ( k = 0 ; k < n ; k++ ) {
	    if ( (k != u) && (k != v) ) {
	      tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) +
		d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) +
		d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) +
		d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	    }
	  }
	  tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	    d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	    d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+
	    d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	  move_values[u][v] = tmp;
	  if (tmp < max_decrease) {
	    max_decrease = tmp;
	    rchosen = u;
	    schosen = v;
	  }
	}
      }
    } else {
      for ( u = 0 ; u < n-1 ; u++) {
	for ( v = u+1 ; v < n ; v++) {
	  if (u == r || v == s || u == s || v == r) {
	    tmp = 0;
	    for ( k = 0 ; k < n ; k++ ) {
	      if ( (k != u) && (k != v) ) {
		tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) +
		  d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) +
		  d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) +
		  d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	      }
	    }
	    tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	      d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	      d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+
	      d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	    move_values[u][v] = tmp;
	    if (tmp < max_decrease) {
	      max_decrease = tmp;
	      rchosen = u;
	      schosen = v;
	      /*   	    printf(" max-decr = %ld\n",tmp); */
	    }
	  } else { /* change derived from move_values */
	    tmp = ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
	      ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] )
	      + ( d[u][r] - d[v][r] + d[v][s] - d[u][s] ) *
	      ( f[q[u]][q[s]] - f[q[v]][q[s]] + f[q[v]][q[r]] - f[q[u]][q[r]] );
 	    tmp += move_values[u][v];
	    move_values[u][v] = tmp;
	  }
	  if (tmp < max_decrease) {
	    max_decrease = tmp;
	    rchosen = u;
	    schosen = v;
	  }
	}
      }
    }
    if ( max_decrease < 0 ) {      /* Obj. function value can be improved */
      assert (rchosen < schosen);
      improvement = TRUE;
      swap(q, rchosen, schosen);
      r = rchosen; /* memorize previously done move */
      s = schosen;/* memorize previously done move */
      trace_print("improvement %ld, exchange %ld and %ld\n",max_decrease, rchosen, schosen);
    }
  }
  free ( move_values );
}

static void
best_2_opt_asymmetric_tabu ( long int * q, int first_it_flag,
                             long int rchosen, long int schosen)
{
/*
      FUNCTION:      best improvement 2-opt local search for asymmetric instances
                     computes the move values for the tabu search
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENTS:      fast evaluation of moves with additional matrix move_values
                     the first local search iteration is slow
                     local search looks for best neighboring solution in each iteration
      NOTES:         obviously, the code could overall be made a bit simpler if using this
                     function also in the iterative improvement algorithm.
*/

  long int   u, v, k;
  long int   r, s;                      /* memorize which is best move in previous iteration */
  long int   tmp;
  long int ** d = instance.distance;
  long int ** f = instance.flow;

  trace_print("best imp, asymmetric case\n");
  r = rchosen;
  s = schosen;

  /* in the first local search iteration the full neighborhood has to be evaluated */
  if (first_it_flag) {
    first_it_flag = FALSE;
    for ( u = 0 ; u < n-1 ; u++) {
      for ( v = u+1 ; v < n ; v++) {
	tmp = 0;
	for ( k = 0 ; k < n ; k++ ) {
	  if ( (k != u) && (k != v) ) {
	    tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) +
	      d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) +
	      d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) +
	      d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	  }
	}
	tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	  d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	  d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+
	  d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	move_values[u][v] = tmp;
      }
    }
  } else {
    for ( u = 0 ; u < n-1 ; u++) {
      for ( v = u+1 ; v < n ; v++) {
	if (u == r || v == s || u == s || v == r) {
	  tmp = 0;
	  for ( k = 0 ; k < n ; k++ ) {
	    if ( (k != u) && (k != v) ) {
	      tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) +
		d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) +
		d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) +
		d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
	    }
	  }
	  tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
	    d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
	    d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+
	    d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
	  move_values[u][v] = tmp;
	} else { /* change derived from move_values */
	  tmp = ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
	    ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] )
	      + ( d[u][r] - d[v][r] + d[v][s] - d[u][s] ) *
	    ( f[q[u]][q[s]] - f[q[v]][q[s]] + f[q[v]][q[r]] - f[q[u]][q[r]] );
	  tmp += move_values[u][v];
	  move_values[u][v] = tmp;
	}
      }
    }
  }
}



static void
best_2_opt_symmetric_tabu ( long int * q, int first_it_flag,
                            int rchosen, int schosen)
{
/*
      FUNCTION:      best improvement 2-opt local search for symmetric instances
                     computes the move values for the tabu search
      INPUT:         pointer to some initial assignment
      OUTPUT:        none
      (SIDE)EFFECTS: initial assignment is locally optimized
      COMMENTS:      faster evaluation of moves with additional matrix move_values
      NOTES:         obviously, the code could overall be made a bit simpler if using this
                     function also in the iterative improvement algorithm.
*/

  long int   u, v, k;
  long int   r, s;
  long int   tmp;
  long int   original_symmetric_factor;
  long int ** d = instance.distance;
  long int ** f = instance.flow;

  trace_print("best imp, symmetric case\n");
  if ( make_symmetric_flag )
    original_symmetric_factor = 1;
  else
    original_symmetric_factor = 2;

  r = rchosen;
  s = schosen;

  /* in the first iteration the full neighborhood has to be evaluated */
  if (first_it_flag) {
    first_it_flag = FALSE;
    for ( u = 0 ; u < n-1 ; u++) {
      for ( v = u+1 ; v < n ; v++) {
	tmp = 0;
	for ( k = 0 ; k < n ; k++ ) {
	  if ( (k != u) && (k != v) ) {
	    tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	  }
	}
	move_values[u][v] = tmp * original_symmetric_factor;
      }
    }
  } else {
    for ( u = 0 ; u < n-1 ; u++) {
      for ( v = u+1 ; v < n ; v++) {
	if (u == r || v == s || u == s || v == r) {
	  tmp = 0.;
	  for (k = 0 ; k < n ; k++) {
	    if ( (k != u) && (k != v) ) {
	      tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
	    }
	  }
	  move_values[u][v] = tmp * original_symmetric_factor;
	} else { /* change derived from prev iteration, u and v differ from rchosen or schosen */
	  tmp = original_symmetric_factor * ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
	    ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] );
	  tmp += move_values[u][v];
	  move_values[u][v] = tmp;
	}
      }
    }
  }
}



static long int ** init_move_values( )
/*
      FUNCTION:      allocate and initialize the move values of speed-up in tabu search
      INPUT:         instance size
      OUTPUT:        pointer to matrix, has to be freed before next tabu search run
      (SIDE)EFFECTS: current assignment is modified
*/
{
  long int i, j;
  long int **matrix;

  trace_print("initialize matrix of move values for tabu search\n");

  matrix = generate_int_matrix (n, n);
  for ( i = 0 ; i < n ; i++ ) {
      for (j = 0 ; j < n ; j++ ) {
          matrix[i][j] = 0;
      }
  }
  return matrix;
}


static long int ** init_tabu(void)
/*
      FUNCTION:      allocate and initialize the tabu values
      INPUT:         instance size
      OUTPUT:        pointer to matrix, has to be freed before next trial
      (SIDE)EFFECTS: current assignment is modified
*/
{
  long int i, j;
  long int     **matrix;

  trace_print("initialze matrix of tabu values\n");

  matrix = generate_int_matrix (n, n);
  for ( i = 0 ; i < n ; i++ ) {
    for (j = 0 ; j < n ; j++ ) {
      matrix[i][j] = 0;
    }
  }
  return matrix;
}

static long int choose_tabu_length(void)
/*
      FUNCTION:      choose a new tabu list length
      INPUT:         instance size
      OUTPUT:        none
*/
{
    double
        min      = 0.1  * n,
        max      = 1.1  * n,
        help     = min + ran01( &seed ) * (max - min);
  if (help < 2 )
      return 2.0;
  return help;
}

static void
make_tabu( long int * q, long int iteration, long int r, long int s )
/*
      FUNCTION:      make an invers move tabu (forbids a location for an item)
      INPUT:         pointer to some assignment, iteration number, two locations involved in move
      OUTPUT:        none
      (SIDE)EFFECTS: matrix of tabu_values is modified
*/
{
  tabu_values[r][q[r]] = iteration + tabu_list_length;
  tabu_values[s][q[s]] = iteration + tabu_list_length;
}

static long int  aspirating_iteration;

static void select_move( long int *q, long int current, long int iteration, long int best_so_far,
                         long int *rchosen, long int *schosen)
/*
      FUNCTION:      select the best move which is not tabu or is aspired
      INPUT:         pointer to some assignment, obj_f_val of current solution, iteration number
      OUTPUT:        none
      (SIDE)EFFECTS: global variables rchosen and schosen are assigned the locations involved in the move
*/
{
  long int   i, j;
  long int   max_decrease;
  long int   taboo, aspired;

  max_decrease = LONG_MAX;
  *rchosen = n; *schosen = n;
  for ( i = 0 ; i < n - 1 ; i++) {
    for (j = (i+1); j < n ; j++) {
      if ( (tabu_values[i][q[j]] > iteration ) && (tabu_values[j][q[i]] > iteration ) )
	taboo = TRUE;
      else
	taboo = FALSE;
      if ( ((current + move_values[i][j]) < best_so_far ) && taboo )
	aspired = TRUE;
      else if ( tabu_values[i][q[j]] < aspirating_iteration &&
	      tabu_values[j][q[i]] < aspirating_iteration )
	aspired = TRUE;
      else
	aspired = FALSE;
      if ( (move_values[i][j] < max_decrease && !taboo) || aspired ) {
          *rchosen = i;
          *schosen = j;
          if ( aspired && !taboo )
              max_decrease = LONG_MIN;
          else
              max_decrease = move_values[i][j];
      }
    }
  }
  assert ( *rchosen >= 0 && *rchosen < n);
  assert ( *schosen >= 0 && *schosen < n);
}



void tabu_search(long int *s, int tabu_search_length)
/*
      FUNCTION:     perform robust tabu search
      INPUT:        s: pointer to initial solution
                    tabu_search_length: how many iterations the tabu search is
                    run, will be multiplied by the instance size
      OUTPUT:       pointer to best solution
      (SIDE)EFFECTS:
*/
{
  long int  i;
  long int  iterations, obj_f_value, *b, *q;
  long int  best_so_far = 0;
  int first_it_flag;
  long int  iterations_before_aspiration;  /* robust tabu search specific, multiple of instance size */
  long int rchosen, schosen;

  trace_print("tabu search\n");
  move_values = init_move_values();
  tabu_values = init_tabu();
  b = malloc( n * sizeof(long int) );
  q = malloc( n * sizeof(long int) );
  for ( i = 0 ; i < n ; i++ ) {
    q[i] = s[i];
    b[i] = s[i];
  }
  obj_f_value       = compute_evaluation_function( b );
  best_so_far       = obj_f_value;
  first_it_flag     = TRUE;
  rchosen           = schosen      = 0;
  iterations        = 1;
  tabu_list_length  = choose_tabu_length();
  tabu_search_length *= n;

  /* the following was used in original Robust Tabu Search; probably never applies here */
  iterations_before_aspiration = (long int)(3 * n * n);

/*    printf("tabu_search_length %ld\n",tabu_search_length); */

  while ( iterations < tabu_search_length ) {

    aspirating_iteration = iterations - iterations_before_aspiration;

    if (is_symmetric())
        best_2_opt_symmetric_tabu ( q, first_it_flag, rchosen, schosen);
    else
        best_2_opt_asymmetric_tabu ( q, first_it_flag, rchosen, schosen);

    first_it_flag = FALSE;

    select_move(q , obj_f_value, iterations, best_so_far,
                &rchosen, &schosen);

    make_tabu( q, iterations, rchosen, schosen ); /* make_tabu has to go before swap */

    swap(q, rchosen, schosen);

    obj_f_value += move_values[rchosen][schosen];

    if ( obj_f_value < best_so_far ) {
      best_so_far = obj_f_value;
/*        printf("best %ld\t time %f\n",best_so_far,elapsed_time( VIRTUAL )); */
      for ( i = 0 ; i < n ; i++ )
	b[i] = q[i];
    }

    if ( !(iterations % ( long int )(2.2 * n + 4))) {
      tabu_list_length = choose_tabu_length();
      trace_print(" tabu_length = %ld, iteration %ld\n", tabu_list_length, iterations);
    }
    iterations++;

  }

  DEBUG (if ( compute_evaluation_function( b ) != best_so_far )
             fprintf(stderr, "Some error must have occurred in local search routine,\n values do not match\n");
      );
  for ( i = 0 ; i < n ; i++ )
    s[i] = b[i];
  free ( b );
  free ( q );
  free ( move_values );
  free ( tabu_values );
}
