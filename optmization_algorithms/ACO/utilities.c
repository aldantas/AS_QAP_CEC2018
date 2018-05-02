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
      File:    utilities.c
      Author:  Thomas Stuetzle
      Purpose: some additional useful procedures
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "InOut.h"
#include "utilities.h"
#include "ants.h"
#include "timer.h"


double mean( long int *values, long int max )
/*
      FUNCTION:       compute the average value of an integer array of length max
      INPUT:          pointer to array, length of array
      OUTPUT:         average
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   m;

  m = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    m += (double)values[j];
  }
  m = m / (double)max;
  return m;
}



double meanr( double *values, long int max )
/*
      FUNCTION:       compute the average value of a floating number array of length max
      INPUT:          pointer to array, length of array
      OUTPUT:         average
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   m;

  m = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    m += values[j];
  }
  m = m / (double)max;
  return m;
}



double std_deviation( long int *values, long int max, double mean )
/*
      FUNCTION:       compute the standard deviation of an integer array
      INPUT:          pointer to array, length of array, mean
      OUTPUT:         standard deviation
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   dev = 0.;

  if (max <= 1)
    return 0.;
  for ( j = 0 ; j < max; j++ ) {
    dev += ((double)values[j] - mean) * ((double)values[j] - mean);
  }
  return sqrt(dev/(double)(max - 1));
}



double std_deviationr( double *values, long int max, double mean )
/*
      FUNCTION:       compute the standard deviation of a floating number array
      INPUT:          pointer to array, length of array, mean
      OUTPUT:         standard deviation
      (SIDE)EFFECTS:  none
*/
{
  long int j;
  double   dev;

  if (max <= 1)
    return 0.;
  dev = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    dev += ((double)values[j] - mean) * ((double)values[j] - mean);
  }
  return sqrt(dev/(double)(max - 1));
}



long int best_of_vector( long int *values, long int l )
/*
      FUNCTION:       return the minimum value in an integer value
      INPUT:          pointer to array, length of array
      OUTPUT:         smallest number in the array
      (SIDE)EFFECTS:  none
*/
{
  long int min, k;

  k = 0;
  min = values[k];
  for( k = 1 ; k < l ; k++ ) {
    if( values[k] < min ) {
      min = values[k];
    }
  }
  return min;
}



long int worst_of_vector( long int *values, long int l )
/*
      FUNCTION:       return the maximum value in an integer value
      INPUT:          pointer to array, length of array
      OUTPUT:         largest number in the array
      (SIDE)EFFECTS:  none
*/
{
  long int max, k;

  k = 0;
  max = values[k];
  for( k = 1 ; k < l ; k++ ) {
    if( values[k] > max ){
      max = values[k];
    }
  }
  return max;
}



double quantil(long int v[], double q, long int l)
/*
      FUNCTION:       return the q-quantil of an ordered integer array
      INPUT:          one array, desired quantil q, length of array
      OUTPUT:         q-quantil of array
      (SIDE)EFFECTS:  none
*/
{
  long int i,j;
  double tmp;

  tmp = q * (double)l;
  if ((double)((long int)tmp) == tmp) {
    i = (long int)tmp;
    j = (long int)(tmp + 1.);
    return ((double)v[i-1] + (double)v[j-1]) / 2.;
  } else {
    i = (long int)(tmp +1.);
    return v[i-1];
  }
}



void swap(long int v[], long int i, long int j)
/*
      FUNCTION:       auxiliary routine for sorting an integer array
      INPUT:          array, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of array are swapped
*/
{
  long int tmp;

  tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
}




void sort(long int v[], long int left, long int right)
/*
      FUNCTION:       recursive routine (quicksort) for sorting an array
      INPUT:          one array, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
  long int k, last;

  if (left >= right)
    return;
  swap(v, left, (left + right)/2);
  last = left;
  for (k=left+1; k <= right; k++)
    if (v[k] < v[left])
      swap(v, ++last, k);
  swap(v, left, last);
  sort(v, left, last);
  sort(v, last+1, right);
}



void swap2(long int v[], long int v2[], long int i, long int j)
/*
      FUNCTION:       auxiliary routine for sorting an integer array
      INPUT:          two arraya, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
  long int tmp;

  tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
  tmp = v2[i];
  v2[i] = v2[j];
  v2[j] = tmp;
}



void sort2(long int v[], long int v2[], long int left, long int right)
/*
      FUNCTION:       recursive routine (quicksort) for sorting one array; second
                      arrays does the same sequence of swaps
      INPUT:          two arrays, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
  long int k, last;

  if (left >= right)
    return;
  swap2(v, v2, left, (left + right)/2);
  last = left;
  for (k=left+1; k <= right; k++)
    if (v[k] < v[left])
      swap2(v, v2, ++last, k);
  swap2(v, v2, left, last);
  sort2(v, v2, left, last);
  sort2(v, v2, last+1, right);
}


/* constants for a random number generator, for details see numerical recipes in C */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran01( long *idum )
/*
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  /* long k; */
  /* double ans; */

/*  trace_print("seed: %ld\n", *idum);*/

  /* k =(*idum)/IQ; */
  /* *idum = IA * (*idum - k * IQ) - IR * k; */
  /* if (*idum < 0 ) *idum += IM; */
  /* ans = AM * (*idum); */
  /* return ans; */
  return (double)rand()/(double)RAND_MAX;

}

long int random_number( long *idum )
/*
      FUNCTION:       generate an integer random number
      INPUT:          pointer to variable containing random number seed
      OUTPUT:         integer random number uniformly distributed in {0,2147483647}
      (SIDE)EFFECTS:  random number seed is modified (important, has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  /* long k; */

  /* k =(*idum)/IQ; */
  /* *idum = IA * (*idum - k * IQ) - IR * k; */
  /* if (*idum < 0 ) *idum += IM; */
  return rand();
}

long int ** generate_int_matrix( long int n, long int m)
/*
      FUNCTION:       malloc a matrix and return pointer to it
      INPUT:          size of matrix as n x m
      OUTPUT:         pointer to matrix
      (SIDE)EFFECTS:
*/
{
  long int i;
  long int **matrix;

  if((matrix = malloc(sizeof(long int) * n * m +
		      sizeof(long int *) * n	 )) == NULL){
    printf("Out of memory, exit.");
    exit(1);
  }
  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (long int *)(matrix + n) + i*m;
  }

  return matrix;
}

double ** generate_double_matrix( long int n, long int m)
/*
      FUNCTION:       malloc a matrix and return pointer to it
      INPUT:          size of matrix as n x m
      OUTPUT:         pointer to matrix
      (SIDE)EFFECTS:
*/
{

  long int i;
  double **matrix;

  if((matrix = malloc(sizeof(double) * n * m +
		      sizeof(double *) * n	 )) == NULL){
    printf("Out of memory, exit.");
    exit(1);
  }
  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (double *)(matrix + n) + i*m;
  }
  return matrix;
}


long int * generate_random_permutation( long int vector_size )
/*
      FUNCTION:       generate a random permutation of the integers 0 .. n-1
      INPUT:          length of the array
      OUTPUT:         pointer to the random permutation
      (SIDE)EFFECTS:  the array holding the random permutation is allocated in this
                      function. Don't forget to free again the memory!
*/
{
   /* http://en.wikipedia.org/wiki/Fisher-Yates_shuffle */
   long int * vector = malloc(vector_size * sizeof(long int));
   int  i;

   for (i = 0; i < vector_size; i++)
       vector[i] = i;

   for (i = vector_size - 1; i > 0 ; i--) {
       /* 0 <= j <= i is the correct range for unbiased shuffling.  */
       int j = (int) (ran01 ( &seed ) * (i + 1));
       assert(j >= 0 && j <= i);
       /* Testing if (i == j) decreases performance for large vectors.  */
       long int tmp = vector[i];
       vector[i] = vector[j];
       vector[j] = tmp;
   }
   return vector;
}

long int check_null_diagonal ( long int **matrix, long int size )
{
/*
      FUNCTION:      check whether the Matrix matrix has a zero diagonal
      INPUT:         pointer to the matrix
      OUTPUT:        TRUE if null diagonal, otherwise FALSE
*/
  long int   i;

  for ( i = 0 ; i < size ; i++ ) {
    if( matrix[i][i] != 0 ) {
      return FALSE;
    }
  }
  return TRUE;
}

bool check_permutation(const long int *t, int size)
{
    bool ok = true;
    int i;
    int * used = calloc (size, sizeof(int));

    if (t == NULL) {
        fprintf (stderr,"\n%s:error: permutation is not initialized!", __FUNCTION__);
        ok = false;
        goto error;
    }

    for (i = 0; i < size; i++) {
        if (used[t[i]]) {
            fprintf(stderr,"\n%s:error: solution vector has two times the value %ld (last position: %d)", __FUNCTION__, t[i], i);
            ok = false;
            goto error;
        }
        else
            used[t[i]] = TRUE;
    }

    for (i = 0; i < size; i++) {
        if (!used[i]) {
            fprintf(stderr,"\n%s:error: vector position %d not occupied", __FUNCTION__, i);
            ok = false;
            goto error;
        }
    }

error:
    free(used);
    return ok;
}

long int check_symmetry ( long int **matrix, long int size )
/*
      FUNCTION:      check whether the Matrix matrix is symmetric
      INPUT:         pointer to the matrix
      OUTPUT:        TRUE if symmetric, otherwise FALSE
*/
{
  long int   i, j;

  for ( i = 0 ; i < size - 1 ; i++ ) {
    for ( j = i + 1 ; j < size ; j++ ) {
      if( matrix[i][j] != matrix[j][i] )
	return FALSE;
    }
  }
  return TRUE;
}


long int num_different_edges (const long int *p1, const long int *p2, int n)
{
    int i, j, h, pos, pred;
    long int distance = 0;
    long int * pos2 = malloc(n * sizeof(long int));

    for ( i = 0 ; i < n ; i++ ) {
	pos2[p2[i]] = i;
    }

    for ( i = 0 ; i < n ; i++ ) {
	j = p1[i];
	h = p1[i+1];
	pos = pos2[j];
	if (pos - 1 < 0)
	    pred = n - 1;
	else
	    pred = pos - 1;
	if (p2[pos+1] == h || p2[pred] == h)
	    ; /* do nothing, edge is common */
	else /* edge (j,h) does not occur in p2 */
	    distance++;
    }

    free ( pos2 );
    return distance;
}

long int num_different_positions (const long int *p1, const long int *p2, int n)
{
    int i;
    long int distance = 0;

    for ( i = 0 ; i < n ; i++ ) {
	if (p1[i] != p2[i])
            distance++;
    }
    return distance;
}

void matrix_long_print (long int **matrix, int n, int m)
{
    int i, j;
    printf("\n");
    for ( i = 0 ; i < n ; i++ ) {
        for ( j = 0 ; j < m ; j++ ) {
            printf(" %ld ", matrix[i][j]);
        }
        printf("\n");
    }
}

void vector_long_fprint (FILE * stream, const long int *v, size_t size)
{
    size_t i;
    for ( i = 0 ; i < size ; i++ ) {
        fprintf(stream, " %ld", v[i]);
    }
}

void vector_long_print (const long int *v, size_t size)
{
    vector_long_fprint (stdout, v, size);
}

