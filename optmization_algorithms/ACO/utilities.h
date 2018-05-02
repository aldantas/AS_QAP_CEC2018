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
      File:    utilities.h
      Author:  Thomas Stuetzle
      Purpose: some additional useful procedures
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
#include <stdbool.h>
#include <math.h>
#define INFTY                 LONG_MAX
#define MAXIMUM_NO_TRIES      100

#define TRUE  1
#define FALSE 0

/* general macros */

#define MAX(x,y)        ((x)>=(y)?(x):(y))
#define MIN(x,y)        ((x)<=(y)?(x):(y))

#ifndef DODEBUG
#define DODEBUG 0
#endif
#if DODEBUG >= 1
#define DEBUG(...) __VA_ARGS__;(void)0
#else
#define DEBUG(...) (void)0
#endif

#ifndef TRACE
#define TRACE 0
#endif
#if TRACE
#define trace_print(...) fprintf(stderr, "trace: " __VA_ARGS__)
#else
#define trace_print(...) (void)0
#endif

#define XOR(x,y) ((x && !y) || (!x && y))

extern long int seed;

static inline double powx(double x, double y)
{
    /* On 64-bits, pow() may be more than 10,000 times slower for some
       inputs than for other nearby inputs.  This affects only pow(), and not
       powf() nor powl().  We use powl to work-around this bug in GNU libc:
       https://sourceware.org/bugzilla/show_bug.cgi?id=13932 */
#if defined(__x86_64__) || defined(__ppc64__)
    return powl(x, y);
#else
    return pow(x, y);
#endif
}

double mean ( long int *values, long int max);

double meanr ( double *values, long int max );

double std_deviation ( long int *values, long int i, double mean );

double std_deviationr ( double *values, long int i, double mean );

long int best_of_vector ( long int *values, long int i );

long int worst_of_vector ( long int *values, long int i );

void swap ( long int v[], long int i, long int j );

void sort ( long int v[], long int left, long int right );

double quantil ( long int vector[], double q, long int numbers );

void swap2(long int v[], long int v2[], long int i, long int j);

void sort2(long int v[], long int v2[], long int left, long int right);

double ran01 ( long *idum );

long int random_number ( long *idum );

void matrix_long_print (long int **matrix, int n, int m);
void vector_long_print (const long int *v, size_t size);
void vector_long_fprint (FILE *stream, const long int *v, size_t size);

long int ** generate_int_matrix( long int n, long int m);

double ** generate_double_matrix( long int n, long int m);

long int * generate_random_permutation( long int n );

bool check_permutation(const long int *t, int size);

long int check_symmetry ( long int **matrix, long int size);

long int check_null_diagonal ( long int **matrix, long int size );

long int num_different_edges (const long int *p1, const long int *p2, int n);

long int num_different_positions (const long int *p1, const long int *p2, int n);
