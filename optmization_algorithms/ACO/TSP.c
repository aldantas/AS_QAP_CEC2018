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
      File:    TSP.c
      Author:  Thomas Stuetzle
      Purpose: TSP related procedures, distance computation, neighbour lists
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


#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "InOut.h"
#include "TSP.h"
#include "ants.h"
#include "ls.h"
#include "utilities.h"

#define M_PI 3.14159265358979323846264

long int n;          /* number of cities in the instance to be solved */

struct problem instance;

static double dtrunc (double x)
{
    int k;

    k = (int) x;
    x = (double) k;
    return x;
}

long int  (*distance)(long int, long int);  /* function pointer */

/*    
      FUNCTION: the following four functions implement different ways of 
                computing distances for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
*/

long int round_distance (long int i, long int j) 
/*    
      FUNCTION: compute Euclidean distances between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd) + 0.5;

    return (long int) r;
}

long int ceil_distance (long int i, long int j) 
/*    
      FUNCTION: compute ceiling distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd);

    return (long int)(ceil (r));
}

long int geo_distance (long int i, long int j) 
/*    
      FUNCTION: compute geometric distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: adapted from concorde code
                for the definition of how to compute this distance see TSPLIB
*/
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    long int dd;
    double x1 = instance.nodeptr[i].x, x2 = instance.nodeptr[j].x, 
	y1 = instance.nodeptr[i].y, y2 = instance.nodeptr[j].y;

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (y1);
    min = y1 - deg;
    longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (y2);
    min = y2 - deg;
    longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    return dd;

}

long int att_distance (long int i, long int j) 
/*    
      FUNCTION: compute ATT distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    long int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}



long int ** compute_distances(void)
/*    
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
*/
{
    long int     i, j;
    long int     **matrix;

    if((matrix = malloc(sizeof(long int) * n * n +
			sizeof(long int *) * n	 )) == NULL){
	fprintf(stderr,"Out of memory, exit.");
	exit(1);
    }
    for ( i = 0 ; i < n ; i++ ) {
	matrix[i] = (long int *)(matrix + n) + i*n;
	for ( j = 0  ; j < n ; j++ ) {
	    matrix[i][j] = distance(i, j);
	}
    }
    return matrix;
}



long int ** compute_nn_lists(struct problem *instance)
/*    
      FUNCTION: computes nearest neighbor lists of depth nn for each city
      INPUT:    none
      OUTPUT:   pointer to the nearest neighbor lists
*/
{
    long int i, node, nn;
    long int *distance_vector;
    long int *help_vector;
    long int **m_nnear;
 
    trace_print("\n computing nearest neighbor lists, ");

    nn = MAX(nn_ls, nn_ants);
    if ( nn >= n ) 
	nn = n - 1;
    DEBUG ( assert( n > nn ) );
    
    trace_print("nn = %ld ... \n",nn); 

    if ((m_nnear = malloc(sizeof(long int) * n * nn
                         + n * sizeof(long int *))) == NULL){
	exit(EXIT_FAILURE);
    }
    distance_vector = calloc(n, sizeof(long int));
    help_vector = calloc(n, sizeof(long int));
 
    for ( node = 0 ; node < n ; node++ ) {  /* compute cnd-sets for all node */
	m_nnear[node] = (long int *)(m_nnear + n) + node * nn;

	for ( i = 0 ; i < n ; i++ ) {  /* Copy distances from nodes to the others */
	    distance_vector[i] = instance->distance[node][i];
	    help_vector[i] = i;
	}
	distance_vector[node] = LONG_MAX;  /* city is not nearest neighbour */
	sort2(distance_vector, help_vector, 0, n-1);
	for ( i = 0 ; i < nn ; i++ ) {
	    m_nnear[node][i] = help_vector[i];
	}
    }
    free(distance_vector);
    free(help_vector);
    trace_print("\n    .. done\n");
    instance->nn_list = m_nnear;
    return instance->nn_list;
}



long int compute_tour_length(const long int *t ) 
/*    
      FUNCTION: compute the tour length of tour t
      INPUT:    pointer to tour t
      OUTPUT:   tour length of tour t
*/
{
    int      i;
    long int tour_length = 0;
  
    for ( i = 0 ; i < n ; i++ ) {
	tour_length += instance.distance[t[i]][t[i+1]];
    }
    return tour_length;
}


int check_solution(const long int *t)
{
    int i;
    const int size = n;

    if (!check_permutation (t, size))
        goto error;

    if (t[0] != t[size]) {
        fprintf(stderr,"\n%s:error: permutation is not a closed tour.", __FUNCTION__);
        goto error;
    }
    return TRUE;

error:
    fprintf(stderr,"\n%s:error: solution_vector:", __FUNCTION__);
    for (i = 0; i < size; i++)
        fprintf(stderr, " %ld", t[i]);
    fprintf(stderr,"\n");
    return FALSE;
}

struct point * read_etsp(const char *tsp_file_name) 
/*    
      FUNCTION: parse and read instance file
      INPUT:    instance name
      OUTPUT:   list of coordinates for all nodes
      COMMENTS: Instance files have to be in TSPLIB format, otherwise procedure fails
*/
{
    FILE         *tsp_file;
    char         buf[LINE_BUF_LEN];
    long int     i, j;
    struct point *nodeptr;

    tsp_file = fopen(tsp_file_name, "r");
    if ( tsp_file == NULL ) {
	fprintf(stderr,"No instance file specified, abort\n");
	exit(1);
    }
    assert(tsp_file != NULL);
    printf("\nreading tsp-file %s ... \n\n", tsp_file_name);

    fscanf(tsp_file,"%s", buf);
    while ( strcmp("NODE_COORD_SECTION", buf) != 0 ) {
	if ( strcmp("NAME", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    trace_print("%s ", buf);
	    fscanf(tsp_file, "%s", buf);
	    strcpy(instance.name, buf);
	    trace_print("%s \n", instance.name);
	    buf[0]=0;
	}
	else if ( strcmp("NAME:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    strcpy(instance.name, buf);
	    trace_print("%s \n", instance.name);
	    buf[0]=0;
	}
	else if ( strcmp("COMMENT", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    trace_print("%s", buf);
	    buf[0]=0;
	}
	else if ( strcmp("COMMENT:", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    trace_print("%s", buf);
	    buf[0]=0;
	}
	else if ( strcmp("TYPE", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    trace_print("%s ", buf);
	    fscanf(tsp_file, "%s", buf);
	    trace_print("%s\n", buf);
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	    buf[0]=0;
	}
	else if ( strcmp("TYPE:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    trace_print("%s\n", buf);
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	    buf[0]=0;
	}
	else if( strcmp("DIMENSION", buf) == 0 ){
	    fscanf(tsp_file, "%s", buf);
	    trace_print("%s ", buf);
	    fscanf(tsp_file, "%ld", &n);
	    instance.n = n;
	    trace_print("%ld\n", n);
	    assert ( n > 2 && n < 6000);
	    buf[0]=0;
	}
	else if ( strcmp("DIMENSION:", buf) == 0 ) {
	    fscanf(tsp_file, "%ld", &n);
	    instance.n = n;
	    trace_print("%ld\n", n);
	    assert ( n > 2 && n < 6000);
	    buf[0]=0;
	}
	else if( strcmp("DISPLAY_DATA_TYPE", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    trace_print("%s", buf);
	    buf[0]=0;
	}
	else if ( strcmp("DISPLAY_DATA_TYPE:", buf) == 0 ) {
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    trace_print("%s", buf);
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE", buf) == 0 ){
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    trace_print("%s ", buf);
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    trace_print("%s\n", buf);
	    if ( strcmp("EUC_2D", buf) == 0 ) {
		distance = round_distance;
	    }
	    else if ( strcmp("CEIL_2D", buf) == 0 ) {
		distance = ceil_distance;
	    }
	    else if ( strcmp("GEO", buf) == 0 ) {
		distance = geo_distance;
	    }
	    else if ( strcmp("ATT", buf) == 0 ) {
		distance = att_distance;
	    }
	    else
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
	    strcpy(instance.edge_weight_type, buf);
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 ){
	    /* set pointer to appropriate distance function; has to be one of 
	       EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    trace_print("%s\n", buf);
		printf("%s\n", buf);
	    printf("%s\n", buf);
	    if ( strcmp("EUC_2D", buf) == 0 ) {
		distance = round_distance;
	    }
	    else if ( strcmp("CEIL_2D", buf) == 0 ) {
		distance = ceil_distance;
	    }
	    else if ( strcmp("GEO", buf) == 0 ) {
		distance = geo_distance;
	    }
	    else if ( strcmp("ATT", buf) == 0 ) {
		distance = att_distance;
	    }
	    else {
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
		exit(1);
	    }
	    strcpy(instance.edge_weight_type, buf);
	    buf[0]=0;
	}
	buf[0]=0;
	fscanf(tsp_file,"%s", buf);
    }


    if( strcmp("NODE_COORD_SECTION", buf) == 0 ){
	trace_print("found section contaning the node coordinates\n");
	    }
    else{
	fprintf(stderr,"\n\nSome error ocurred finding start of coordinates from tsp file !!\n");
	exit(1);
    }

    if( (nodeptr = malloc(sizeof(struct point) * n)) == NULL )
	exit(EXIT_FAILURE);
    else {
	for ( i = 0 ; i < n ; i++ ) {
	    fscanf(tsp_file,"%ld %lf %lf", &j, &nodeptr[i].x, &nodeptr[i].y );
	}
    }
    trace_print("number of cities is %ld\n",n);
    trace_print("\n... done\n");
	return (nodeptr);
}

void
read_instance (const char* filename, struct problem *instance)
{
    instance->nodeptr = read_etsp(filename);
    printf("calculating distance matrix ..\n\n");
    instance->distance = compute_distances();
    printf(" .. done\n");
}

void free_instance (struct problem *instance)
{
    free( instance->distance );
    free( instance->nn_list );
}

const char * get_instance_name(const struct problem *instance)
{
    return instance->name;
}

void printHeur(void)
/*    
      FUNCTION:       print heuristic information 
      INPUT:          none
      OUTPUT:         none
*/
{
  long int i, j;

  printf("Heuristic information:\n");
  for ( i = 0 ; i < n ; i++) {
    printf("From %ld:  ",i);
    for ( j = 0 ; j < n - 1 ; j++ ) {
      printf(" %.3f ", HEURISTIC(i,j));
    }
    printf(" %.3f\n", HEURISTIC(i,j));
    printf("\n");
  }
  printf("\n");
}

void printDist(void) 
/*    
      FUNCTION:       print distance matrix 
      INPUT:          none
      OUTPUT:         none
*/
{
  long int i,j;

  printf("Distance Matrix:\n");
  for ( i = 0 ; i < n ; i++) {
    printf("From %ld:  ",i);
    for ( j = 0 ; j < n - 1 ; j++ ) {
      printf(" %ld", instance.distance[i][j]);
    }
    printf(" %ld\n", instance.distance[i][n-1]);
    printf("\n");
  }
  printf("\n");
}



