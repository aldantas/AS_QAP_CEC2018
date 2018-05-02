#ifndef _PROBLEM_H_
#define _PROBLEM_H_
extern struct problem instance;
extern long int n; /* number of cities in the instance to be solved */

long int compute_tour_length(const long int *p );
long int compute_evaluation_function (long int *p);
void read_instance (const char* filename, struct problem *instance);
void free_instance (struct problem *instance);
const char * get_instance_name(const struct problem *instance);
long int ** compute_nn_lists (struct problem *instance);
int check_solution(const long int *t);

extern const char const * PROG_ID_STR;


#endif
