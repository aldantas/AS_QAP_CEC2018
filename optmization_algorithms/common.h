#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef int*   type_vector;
typedef long** type_matrix;

char *input_file = NULL;
char *output_file = NULL;
int seed = 0;
long max_iterations;
bool is_iteration_set = false;
bool has_zero_matrix;


void copy_vector(type_vector a, type_vector b, long n)
{
	for(int i = 0; i < n; i++) {
		a[i] = b[i];
	}
}

void copy_matrix(type_matrix a, type_matrix b, long n)
{
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			a[i][j] = b[i][j];
		}
	}
}

/*Reads the parameter for the problem*/
void readParameter(int argc, char **argv)
{
	if(argc < 3){
		printf("Usage:\n %s -f [input file path]\n", argv[0]);
		exit(-1);
	}

	int i=1;
	if(argv[i][0] != '-'){
		printf("Parameter %s not defined.\n", argv[i]);
		exit(-1);
	} else if(argv[i][1] != 'f'){
		printf("Expected first parameter -f [input file path], found %s .\n", argv[i]);
		exit(-1);
	}

	i++;
	input_file = argv[i]; //input file path

	for(i=3;i<argc;i++){
		if(argv[i][0] == '-'){
			switch(argv[i][1]){
				case 's'://random seed
					i++;
					seed = atoi(argv[i]);
					break;
				case 'i'://max iterations
					i++;
					max_iterations = atoi(argv[i]);
					is_iteration_set = true;
					break;
				case 'o'://output file path
					i++;
					output_file = argv[i];
					break;
				default:
					printf("Parameter %s not defined.\n", argv[i]);
					exit(-1);
			}
		}else{
			printf("Parameter %s not defined.\n", argv[i]);
			exit(-1);
		}
	}
	if(!seed) seed = time(NULL);
}

int load_problem(int &n, type_matrix &a, type_matrix &b)
{
	FILE *file = fopen(input_file, "r");

	if(file == NULL)
		return -1;

	int temp;
	//reads the number of facilities
	temp = fscanf(file, "%i\n\n", &n);

	a = new long* [n];
	b = new long* [n];
	for (int i = 0; i < n; i++) {
		a[i] = new long[n];
		b[i] = new long[n];
	}

	has_zero_matrix = true;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			temp = fscanf(file, "%ld", &a[i][j]);
			if (a[i][j] != 0) {
				has_zero_matrix = false;
			}
		}
	}
	if(!has_zero_matrix) has_zero_matrix = true;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			temp = fscanf(file, "%ld", &b[i][j]);
			if (has_zero_matrix && a[i][j] != 0) {
				has_zero_matrix = false;
			}
		}
	}

	/* Just to avoid compilation warning */
	temp+=1;

	fclose(file);
	return 0;
}

int write_results(long cost, clock_t best_time, clock_t total_time, long
		best_evals, long total_evals)
{
	FILE *file = fopen(output_file, "wt");

	if(file == NULL)
		return -1;

	fprintf(file, "%ld\n", cost);
	fprintf(file, "%.2f\n", best_time/static_cast<double>(CLOCKS_PER_SEC));
	fprintf(file, "%ld\n", best_evals);
	fprintf(file, "%.2f\n", total_time/static_cast<double>(CLOCKS_PER_SEC));
	fprintf(file, "%ld\n", total_evals);

	fclose(file);
	return 0;
}
