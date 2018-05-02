#include <stdio.h>
#include <stdlib.h>
#include <time.h>

char *input_file = NULL;
char *output_file = NULL;
char *seed = "0";
char *ls_choice = "2";
long max_iterations;
int is_iteration_set = 0;

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
	}else if(argv[i][1] != 'f'){
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
					seed = argv[i];
					break;
				case 'i'://max iterations
					i++;
					max_iterations = atoi(argv[i]);
					is_iteration_set = 1;
					break;
				case 'l'://local choice
					i++;
					ls_choice = argv[i];
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

int write_results(long cost, double best_time, double total_time)
{
	FILE *file = fopen(output_file, "wt");

	if(file == NULL)
		return -1;

	fprintf(file, "%ld\n", cost);
	fprintf(file, "%.2f\n", best_time);
	fprintf(file, "%ld\n", 0);
	fprintf(file, "%.2f\n", total_time);
	fprintf(file, "%ld\n", 0);

	fclose(file);
	return 0;
}
