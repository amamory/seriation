#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

typedef struct graph_ptr{
	short cur;
	short ini;
} graph_ptr;

typedef struct args{
	short plot;
	short video;
	float alpha;
	char *file_name;
	char *order_file;
	float cooling_factor;
	short cooling_interval;
	float percentual_energy;
	unsigned short mc_steps;
	double seed;
} args;

struct timespec diff(struct timespec start, struct timespec end){
	struct timespec diff_time;

	if((end.tv_nsec - start.tv_nsec) < 0){
		diff_time.tv_sec = end.tv_sec - start.tv_sec - 1;
		diff_time.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
	}else{
		diff_time.tv_sec = end.tv_sec - start.tv_sec;
		diff_time.tv_nsec = end.tv_nsec - start.tv_nsec;
	}
	return diff_time;
}

void print_help(){
	printf("\nUsage: ./cfm_seriation [OPTION...]\n");
	printf("\n");
	printf(" Seriation Parameters:\n");
	printf("   f=[NETWORK FILE]           Network file path name\n");
	printf("   o=[ORDER FILE]             Apply initial order\n");
	printf("   i=[INTERVAL]               Number of isothermal steps\n");
	printf("   m=[STEPS]                  Number of steps\n");
	printf("   c=[FACTOR]                 Cooling factor\n");
	printf("   a=[ALPHA]                  Alpha value\n");
	printf("   p=[PERCENTUAL]             Percentual energy for initial temperature\n");
	printf("   s=[SEDD]                   Random seed\n");
	printf("   P                          Plot graphs\n");
	printf("   v                          Generate video\n");
}

args args_parser(int argc, char *argv[]){
	args options;
	short x;

	options.plot = 0;
	options.video = 0;
	options.alpha = 1.0;
	options.mc_steps = 100;
	options.file_name = NULL;
	options.order_file = NULL;
	options.seed = time(NULL);
	options.cooling_factor = 0.4;
	options.cooling_interval = 100;
	options.percentual_energy = 0.01;

	for(x=1;x<argc;x++){
		switch(argv[x][0]){
			case 'f':
				options.file_name = &argv[x][2];
				break;
			case 'o':
				options.order_file = &argv[x][2];
				break;
			case 'i':
				options.cooling_interval = atoi(&argv[x][2]);
				break;
			case 'm':
				options.mc_steps = atoi(&argv[x][2]);
				break;
			case 'c':
				options.cooling_factor = atof(&argv[x][2]);
				break;
			case 'a':
				options.alpha = atof(&argv[x][2]);
				break;
			case 'p':
				options.percentual_energy = atof(&argv[x][2]);
				break;
			case 's':
				options.seed = atoi(&argv[x][2]);
				break;
			case 'P':
				options.plot = 1;
				break;
			case 'v':
				options.video = 1;
				break;
			case 'h':
				print_help();
				exit(1);
		}
	}
	return options;
}


void mat2str(char *matrix[], graph_ptr *order, unsigned short n_nodes){
	short x, y;
	short col, row;
	graph_ptr pos_i, pos_j;

	printf(" ");
	for(col=0;col<n_nodes;col++){
		if(col == (n_nodes-1))
			printf("%d\n",col);
		else
			printf("%d, ",col);
	}

	for(row=(n_nodes-1);row>=0;row--){
		printf("[");
		pos_i = order[row];
		for(col=0;col<n_nodes;col++){
			pos_j = order[col];
			if(col == (n_nodes-1))
				printf("%d] : %d\n",matrix[pos_i.cur][pos_j.cur], row);
			else
				printf("%d, ",matrix[pos_i.cur][pos_j.cur]);
		}
	}
}

void swap(graph_ptr *order, short a, short b){
	short aux_swap;

	aux_swap = order[order[a].cur].ini;
	order[order[a].cur].ini = order[order[b].cur].ini;
	order[order[b].cur].ini = aux_swap;

	aux_swap = order[a].cur;
	order[a].cur = order[b].cur;
	order[b].cur = aux_swap;
}

void get_size_column(char *matrix[], short *size_column, unsigned short n_nodes){
	short col, row;

	for(row=(n_nodes-1);row>=0;row--){
		size_column[row] = 0;
		for(col=0;col<n_nodes;col++){
			if(matrix[row][col] == 1){
				size_column[row]++;
			}
		}
	}
}

void compress_matrix(char *matrix[], short *ones_list[], unsigned short n_nodes){
	short col, row, aux_col;

	for(row=0;row<n_nodes;row++){
		aux_col = 0;
		for(col=0;col<n_nodes;col++){
			if(matrix[row][col] == 1){
				ones_list[row][aux_col] = col;
				aux_col++;
			}
		}
	}
}

double getPartEnergy(char *matrix[], graph_ptr *order, short a, short b, short *size_column, short *ones_list[], unsigned short n_nodes, float alpha){
	short index_a = a, index_b = b, init_a, init_b, row, col, index_cnt;
	double col_a_energy = 0, col_b_energy = 0, row_a_energy = 0, row_b_energy = 0, part_energy, abs;
	unsigned char index_a_cnt, index_b_cnt, neighbors;
	graph_ptr pos_i, pos_j;

	if(a > b){
		index_a = b;
		index_b = a;
	}

	init_a = index_a;
	init_b = index_b;

	if((index_a+1) == index_b){
		if(index_a == 0){
			index_a_cnt = 2;
		}else{
			index_a_cnt = 3;
			index_a = index_a-1;
		}

		if(index_b == (n_nodes-1)){
			index_b_cnt = 0;
		}else{
			index_b_cnt = 1;
			index_b = index_b+1;
		}

	}else if((index_a+1) == (index_b-1)){
		if(index_a == 0){
			index_a_cnt = 2;
		}else{
			index_a_cnt = 3;
			index_a = index_a-1;
		}

		if(index_b == (n_nodes-1))
			index_b_cnt = 1;
		else
			index_b_cnt = 2;

	}else{
		if(index_a == 0){
			index_a_cnt = 2;
		}else{
			index_a_cnt = 3;
			index_a = index_a-1;
		}

		if(index_b == (n_nodes-1))
			index_b_cnt = 2;
		else
			index_b_cnt = 3;

		index_b = index_b-1;
	}

	index_cnt = index_a + index_a_cnt;
	/* column a */
	for(col=index_a;col<index_cnt;col++){
		part_energy = 0;
		pos_i = order[col];
		for(row=0;row<size_column[pos_i.cur];row++){
			pos_j = order[ones_list[pos_i.cur][row]];
			if(col < pos_j.ini){
				neighbors = 0;
				if(pos_j.ini != 0)
					neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
				if(pos_j.ini != (n_nodes-1))
					neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
				if(col != 0)
					neighbors = neighbors + (1 - matrix[order[col-1].cur][order[pos_j.ini].cur]);
				if(col != (n_nodes-1))
					neighbors = neighbors + (1 - matrix[order[col+1].cur][order[pos_j.ini].cur]);

				abs = pow((pos_j.ini - col), alpha);
				part_energy = part_energy + (abs * neighbors);
			}
		}
		col_a_energy = col_a_energy + part_energy;
	}

	/* row a */
	for(row=(index_cnt-1);row>=index_a;row--){
		part_energy = 0;
		pos_i = order[row];
		for(col=0;col<size_column[pos_i.cur];col++){
			pos_j = order[ones_list[pos_i.cur][col]];
			if(row > pos_j.ini){
				if((pos_j.ini!= init_a-1) && (pos_j.ini != init_a) && (pos_j.ini != init_a+1) && (pos_j.ini != init_b-1) && (pos_j.ini != init_b) && (pos_j.ini != init_b+1)){
					neighbors = 0;
					if(pos_j.ini != 0)
						neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
					if(pos_j.ini != (n_nodes-1))
						neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
					if(row != 0)
						neighbors = neighbors + (1 - matrix[order[row-1].cur][order[pos_j.ini].cur]);
					if(row != (n_nodes-1))
						neighbors = neighbors + (1 - matrix[order[row+1].cur][order[pos_j.ini].cur]);

					abs = pow((row - pos_j.ini), alpha);
					part_energy = part_energy + (abs * neighbors);
				}
			}
		}
		row_a_energy = row_a_energy + part_energy;
	}

	if(index_b_cnt != 0){
		index_cnt = index_b + index_b_cnt;
		/* column b */
		for(col=index_b;col<index_cnt;col++){
			part_energy = 0;
			pos_i = order[col];
			for(row=0;row<size_column[pos_i.cur];row++){
				pos_j = order[ones_list[pos_i.cur][row]];
				if(col < pos_j.ini){
					neighbors = 0;
					if(pos_j.ini != 0)
						neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
					if(pos_j.ini != (n_nodes-1))
						neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
					if(col != 0)
						neighbors = neighbors + (1 - matrix[order[col-1].cur][order[pos_j.ini].cur]);
					if(col != (n_nodes-1))
						neighbors = neighbors + (1 - matrix[order[col+1].cur][order[pos_j.ini].cur]);

					abs = pow((pos_j.ini - col), alpha);
					part_energy = part_energy + (abs * neighbors);
				}
			}
			col_b_energy = col_b_energy + part_energy;
		}

		/* row b */
		for(row=(index_cnt-1);row>=index_b;row--){
			part_energy = 0;
			pos_i = order[row];
			for(col=0;col<size_column[pos_i.cur];col++){
				pos_j = order[ones_list[pos_i.cur][col]];
				if(row > pos_j.ini){
					if((pos_j.ini!= init_b-1) && (pos_j.ini != init_b) && (pos_j.ini != init_b+1) && (pos_j.ini != init_a-1) && (pos_j.ini != init_a) && (pos_j.ini != init_a+1)){
						neighbors = 0;
						if(pos_j.ini != 0)
							neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
						if(pos_j.ini != (n_nodes-1))
							neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
						if(row != 0)
							neighbors = neighbors + (1 - matrix[order[row-1].cur][order[pos_j.ini].cur]);
						if(row != (n_nodes-1))
							neighbors = neighbors + (1 - matrix[order[row+1].cur][order[pos_j.ini].cur]);

						abs = pow((row - pos_j.ini), alpha);
						part_energy = part_energy + (abs * neighbors);
					}
				}
			}
			row_b_energy = row_b_energy + part_energy;
		}
	}
	return (col_a_energy + col_b_energy + row_a_energy + row_b_energy);
}

double getMatEnergy(char *matrix[], graph_ptr *order, short *size_column, short *ones_list[], unsigned short n_nodes, float alpha){
	double energy = 0, row_energy, abs;
	unsigned char neighbors;
	short row, col;
	graph_ptr pos_i, pos_j;

	for(row=(n_nodes-1);row>=0;row--){
		row_energy = 0;
		pos_i = order[row];
		for(col=0;col<size_column[pos_i.cur];col++){
			pos_j = order[ones_list[pos_i.cur][col]];
			if(row < pos_j.ini){
				neighbors = 0;
				if(pos_j.ini != 0)
					neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini-1].cur]);
				if(pos_j.ini != (n_nodes-1))
					neighbors = neighbors + (1 - matrix[pos_i.cur][order[pos_j.ini+1].cur]);
				if(row != 0)
					neighbors = neighbors + (1 - matrix[order[row-1].cur][order[pos_j.ini].cur]);
				if(row != (n_nodes-1))
					neighbors = neighbors + (1 - matrix[order[row+1].cur][order[pos_j.ini].cur]);

				abs = pow((pos_j.ini - row), alpha);
				row_energy = row_energy + (abs * neighbors);
			}
		}
		energy = energy + row_energy;
	}
	return energy;
}

void buildMat(FILE *inputfile, char ***matrix, char **nodes, unsigned short n_nodes){
	char line[80], *line_ptr, *node_left, *node_right;
	short x, y, last_node = 0, node_index_left, node_index_right;
	char **alloc_matrix = (char**) malloc(n_nodes * sizeof(char*));
	for(x=0;x<n_nodes;x++)
		alloc_matrix[x] = (char *) calloc(n_nodes,sizeof(char));

	while(fgets(line, 80, inputfile) != NULL){
		line_ptr = line;
		node_left = strsep(&line_ptr," \t\n");
		node_right = strsep(&line_ptr," \t\n");
		if((node_left != NULL) && (node_right != NULL)){
			node_index_left = -1;
			node_index_right = -1;
			for(x=0;x<n_nodes;x++){
				if(!strcmp(nodes[x],node_left))
					node_index_left = x;
				if(!strcmp(nodes[x],node_right))
					node_index_right = x;
				if((node_index_left>0) && (node_index_right>0))
					break;
			}
		}
		alloc_matrix[node_index_left][node_index_right] = 1;
		alloc_matrix[node_index_right][node_index_left] = 1;
	}
	*matrix = alloc_matrix;
}

void applyOrderFile(FILE *inputfile, char **nodes){
	char line[80], *line_ptr, *protein, *pos_ptr, *end_ptr;
	short protein_length;
	int pos;

	while(fgets(line, 80, inputfile) != NULL){
		line_ptr = line;
		protein = strsep(&line_ptr," \t\n");
		protein_length = strlen(protein) + 1;
		pos_ptr = strsep(&line_ptr," \t\n");

		if((protein != NULL) && (pos_ptr != NULL)){
			if(strcmp(protein,"Protein")){
				pos = (int) strtol(pos_ptr, &end_ptr,10);
				if(pos_ptr == end_ptr){
					printf("Invalid order file format.\n");
					exit(1);
				}else{
					nodes[pos] = (char*) calloc(protein_length,sizeof(char));
					strcpy(nodes[pos],protein);
				}
			}
		}else{
			printf("ERROR while reading order file.\n");
			exit(1);
		}
	}
}

unsigned short readProteinList(FILE *inputfile, char ***protein_list){
	char line[80], *line_ptr, *node_left, *node_right;
	short x, y, last_node = 0, node_left_length, node_right_length;
	unsigned short left_found = 0, right_found = 0, n_nodes = 0;
	char **alloc_nodes = NULL;
	unsigned int n_edges = 0;

	while(fgets(line, 80, inputfile) != NULL){
		line_ptr = line;
		node_left = strsep(&line_ptr," \t\n");
		node_left_length = strlen(node_left) + 1;
		node_right = strsep(&line_ptr," \t\n");
		node_right_length = strlen(node_right) + 1;
		if((node_left != NULL) && (node_right != NULL)){
			if(n_nodes == 0){
				alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
				alloc_nodes[n_nodes] = (char *) calloc(node_left_length,sizeof(char));
				strcpy(alloc_nodes[n_nodes],node_left);
				left_found = 1;
				n_nodes++;
				alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
				alloc_nodes[n_nodes] = (char *) calloc(node_right_length,sizeof(char));
				strcpy(alloc_nodes[n_nodes],node_right);
				right_found = 1;
				n_nodes++;
			}else{
				if(!strcmp(alloc_nodes[last_node],node_left)){
					right_found = 0;
					for(x=0;x<n_nodes;x++){
						if(!strcmp(alloc_nodes[x],node_right)){
							right_found = 1;
							break;
						}else if(x == (n_nodes-1)){
							alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
							alloc_nodes[n_nodes] = (char *) calloc(node_right_length,sizeof(char));
							strcpy(alloc_nodes[n_nodes],node_right);
							right_found = 1;
							n_nodes++;
						}
					}
				}else{
					left_found = 0;
					right_found = 0;
					for(x=0;x<n_nodes;x++){
						if(!strcmp(alloc_nodes[x],node_left) && !left_found)
							left_found = 1;
						if(!strcmp(alloc_nodes[x],node_right) && !right_found)
							right_found = 1;
						if(left_found && right_found)
							break;
						if(x == (n_nodes-1)){
							if(!left_found){
								alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
								alloc_nodes[n_nodes] = (char *) calloc(node_left_length,sizeof(char));
								strcpy(alloc_nodes[n_nodes],node_left);
								last_node = n_nodes;
								n_nodes++;
							}
							if(!right_found){
								alloc_nodes = (char **) realloc(alloc_nodes,(n_nodes+1) * sizeof(char *));
								alloc_nodes[n_nodes] = (char *) calloc(node_right_length,sizeof(char));
								strcpy(alloc_nodes[n_nodes],node_right);
								n_nodes++;
							}
						}
					}
				}
			}
			n_edges++;
		}
	}
	printf("\tProteins: %d\n",n_nodes);
	printf("\tInteractions: %d\n",n_edges/2);

	*protein_list = alloc_nodes;
	return n_nodes;
}

void rand_values(short *random_list, unsigned short n_nodes){
	short x, y, value;
	unsigned short repeat;

	for(x=0;x<n_nodes;x++){
		repeat = 1;
		random_list[x] = -1;
		do{
			value = (int)(rand()/(RAND_MAX + 1.0) * n_nodes);
			for(y=0;y<=x;y++){
				if((random_list[y]) == value){
					break;
				}else if(y == x){
					random_list[y] = value;
					repeat = 0;
				}
			}
		}while(repeat);
	}
}

void progressBar(unsigned short step, unsigned short monte_carlo_steps){
	float ratio;
	unsigned short incomplete, x;

	ratio = (float)step/monte_carlo_steps;
	incomplete = ratio * 100;
	printf("%3d%% [", (int)(ratio*100) );

	for(x=0;x<incomplete;x++)
		printf("=");
	for(x=incomplete;x<100;x++)
		printf(" ");
	printf("]\r");
}

void buildFrame(graph_ptr *order, short *size_column, short *ones_list[], char *frame_path, char *gnuplot_path, short int monte_carlo_steps, unsigned short n_nodes, short frame_cnt){
	short row, col;
	graph_ptr pos_i, pos_j;
	char *frame_name = NULL, *command = NULL;

	FILE *frame_image;

	frame_name = (char*) malloc((strlen(frame_path) + 12)*sizeof(char));
	sprintf(frame_name,"%s/frame%d",frame_path,frame_cnt);
	command = (char*) malloc((strlen(frame_name) + (strlen(gnuplot_path) + strlen(frame_name)) + 143 + 58)*sizeof(char));
	sprintf(command,"gnuplot -p -e \"load '/usr/share/cfm-seriation/etc/matrix.gnu';set output '%s.png';filename='%s';set multiplot;p '%s.dat' w dots lt 0;set xrange [0:%d];load '/usr/share/cfm-seriation/etc/video.gnu'\"",frame_name,gnuplot_path,frame_name,monte_carlo_steps);
	frame_name = (char*) realloc(frame_name,(strlen(frame_path) + 16)*sizeof(char));
	sprintf(frame_name,"%s/frame%d.dat",frame_path,frame_cnt);
	frame_image = fopen(frame_name,"w+");

	for(row=(n_nodes-1);row>=0;row--){
		pos_i = order[row];
		for(col=0;col<size_column[pos_i.cur];col++){
			pos_j = order[ones_list[pos_i.cur][col]];
			if(row < pos_j.ini){
				fprintf(frame_image,"%f\t%f\n",(float)pos_j.ini/n_nodes,(float)row/n_nodes);
				fprintf(frame_image,"%f\t%f\n",(float)row/n_nodes,(float)pos_j.ini/n_nodes);
			}
		}
	}

	fflush(frame_image);
	fclose(frame_image);
	system(command);
	sprintf(command,"rm %s",frame_name);
	system(command);
}

void print_config(){
	printf("\nUsage: cfm-seriation [OPTION...]\n");
	printf("\n");
	printf(" Seriation Parameters:\n");
	printf("   f=[NETWORK FILE].dat       Network file path name\n");
	printf("   o=[ORDER FILE].dat         Apply initial order\n");
	printf("   i=[INTERVAL]               Number of isothermal steps\n");
	printf("   m=[STEPS]                  Number of steps\n");
	printf("   c=[FACTOR]                 Cooling factor\n");
	printf("   a=[ALPHA]                  Alpha value\n");
	printf("   p=[PERCENTUAL]             Percentual energy for initial temperature\n");
	printf("   s=[SEDD]                   Random seed\n");
	printf("   P                          Plot graphs\n");
	printf("   v                          Generate video\n");
}

/* ****************************************************************************************** */
/* -------------------------------------- MAIN PROGRAM -------------------------------------- */
/* ****************************************************************************************** */
int main (int argc, char *argv[]){
	args options;
	unsigned short step, n_nodes;
	graph_ptr *graph_order = NULL, pos_i, pos_j;
	short n, a, b, x, y, *random_list = NULL, **ones_list = NULL, *size_column = NULL;
	double initial_energy, curr_energy, new_energy, final_energy, delta_energy, curr_part_energy, new_part_energy;
	double temperature, expo, swap_cnt = 0;
	long int bar_char, ext_str;
	char **matrix = NULL, **protein_list = NULL, **nodes = NULL;
	FILE *datfile, *order, *gnuplot, *energy, *result;
	char *dir_name = NULL, *gnuplot_path = NULL, *order_path = NULL, *result_path = NULL, *command = NULL, *dat_name = NULL, *path_name = NULL, *frame_path = NULL, time_name[20];
	struct tm *loctime;
	struct timespec init_time, final_time, elapsed_time;
	time_t curtime;

	options = args_parser(argc, argv);
	print_config();

	if((options.file_name != NULL) && (strstr(options.file_name,".dat") != NULL)){
		/* Check if file exists. */
		datfile = fopen(options.file_name,"r");
		if(datfile == NULL){
			printf("ERROR! Invalid file path \"%s\".\n",options.file_name);
			exit(1);
		}
		/* Get the current time. */
		curtime = time(NULL);

		/* Convert it to local time representation. */
		loctime = localtime(&curtime);
		strftime(time_name, 21, "%Y-%m-%d_%Hh%Mm%Ss",loctime);

		/* Extract path and file name from argument. */
		if(strrchr(options.file_name,'/') == NULL){
			ext_str = strstr(options.file_name,".dat") - options.file_name;
			dat_name = (char*) malloc((ext_str + 1)*sizeof(char));
			strncpy(dat_name,options.file_name,ext_str);
			dat_name[ext_str] = '\0';
		}else{
			bar_char = strrchr(options.file_name,'/') - options.file_name;
			ext_str = strstr(options.file_name,".dat") - options.file_name;
			dat_name = (char*) malloc((ext_str - bar_char)*sizeof(char));
			strncpy(dat_name,&options.file_name[bar_char + 1],(ext_str - bar_char - 1));
			dat_name[(ext_str - bar_char - 1)] = '\0';
			// future use (deb package)
			path_name = (char*) malloc((bar_char + 1)*sizeof(char));
			strncpy(path_name,options.file_name,bar_char);
			path_name[bar_char] = '\0';
		}

		/* Create dirs with input file name. */
		dir_name = (char*) malloc((strlen(dat_name) + strlen(time_name) + 2)*sizeof(char));
		sprintf(dir_name,"%s_%s",dat_name,time_name);
		mkdir(dir_name, 0777);

		if(options.video == 1){
			frame_path = (char*) malloc((strlen(dir_name) + 8)*sizeof(char));
			sprintf(frame_path,"%s/frames",dir_name);
			mkdir(frame_path, 0777);
		}

		/* Initialize results file. */
		result_path = (char*) malloc((strlen(dir_name) + 12)*sizeof(char));
		sprintf(result_path,"%s/result.txt",dir_name);
		result = fopen(result_path,"w+");
		if(result == NULL){
			printf("ERROR while creating %s, please check!\n",result_path);
			exit(1);
		}
		fputs("\n**** TEST Command ****\n",result);
		for(x=0;x<argc;x++){
			fputs(argv[x],result);
			fputs(" ",result);
		}
		fputs("\n",result);
		fputs("\n**** TEST Parameters ****\n",result);
		fprintf(result,"Network file path name: %s\n",options.file_name);
		if(options.order_file == NULL)
			fprintf(result,"Random seed: %.lf\n",options.seed);
		else
			fprintf(result,"Apply initial order: %s\n",options.order_file);
		fprintf(result,"Number of isothermal steps: %d\n",options.cooling_interval);
		fprintf(result,"Number of steps: %d\n",options.mc_steps);
		fprintf(result,"Cooling factor: %f\n",options.cooling_factor);
		fprintf(result,"Alpha value: %f\n",options.alpha);
		fprintf(result,"Percentual energy for initial temperature: %f\n",options.percentual_energy);

		/* Apply random order or initial order file and build matrix. */
		datfile = fopen(options.file_name,"r");
		printf("\nReading file...\n");
		n_nodes = readProteinList(datfile, &protein_list);
		rewind(datfile);

		graph_order = (graph_ptr*) malloc(n_nodes*sizeof(graph_ptr));
		nodes = (char**) malloc(n_nodes*sizeof(char*));

		if(options.order_file != NULL){
			order = fopen(options.order_file, "r");
			if(order == NULL){
				printf("ERROR! File %s does not exist.\n",options.order_file);
				exit(1);
			}
			printf("Applying initial order...\n");
			applyOrderFile(order, nodes);
			fclose(order);
			for(n=0;n<n_nodes;n++){
				free(protein_list[n]);
				graph_order[n].cur = n;
				graph_order[n].ini = n;
			}
		}else{
			printf("Applying random order...\n");
			srand(options.seed);
			random_list = (short*) malloc(n_nodes*sizeof(short));
			rand_values(random_list, n_nodes);

			for(n=0;n<n_nodes;n++){
				nodes[random_list[n]] = protein_list[n];
				graph_order[n].cur = n;
				graph_order[n].ini = n;
			}
			free(random_list);
		}
		buildMat(datfile, &matrix, nodes, n_nodes);
		free(protein_list);
		fclose(datfile);

		/* Create a list with number of elements of each column. */
		size_column = (short*) malloc(n_nodes*sizeof(short));
		get_size_column(matrix, size_column, n_nodes);

		/* Create a list representation to compress the matrix. */
		ones_list = (short**) malloc(n_nodes*sizeof(short*));
		for(n=0;n<n_nodes;n++){
			ones_list[n] = (short*) malloc(size_column[n]*sizeof(short));
		}
		compress_matrix(matrix, ones_list, n_nodes);

		/* Generate initial order file. */
		if(options.plot == 1)
			printf("Saving and plotting initial order...\n");
		else
			printf("Saving initial order...\n");
		order_path = (char*) malloc((strlen(dir_name) + 19)*sizeof(char));
		gnuplot_path = (char*) malloc((strlen(dir_name) + 13)*sizeof(char));
		sprintf(order_path,"%s/initial_order.dat",dir_name);
		order = fopen(order_path,"w+");
		if(order == NULL){
			printf("ERROR! File %s does not exist.\n",order_path);
			exit(1);
		}

		if(options.plot == 1){
			sprintf(gnuplot_path,"%s/initial",dir_name);
			command = (char*) malloc((strlen(dir_name) + (strlen(gnuplot_path) *2) + 125 + 29)*sizeof(char));
			sprintf(command,"gnuplot -p -e \"load '/usr/share/cfm-seriation/etc/matrix.gnu'; set output '%s.png';p '%s.dat' w dots lt 0\"",gnuplot_path,gnuplot_path);
			sprintf(gnuplot_path,"%s/initial.dat",dir_name);

			gnuplot = fopen(gnuplot_path,"w+");
			if(gnuplot == NULL){
				printf("ERROR! File %s does not exist.\n",gnuplot_path);
				exit(1);
			}
		}

		fprintf(order,"Protein\tdim1\n");
		for(x=(n_nodes-1);x>=0;x--){
			pos_i = graph_order[x];
			fprintf(order,"%s\t%d\n",nodes[pos_i.cur],pos_i.cur);
			if(options.plot == 1){
				for(y=0;y<size_column[pos_i.cur];y++){
					pos_j = graph_order[ones_list[pos_i.cur][y]];
					if(x < pos_j.ini){
						fprintf(gnuplot,"%f\t%f\n",(float)pos_j.ini/n_nodes,(float)x/n_nodes);
						fprintf(gnuplot,"%f\t%f\n",(float)x/n_nodes,(float)pos_j.ini/n_nodes);
					}
				}
			}
		}

		if(options.plot == 1){
			fflush(gnuplot);
			fclose(gnuplot);
			system(command);
		}
		fflush(order);
		fclose(order);
	}else{
		printf("ERROR! Invalid file name.\n");
		exit(1);
	}

	initial_energy = getMatEnergy(matrix, graph_order, size_column, ones_list, n_nodes, options.alpha);
	printf("INITIAL Energy: %.f\n",initial_energy*2);
	curr_energy = initial_energy;
	temperature = curr_energy * options.percentual_energy;

	/* Initialize energy file. */
	sprintf(gnuplot_path,"%s/energy.dat",dir_name);
	energy = fopen(gnuplot_path,"w+");
	if(energy == NULL){
		printf("ERROR! File %s does not exist.\n",gnuplot_path);
		exit(1);
	}
	srand(options.seed);
	fprintf(energy,"steps\tenergy\ttemperature\tswaps\tH/H0\tT/H0\n");
	fprintf(energy,"0\t%.f\t%.3lE\t%.lf\t%.3f\t%.3lE\n",curr_energy*2, temperature*2, swap_cnt, (curr_energy/initial_energy), (temperature/initial_energy));
	fflush(energy);

	if(options.video == 1)
		buildFrame(graph_order, size_column, ones_list, frame_path, gnuplot_path, options.mc_steps, n_nodes, 0);

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &init_time);

	printf("Ordering...\n");
	for(step=1;step<=options.mc_steps;step++){
		for(n=0;n<n_nodes;n++){
			/* Get random nodes. */
			a = (int)(rand()/(RAND_MAX + 1.0) * n_nodes);
			do{
				b = (int)(rand()/(RAND_MAX + 1.0) * n_nodes);
			}while(a == b);

			/* Calculate current partial energy. */
			curr_part_energy = getPartEnergy(matrix, graph_order, a, b, size_column, ones_list, n_nodes, options.alpha);

			/* Swap chosen nodes. */
			swap(graph_order, a, b);

			/* Calculate new partial energy. */
			new_part_energy = getPartEnergy(matrix, graph_order, a, b, size_column, ones_list, n_nodes, options.alpha);

			/* Calculate new energy. */
			new_energy = curr_energy - curr_part_energy + new_part_energy;

			/* Generate delta. */
			delta_energy = new_energy - curr_energy;

			/* If the new energy is better accept and assign. */
			if(delta_energy <= 0){
				curr_energy = new_energy;
				swap_cnt++;
			}else{
				/* If new distance is worse accept but with a probability level. */
				expo = expl(-delta_energy/temperature);

				if(expo > (rand()/RAND_MAX)){
					curr_energy = new_energy;
					swap_cnt++;
				}else{
					swap(graph_order, a, b);
				}
			}
		}

		/* Cooling process. */
		if(step%options.cooling_interval == 0)
			temperature = temperature * options.cooling_factor;

		fprintf(energy,"%d\t%.f\t%.3lE\t%.lf\t%.3f\t%.3lE\n",step, curr_energy*2, temperature*2, swap_cnt, (curr_energy/initial_energy), (temperature/initial_energy));
		fflush(energy);

		if(options.video == 1)
			buildFrame(graph_order, size_column, ones_list, frame_path, gnuplot_path, options.mc_steps, n_nodes, step);
		progressBar(step,options.mc_steps);
		fflush(stdout);
	}

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &final_time);
	elapsed_time = diff(init_time, final_time);
	final_energy = curr_energy;
	printf("\nFINAL Energy: %.f\n",final_energy*2);

	/* Plot and save final results. */
	if(options.plot == 1)
		printf("Saving and plotting final order...\n");
	else
		printf("Saving final order...\n");
	fprintf(result,"\n**** TEST Results ****\nCPU Time: %lld.%ld\n", (long long)elapsed_time.tv_sec, elapsed_time.tv_nsec);
	fprintf(result,"INITIAL Energy: %.f\nFINAL Energy: %.f",initial_energy*2, final_energy*2);
	fclose(result);
	fflush(energy);
	fclose(energy);
	if(options.plot == 1){
		sprintf(gnuplot_path,"%s/energy",dir_name);
		if(options.mc_steps > 0){
			sprintf(command,"gnuplot -p -e \"filename='%s.dat';set output '%s.png';set xrange [0:%d];load '/usr/share/cfm-seriation/etc/energy.gnu'\"",gnuplot_path,gnuplot_path,options.mc_steps);
			system(command);
		}
		sprintf(gnuplot_path,"%s/final",dir_name);
		sprintf(command,"gnuplot -p -e \"load '/usr/share/cfm-seriation/etc/matrix.gnu'; set output '%s.png';p '%s.dat' w dots lt 0\"",gnuplot_path,gnuplot_path);
		sprintf(gnuplot_path,"%s/final.dat",dir_name);
		gnuplot = fopen(gnuplot_path,"w+");
		if(gnuplot == NULL){
			printf("ERROR! File %s does not exist.\n",gnuplot_path);
			exit(1);
		}
	}
	sprintf(order_path,"%s/final_order.dat",dir_name);
	order = fopen(order_path,"w+");
	if(order == NULL){
		printf("ERROR! File %s does not exist.\n",order_path);
		exit(1);
	}
	fprintf(order,"Protein\tdim1\n");
	for(x=(n_nodes-1);x>=0;x--){
		pos_i = graph_order[x];
		fprintf(order,"%s\t%d\n",nodes[pos_i.cur],pos_i.cur);
		if(options.plot == 1){
			for(y=0;y<size_column[pos_i.cur];y++){
				pos_j = graph_order[ones_list[pos_i.cur][y]];
				if(x < pos_j.ini){
					fprintf(gnuplot,"%f\t%f\n",(float)pos_j.ini/n_nodes,(float)x/n_nodes);
					fprintf(gnuplot,"%f\t%f\n",(float)x/n_nodes,(float)pos_j.ini/n_nodes);
				}
			}
		}
	}

	if(options.video == 1){
		if(options.plot == 0)
			command = (char*) malloc((strlen(frame_path) + (strlen(dir_name) *2) + 74)*sizeof(char));
		sprintf(command,"avconv -v quiet -r 10 -i %s/frame%s.png -r 10 -vcodec libx264 -b 3000k %s/%s.mp4",frame_path,"%d",dir_name,dir_name);
		system(command);
		sprintf(command,"rm -rf %s",frame_path);
		system(command);
	}

	if(options.plot == 1){
		fflush(gnuplot);
		fclose(gnuplot);
	}
	fflush(order);
	fclose(order);
	system(command);

	for(n=0;n<n_nodes;n++){
		free(nodes[n]);
		free(matrix[n]);
		free(ones_list[n]);
	}
	free(nodes);
	free(matrix);
	free(ones_list);
	free(graph_order);

	printf("Done!\n");
	return 0;
}
