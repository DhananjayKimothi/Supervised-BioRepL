//  non overlapping window for  training 
// overlapping window for testing 
// hardcoded -- The path of input fasta file, sequence vec, kmer vec and vocab
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h> // computation time

#define MaxSequences 200000 // ?
#define K 3 // kmer size
#define MAXLINE 10000  //? it was 128
#define SigLength 100
#define Window 1 // context size
#define min_alpha 0.005000 // 
#define n_inliers 2
#define n_outliers 2
#define gamma 0.5
#define r_families 5 // randomly selected families for each kmer except the one to which the sequence belongs.

int file_n;

clock_t start, end1, end2;
double training_time;
double inference_time;


typedef float real;


int* lengths;
int alphabet_size;


int num_sequences, num_kmers, num_families, num_nodes; // counts
char*sequences[MaxSequences]; // ?
int family_name[10000]; // 10
int family_size[10000]; // number of samples in a family
int family_start[10000]; // starting index of family
int alphabet_size; // 
int depth; // depth of the binary tree (herierchical softmax)
int last_family; // variable to store last family


int batch_words = 1000; // kmers in a batch after which step size is changed
int MaxNode = 0, MinNode = 1000000;  // variables for stor
int point_size = 100; // H-softmax parameter

double** seq_vector, **kmer_vector, **vocab_vector; // store vectors
int** codes, **points; // To store the binary code and the path for herierchical tree

int reverse[256]; // mapping kmer to dictionary index



char* path_save_model;
char* path_file_to_be_processed; // fasta file to be processed
char* path_save_vectors; 
char* vocab_path; // python generated file, contains the initialization vectors, binary code, path (from herierchical-softmax) for kmers 
char* doctag_path; // python generated file, contains the initialization vectors for the sequences
char* path; 

void ReadFASTA(char *filename)
{

	printf("\treading %s...\n", filename);
	num_families = 0;
	num_sequences = 0;
	last_family = -1;
	char  line[MAXLINE];
	FILE *fp = fopen(filename, "r");
	fgets(line, MAXLINE, fp);


	while (1)
	{
		if (line[0] != '>')
			break;

		int family = atoi(strrchr(line, '_') + 1);

		char* seq_start = NULL;
		int seq_length = 0;

		while (fgets(line, MAXLINE, fp))
		{
			if (line[0] == '>')
				break;

			int line_length = strlen(line) - 1;
			line[line_length] = '\0';

			seq_start = realloc(seq_start, seq_length + line_length + 1);
			strcpy(seq_start + seq_length, line);
			seq_length += line_length;
		}

		if (family != last_family)
		{
			family_name[num_families] = family;
			family_start[num_families] = num_sequences;
			family_size[num_families] = 1;
			last_family = family;
			num_families++;
		}
		else
			family_size[num_families - 1]++;

		sequences[num_sequences++] = seq_start;

	}
	fclose(fp);
}
void readfile_doctag(char* filename)
{
	double number;
	char  line[MAXLINE];
	FILE *fp = fopen(filename, "r");
	char *token;
	double val;
	int j = 0;
	seq_vector = (double**)malloc(sizeof(double*) * num_sequences);
	fgets(line, MAXLINE, fp);

	while (fgets(line, MAXLINE, fp))
	{
		seq_vector[j] = (double*)malloc(sizeof(double) * SigLength);
		if (line[0] == '>')
			continue;
		else
		{


			token = strtok(line, " ");
			int i = 0;

			while (token != NULL)
			{

				val = atof(token); // char to float
								   //printf("\n%f", val);

								   //printf("i and j (%d, %d)", j, i);
				seq_vector[j][i] = val;
				i = i++;

				token = strtok(NULL, " ");

			}

			j++;
		}
	}
	fclose(fp);

}



int RandomInt(int max)
{
	return max * ((double)rand() / RAND_MAX);
}
double RandomDouble()
{
	return ((double)rand() / RAND_MAX) - 0.5;
}

void readfile_vocab(char* filename)
{
	FILE *fp = fopen(filename, "r");

	char  kmerLine[MAXLINE];
	char  codeLine[MAXLINE];
	char  codeLineCopy[MAXLINE];
	char  pointsLine[MAXLINE];
	char  vectorLine[MAXLINE];


	num_kmers = pow(alphabet_size, K);
	codes = (int**)malloc(sizeof(int*) * num_kmers);
	points = (int**)malloc(sizeof(int*) * num_kmers);
	kmer_vector = (double**)malloc(sizeof(double*) * num_kmers);
	lengths = (int*)malloc(sizeof(int) * num_kmers); // length of path for kmers

													 // fill with random numbers
	for (int kmer = 0; kmer < num_kmers; kmer++)
	{
		kmer_vector[kmer] = (double*)malloc(sizeof(double)*SigLength);
		for (int i = 0; i < SigLength; i++)
		{
			kmer_vector[kmer][i] = 0.0;
			//kmer_vector[kmer][i] = RandomDouble() / SigLength;
		}
	}


	// fill with random numbers
	for (int kmer = 0; kmer < num_kmers; kmer++)
	{
		points[kmer] = (int*)malloc(sizeof(int) * 100);
		for (int i = 0; i < point_size; i++)
		{
			points[kmer][i] = 10;
		}
	}

	//fill with random numebrs

	for (int kmer = 0; kmer < num_kmers; kmer++)
	{
		codes[kmer] = (int*)malloc(sizeof(int) * 100);
		for (int i = 0; i < point_size; i++)
		{
			codes[kmer][i] = 0;
		}
	}

	// fill with random numebrs
	for (int kmer = 0; kmer < num_kmers; kmer++)
	{
		lengths[kmer] = 100;
	}


	while (fgets(kmerLine, MAXLINE, fp))
	{
		int index = KmerIndex(kmerLine + 1);

		fgets(codeLine, MAXLINE, fp);
		strcpy(codeLineCopy, codeLine);

		char* token = strtok(codeLineCopy, " ");
		int Length = 0;
		while (token)
		{
			Length++;
			token = strtok(NULL, " ");
		}

		lengths[index] = Length;
		codes[index] = (int*)malloc(sizeof(int) * Length);
		points[index] = (int*)malloc(sizeof(int) * Length);
		kmer_vector[index] = (double*)malloc(sizeof(double) * SigLength);

		token = strtok(codeLine, " ");
		for (int i = 0; i < Length; i++)
		{
			codes[index][i] = atoi(token);
			token = strtok(NULL, " ");
		}

		fgets(pointsLine, MAXLINE, fp);
		token = strtok(pointsLine, " ");
		for (int i = 0; i < Length; i++)
		{
			points[index][i] = atoi(token);
			if (points[index][i] > MaxNode)
				MaxNode = points[index][i];
			if (points[index][i] < MinNode)
				MinNode = points[index][i];
			token = strtok(NULL, " ");
		}

		fgets(vectorLine, MAXLINE, fp);
		token = strtok(vectorLine, " ");
		for (int i = 0; i < SigLength; i++)
		{
			kmer_vector[index][i] = atof(token);
			token = strtok(NULL, " ");
		}
	}
	num_nodes = MaxNode - MinNode + 1;
	vocab_vector = (double**)malloc(sizeof(double*) * num_nodes);
	for (int j = 0; j < num_nodes; j++)

	{
		vocab_vector[j] = (double*)malloc(sizeof(double)*SigLength);
		for (int i = 0; i < SigLength; i++)
			vocab_vector[j][i] = 0.0;
	}

	//printf("%d - %d\n", MinNode, MaxNode);
	fclose(fp);
}
void SaveModel(char* path)
{
	char filename[256];
	sprintf(filename, "%s", path);
	int Leng;
	//sprintf(filename, "%s\\Model_SuV%d.binary", path, file_n);
	//sprintf(filename, "%s\\Model1%d.binary", path, test);
	printf("\twriting %s...\n", filename);
	FILE* file = fopen(filename, "wb");

	fwrite(&num_kmers, sizeof(int), 1, file);
	fwrite(&num_nodes, sizeof(int), 1, file);
	fwrite(lengths, sizeof(int), num_kmers, file);

	for (int i = 0; i < num_kmers; i++)
		fwrite(kmer_vector[i], sizeof(double), SigLength, file);

	for (int i = 0; i < num_nodes; i++)
		fwrite(vocab_vector[i], sizeof(double), SigLength, file);

	for (int i = 0; i < num_kmers; i++)
		//Leng = Lengths[i] // check 
		fwrite(codes[i], sizeof(int), lengths[i], file);

	for (int i = 0; i < num_kmers; i++)
		//Leng = Lengths[i] // check 
		fwrite(points[i], sizeof(int), lengths[i], file);
	fclose(file);
}
void SaveSequences(char* path)
{
	char filename[256];
	sprintf(filename, "%s", path);
	
	
	//if (Training)
		//sprintf(filename, "%s\\t_SuV%d.binary", path, file_n);
	//sprintf(filename, "%s\\Database1%d.binary", path, test);
	//else
		//sprintf(filename, "%s\\DB_SuV%d.binary", path, file_n);
	//sprintf(filename, "%s\\Test1%d.binary", path, test);

	printf("\twriting %s...\n", filename);
	FILE* file = fopen(filename, "wb");

	fwrite(&num_families, sizeof(int), 1, file);
	fwrite(&family_name, sizeof(int), num_families, file);
	fwrite(&family_size, sizeof(int), num_families, file);

	fwrite(&num_sequences, sizeof(int), 1, file);
	for (int i = 0; i < num_sequences; i++)
		fwrite(seq_vector[i], sizeof(double), SigLength, file);

	fclose(file);
}
int KmerIndex(char* kmer)
{
	int index = 0;
	for (int i = 0; i < K; i++)
		index = index * alphabet_size + reverse[kmer[i]];
	return index;
}


void Process(int testNum, int Training)
{

	char filename[256];
	sprintf(filename, "%s\\out.txt", path);
	
	printf("test = %d, training = %d\n", testNum, Training);
	char fastaFilename[256]; // ?
	sprintf(fastaFilename, "%s", path_file_to_be_processed);

	int Iterations = 5;
	int Steps = 1;


	char* alphabet = "ABCDEFGHIKLMNOPQRSTUVWXYZ"; // 25 alphabets --- not so sure alphabets are also considered
												  //char* alphabet = "ACDEFGHIKLMNPQRSTUVWY";

	srand(42); // ? probably its like seed
	alphabet_size = strlen(alphabet); // string size : 25 in this case
	for (int i = 0; i < alphabet_size; i++)
		reverse[alphabet[i]] = i;

	ReadFASTA(fastaFilename); // reading fasta file .. see function

	printf("\t%s: %d Families, %d Sequences\n", fastaFilename, num_families, num_sequences);


	double next_alpha; // step size
	double next_alpha1;
	double alpha;
	int pushed_words = 0; // to track the progress the alpha depends on the progress
	int total_words = 0; // total number of words
	int total_pushed_words = 0;
	int pushed_examples = 0;
	int total_examples = num_sequences * Iterations;
	int document_based = 1;
	int word_based = 0;
	

	// read vocab file: includes kmer intialization, code and nodes traversed in the herierchical tree; read doctag file: initialization for sequence tags
	readfile_vocab(vocab_path);
	readfile_doctag(doctag_path);
	alpha = 0.025; // initialize alpha, decreases with the progress 
	
	for (int seq_idx = 0; seq_idx < num_sequences; seq_idx++)
		total_words += strlen(sequences[seq_idx]);
	total_words *= Iterations;
	next_alpha1 = alpha;
	
	//Inference stage
	double min_alphai = 0.0001;

	// Algorithm
	double neu2e_inlier[SigLength];
	double neu2e_outlier[SigLength];
	double neu2e[SigLength];
	double neu1e[SigLength];
	double Sig_seq_inlier = 0;
	int batch_size = 0;
	next_alpha = next_alpha1;
	int seq_idx_tst = 79;

	for (int iterations = 0; iterations < Iterations; iterations++) // iteration over multiple epochs
	{
		printf("iteration no. %d", iterations + 1);
		for (int family = 0; family < num_families; family++)
		{
			printf("Family processed %d", family);
			int start = family_start[family]; // starting index of family
			for (int family_idx = 0; family_idx < family_size[family]; family_idx++) // over sequences in the family
			{
				int seq_len = 0; // kmers in the sequence
				double alphai = 0.1; // for each sequence
				int seq_idx = start + family_idx; // global seq id
				//printf(" %d", seq_idx);

				char* seq = sequences[seq_idx];  // seq is a pointer which is pointing to the starting location of sequence			
												 // Non Overlap
				int seqlen = strlen(seq);
				for (int i = 0; i < K; i++)
					seq_len += floor((seqlen - i) / K);

				if (batch_size + seq_len <= batch_words)
					batch_size += seq_len;
				else
				{
					if (min_alpha < next_alpha && word_based)
					{
						double progress = (double)total_pushed_words / total_words;
						next_alpha1 = max(min_alpha, alpha - (alpha - min_alpha) * progress);
					}

					if (min_alpha < next_alpha && document_based)
					{
						double progress = (double)pushed_examples / total_examples;
						next_alpha1 = max(min_alpha, alpha - (alpha - min_alpha) * progress);
					}
					next_alpha = next_alpha1;
					pushed_words = 0;
					batch_size = 0;
					batch_size += seq_len;
					//printf("\n%f", next_alpha);
				}
				pushed_examples += 1;

				for (int slide = 0; slide < K; slide++)
				{
					pushed_words += floor((seqlen - slide) / K);
					total_pushed_words += floor((seqlen - slide) / K);
					int pos = 0;
					for (pos = pos + slide; pos < seqlen - K + 1; pos = pos + K) // iteration over a squence
					{
						char* kmer = seq + pos; // kmer is a pointing to the starting location of "to be predicted kmer"										
						double l1[SigLength]; // input signal
						for (int i = 0; i < SigLength; i++)
						{
							l1[i] = seq_vector[seq_idx][i];
						}

						int context[Window * 2]; // 
						int context_size = 0; // temp variable
						int reduced_window = RandomInt(Window - 1);
						for (int offset = K; offset <= (Window*K) - reduced_window; offset = offset + K) // predicting kmer given its context
						{
							if (pos + offset < seqlen - K + 1)
							{
								char* kmer2 = kmer + offset; // collecting kmers from right
								int kmer2_idx = KmerIndex(kmer2);
								context[context_size++] = kmer2_idx; // storing the indices in the array (is it array?) 
								for (int i = 0; i < SigLength; i++)
								{
									l1[i] += kmer_vector[kmer2_idx][i];
								}

							}
							if (0 <= pos - offset)
							{
								char* kmer2 = kmer - offset; // collecting kmers from the left 				
								int kmer2_idx = KmerIndex(kmer2);
								context[context_size++] = kmer2_idx;
								for (int i = 0; i < SigLength; i++)
								{
									l1[i] += kmer_vector[kmer2_idx][i];
								}
							}
						}

						// Imp: We are using average of the context for feeding in the NN, addition or concatenation are other two important aspects of it.
						int count = 1 + context_size;
						if (count > 1)
						{
							for (int i = 0; i < SigLength; i++)
								l1[i] /= count;
						}
						// Train CBOW ...
						int kmer_index = KmerIndex(kmer); //index of the kmer 
						memset(neu1e, 0, sizeof(neu1e)); // initializing with zeros
						int* path = points[kmer_index]; // path is a pointer to the 
						depth = lengths[kmer_index];

						for (int h = 0; h < depth; h++)
						{
							int node_idx = path[h];
							double dot = 0;
							for (int i = 0; i < SigLength; i++)
							{
								dot += vocab_vector[node_idx][i] * l1[i];
							}
							double fa = 1.0 / (1.0 + exp(-dot));
							double ga = (1.0 - codes[kmer_index][h] - fa) * next_alpha;
							for (int i = 0; i < SigLength; i++)
							{
								neu1e[i] += ga * vocab_vector[node_idx][i];
							}
							for (int i = 0; i < SigLength; i++)
							{
								vocab_vector[node_idx][i] += ga * l1[i];
							}
						}
						for (int j = 0; j < context_size; j++)
						{
							int kmer_idx = context[j];
							for (int i = 0; i < SigLength; i++)
							{
								kmer_vector[kmer_idx][i] += neu1e[i];
							}
						}
						memset(neu2e, 0, sizeof(neu2e));
						for (int n = 0; n < n_inliers; n++)
						{

							int global_inlier_idx = family_start[family] + RandomInt(family_size[family] - 1); // Randomness
							if (global_inlier_idx != seq_idx)
							{
								memset(neu2e_inlier, 0, sizeof(neu2e_inlier));
								double dot = 0;
								for (int i = 0; i < SigLength; i++)
									dot += seq_vector[seq_idx][i] * seq_vector[global_inlier_idx][i];

								Sig_seq_inlier = 1.0 / (1.0 + exp(dot));
								for (int i = 0; i < SigLength; i++)
									neu2e_inlier[i] += Sig_seq_inlier * seq_vector[global_inlier_idx][i];

								memset(neu2e_outlier, 0, sizeof(neu2e_outlier));

								for (int otherfamily = 0; otherfamily < num_families; otherfamily++) // for outliers from each of the remaining classes. 
																									 //int otherfamily;
																									 // randomly selecting the classes for iteration through each kmer
																									 //for (int i = 0; i<r_families; i++)
								{
									//otherfamily = RandomInt(num_families - 1);
									// generate random number five times between 0 and 633
									if (otherfamily != family)
										for (int m = 0; m < n_outliers; m++)
										{
											int global_outlier_idx = family_start[otherfamily] + RandomInt(family_size[otherfamily] - 1); // randomness
											double dot2 = 0;
											for (int i = 0; i < SigLength; i++)
												dot2 += seq_vector[seq_idx][i] * seq_vector[global_outlier_idx][i];
											double Sig_seq_outlr = 1.0 / (1.0 + (exp(-dot2)));
											for (int i = 0; i < SigLength; i++)
												neu2e_outlier[i] += -Sig_seq_outlr * seq_vector[global_outlier_idx][i];
										}
								}

							}
							for (int i = 0; i < SigLength; i++)
								neu2e[i] = neu2e_inlier[i] + neu2e_outlier[i];

						}
						for (int i = 0; i < SigLength; i++)
							neu2e[i] *= next_alpha * gamma;

						// updation								
						for (int i = 0; i < SigLength; i++)
							seq_vector[seq_idx][i] += (neu2e[i] + neu1e[i]);
						for (int j = 0; j < context_size; j++)
						{
							int kmer_idx = context[j];
						}
					}
				}
							}
		}
	}

	SaveSequences(path_save_vectors);
	SaveModel(path_save_model);

	for (int i = 0; i < num_kmers; i++)
	{
		free(codes[i]);
		free(points[i]);
	}

	for (int j = 0; j < num_nodes; j++)
		free(vocab_vector[j]);

	for (int kmer = 0; kmer < num_kmers; kmer++)
		free(kmer_vector[kmer]);

	for (int j = 0; j < num_sequences; j++)
		free(seq_vector[j]);

	free(vocab_vector);
	free(kmer_vector);
	free(seq_vector);
	free(codes);
	free(points);
}

int main(int argc, char *argv[])
{
	
	start = clock();
	path_save_model = "C:\\Users\\Dhananjay\\Desktop\\Review_PlosOne\\DATA\\uniref_exp3\\Herierchical\\supervec\\M_suv_2.binary";
	path_file_to_be_processed = "C:\\Users\\Dhananjay\\Desktop\\Review_PlosOne\\DATA\\uniref_exp3\\Herierchical\\split2\\T_100_2.fasta";
	path_save_vectors = "C:\\Users\\Dhananjay\\Desktop\\Review_PlosOne\\DATA\\uniref_exp3\\Herierchical\\supervec\\T_suv_2.binary";
	doctag_path = "C:\\Users\\Dhananjay\\Desktop\\Review_PlosOne\\DATA\\uniref_exp3\\Herierchical\\split2\\doctag_txt2.txt";
	vocab_path = "C:\\Users\\Dhananjay\\Desktop\\Review_PlosOne\\DATA\\uniref_exp3\\Herierchical\\split2\\vocab_txt2.txt";

	
	/*
	path_save_model = argv[1];
	path_file_to_be_processed = argv[2];
	path_save_vectors = argv[3];
	doctag_path = argv[4];
	vocab_path = argv[5];
	Process(1, 1);
	*/
	Process(1, 1);

	
	end1 = clock();
	training_time = ((double)(end1 - start)) / CLOCKS_PER_SEC;
	printf("training in %f sec ", training_time);

	int a;

	printf("Please input an integer value: ");
	scanf("%d", &a);
	
}