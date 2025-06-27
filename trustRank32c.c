
/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2024/25
* 
* Progetto dell'algoritmo TrustRank per il calcolo della fiducia delle pagine web
* in linguaggio assembly x86-32 + SSE
* 
* F. Angiulli F. Fassetti S. Nisticò, marzo 2025
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf32 pst32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o pst32.o pst32c.c -o pst32c -lm && ./pst32c $pars
* 
* oppure
* 
* ./runpst32
* 
*/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <locale.h>
#include <stdbool.h>
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

float alfaI;
float alfaB;
float unoAlfaB;

typedef struct {
	int* graph;
	MATRIX tranMat; // Transition Matrix
	MATRIX tranMatInv; // Reverse Transition Matrix
	int numPages; // Numero Pagine
	int limitOracle;
	type alfaB;
	int maxBias;
	type alfaI;
	VECTOR valoriOracolo;
	int silent;
	int display;
	VECTOR results
	;
} params;


/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (type*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (type**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16);
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

VECTOR alloc_vector(int n) {
	return (VECTOR) get_block(sizeof(type),n);
}

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

void dealloc_vector(void* mat) {
	free_block(mat);
}



/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri typeing-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/

int* load_data_int(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status;
	fp = fopen(filename, "rb");

	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	status = (int) fread(&cols, sizeof(int), 1, fp);
	status = (int) fread(&rows, sizeof(int), 1, fp);

	int* data = alloc_int_matrix(rows,cols);

	status = (int) fread(data, sizeof(int), rows*cols, fp);
	fclose(fp);

	*n = rows;
	*k = cols;

	return data;
}

MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status;
	fp = fopen(filename, "rb");

	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = (int) fread(&rows, sizeof(int), 1, fp);
	status = (int) fread(&cols, sizeof(int), 1, fp);

	MATRIX data = alloc_matrix(rows,cols);

	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 
* 	load_seq
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*1 byte: matrix data in row-major order --> charatteri che compongono la stringa
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
// char* load_seq(char* filename, int *n, int *k) {
// 	FILE* fp;
// 	int rows, cols, status, i;
//
// 	fp = fopen(filename, "rb");
//
// 	if (fp == NULL){
// 		printf("'%s': bad data file name!\n", filename);
// 		exit(0);
// 	}
//
// 	status = fread(&cols, sizeof(int), 1, fp);
// 	status = fread(&rows, sizeof(int), 1, fp);
//
//
// 	char* data = alloc_char_matrix(rows,cols);
// 	status = fread(data, sizeof(char), rows*cols, fp);
// 	fclose(fp);
//
// 	*n = rows;
// 	*k = cols;
//
// 	return data;
// }

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o typeing-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti l'intero 1 		--> numero intero a 32 bit
* 	successivi 4 byte: numero di elementi k     --> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> k numero typeing-point a precisione singola
*/
void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	//int i;
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

/*
 * Procedure Assembly

Per generare l' eseguibile:
*
* nasm -f elf32 trustRank32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o trustRank32.o trustRank32c.c -o trustRank32c -lm && ./trustRank32c $pars
* pars: -tm graph_50.ds -or t0_50_32.ds -re results_50_0.7_ds.2 -np 50 -lo 50 -ab 0.85 -mb 50 -ai 0.85
*/

extern VECTOR funzione_unica(MATRIX tranMatInv, int numPages, type decay, int max_outer_iterations, int* indici, VECTOR d, VECTOR ret, VECTOR somma, bool funz1, MATRIX tranMatParam);

/*
 * Procedure C
 */

VECTOR selectSeed(MATRIX tranMatInv, int numPages, type alfaI, int mI, int* indici, VECTOR d)
{
	VECTOR s = alloc_vector(numPages);

	for (int i = 0; i < numPages; i++)
	{
		s[i] = 1.0;
	}

	// int iterazione = 0;

	// s = alfaI * tranMatInv * s + (1 - alfa) * (1/n) * unoN
	// vettore = scalare * Matrice * vettore + scalare * scalare * vettore
	type somma = (1 - alfaI) / (type) numPages;


	return funzione_unica(tranMatInv, numPages, alfaI, mI, indici, d, s, &somma, true, tranMatInv);

	// for (int m = 0; m < mI; m++)
	// {
	// 	for (int i = 0; i < numPages; i++)
	// 	{
	// 		if (m == 0)
	// 		{
	// 			indici[i] = i; // inizializzo il vettore di indici, risparmio un'iterazione in n inizializzandolo qua dentro
	// 			d[i] = 0;
	// 		}
	// 		type riga = 0;
	// 		for (int j = 0; j < numPages; j++)
	// 		{
	// 			if (m == 0 && i == 0) // Così risparmio un ciclo in n in cui inizializzo tutti gli elementi di s a 1
	// 			{
	// 				s[j] = 1;
	// 			}
	//
	//
	// 			int x = i * numPages + j;
	//
	// 			riga = riga + alfaI * tranMatInv[x] * s[j];
	// 		}
	// 		s[i] = riga + somma;
	// 	}
	// }
	//
	// return s;
}


void merge(VECTOR arr, int index[], int left, int mid, int right)
{
	int i, j, k;
	int n1 = mid - left + 1;
	int n2 = right - mid;

	// Create temporary arrays
	VECTOR leftArr[n1], rightArr[n2];
	int leftIndex[n1], rightIndex[n2];

	// Copy data to temporary arrays
	for (i = 0; i < n1; i++)
	{
		leftArr[i] = &arr[left + i];
		leftIndex[i] = index[left + i];
	}
	for (j = 0; j < n2; j++)
	{
		rightArr[j] = arr[mid + 1 + j];
		rightIndex[j] = index[mid + 1 + j];
	}
	// Merge the temporary arrays back into arr[left..right]
	i = 0;
	j = 0;
	k = left;
	while (i < n1 && j < n2)
	{
		if (leftArr[i] <= &rightArr[j])
		{
			arr[k] = *leftArr[i];
			index[k] = leftIndex[i];
			i++;
		} else
		{
			arr[k] = rightArr[j];
			index[k] = rightIndex[j];
			j++;
		}
		k++;
	}

	// Copy the remaining elements of leftArr[], if any
	while (i < n1)
	{
		arr[k] = *leftArr[i];
		index[k] = leftIndex[i];
		i++;
		k++;
	}

	// Copy the remaining elements of rightArr[], if any
	while (j < n2)
	{
		arr[k] = rightArr[j];
		index[k] = rightIndex[j];
		j++;
		k++;
	}
}


// The subarray to be sorted is in the index range [left-right]
void mergeSort(VECTOR arr, int index[], int left, int right)
{
	if (left < right)
	{

		// Calculate the midpoint
		int mid = left + (right - left) / 2;

		// Sort first and second halves
		mergeSort(arr, index, left, mid);
		mergeSort(arr, index, mid + 1, right);

		// Merge the sorted halves
		merge(arr, index, left, mid, right);
	}
}


int* rank(int* index, VECTOR s, int numPages) // è praticamente un sort, ma funziona su due liste
{
	mergeSort(s, index, 0, numPages - 1);
	return index;
}

VECTOR copy_vector(VECTOR src, int size) {
	VECTOR dest = alloc_vector(size);
	memcpy(dest, src, sizeof(type) * size);
	return dest;
}


VECTOR computeScores(MATRIX tranMat, type alfaB, int maxBias, VECTOR d, int numPages)
{
	VECTOR ret = copy_vector(d, numPages);
	type unoAlfaB = (type) 1 - alfaB;
	VECTOR somma = alloc_vector(numPages);
	memset(somma, 0, numPages * 1 * sizeof(type));
	return funzione_unica(tranMat, numPages, alfaB, maxBias, NULL, d, ret, somma, false, tranMat);

	/*
	 *				 | A B C |	 | 1 |   | X |   | (alfaB*1*A + alfaB*2*B + alfaB*3*C) + X |
	 * ret = alfaB * | D E F | * | 2 | + | Y | = | (alfaB*1*D + alfaB*2*E + alfaB*3*F) + Y |
	 *				 | G H I |	 | 3 |   | Z |   | (alfaB*1*G + alfaB*2*H + alfaB*3*I) + Z |
	 */
	//
	// for (int b = 0; b < maxBias; b++)
	// {
	// 	// ret = alfaB * tranMat * ret + (1 - alfaB) * d
	// 	// vettore = scalare * matrice * vettore + scalare * vettore
	//
	// 	for (int i = 0; i < numPages; i++)
	// 	{
	// 		somma[i] = unoAlfaB * d[i];
	// 		type riga = 0;
	// 		for (int j = 0; j < numPages; j++)
	// 		{
	// 			riga = riga + alfaB * ret[j] * tranMat[i * numPages + j];
	// 		}
	//
	// 		ret[i] = riga + somma[i];
	// 	}
	// }
	// return ret;
}


MATRIX trustRank(MATRIX tranMat, MATRIX tranMatInv, int numPages, int limitOracle, type alfaB, int maxBias, type alfaI, VECTOR valoriOracolo)
{
	int* indici = alloc_int_matrix(numPages, 1);
	VECTOR d = alloc_vector(numPages);
	VECTOR s = selectSeed(tranMatInv, numPages, alfaI, maxBias, indici, d);

	int* sigma = rank(indici, s, numPages); //rank restituisce una lista ordinata per l'affidabilità delle pagine (CONTIENE INDICI PAG)

	for (int i = 0; i < limitOracle; i++) //Singolo FOR
	{
		if (valoriOracolo[sigma[i]] == 1) // Al posto della chiamata a funzione
		{
			d[sigma[i]] = (type) 1 / (type) numPages; // MEMORIZZO GIà NORMALIZZATO SULLA LUNGHEZZA
		}
	}

	// Nel file trustRank32c.c, dentro la funzione trustRank, dopo il loop dell'oracolo:

	// Aggiungi queste righe per debug (puoi commentarle o rimuoverle dopo aver verificato)

	printf("Vector d after oracle step:\n");
	for (int i = 0; i < numPages; i++) {
		printf("%f ", d[i]);
	}
	printf("\n");


	return computeScores(tranMat, alfaB, maxBias, d, numPages);

	return computeScores(tranMat, alfaB, maxBias, d, numPages);
}

void exec(params* input)
{
	alfaI = input->alfaI;
	alfaB = input->alfaB;
	unoAlfaB = (float) 1 - input->alfaB;
	input->results= trustRank(input->tranMat, input->tranMatInv, input->numPages, input->limitOracle, input->alfaB, input->maxBias, input->alfaI, input->valoriOracolo);

}

void loadTranMat(params* input, int archi)
{
	MATRIX tranMat = alloc_matrix(input->numPages, input->numPages);
	MATRIX tranMatInv = alloc_matrix(input->numPages, input->numPages);

	int* tempIng = alloc_int_matrix(input->numPages, 1);
	int* tempUsc = alloc_int_matrix(input->numPages, 1);
	int* tempArco = alloc_int_matrix(input->numPages, input->numPages);


	int* graph = input->graph;

	for (int i = 0; i < input->numPages; i++)
	{
		tempIng[i] = 0;
		tempUsc[i] = 0;
	}
	for (int i = 0; i < archi * 2; i += 2)
	{
		int p = graph[i];
		int q = graph[i + 1];

		tempUsc[p] = tempUsc[p] + 1;
		tempIng[q] = tempIng[q] + 1;
		tempArco[q * input->numPages + p] = 1;
	}

	for (int i = 0; i < input->numPages; i++)
	{
		for (int j = 0; j < input->numPages; j++)
		{
			if (tempArco[i * input->numPages + j] == 1)
			{
				tranMat[i * input->numPages + j] = (type) 1 / (type) tempUsc[j];
				tranMatInv[j * input->numPages + i] = (type) 1 / (type) tempIng[j];
			}
		}
	}

	input->tranMat = tranMat;
	input->tranMatInv = tranMatInv;

	dealloc_matrix(tempIng);
	dealloc_matrix(tempUsc);
	dealloc_matrix(tempArco);
}


int main(int argc, char** argv)
{

	char* fname_graph = NULL;
	char* fname_oracle = NULL;
	char* fname_result = NULL;

	clock_t t;
	int d;
	int righe;
	int archi;
	type time;

	//
	// Imposta i valori di default dei parametri
	//
	params* input = malloc(sizeof(params));
	input->graph = NULL;
	input->tranMat = NULL;
	input->tranMatInv = NULL;
	input->results = NULL;
	input->numPages = 0;
	input->limitOracle = 0;
	input->alfaB = 0;
	input->maxBias = 0;
	input->alfaI = 0;
	input->valoriOracolo = NULL;
	input->silent= -1;
	input->display= -1;

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if(argc <= 7){
		printf("%s -tm <TM> -or <OR> -re <RE> -np <NP> -lo <LO> -ab <AB> -mb <MB> -ai <AI> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tTM: nome del file contenente la transition matrix delle pagine web\n");
		printf("\tOR: nome del file contenente le risposte dell'oracolo umano\n");
		printf("\tRE: nome del file in cui verranno salvati i risultati dell'algoritmo\n");
		printf("\tNP: numero delle pagine web\n");
		printf("\tLO: limite di iterazioni con l'oracolo\n");
		printf("\tAB: decay factor bias\n");
		printf("\tMB: numero massimo di iterazioni con bias\n");
		printf("\tAI: decay factor\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc)
	{
		if (strcmp(argv[par],"-s") == 0)
		{
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0)
		{
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-tm") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing graph file name!\n");
				exit(1);
			}
			fname_graph = argv[par];
			par++;
		} else if (strcmp(argv[par],"-or") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing oracle file name!\n");
				exit(1);
			}
			fname_oracle = argv[par];
			par++;
		} else if (strcmp(argv[par],"-re") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing result file name!\n");
				exit(1);
			}
			fname_result= argv[par];
			par++;
		} else if (strcmp(argv[par],"-np") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing np value!\n");
				exit(1);
			}
			input->numPages = atoi(argv[par]); //atof --> sto assegnando al parametro to atof dell'oggetto param
			par++;
		} else if (strcmp(argv[par],"-lo") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing lo value!\n");
				exit(1);
			}
			input->limitOracle = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-ab") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing ab value!\n");
				exit(1);
			}
			input->alfaB = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-mb") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing seed value!\n");
				exit(1);
			}
			input->maxBias = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-ai") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing seed value!\n");
				exit(1);
			}
			input->alfaI = atof(argv[par]);
			par++;
		}else
		{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}


	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(fname_graph == NULL || strlen(fname_graph) == 0)
	{
		printf("Missing graph file name!\n");
		exit(1);
	}
	if(fname_oracle == NULL || strlen(fname_oracle) == 0)
	{
		printf("Missing oracle file name!\n");
		exit(1);
	}

	input->graph = load_data_int(fname_graph, &archi, &d);

	// printf("\n[");
	// for (int i = 0; i < archi * 2; i += 2)
	// {
	// 	printf("\t[%d, %d]\n", input->graph[i], input->graph[i + 1]);
	//
	// }
	// printf("]\n");

	if(d != 2){
		printf("Invalid size of graph file, should be (NumeroArchi)x2!\n");
		exit(1);
	}

	loadTranMat(input, archi);


	righe = input->numPages;
	d = 1;

	input->valoriOracolo= load_data(fname_oracle, &righe, &d);
	// dealloc_matrix(fname_oracle);

	// printf("\n[");
	// for (int i = 0; i < input->numPages; i++)
	// {
	// 	printf(" %f, ", input->valoriOracolo[i]);
	//
	// }
	// printf("]\n");

	if(d != 1)
	{
		printf("Invalid size of oracle file, should be %ix1!\n", input->numPages);
		exit(1);
	}

	if(input->alfaB <= 0)
	{
		printf("Invalid value of ab parameter!\n");
		exit(1);
	}

	if(input->maxBias <= 0)
	{
		printf("Invalid value of mb parameter!\n");
		exit(1);
	}

	if(input->alfaI <= 0)
	{
		printf("Invalid value of ai parameter!\n");
		exit(1);
	}

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent)
	{
		printf("Graph file name: '%s'\n", fname_graph);
		printf("Oracle file name: '%s'\n", fname_oracle);
		printf("Number of web pages: %d\n", input->numPages);
	}

	// COMMENTARE QUESTA RIGA!
	//prova(input);
	//

	//
	// Calcolo dei valori di Trust Rank
	//
	t = clock();
	exec(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
	{
		printf("Execution time = %.3f secs\n", time);
    }
	else
	{
		printf("\n%.3f\n", time);
	}

	//
	// Salva il risultato
	//
	sprintf(fname_result, "out32_%d_.ds2", input->numPages);
	save_out(fname_result, input->results, input->numPages);
	//dealloc_matrix(fname_result);

	if(input->display)
	{
		if(input->results == NULL)
			printf("out: NULL\n");
		else
		{
			int i,j;
			printf("results: [");
			for(i=0; i<input->numPages; i++)
			{
				printf("%f,", input->results[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	free(input);

	return 0;

}
