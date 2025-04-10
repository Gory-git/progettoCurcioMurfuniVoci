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
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

typedef struct {
	MATRIX tranMat; // Transition Matrix
	int numP; // Numero Pagine
	int lim;
	type alfaB;
	int mB;
	type alfaI;
	VECTOR scoreDistVect;
	MATRIX tranMatInv;
	VECTOR s;
	VECTOR index;
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
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
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
char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	
	char* data = alloc_char_matrix(rows,cols);
	status = fread(data, sizeof(char), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

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
	int i;
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
 * Funzioni ad-hoc
 */

VECTOR selectSeed(MATRIX tranMatInv, int numPages, type alfaI, int mI, int* indici)
{
	VECTOR s = alloc_vector(numPages);
	int iterazione = 0;

	// s = alfaI * tranMatInv * s + (1 - alfa) * (1/n) * unoN
	// vettore = scalare * Matrice * vettore + scalare * scalare * vettore
	type somma = (1 - alfaI) / (type) numPages;

	for (int m = 0; m < mI; m++)
	{
		for (int i = 0; i < numPages; i++)
		{
			indici[i] = i; // inizializzo il vettore di indici, risparmio un'iterazione in n inizializzandolo qua dentro
			type riga = 0;
			for (int j = 0; j < numPages; j++)
			{
				if (i == 0) // Così risparmio un ciclo in n in cui inizializzo tutti gli elementi di s a 1
				{
					s[j] = 1;
				}

				riga = riga + alfaI * tranMatInv[i * numPages + j] * s[j];
			}
			s[i] = riga + somma;
		}
	}

	return s;
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
		leftIndex[i] = &index[left + i];
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


// int oracle(int index, int* valoriOracolo)
// {
// 	return valoriOracolo[index];
// }

// NON SERVE
// VECTOR normalize(int* scoreDistVect, int numPages)
// {
// 	VECTOR ret = alloc_vector(numPages);
// 	for (int i = 0; i < numPages; i++)
// 	{
// 		ret[i] = (type) scoreDistVect[i] / (type) numPages; // TODO RIGUARDARE quando conosciamo megio le matrici e i vettori
// 	}
// 	return ret;
// }

VECTOR computeScores(MATRIX tranMat, type alfaB, int maxBias, VECTOR d, int numPages)
{
	VECTOR ret = d; // ATTUALMENTE SHALLOW COPY, TODO DEEP COPY
	type unoAlfaB = (type) 1 - alfaB;
	VECTOR somma = alloc_vector(numPages);

	/*
	 *				 | A B C |	 | 1 |   | X |   | (alfaB*1*A + alfaB*2*B + alfaB*3*C) + X |
	 * ret = alfaB * | D E F | * | 2 | + | Y | = | (alfaB*1*D + alfaB*2*E + alfaB*3*F) + Y |
	 *				 | G H I |	 | 3 |   | Z |   | (alfaB*1*G + alfaB*2*H + alfaB*3*I) + Z |
	 */

	for (int b = 0; b < maxBias; b++)
	{
		// ret = alfaB * tranMat * ret + (1 - alfaB) * d
		// vettore = scalare * matrice * vettore + scalare * vettore

		for (int i = 0; i < numPages; i++)
		{
			somma[i] = unoAlfaB * d[i];
			type riga = 0;
			for (int j = 0; j < numPages; j++)
			{
				riga = riga + alfaB * ret[j] * tranMat[i * numPages + j];
			}

			ret[i] = riga + somma[i];
		}
	}
	return ret;
}

MATRIX reverseMat(MATRIX mat, int numPages)
{
	/* T
	 * |1 2 3|
	 * |4 5 6|
	 * |7 8 9|
	 */
	/*
	 * U
	 * |9 8 7|
	 * |6 5 4|
	 * |3 2 1|
	 */
	// [i, j] -> [n - i - 1, m - j - 1]; n = num righe, m = num col; n = m = numPages

	MATRIX ret = alloc_matrix(numPages,numPages);
	for (int x = 0; x <= numPages * numPages; x++)
	{
		int y = numPages * numPages - x - 1;
		ret[y] = mat[x];
	}
	return ret;
}

MATRIX trustRank(MATRIX tranMat, int numPages, int limitOracle, type alfaB, int maxBias, type alfaI, int* valoriOracolo)
{
    MATRIX tranMatInv = reverseMat(tranMat, numPages);
	int mI = 100; // numero massimo di iterazioni, deciso empiricamente AC-DC
	int* indici = alloc_int_matrix(numPages, 1);
	VECTOR s = selectSeed(tranMatInv, numPages, alfaI, mI, indici);
	int* sigma = rank(indici, s, numPages); //rank restituisce una lista ordinata per l'affidabilità delle pagine (CONTIENE INDICI PAG)

	VECTOR d = alloc_vector(numPages);
	for (int i = 0; i < limitOracle; i++) //Singolo FOR
	{
		// if (i <= numPages - 1)
		// {
		// 	int j = i;
		// }
		// else
		// {
		// 	int j = i % numPages;
		// } Non ho capito il senso di questi if

		if (valoriOracolo(sigma[i]) == 1) // Al posto della chiamata a funzione
		{
			d[sigma[i]] = (type) 1 / (type) numPages; // MEMORIZZO GIà NORMALIZZATO SULLA LUNGHEZZA
		}
	}

	// VECTOR dNormalized = normalize(d, numPages);						//somma elementi = 1 NON SERVE, POSSIAMO FARLO
																		// DIRETTAMENTE NEL CICLO DI SOPRA, SICCOME
																		// SAPPIAMO GIà LA LUNGHEZZA DELLA LISTA
	return computeScores(tranMat, alfaB, maxBias, d /* dNormalized INUTILE*/, numPages);
	//return d/*Normalized*/;
}

int main(int argc, char** argv)
{

}
