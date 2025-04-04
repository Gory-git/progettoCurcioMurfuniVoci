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
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};		// hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};		// volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};		// charge

typedef struct {
	MATRIX tranMat; // Transition Matrix
	int numP; // Numero Pagine
	int lim;
	float alfaB;
	int mB;
	float alfaI;
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

int* selectSeed(MATRIX tranMatInv, int numP, float alfaI, int mI)
{

}

int* rank(int* index, int* s)
{

}

int oracle(int index)
{
	// TODO CAPIRE come implementare
}

VECTOR normalize(int* scoreDistVect, int numPages)
{
	VECTOR ret = alloc_vector(numPages);
	for (int i = 0; i < numPages; i++)
	{
		ret[i] = (type) scoreDistVect[i] / numPages; // TODO RIGUARDARE quando conosciamo megio le matrici e i vettori
	}
	return ret;
}

void computeScores(MATRIX tranMat, float alfaB, int maxBias, int* d, VECTOR trustScores){
	for (int i = 0; i < maxBias; i++)
	{
		// TODO IMPOSTARE trustScores[i] = alfaB * tranMat * trustScores + (1 - alfaB) * d
		// TODO CAPIRE
	}
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
	for (int i = 0; i < numPages; i++)
	{
		for (int j = 0; j < numPages; j++)
		{
			// TODO SCRIVERE BENE
			// ret[x] = mat[y]; dove x = i + ((numpages - 1) * i) + j
			//					dove y = (numpages - 1 - i) + ((numpages - 1) * (numpages - 1 - i) + (numpages - 1 - j)
			// TODO FORSE MEGLIO un solo for e reverso la lista, più leggibile
		}
	}
	return ret;
}



MATRIX trustRank(MATRIX tranMat, int numPages, int limitOracle, float alfaB, int maxBias, float alfaI)
{
    MATRIX tranMatInv = reverseMat(tranMat, numPages);
	int mI = 100; // numero massimo di iterazioni
	int* s = selectSeed(tranMatInv, numPages, alfaI, mI);
	int* indici = alloc_int_matrix(numPages, 1);
	int* sigma = rank(indici, s); //rank restituisce una lista ordinats per l'affidabilità delle pagine CONTIENE INDICI PAG
	int* d = alloc_int_matrix(numPages, 1);
	for (int i = 0; i < limitOracle; i++)
	{
		for (int j = 0; j < numPages; j++)
		{
			if (oracle(sigma[j]) == 1)
			{
				d[sigma[j]]=1;
			}
		}
	}
	VECTOR dNormalized = normalize(d, numPages); //somma elementi = 1
	VECTOR trustScores = dNormalized; // TODO controllare se shallow o deep copy
	computeScores(tranMat, alfaB, maxBias, d, trustScores);
	return trustScores;
}

int main(int argc, char** argv){

}