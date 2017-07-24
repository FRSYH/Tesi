#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

//Scambio di due righe della matrice m.
void swap_rows(long long **m, int row, int col, int j, int i){
	
	int k;
	long long tmp;
	if( j!=i ){
		#pragma omp parallel for private (k,tmp) shared (m)
		for(k=0;k<col;k++){
			tmp = m[i][k];
			m[i][k] = m[j][k];
			m[j][k] = tmp;
		}
	}
}

//Stampa formattata della matrice
void print_matrix (long long **m, int row, int col){

	int i,j;	
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++){
			printf("%lli ",m[i][j]);						
		}
		printf("\n\n\n");
	}
	printf("\n");
}

//copia il vettore vet2 in vet1, entrambi di lunghezza len
void vctcpy(int *vet1, const int *vet2, int len) {
	for (int i = 0; i < len; i++)
		vet1[i] = vet2[i];
	return;
}

int null_rows(long long **m, int row, int col){
//calcola il numero di righe nulle presenti nella matrice m.

	int i,j,last,null_rows;
	null_rows = 0;
	for(i=0; i<row; i++){
		last = -1;
		for(j=col-1; j>-1; j--){
			if(m[i][j] != 0 ){
				last = j;
				break;
			}
		}
		if( last == -1 )
			null_rows++;
	}
	return null_rows;
}

void eliminate_null_rows(long long ***m, int *row, int col){
//Elimina dalla matrice m le righe nulle.
//N.B. questa procedura elimina le ultime righe nulle della matrice.
//Questa funzione DEVE essere utilizzata dopo la riduzione di Gauss.
//La riduzione di Gauss sposta nelle ultime posizioni tutte le righe nulle.
//Se non si esegue questa funzione dopo Gauss si possono eliminare righe non nulle.	

	int null_row = null_rows(*m,*row,col);
	if(null_row != 0){
		*m = realloc( *m , (*row - null_row ) * sizeof (long long *));
		*row = *row - null_row;
	}
}

void matrix_free_long(long long ***m, int row, int col){
//Deallocazione di una matrice di tipo long long con dimensioni indicate.	
	for (int i=0; i<row; i++)      
		free((*m)[i]);
	free(*m);	
}

void matrix_alloc_int(int ***m, int row, int col){
//Allocazione di una matrice di tipo int con dimensioni indicate.	
	*m = malloc(row * sizeof (int *) );
	if( *m != NULL )
		for (int i=0; i<row; i++)
			(*m)[i] = calloc(col , sizeof (int) );	
}

void matrix_free_int(int ***m, int row, int col){
//Deallocazione di una matrice di tipo int con dimensioni indicate.	
	for (int i=0; i<row; i++)      
		free((*m)[i]);
	free(*m);	
}

//copia la matrice m1 nella matrice m2
void matrix_cpy(long long **m1, int row, int col, long long **m2){
	
	int i,j;
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			m2[i][j] = m1[i][j];
		}
	}

}


void matrix_alloc_long(long long ***m, int row, int col){
//Allocazione di una matrice di tipo int con dimensioni indicate.	
	*m = malloc(row * sizeof (long long *) );
	if( *m != NULL )
		for (int i=0; i<row; i++){
			(*m)[i] = calloc(col , sizeof (long) );	
			if ((*m)[i]==NULL)
			{
				break;
			}
		}
}

void matrix_realloc_long(long long ***m, int new_row, int new_col){
	int i;
	*m = realloc( *m , (new_row) * sizeof (long long *));
	for(i=0; i<new_row; i++){
		(*m)[i] = realloc((*m)[i] , new_col * sizeof (long long) );
	}
}

//aggiunge la riga r alla matrice m, r deve avere linghezza uguale al numero delle colonne di m
void add_row_to_matrix(long long ***m, int *row, int col, long long *r){
	
	int i;
	*m = realloc( *m , (*row+1) * sizeof (long long *));
	(*m)[*row] = malloc(col * sizeof (long long) );
	for(i=0; i<col; i++){
		(*m)[*row][i] = r[i];
	}
	*row = *row + 1;	

}

void append_matrix(long long ***m1, int *row1, int col1, long long **m2, int row2, int col2){
	int i=0;
	if( col1 == col2 ){ //se le matrici hanno lo stesso numero di colonne
		for(i=0; i<row2; i++){
			add_row_to_matrix(m1,row1,col1,m2[i]);
		}
	}
}


//funzione di confronto tra la righa rowA con rowB, scorrendo le colonne da destra a sinistra
//restituisce 1 se rowA > rowB, -1 se rowB > rowA, 0 altrimenti. Compatibile con qsort_r
int compare_rows(const void *rowA, const void *rowB, void *columns) {

	long long *row1, *row2;
	int col;
	
	col = *((int *) columns);
	row1 = *((long long **) rowA);
	row2 = *((long long **) rowB);
	
	for (int i = col-1; i >= 0; i--) {
		if (row1[i] > row2[i])
			return 1;
		else if (row1[i] < row2[i])
			return -1;
	}
	
	return 0;
}


//tolgo da m2 le righe iniziali uguali a m1
void eliminate_equal_rows(long long **m1, int row1, long long ***m2, int *row2, int col) {
		
		int eliminated_rows = 0;
		
		//per ogni riga consecuitiva in m1 uguale a quella in m2
		while (eliminated_rows < row1 && !compare_rows(&m1[eliminated_rows], &(*m2)[eliminated_rows], &col))
			eliminated_rows++;
			
		
		//elimino le righe iniziali uguali trovate
		for (int i = 0; i < eliminated_rows; i++)
			free((*m2)[i]);
		
		//faccio iniziare m2 da dopo le righe eliminate e aggiorno il numero di righe
		*m2 = &((*m2)[eliminated_rows]);
		*row2 = *row2 - eliminated_rows;
		return;
}

