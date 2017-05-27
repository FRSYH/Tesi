#ifndef MATRIX_H_  /* Include guard */
#define MATRIX_H_

void swap_rows(long long **m, int row, int col, int j, int i);  // scambia tra di loro due righe della matrice

void print_matrix (long long **m, int row, int col); // stampa la matrice

//copia il vettore vet2 in vet1, entrambi di lunghezza len
void vctcpy(int *vet1, int const *vet2, int len);

void matrix_free_long(long long ***m, int row, int col);

void matrix_alloc_int(int ***m, int row, int col);

void matrix_free_int(int ***m, int row, int col);

void matrix_cpy(long long **m1, int row, int col, long long **m2);

void matrix_alloc_long(long long ***m, int row, int col);

void matrix_realloc_long(long long ***m, int new_row, int new_col);

void add_row_to_matrix(long long ***m, int *row, int col, long long *r);

//compute the number of null rows (rows full of 0)
int null_rows(long long **m, int row, int col);

//eliminate the matrix null rows (reallocation - resize)
void eliminate_null_rows(long long ***m, int *row, int col);

void append_matrix(long long ***m1, int *row1, int col1, long long **m2, int row2, int col2);

#endif //MATRIX_H_
