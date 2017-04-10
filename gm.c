#define _GNU_SOURCE 
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>

int max_degree = 0;
int module = 0;

int mod(int n, int p); //n mod p

int add_mod(int a, int b, int p); // a + b mod p

int sub_mod(int a, int b, int p); // a - b mod p

int mul_mod(int a, int b, int p); // a * b mod p

int invers(int n, int p);  // n^-1 mod p

void gauss(int **m, int row, int col); // riduzione di gauss della matrice m

void riduzione(int **m, int row, int col, int riga_pivot, int j);

void swap_rows(int **m, int row, int col, int j, int i);  // scambia le due righe indicate della matrice

void print_matrix (int **m, int row, int col); // stampa la matrice

void init_matrix(int **m, int row, int col); //inizializza la matrice dei coefficienti

//moltiplica la riga indicata per tutti i possibili monomi
void moltiplica_riga(int **m, int *row, int col, int riga, int **map,int * degree, int * degree_position); 

//returns the number of all possible monomial with n variables and degree <= m
int monomial_combinations(int n, int m); 

//initialize the vector that keeps the number of monomial with the same grade and their position
void init_degree_vector(int * degree, int * degree_position, int num_var);

//return the grade of a monomial
int grado_monomio(int posizione, int *degree_position);

int combination(int n, int k);

int main (void){


	int row, col, i, num_var, degree[max_degree+1], degree_position[max_degree+1];
	int **m, *d_row, row_max;
	
	row = 1;
	row_max = 100;
	col = 120;
	module = 773;	
	max_degree = 7;
	num_var = 3;
	d_row = &row;


//###########################################################################
	m = malloc(row * sizeof (int *) );            // allocazione della matrice
	if( m != NULL ){
		for (i=0; i<row; i++)
		{
			m[i] = calloc(col , sizeof (int) );
		}

	}
//##########################################################################
	
	init_matrix(m,row,col); 

	//gauss(m, row, col);

	moltiplica_riga(m,d_row,col,0,map,degree,degree_position);

	print_matrix(m, *d_row, col);	
	
	init_degree_vector(degree,degree_position,num_var);


//##################################################################  	
	for (i=0; i<row; i++)      // deallocazione della matrice
	{
		free(m[i]);
	}
	free(m);	
//#################################################################

	return 0;
}


int mod(int n, int p){
	int v = n;
	if( v >= p ){
		do{
			v -= p;
		}while( v >= p );
	}else{
		if( v < 0 ){
			do{
				v += p;
			}while( v < 0 );
		}
	}	
	return v;
}


int invers(int n, int p){
	int b0 = p, t, q;
	int x0 = 0, x1 = 1;
	if (p == 1) return 1;
	while (n > 1) {
		q = n / p;
		t = p, p = (n % p), n = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;	
}


int add_mod(int a, int b, int p){
	return mod((a+b),p);
}


int sub_mod(int a, int b, int p){
	return mod((a-b),p);
}


int mul_mod(int a, int b, int p){
	return mod((a*b),p);
}


void gauss(int **m, int row, int col){
	
	int pivot_riga, pivot_colonna,righe_trovate,j;
	
	righe_trovate = -1;
	printf("entro gauss\n");
	for( j=0; j<col; j++){
		pivot_colonna = j;
		for( pivot_riga=(righe_trovate+1); pivot_riga<row; pivot_riga++ ){
			if( m[pivot_riga][pivot_colonna] != 0.0 ){
				riduzione(m,row,col,pivot_riga,pivot_colonna);
				righe_trovate++;
				swap_rows(m,row,col,righe_trovate,pivot_riga);
				break;
			}
		}
	}

}


void riduzione(int **m, int row, int col, int riga_pivot, int j){
	
	int r,k,s,inv,a;

	printf("entro riduzione\n");
	printf("riga_pivot %d\n", riga_pivot);
	for( r=riga_pivot+1; r<row; r++ ){
		if( m[r][j] != 0.0 ){
			inv = invers(m[riga_pivot][j],module);			//calcola l'inverso moltiplicativo di m[riga_pivot][j] nel campo indicato
			s = mul_mod(inv,m[r][j],module);
			for( k=0; k<col; k++ ){
				a = mul_mod(s,m[riga_pivot][k],module);
				m[r][k] = sub_mod(m[r][k],a,module);
			}
		}
	}
	
}


void swap_rows(int **m, int row, int col, int j, int i){
	
	int k,tmp;
	for(k=0;k<row;k++){
		tmp = m[i][k];
		m[i][k] = m[j][k];
		m[j][k] = tmp;
	}
}


void print_matrix (int **m, int row, int col){
	

	int i,j;	
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++)
		{
			printf("%d ",m[i][j]);						
		}
		printf("\n\n\n");
	}
	printf("\n");
}


void init_matrix(int **m, int row, int col){

	m[0][0] = 640;
	m[0][2] = 640;
	m[0][3] = 640;
	m[0][5] = 640;
	m[0][6] = 640;
	m[0][8] = 640;
	m[0][23] = 640;
	m[0][24] = 640;
	m[0][25] = 640;
	m[0][27] = 640;
	m[0][28] = 640;
	m[0][32] = 640;
}

//returns the number of all possible monomial with n variables and degree <= m
int monomial_combinations(int n, int m) {

	double result = 0., num, den;
	
	// result = Summation (j = 1 -> m) {(j+n-1)! / j!*(n-1)!}
	// simplified to Summation (j = 1 -> m) {(j+n-1)*(j+n-2)* ... *(n) / j!}
	for (int j = 1; j <= m; j++) {
		num = 1.;
		for (int k = j; k > 0 ; k--)
			num = num * (n+k-1);
		den = gsl_sf_fact (j);
		result += (num / den);
	}
	
	return (int) result;

}


void moltiplica_riga(int **m, int * row, int col, int riga, int **map,int * degree, int * degree_position){
//moltiplica la riga indicata per ogni monomio in modo tale che il prodotto abbia grado <= del grado massimo

	int grado_massimo_riga, grado_massimo_monomio, grado_a, grado_b, grado_prodotto;
	int i,j,k,last,sum, posizione_nuovo_monomio, offset, offset_a, offset_b,v,value;

	
    v = 0;
	last = -1;
	//cerco la posizione dell'ultimo coefficiente non nullo del polinomio rappresentato nella riga, grado più alto
	for(i=col-1; i>0; i--){
		v = m[riga][i];		
		if( v != 0.0 ){
			last = i;
			break;
		}
	}
	printf("last %d\n",last);
	//risalgo al grado del monomio appena trovato
	//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado
	if( last != -1 ){

		grado_massimo_riga = grado_monomio(last,degree_position);
		
		//calcolo il grado massimo che deve avere il monomio per cui moltiplicare		
		grado_massimo_monomio = max_degree - grado_massimo_riga;		
		
		printf("grado massimo riga %d grado massimo monomio %d\n", grado_massimo_riga,grado_massimo_monomio);
		//moltiplico la riga per ogni monomio possibile
		sum = 1;
		for(i=1; i<degree_position[grado_massimo_monomio+1]; i++){     //scorre tutti i gradi per i quali posso moltiplicare
			
			for(j=0; j<last; j++){     //scorre fino all'ultimo elemento della riga
				m[riga][ map[j][i] ] = m[riga][j];
			}			
				
		}
		printf("moltiplicazione eseguita con successo !!\n");

	}

}


int grado_monomio(int posizione, int *degree_position){
	
	int i,grado;
	for(i=0; i<max_degree+1; i++){
		if( degree_position[i] > posizione ){
			grado = i-1;
			break;	
		}			
	}	
	return grado;
}

void init_degree_vector(int * degree, int * degree_position, int num_var){
//inizializza degree con il numero di monomi di grado i-esimo <= del grado massimo
//inizializza degree_position con la posizione di inizio del gruppo di monomi dello stesso grado

	int i,j,c;
	for(i=0; i<max_degree+1; i++){
		c = combination(num_var,i);
		degree[i] = c;
	}

	degree_position[0] = 0;
	for(i=1; i<max_degree+1; i++){
		c = 0;		
		for(j=0; j<i; j++){
			c += degree[j];
		}
		degree_position[i] = c;
	}
}


int combination(int n, int k){
	
	int a,b,c;
	a = gsl_sf_fact(n+k-1);
	b = gsl_sf_fact(k);
	c = gsl_sf_fact((n+k-1)-k);
	return  a/(c*b);
}


