#include <stdlib.h>
#include "matrix.h"
#include "linalg.h"
#include <stdio.h>
//n mod p 
//Riduzione di n in modulo p.
long long mod(long long n, long long p){
	long long v = n,x =0;

	if( v >= p ){
		v = n%p;
	}else{
		if( v < 0 ){
			x = n/p;
			v = n-(x*p);
			v += p;
		}
	}
	return v;
}

//inverso moltiplicativo di n in modulo p (con p primo).
long long invers(long long n, long long p){
	long long b0 = p, t, q;
	long long x0 = 0, x1 = 1;
	if (p == 1) return 1;
	while (n > 1) {
		q = n / p;
		t = p, p = (n % p), n = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;	
}


// a + b mod p
//sommatoria di a e b in modulo p
long long add_mod(long long a, long long b, long long p){
	return mod((a+b),p);
}

// a - b mod p
//sottrazione di a e b in modulo p
long long sub_mod(long long a, long long b, long long p){
	return mod((a-b),p);
}

// a * b mod p
//prodotto di a e b in modulo p
long long mul_mod(long long a, long long b, long long p){
	return mod((a*b),p);
}


//restituisce il numero di tutti i possibili monomi con n variabili e grado <= m
int monomial_combinations(int n, int m) {

	// dichiarato double per compatibilità con il fattoriale della libreria gsl
	int result = 0, num, den;
	// result = Sommatoria (per j da 1 a m) {(j+n-1)! / j!*(n-1)!}
	// semplificato a: Sommatoria (per j da 1 a m) {(j+n-1)*(j+n-2)* ... *(n) / j!}
	for (int j = 1; j <= m; j++) {
		num = 1;
		for (int k = j; k > 0 ; k--)
			num = num * (n+k-1);
		den = factorial(j);
		result += (num / den);
	}
	return  result;
}

//Calcola la riduzione di Gauss della matrice m (matrice di grandezza row e col).
//La riduzione è calcolata sulla triangolare superiore sinistra.
void gauss(long long **m, int row, int col, int modulo){
	
	int pivot_riga, pivot_colonna,righe_trovate,j;

	righe_trovate = -1;
	for( j=0; j<col; j++){
		pivot_colonna = col-(j+1);
		for( pivot_riga=(righe_trovate+1); pivot_riga<row; pivot_riga++ ){
			if( m[pivot_riga][pivot_colonna] != 0 ){
				riduzione(m,row,col,pivot_riga,pivot_colonna, modulo);
				righe_trovate++;
				swap_rows(m,row,col,righe_trovate,pivot_riga);
				break;
			}
		}
	}
}

void gauss2(long long **m, int row, int col, int module){
	
	int pivot_riga = 0,r = 0,righe_trovate = 0,i,k;
	long long s,inv,a;

	for(int pivot_colonna = col-1; pivot_colonna >= 0; pivot_colonna-- ){
		r = righe_trovate;

		while( r < row && m[r][pivot_colonna] == 0 ){
			r++;
		}
		// ho trovato la prima riga con elemento non nullo in posizione r e pivot_colonna oppure non esiste nessuna riga con elemento non nullo in posizione pivot_colonna

		if( r < row ){ //significa che ho trovato un valore non nullo
			swap_rows(m,row,col,righe_trovate,r); //sposto la riga appena trovata nella posizone corretta
			pivot_riga = righe_trovate;
			righe_trovate++;
			#pragma omp parallel for private(i,inv,s,k,a)		
			for( i = righe_trovate; i < row; i++ ){
				if( m[i][pivot_colonna] != 0 ){
					inv = invers(m[pivot_riga][pivot_colonna],module);		//inverso dell´ elemento in m[r][pivot_colonna]
					s = mul_mod(inv,m[i][pivot_colonna],module);						
					//#pragma omp parallel for private (k,a) shared (m)	
					for( k = 0; k < pivot_colonna+1; k++ ){
						a = mul_mod(s,m[pivot_riga][k],module);
						m[i][k] = sub_mod(m[i][k],a,module);

					}
				}
			}
		}
	}
}


//Calcola la riduzione di Gauss di una singola riga della matrice m.
void riduzione(long long **m, int row, int col, int riga_pivot, int j, int module){
	
	int r,k;
	long long s,inv,a;
	//#pragma omp parallel for private(inv,s,k,a) shared (r,m)
	for( r=riga_pivot+1; r<row; r++ ){

		if( m[r][j] != 0 ){
			inv = invers(m[riga_pivot][j],module);			//calcola l'inverso moltiplicativo di m[riga_pivot][j] nel campo indicato
			s = mul_mod(inv,m[r][j],module);
			//#pragma omp parallel for private(a)
			for( k=0; k<col; k++ ){
				a = mul_mod(s,m[riga_pivot][k],module);
				m[r][k] = sub_mod(m[r][k],a,module);

			}
			
		}
	}
}


//confronta due monomi di *arg variabili secondo l'ordinamento grevlex
//restituisce un intero positivo se monom1 > monom2, zero se sono uguali, uno negativo altrimenti
//i monomi sono sempre rappresentati come array di lunghezza pari al numero delle variabili
//sono fatti diversi cast perchè il tipo degli argomenti è compatibile con qsort_r
int grevlex_comparison(const void *monom1, const void *monom2, void *arg) {

	int degree1 = 0, degree2 = 0, n, *mon1, *mon2;
	n = *((int *) arg);
	mon1 = *((int **) monom1);
	mon2 = *((int **) monom2);
	
	//calcolo i gradi dei monomi
	for (int v = 0; v < n; v++) {
		degree1 += mon1[v];
		degree2 += mon2[v];
	}
	if (degree1 > degree2)
		return 1;
	else if (degree1 < degree2)
		return -1;
	//se il grado è uguale guardo l'utlima cifra non nulla
	//del array risultante dalla sottrazione dei monomi
	else {
		int *temp = malloc(n * sizeof(int));
		for (int v = 0; v < n; v++)
			temp[v] = mon1[v] - mon2[v];
		for (int v = (n-1); v >= 0; v--) {
			if (temp[v] != 0) {
				return -temp[v];
				free(temp);
				//per evitare di fare free due volte sul  puntatore lo setto a NULL dopo la free
				temp = NULL;
			}
		}
		free(temp);
	}	
	return 0;
}


//calcola il fattoriale di n (se n è negativo return -1)
int factorial(int n){
	int k;

	if (n<0) //se n è negativo non esiste il fattoriale
	{
		return -1; //codice di errore
	}else{ //altrimenti calcolo il fattoriale

		if( n==0 || n==1 ){
			return 1;
		}else{
			k=1;
			for (int i = 2; i <= n; i++){
				k *= i;	
			}
			return k;
		}
	}
}


int combination(int n, int k){
	int a,b,c;
	a = factorial(n+k-1);
	b = factorial(k);
	c = factorial((n+k-1)-k);
	return  a/(c*b);
}

//https://git.devuan.org/jaretcantu/eudev/commit/a9e12476ed32256690eb801099c41526834b6390
//mancante nella stdlib, controparte di qsort_r
//effettua una ricerca binaria di key nell'array base di lunghezza nmemb i cui elementi
//hanno dimensione size, e restituisce un puntatore all'elemento uguale a key se c'è, altrimenti NULL.
//compar è la funzione di ordinamento con cui viene confrontato key con base
//arg è il terzo argomento di compar
void *bsearch_r(const void *key, const void *base, size_t nmemb, size_t size,
                 int (*compar) (const void *, const void *, void *),
                 void *arg) {
	size_t l, u, idx;
	const void *p;
	int comparison;

	l = 0;
	u = nmemb;
	while (l < u) {
		idx = (l + u) / 2;
		p = (void *)(((const char *) base) + (idx * size));
		comparison = compar(key, p, arg);
		if (comparison < 0)
			u = idx;
		else if (comparison > 0)
			l = idx + 1;
		else
			return (void *)p;
	}
	return NULL;
}

