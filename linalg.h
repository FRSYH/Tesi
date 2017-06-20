#ifndef LINALG_H_  /* Include guard */
#define LINALG_H_

long long mod(long long n, long long p); //n mod p

long long add_mod(long long a, long long b, long long p); // a + b mod p

long long sub_mod(long long a, long long b, long long p); // a - b mod p

long long mul_mod(long long a, long long b, long long p); // a * b mod p

long long invers(long long n, long long p);  // inverso moltiplicativo in mod p

void gauss(long long **m, int row, int col, int modulo); // riduzione di gauss della matrice m

void gauss2(long long **m, int row, int col, int modulo, int start); // riduzione di gauss della matrice m

void riduzione(long long **m, int row, int col, int riga_pivot, int j, int module);

int combination(int n, int k);

//restituisce il numero di tutti i possibili monomi con n variabili e grado <=m
int monomial_combinations(int n, int m); 

//confronta due monomi di *arg variabili secondo l'ordinamento grevlex
//restituisce un intero positivo se monom1 > monom2, zero se sono uguali, uno negativo altrimenti
int grevlex_comparison(const void *mon1, const void *mon2, void *arg);

//calcola il fattoriale di n
int factorial(int n);

//mancante nella stdlib, controparte di qsort_r
void *bsearch_r(const void *key, const void *base, size_t nmemb, size_t size,
                 int (*compar) (const void *, const void *, void *),
                 void *arg);

#endif //LINALG_H_
