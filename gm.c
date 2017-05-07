#define _GNU_SOURCE 
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>
#include <ctype.h>


int max_degree = 0;
long long module = 0;

long long mod(long long n, long long p); //n mod p

long long add_mod(long long a, long long b, long long p); // a + b mod p

long long sub_mod(long long a, long long b, long long p); // a - b mod p

long long mul_mod(long long a, long long b, long long p); // a * b mod p

long long invers(long long n, long long p);  // inverso moltiplicativo in mod p

void gauss(long long **m, int row, int col); // riduzione di gauss della matrice m

void riduzione(long long **m, int row, int col, int riga_pivot, int j);

void swap_rows(long long **m, int row, int col, int j, int i);  // scambia tra di loro due righe della matrice

void print_matrix (long long **m, int row, int col); // stampa la matrice

int init_matrix(long long **m, int row, int col,int **vet_grd, char *v, int num_var); //inizializza la matrice dei coefficienti

//moltiplica la riga indicata per tutti i possibili monomi
void moltiplica_riga(long long ***m, int *row, int col, int riga, int **map,int * degree, int **vet, int num_var); 

//initialize the vector that keeps the number of monomial with the same grade and their position
void init_degree_vector(int * degree, int num_var);

//return the grade of a monomial
int grado_monomio(int posizione, int **vet, int num_var);

int combination(int n, int k);

//returns the number of all possible monomial with n variables and degree <= m
int monomial_combinations(int n, int m); 

//restituisce un array contenente tutti i len monomi con n variabili e grado <= m
//len è il numero di possibili monomi con n variabili e grado <= m
int **monomial_computation(int n, int m, int len);

//computes all possible monomials with n variables and degree <= m
//saves them into array vet, see the definition for more informations
void monomial_computation_rec(int n, int m,  int **vet, int turn, int *monomial, int *pos);

//copies vector vet2 into vet1 of length len
void vctcpy(int *vet1, int const *vet2, int len);

//compares two monomials of *(arg) variables following the grevlex order
//returns a positive number if mon1 > mon2, 0 if they are equal, negative otherwise
int grevlex_comparison(const void *mon1, const void *mon2, void *arg);

//maps all the possible multiplications of the monomials of n variables in the array
//vet of length len into the matrix map[len][len]
void setup_map(int **map, int **vet, int len, int n, int m);


//missing in stdlib, counterpart of qsort_r
void *bsearch_r(const void *key, const void *base, size_t nmemb, size_t size,
                 int (*compar) (const void *, const void *, void *),
                 void *arg);


//multiply the entire matrix for all possible correct monomial
void moltiplica_matrice(long long ***m, int *row, int col, int **map, int *degree, int **vet, int num_var);

//compute the different degree of the matrix rows
void matrix_degree(long long **m, int row, int col, int *m_deg, int **vet, int num_var);

//formatted print of matrix degree
void print_matrix_degree(int *m_deg);

//compute the number of null rows (rows full of 0)
int null_rows(long long **m, int row, int col);

//eliminate the matrix null rows (reallocation - resize)
void eliminate_null_rows(long long ***m, int *row, int col);

//execute the multiplication, gauss reduction, elimination null_rows and compare the result to target
void execute(long long ***m, int *d_row, int col, int **map, int *degree, int **vet, int num_var);

//compare the matrix degree with the target degree wich are {0,1,2,3,4,5...max_degree} return 0 if equal -1 if not
int target_degree(int *v);

void allocation(long long ***m, int *row, int *col, int *num_var, char **v);

void matrix_free_long(long long ***m, int row, int col);

void matrix_alloc_int(int ***m, int row, int col);

void matrix_free_int(int ***m, int row, int col);

int parse(int num_var, char *vet, long long **m, int **vet_grd, int len);

int parse_mon(char * mon, int len,long long * val, int num_var, char *vet, int *grade, int pos_pol);

int main (void){
	

	int row, col, num_var;
	int *d_row,**vet, len, **map;
	long long **m;
	char *v;
	row = col = num_var = 0;
	allocation(&m,&row,&col,&num_var,&v);  //predispone la matrice dei coefficienti
	d_row = &row;

	int degree[max_degree+1];
	len = col;

	vet = monomial_computation(num_var, max_degree, len);


	qsort_r(vet, len, sizeof(int*), grevlex_comparison, &num_var);

	if( init_matrix(m,row,col,vet,v,num_var) == -1 ){ //inizializzazione matrice (lettura dati input)
		printf("Errore di input !!!\n TERMINAZIONE PROGRAMMA"); //se l'input è in formato scorrettro abort del programma
		return 0;
	}
	 

	matrix_alloc_int(&map,len,len);
	setup_map(map, vet, len, num_var, max_degree);     // a questo punto posso utilizzare la mappa

	//RISOLUZIONE PROBLEMA
	init_degree_vector(degree,num_var); 
	execute(&m,d_row,col,map,degree,vet,num_var);  //soluzione trovata
	print_matrix(m, row, col);	 //stampa la matrice soluzione

	matrix_free_long(&m,row,col);
	matrix_free_int(&map,len,len);
	matrix_free_int(&vet,len,num_var);

	return 0;
}


//n mod p 
//Riduzione di n in modulo p.
long long mod(long long n, long long p){
	long long v = n,x =0;
	if( v >= p ){
		x = n/p;
		v = n-(x*p);
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

//Calcola la riduzione di Gauss della matrice m (matrice di grandezza row e col).
//La riduzione è calcolata sulla triangolare superiore sinistra.
void gauss(long long **m, int row, int col){
	
	int pivot_riga, pivot_colonna,righe_trovate,j;
	
	righe_trovate = -1;
	for( j=0; j<col; j++){
		pivot_colonna = col-(j+1);
		for( pivot_riga=(righe_trovate+1); pivot_riga<row; pivot_riga++ ){
			if( m[pivot_riga][pivot_colonna] != 0 ){
				riduzione(m,row,col,pivot_riga,pivot_colonna);
				righe_trovate++;
				swap_rows(m,row,col,righe_trovate,pivot_riga);
				break;
			}
		}
	}
}

//Calcola la riduzione di Gauss di una singola riga della matrice m.
void riduzione(long long **m, int row, int col, int riga_pivot, int j){
	
	int r,k;
	long long s,inv,a;
	for( r=riga_pivot+1; r<row; r++ ){
		if( m[r][j] != 0 ){
			inv = invers(m[riga_pivot][j],module);			//calcola l'inverso moltiplicativo di m[riga_pivot][j] nel campo indicato
			s = mul_mod(inv,m[r][j],module);
			for( k=0; k<col; k++ ){
				a = mul_mod(s,m[riga_pivot][k],module);
				m[r][k] = sub_mod(m[r][k],a,module);
			}
			
		}
	}
}

//Scambio di due righe della matrice m.
void swap_rows(long long **m, int row, int col, int j, int i){
	
	int k;
	long long tmp;
	for(k=0;k<col;k++){
		tmp = m[i][k];
		m[i][k] = m[j][k];
		m[j][k] = tmp;
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


int init_matrix(long long **m, int row, int col, int **vet_grd, char *v, int num_var){
//Inizializza la matrice principale (dei coefficienti) con i coefficienti dei polinomi forniti come input.
	return parse(num_var,v,m,vet_grd,col);
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


/*
Moltiplica la riga indicata (riga) della matrice m per ogni monomio in modo tale che il risultato abbia grado <= del grado massimo.
Il prodotto avviene su monomi con coefficiente sempre uguale a 1.
Il prodotto consiste quindi in uno shift della posizione del monomio in esame.
Il parametro map fornisce una mappa delle posizioni in cui inserire il prodotto di due monomi.
La matrice aumenta il numero di righe in base a quanti prodotti devo eseguire.
Gli ultimi due parametri servono per il calcolo del grado di un monomio.
*/
void moltiplica_riga(long long ***m, int * row, int col, int riga, int **map,int * degree, int **vet, int num_var){


	int grado_massimo_riga, grado_massimo_monomio,i,j,last,new_row;
	last = -1;
	//cerco la posizione dell'ultimo coefficiente non nullo del polinomio rappresentato nella riga.
	for(i=col-1; i>0; i--){
		if( (*m)[riga][i] != 0 ){
			last = i;
			break;
		}
	}
	//risalgo al grado del monomio appena trovato
	//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado
	if( last != -1 ){

		grado_massimo_riga = grado_monomio(last,vet,num_var);		

		//calcolo il grado massimo che deve avere il monomio per cui moltiplicare		
		grado_massimo_monomio = max_degree - grado_massimo_riga;		
		// a questo punto conosco per quanti monomi devo moltiplicare e quindi 
		// conosco il numero di righe che devo aggiungere alla matrice
		new_row = 0;
		for(i=1; i<(grado_massimo_monomio+1); i++){
			new_row += degree[i];
		}

		*m = realloc( *m , (*row + new_row ) * sizeof (long long *));

		for (i=(*row); i< (*row + new_row ); i++)
			(*m)[i] = calloc(col , sizeof (long long) );

		for(i=1; i<(new_row+1); i++){     								//scorre tutti i monomi per i quali posso moltiplicare
			for(j=0; j<(last+1); j++)     								//scorre fino all'ultimo elemento della riga
				(*m)[*row][ map[i][j] ] = (*m)[riga][j];  				//shift nella posizione corretta indicata dalla mappa
			*row = *row + 1;											//aumento del conteggio delle righe
		}
	}
}



int grado_monomio(int posizione, int **vet, int num_var){
//Calcola il grado del monomio a partire dalla posizione occupata nel vettore (ordinato) delle posizioni rispetto l'ordinamento scelto.
//(la posizione occupata deve essere corretta).	
	int i,grado;
	grado = 0;
	for(i=0; i<num_var; i++){
		grado += vet[posizione][i];
	}
	return grado;
}


void init_degree_vector(int * degree, int num_var){
//inizializza il vettore degree con il numero di monomi di grado i-esimo <= del grado massimo
	int i,j,c;
	for(i=0; i<max_degree+1; i++){
		c = combination(num_var,i);
		degree[i] = c;
	}
}


int combination(int n, int k){
	int a,b,c;
	a = gsl_sf_fact(n+k-1);
	b = gsl_sf_fact(k);
	c = gsl_sf_fact((n+k-1)-k);
	return  a/(c*b);
}

/*
mappa tutte le possibili moltiplicazioni dei monomi di n variabili e grado <= m
dell'array vet di lunghezza len, nella matrice map[len][len].
Al termine map[x][y] contiene la posizione all'interno di vet del
monomio risultato dal prodotto di vet[x]*vet[y]
Esempio: vet[4] * vet[10] = vet [map[4][10]] 
se il grado del prodotto supera m viene messo il valore -1 nella matrice
la matrice map deve essere già correttamente allocata
l'arrey vet deve essere ordinato secondo grevlex
*/
void setup_map(int **map, int **vet, int len, int n, int m) {

	int sum, *temp = malloc(n * sizeof(int));
	
	//per ogni monomio in vet
	for (int row = 0; row < len; row++)
		//provo a moltiplicarlo con ogni monomio in vet
		for (int col = 0; col < len; col++) {
			sum = 0;
			//eseguo il prodotto (sum è la somma dei gradi)
			for (int v = 0; v < n; v++) {
				temp[v] = vet[row][v] + vet[col][v];
				sum += temp[v];
			}
			//se il grado del prodotto > grado massimo tutti i restanti prodotti
			//su quella riga sono > grado massimo, setto a -1 il resto della riga
			if (sum > m) {
				for (int i = col; i < len; i++)
					map[row][i] = -1;
				break;	
			}
			//altrimenti cerco il prodotto in vet e metto l'indice in map
			else
				map[row][col] = (int **)(bsearch_r((void *) &temp, (void *) vet, len, (sizeof(int*)), grevlex_comparison, &n)) - vet;
		}
	free(temp);
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


/*restituisce un array contenente tutti i len monomi con n variabili e grado <= m
len è il numero di possibili monomi con n variabili e grado <= m
i monomi sono array di interi di lunghezza n dove il valore di ogni posizione rappresenta
il grado della variabile in quella posizione. Esempio: n=3, x^2*y*z = [2,1,1]
len viene passato come argomento per evitare di ricalcolarlo internamente
*/
int **monomial_computation(int n, int m, int len) {

	int **vet, *monomial;
	
	//alloco la memoria per l'array
	matrix_alloc_int(&vet,len,n);
	
	//strutture di supporto necessarie per il calcolo
	monomial = malloc(n * sizeof(int));
	int pos = 0;
	
	//il calcolo è fatto dalla funzione ricorsiva correttemente parametrizzata
	monomial_computation_rec(n, m, vet, 0, monomial, &pos);
	
	free(monomial);
	
	return vet;
}



/*funzione ricorsiva che calcola tutti i possibili monomi con n variabili e grado <= m
e li inserisce nell'array vet. I monomi sono rappresentati come array di interi dove
il valore di ogni posizione rappresenta il grado della variabile in quella posizione.
Esempio: n=3, x^2*y*z = [2,1,1].
L'array vet deve essere già allocato correttamente. Gli altri parametri sono necessari
per la struttura ricorsiva della funzione e alla prima chiamata devono essere:
- turn = 0, rappresenta la posizione della variabile nel monomio
- monomial = array di interi di lunghezza n già allocato e usato per calcolare i vari monomi
- *pos = 0 puntatore ad intero, rappresenta la prima posizione libera nell'array vet
*/
void monomial_computation_rec(int n, int m, int **vet, int turn, int *monomial, int *pos) {

	//per ogni variabile provo tutti i gradi da 0 a m
	for (int degree = 0; degree <= m; degree++) {
		//se questa è la prima variabile azzero il monomio
		if (turn == 0) {
			//azzero il monomio lasciando solo il grado della prima variabile
			monomial[0] = degree;
			for (int v = 1; v < n; v++)
				monomial[v] = 0;		
		}
		//altrimenti le altre variabili aggiungo il proprio grado al monomio
		else
			monomial[turn] = degree;
		
		
		//ottengo il grado del monomio sommando i gradi delle variabili
		int sum = 0;
		for (int v = 0; v <= turn; v++)
			sum += monomial[v];
		//se il grado del monomio supera quello massimo non ha senso continuare a cercare
		//altri monomi partendo da questo, perchè tutti avranno grado maggiore o uguale
		if (sum > m)
			break;

		//se questa è l'ultima variabile copia il monomio nell'array vet and incrementa l'indice pos
		if (turn == (n-1)) {
			vctcpy(vet[(*pos)], monomial, n);
			(*pos)++;
		}
		//altrimenti richiama se stessa cambiando la variabile (turn)
		else
			monomial_computation_rec(n, m, vet, turn+1, monomial, pos);
	}
	
	return;
}


//copia il vettore vet2 in vet1, entrambi di lunghezza len
void vctcpy(int *vet1, const int *vet2, int len) {
	for (int i = 0; i < len; i++)
		vet1[i] = vet2[i];
	return;
}


void moltiplica_matrice(long long ***m, int *row, int col, int **map, int *degree, int **vet, int num_var){
//Moltiplica tutte le righe della matrice m per tutti i monomi possibili che forniscono un risultato che ha grado <= grado massimo.	
	int n,i;
	n = *row;    //n conta il numero di righe della matrice di partenza che devo moltiplicare
	for(i=0; i<n; i++){
		moltiplica_riga(m,row,col,i,map,degree,vet,num_var);
	}
}


void matrix_degree(long long **m, int row, int col, int *m_deg, int **vet, int num_var){
//m_deg è un vettore che ha lunghezza pari al grado massimo.
//la funzione calcola i gradi dei polinomi presenti nella matrice.
//Ogni cella del vettore m_deg rappresenta un grado, se esso compare nella matrice allora viene impostato a 1 o altrimenti.	

	int i,j,last,grado;
	for(i=0; i<row; i++){
		for(j=col-1; j>0; j--){
			if( m[i][j] != 0 ){
				last = j;           //posizione dell'ultimo coefficiente della riga
				break;
			}
		}
		grado = grado_monomio(last,vet,num_var);	
		m_deg[grado] = 1;		
	}
}



void print_matrix_degree(int *m_deg){
//stampa il vettore dei gradi della matrice.	
	int i;	
	printf("Gradi della matrice = {");
	for(i=0; i<max_degree+1; i++)
		if( m_deg[i] != 0 )	printf(" %d ",i);		
	printf("}\n");
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
	*m = realloc( *m , (*row - null_row ) * sizeof (long long *));
	*row = *row - null_row;
}


void execute(long long ***m, int *d_row, int col, int **map, int *degree, int **vet, int num_var){
/*
Questa funzione itera la procedura di moltiplicazione della matrice e la riduzione di Gauss fino a che non si raggiunge la terminazione.
La terminazione è data da:
	- raggiunta dei gradi {1,2,3,4,...,max_degree}
	- iterazione infinita, non c'è soluzione con gradi completi

*/
	int *m_deg = calloc(max_degree+1, sizeof(int));
	printf("Inizio procedura\n");
	matrix_degree(*m,*d_row,col,m_deg,vet,num_var);

	int flag,old,new;
	flag = old = new = 0;
	old = *d_row;
	while( flag != 1 ){
		printf("\n -Eseguo moltiplicazione, ");
		fflush(stdout);
		moltiplica_matrice(m,d_row,col,map,degree,vet,num_var);  //moltiplico la matrice per tutti i monomi possibili
		printf("numero righe: %d", *d_row);
		printf("\n -Eseguo Gauss, ");
		fflush(stdout);	
		gauss(*m, *d_row, col);                                     //applico la riduzione di Gauss
		eliminate_null_rows(m,d_row,col);							//elimino le righe nulle della matrice
  		matrix_degree(*m,*d_row,col,m_deg,vet,num_var);
		printf("numero righe: %d\n", *d_row);
		print_matrix_degree(m_deg);

		new = *d_row;

		if( old == new  ){ //se per due volte trovo una matrice con le stesso numero di righe mi fermo
			flag = 1;
		}else{
			if( target_degree(m_deg) == 0 ){  //se trovo una matrice con gradi [1,2,3...,max_degree] mi fermo
				flag = 1;
			}else{
				old = new; 
			}			
		}
	}
	printf("\nFine procedura, target raggiunto\n\n");
	free(m_deg);
}


int target_degree(int *v){
//Controlla se il vettore v rappresenta la condizione di terminazione con gradi completi {1,2,3,...,max_degree}
//Se la condizione è soddisfatta return 0 altrimenti -1	
	
	int i,flag;
	flag = 0;
	for(i=1; i<max_degree+1; i++){
		if( v[i] != 1 ){
			flag = -1;
			break;
		}
	}
	return flag;	
}


void allocation(long long ***m, int *row, int *col, int *num_var, char **v){
/*
Legge da input le seguenti informazioni:
	- modulo dei coefficienti
	- grado massimo
	- numero dei polinomi di partenza
	- variabili utilizzate nei polinomi

con queste informazioni alloca la matrice principale (matrice che conterrà i polinomi) e stabilisce il numero di variabili utilizzate.
*/
	scanf("%lli",&module); //leggo il modulo
	getchar();
	scanf("%d",&max_degree); //leggo il grado massimo
	getchar();
	scanf("%d",row);  //leggo numero dei polinomi di partenza
	getchar();

	int i,j,k,pos_pol,num_pol;
	char c;

	i=0;
	pos_pol = 0;
	*v = malloc(sizeof(char));
	c = getchar();
	while( c != '\n' ){
		(*v)[i] = c;
		i++;
		(*num_var)++;
		*v = realloc(*v, (i+1)*sizeof(char) );
		c = getchar();
	}

	*col = 1+monomial_combinations(*num_var, max_degree);

	*m = malloc((*row) * sizeof (long long *) );            // allocazione della matrice dei coefficienti
	if( *m != NULL )
		for (int i=0; i<(*row); i++)
			(*m)[i] = calloc((*col) , sizeof (long long) );	


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

int parse(int num_var, char *vet, long long **m, int **vet_grd, int len){
/*
Esegue la lettura (parse) dei polinomi di partenza nel seguente modo.
Si legge un monomio alla volta. 
Il monomio viene scomposta dalla funzione parse_mon.
Si inserisce il coefficiente del monomio nella matrice principale (matrice dei coefficienti) nella posizione corretta.
La posizione corretta è indicata da vet_grd.
Si leggono tutti i monomi di tutti i polinomi di partenza.
In caso di errore di formato nell'input la funzione si interrompe restituendo segnale di errore -1.
*/
	int pos_pol = 0,i,col;
	char c,* mon;
	long long cof = 0;
	c = getchar();

	int *grade;

	grade = calloc(3,sizeof(int));

	while( c != EOF ){

		while( c != '\n' && c != EOF){
			mon = malloc( sizeof(char) );
			i = 0;	
			while( c != '+' && c != EOF && c != '\n'){
				mon = realloc(mon, (i+1)* sizeof(char));
				mon[i] = c;
				i++;
				c = getchar();
			}
			if( parse_mon(mon,i,&cof,num_var,vet,grade,pos_pol) == -1 ){
				return -1;
			}
			//inserire monomio in posizione corretta
			col = (int **)(bsearch_r((void *) &grade, (void *) vet_grd, len, (sizeof(int*)), grevlex_comparison, &num_var)) - vet_grd;
			m[pos_pol][col] = cof;
			if(c=='\n'){
				pos_pol++;
			}
			mon = calloc(i+1,sizeof(char));
			grade = calloc(3,sizeof(int));
			c = getchar();
		}
		c = getchar();	
	}
	free(mon);
}

int parse_mon(char * mon, int len,long long * val, int num_var, char *vet, int *grade, int pos_pol){
/*
La funzione esegue il parse di un singolo monomio.
In particolare si divide il monomio in:
	- coefficiente (se esiste, altrimenti = 1)
	- vettore rappresentante i gradi delle singole variabili

es. 3 variabili xyz, Monomio: 345x^2*y 
	- coefficiente 345
	- vettore dei gradi [2,1,0] ogni posizione rappresenta il grado della variabile corrispondente (le variabili seono inserite in ordine alfabetico).

Se il coefficiente non è presente allora = 1
Se l'esponente di una variabile non è presente allora = 1
se una variabile non compare nel monomio allora grado = 0
*/
	int i,k,pos_var;
	char c,* cof,*exp;
	cof = malloc( sizeof(char) );
	exp = malloc( sizeof(char) );
	i = 0;
	pos_var = 0;
	if( isdigit(mon[i]) != 0 ){  // se c è un numero
		i = 0;
		while( isdigit(mon[i]) && i<len){
			cof = realloc(cof, (i+1) * sizeof(char));
			cof[i] = mon[i];
			i++;                   
		}
		*val = mod(atoll(cof),module);

		//printf("val: %lli\n", *val);
	}else{
		if( isalpha(mon[i]) != 0 ){
			*val = 1;
			//printf("val: %lli\n", *val);
		}else{
			//errore
			return -1;
		}	
	}
	if( i < len ){
		while( i < len ){
			if( mon[i] == '*' ){ 
				i++;    
			}
			if( i<len && (isalpha(mon[i]) != 0) ){
				for(k=0; k<num_var; k++){
					if( vet[k] == mon[i] ){
						pos_var = k;
						break;
					}
				}
				i++;
				if( i<len ){
					if( mon[i] == '^' ){ //ho trovato il grado della variabile
						i++;
						if( isdigit(mon[i]) != 0 ){
							k = 0;
							while( isdigit(mon[i]) != 0 && i < len ){
								exp = realloc(exp, (k+1) * sizeof(char));
								exp[k] = mon[i];
								i++;
								k++;
							}
						grade[pos_var] = atoi(exp);
						}else{
							//errore
							return -1;									
						}
					}else{
						grade[pos_var] = 1;	
					}
				}else{
					grade[pos_var] = 1;
				}								
			}else{
				//errore
				return -1;				
			}
			i++;
		}		
	}
	free(exp);
	free(cof);
}
