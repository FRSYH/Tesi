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
void moltiplica_riga(int ***m, int *row, int col, int riga, int **map,int * degree, int * degree_position); 

//initialize the vector that keeps the number of monomial with the same grade and their position
void init_degree_vector(int * degree, int * degree_position, int num_var);

//return the grade of a monomial
int grado_monomio(int posizione, int *degree_position);

int combination(int n, int k);

//returns the number of all possible monomial with n variables and degree <= m
int monomial_combinations(int n, int m); 

//computes all possible monomials with n variables and degree <= m
//saves them into array vet, see the definition for more informations
void monomial_computation(int n, int m,  int **vet, int turn, int *monomial);

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
void moltiplica_matrice(int ***m, int *row, int col, int **map, int *degree, int *degree_position);

//compute the different degree of the matrix rows
void matrix_degree(int **m, int row, int col, int *m_deg, int *degree_position);

//formatted print of matrix degree
void print_matrix_degree(int *m_deg);

//compute the number of null rows (rows of full 0)
int null_rows(int **m, int row, int col);

//eliminate the null rows of the matrix (reallocation)
void eliminate_null_rows(int ***m, int *row, int col);

//execute the multiplication, gauss reduction, elimination null_rows and compare the result to target
void execute(int ***m, int *d_row, int col, int **map, int *degree, int *degree_position);

//compare the matrix degree with the target degree wich are {0,1,2,3,4,5...max_degree} return 0 if equal -1 if not
int target_degree(int *v);

int main (void){

	max_degree = 7;

	int row, col, i, num_var, degree[max_degree+1], degree_position[max_degree+1];
	int **m, *d_row;
	
	row = 4;
	col = 120;
	module = 773;	
	num_var = 3;
	d_row = &row;



//########################################################################
	//map allocation and creation

	int **vet, len, **map;

	len = 1+monomial_combinations(num_var, max_degree);

	vet = malloc(sizeof(int*) * len);

	for (int i = 0; i < len; i++)
		vet[i] = malloc(num_var * sizeof(int));
	
	int *mon = malloc(num_var*sizeof(int));

	monomial_computation(num_var, max_degree, vet, 0, mon);

	free(mon);

	qsort_r(vet, len, sizeof(int*), grevlex_comparison, &num_var);

	map = malloc(len * sizeof(int *));
	for (int i = 0; i < len; i++)
		map[i] = malloc(len * sizeof(int));

	setup_map(map, vet, len, num_var, max_degree);           // a questo punto posso utilizzare la mappa


//###########################################################################
	m = malloc(row * sizeof (int *) );            // allocazione della matrice
	if( m != NULL ){
		for (i=0; i<row; i++)
		{
			m[i] = calloc(col , sizeof (int) );
		}

	}

//##########################################################################
	
	//RISOLUZIONE PROBLEMA

	init_matrix(m,row,col); 

	init_degree_vector(degree,degree_position,num_var);

	execute(&m,d_row,col,map,degree,degree_position);

	print_matrix(m, row, col);	

	
	

//##################################################################  	

	for (i=0; i<row; i++)      // deallocazione matrice
	{
		free(m[i]);
	}
	free(m);	

	for (int i = 0; i < len; i++) {
		free(vet[i]);
	}
	free(vet);	

	for (int i = 0; i < len; i++) {
		free(map[i]);
	}
	free(map);	

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
	
	int pivot_riga, pivot_colonna,righe_trovate,j,n_col;
	
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


void riduzione(int **m, int row, int col, int riga_pivot, int j){
	
	int r,k,s,inv,a,b;

	for( r=riga_pivot+1; r<row; r++ ){
		if( m[r][j] != 0 ){
			inv = invers(m[riga_pivot][j],module);			//calcola l'inverso moltiplicativo di m[riga_pivot][j] nel campo indicato
			s = mul_mod(inv,m[r][j],module);
			for( k=0; k<col; k++ ){
				a = mul_mod(s,m[riga_pivot][k],module);
				b = m[r][k];
				m[r][k] = sub_mod(m[r][k],a,module);
			}
		}
	}
	
}


void swap_rows(int **m, int row, int col, int j, int i){
	
	int k,tmp;
	for(k=0;k<col;k++){
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

	//prima riga
	m[0][0] = 640;
	m[0][1] = 328;
	m[0][2] = 328;
	m[0][3] = 328;
	m[0][5] = 431;
	m[0][6] = 431;
	m[0][8] = 431;
	m[0][23] = 1;
	m[0][24] = 771;
	m[0][25] = 1;
	m[0][27] = 771;
	m[0][28] = 771;
	m[0][32] = 1;
	//seconda riga
	m[1][0] = 131;
	m[1][3] = 104;
	m[1][9] = 64;
	m[1][19] = 413;
	m[1][34] = 1;
	//terza riga
	m[2][0] = 356;
	m[2][2] = 486;
	m[2][7] = 294;
	m[2][16] = 657;
	m[2][30] = 1;
	//quarta riga
	m[3][0] = 691;
	m[3][1] = 52;
	m[3][4] = 393;
	m[3][10] = 1;
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


void moltiplica_riga(int ***m, int * row, int col, int riga, int **map,int * degree, int * degree_position){
//moltiplica la riga indicata per ogni monomio in modo tale che il prodotto abbia grado <= del grado massimo

	int grado_massimo_riga, grado_massimo_monomio, grado_a, grado_b, grado_prodotto;
	int i,j,k,last, posizione_nuovo_monomio, offset, offset_a, offset_b,v,value;
	int new_row;
	
    v = 0;
	last = -1;
	//cerco la posizione dell'ultimo coefficiente non nullo del polinomio rappresentato nella riga, grado più alto
	for(i=col-1; i>0; i--){
		v = (*m)[riga][i];		
		if( v != 0 ){
			last = i;
			break;
		}
	}

	//risalgo al grado del monomio appena trovato
	//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado
	if( last != -1 ){

		grado_massimo_riga = grado_monomio(last,degree_position);
		
		//calcolo il grado massimo che deve avere il monomio per cui moltiplicare		
		grado_massimo_monomio = max_degree - grado_massimo_riga;		
	
		// a questo punto conosco per quanti monomi devo moltiplicare e quindi 
		// conosco il numero di righe che devo aggiungere alla matrice

		new_row = degree_position[grado_massimo_monomio+1]-1; //righe da aggiungere 

		*m = realloc( *m , (*row + new_row ) * sizeof (int *));
		for (i=(*row); i< (*row + new_row ); i++)
		{
			(*m)[i] = calloc(col , sizeof (int) );
		}
	
		for(i=1; i<degree_position[grado_massimo_monomio+1]; i++){     //scorre tutti i gradi per i quali posso moltiplicare
			
			for(j=0; j<(last+1); j++){     //scorre fino all'ultimo elemento della riga
				(*m)[*row][ map[i][j] ] = (*m)[riga][j];  //shift nella posizione corretta indicata dalla mappa
			}
			*row = *row + 1;			//aumento del conteggio delle righe
				
		}

	}

}


int grado_monomio(int posizione, int *degree_position){
	
	int i,grado;
	grado = 7;
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



//maps all the possible multiplications of the monomials of n variables and
//degree <= m in the array vet of length len into the matrix map[len][len]
//sets the value to -1 if the multiplication exceeds vet
void setup_map(int **map, int **vet, int len, int n, int m) {

	int sum, *temp = malloc(n * sizeof(int));

	for (int row = 0; row < len; row++)
		for (int col = 0; col < len; col++) {
			sum = 0;
			for (int v = 0; v < n; v++) {
				temp[v] = vet[row][v] + vet[col][v];
				sum += temp[v];
			}
			if (sum > m) {
				for (int i = col; i < len; i++)
					map[row][i] = -1;
				break;	
			}
			else
				map[row][col] = (int **)(bsearch_r((void *) &temp, (void *) vet, len, (sizeof(int*)), grevlex_comparison, &n)) - vet;
		}
		
	free(temp);
}


//https://git.devuan.org/jaretcantu/eudev/commit/a9e12476ed32256690eb801099c41526834b6390
//missing in stdlib, counterpart of qsort_r
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

//compares two monomials of n variables following the grevlex order
//returns a positive number if mon1 > mon2, 0 if they are equal, negative otherwise
//various casts are made to achive compatibility with qsort_r
int grevlex_comparison(const void *monom1, const void *monom2, void *arg) {

	int degree1 = 0, degree2 = 0, n, *mon1, *mon2;
	n = *((int *) arg);
	mon1 = *((int **) monom1);
	mon2 = *((int **) monom2);
	
	for (int v = 0; v < n; v++) {
		degree1 += mon1[v];
		degree2 += mon2[v];
	}
	
	if (degree1 > degree2)
		return 1;
	else if (degree1 < degree2)
		return -1;
	else {
		int *temp = malloc(n * sizeof(int));
		for (int v = 0; v < n; v++)
			temp[v] = mon1[v] - mon2[v];
		for (int v = (n-1); v >= 0; v--) {
			if (temp[v] != 0) {
				return -temp[v];
				free(temp);
				//to avoid freeing the same pointer twice
				temp = NULL;
			}
		}
				free(temp);
	}
	
	return 0;
	
}


/*recursive function that computes all possible monomials with n variables
and degree <= m, saves them into array vet. A monomial is represented as an
array of int, where each positions describes the degree of the variable.
The array has to be initialized with the correct length.
Other parameters are need for recursive structure and should always be
turn = 0 represents the position of the variable in the monomial
monomial = array of int of length n used to compute the various permutations
*/
void monomial_computation(int n, int m, int **vet, int turn, int *monomial) {

	//s keeps track of the monomial already saved in the array vet
	//static becasue it must keep values through the recursive calls
	//no thread safe
	static int s = 0;

	//for every variable try all degrees from 0 to m
	for (int degree = 0; degree <= m; degree++) {
		//if this is the first variable all other variables have looped
		if (turn == 0) {
			//reset the monomial with only your degree
			monomial[0] = degree;
			for (int v = 1; v < n; v++)
				monomial[v] = 0;		
		}
		//other variables add their values to the monomial
		else
			monomial[turn] = degree;
		
		
		//get total degree of monomial by adding the variables degrees together
		int sum = 0;
		for (int v = 0; v <= turn; v++)
			sum += monomial[v];
		//if the degree of the monomial exceeds the maximum it's pointless to continue
		//looping searching for other monomials since they all have degree >=
		if (sum > m)
			break;


		//if this is the last variable and we haven't broke out with the previous condition
		//copy this monomial to the vecor vet and increase the index to the vector
		if (turn == (n-1)) {
			vctcpy(vet[s], monomial, n);
			s+=1;
		}
		//otherwise recursively call itself changing the variable (turn)
		else
			monomial_computation(n, m, vet, turn+1, monomial);
	}
	
	//before finishing the function (the first variable has finished looping)
	//set the array index to zero so that the function can be called again
	//because it wouldn't reset if called again since it declared static 
	if (turn == 0)
		s = 0;
	return;
}


//copies vector vet2 into vet1 of length len
void vctcpy(int *vet1, const int *vet2, int len) {
	for (int i = 0; i < len; i++)
		vet1[i] = vet2[i];
	return;
}


void moltiplica_matrice(int ***m, int *row, int col, int **map, int *degree, int *degree_position){
	
	int n,i;
	n = *row;    //n conta il numero di righe della matrice di partenza che devo moltiplicare
	for(i=0; i<n; i++){
		moltiplica_riga(m,row,col,i,map,degree,degree_position);
	}

}

void matrix_degree(int **m, int row, int col, int *m_deg, int *degree_position){
	
	int i,j,last,grado;
	for(i=0; i<row; i++){
		for(j=col-1; j>0; j--){
			if( m[i][j] != 0 ){
				last = j;           //posizione dell'ultimo coefficiente della riga
				break;
			}
		}

		grado = grado_monomio(last,degree_position);	
		m_deg[grado] = 1;		
	}
}



void print_matrix_degree(int *m_deg){
	int i;	
	printf("Gradi della matrice = {");
	for(i=0; i<max_degree+1; i++){
		if( m_deg[i] != 0 ){
			printf(" %d ",i);
		
		}
	}
	printf("}\n");
}



int null_rows(int **m, int row, int col){
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
		if( last == -1 ){
			null_rows++;
		}
	}
	return null_rows;
}


void eliminate_null_rows(int ***m, int *row, int col){
	
	int i,null_row;
	null_row = null_rows(*m,*row,col);
	*m = realloc( *m , (*row - null_row ) * sizeof (int *));
	*row = *row - null_row;

}


void execute(int ***m, int *d_row, int col, int **map, int *degree, int *degree_position){

	int *m_deg;

	m_deg = calloc(max_degree+1, sizeof(int));

	matrix_degree(*m,*d_row,col,m_deg,degree_position);

	while( target_degree(m_deg) != 0 ){
		
		moltiplica_matrice(m,d_row,col,map,degree,degree_position);
	
		gauss(*m, *d_row, col);

		eliminate_null_rows(m,d_row,col);

  		matrix_degree(*m,*d_row,col,m_deg,degree_position);

		printf("Numero righe: %d\n", *d_row);

		print_matrix_degree(m_deg);
		
	}

//	printf("numero righe %d righe nulle %d\n",*d_row,null_rows(*m,*d_row,col));

	free(m_deg);
}


int target_degree(int *v){
	
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











