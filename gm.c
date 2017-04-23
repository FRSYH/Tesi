#define _GNU_SOURCE 
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>
#include <ctype.h>


int max_degree = 0;
int module = 0;

long long mod(long long n, long long p); //n mod p

long long add_mod(long long a, long long b, long long p); // a + b mod p

long long sub_mod(long long a, long long b, long long p); // a - b mod p

long long mul_mod(long long a, long long b, long long p); // a * b mod p

long long invers(long long n, long long p);  // inverso moltiplicativo in mod p

void gauss(long long **m, int row, int col); // riduzione di gauss della matrice m

void riduzione(long long **m, int row, int col, int riga_pivot, int j);

void swap_rows(long long **m, int row, int col, int j, int i);  // scambia tra di loro due righe della matrice

void print_matrix (long long **m, int row, int col); // stampa la matrice

void init_matrix(long long **m, int row, int col,int **vet_grd, char *v, int num_var); //inizializza la matrice dei coefficienti

//moltiplica la riga indicata per tutti i possibili monomi
void moltiplica_riga(long long ***m, int *row, int col, int riga, int **map,int * degree, int * degree_position); 

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
void moltiplica_matrice(long long ***m, int *row, int col, int **map, int *degree, int *degree_position);

//compute the different degree of the matrix rows
void matrix_degree(long long **m, int row, int col, int *m_deg, int *degree_position);

//formatted print of matrix degree
void print_matrix_degree(int *m_deg);

//compute the number of null rows (rows full of 0)
int null_rows(long long **m, int row, int col);

//eliminate the matrix null rows (reallocation - resize)
void eliminate_null_rows(long long ***m, int *row, int col);

//execute the multiplication, gauss reduction, elimination null_rows and compare the result to target
void execute(long long ***m, int *d_row, int col, int **map, int *degree, int *degree_position);

//compare the matrix degree with the target degree wich are {0,1,2,3,4,5...max_degree} return 0 if equal -1 if not
int target_degree(int *v);

void allocation(long long ***m, int *row, int *col, int *num_var, char **v);

void matrix_free_long(long long ***m, int row, int col);

void matrix_alloc_int(int ***m, int row, int col);

void matrix_free_int(int ***m, int row, int col);

void parse(int num_var, char *vet, long long **m, int **vet_grd, int len);

void parse_mon(char * mon, int len,int * val, int num_var, char *vet, int *grade, int pos_pol);

int main (void){
	

	int row, col, num_var;
	int *d_row,**vet, len, **map;
	long long **m;
	char *v;
	row = col = num_var = 0;
	allocation(&m,&row,&col,&num_var,&v);
	d_row = &row;

	int degree[max_degree+1], degree_position[max_degree+1];
	len = col;
	matrix_alloc_int(&vet,len,num_var);

	int *mon = malloc(num_var*sizeof(int));

	monomial_computation(num_var, max_degree, vet, 0, mon);

	free(mon);

	qsort_r(vet, len, sizeof(int*), grevlex_comparison, &num_var);

	init_matrix(m,row,col,vet,v,num_var); //inizializzazione matrice

	matrix_alloc_int(&map,len,len);
	setup_map(map, vet, len, num_var, max_degree);     // a questo punto posso utilizzare la mappa

	//RISOLUZIONE PROBLEMA
	init_degree_vector(degree,degree_position,num_var); 
	execute(&m,d_row,col,map,degree,degree_position);  //soluzione trovata
	print_matrix(m, row, col);	 //stampa la matrice soluzione

	matrix_free_long(&m,row,col);
	matrix_free_int(&map,len,len);
	matrix_free_int(&vet,len,num_var);

	return 0;
}

long long mod(long long n, long long p){
	long long v = n;
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


long long add_mod(long long a, long long b, long long p){
	return mod((a+b),p);
}


long long sub_mod(long long a, long long b, long long p){
	return mod((a-b),p);
}


long long mul_mod(long long a, long long b, long long p){
	return mod((a*b),p);
}

//compute the Gauss reduction over the matrix m with specified row and column size.
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

//compute the reduction of a single row of the matrix int the gaussian reduction
void riduzione(long long **m, int row, int col, int riga_pivot, int j){
	
	int r,k;
	long s,inv,a;
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


void swap_rows(long long **m, int row, int col, int j, int i){
	
	int k;
	long long tmp;
	for(k=0;k<col;k++){
		tmp = m[i][k];
		m[i][k] = m[j][k];
		m[j][k] = tmp;
	}
}


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


void init_matrix(long long **m, int row, int col, int **vet_grd, char *v, int num_var){

	parse(num_var,v,m,vet_grd,col);

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


void moltiplica_riga(long long ***m, int * row, int col, int riga, int **map,int * degree, int * degree_position){
//moltiplica la riga indicata per ogni monomio in modo tale che il prodotto abbia grado <= del grado massimo

	int grado_massimo_riga, grado_massimo_monomio,i,j,last,new_row;
	last = -1;
	//cerco la posizione dell'ultimo coefficiente non nullo del polinomio rappresentato nella riga, grado più alto
	for(i=col-1; i>0; i--){
		if( (*m)[riga][i] != 0 ){
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
		*m = realloc( *m , (*row + new_row ) * sizeof (long long *));
		for (i=(*row); i< (*row + new_row ); i++)
			(*m)[i] = calloc(col , sizeof (long long) );
		for(i=1; i<degree_position[grado_massimo_monomio+1]; i++){     //scorre tutti i gradi per i quali posso moltiplicare
			for(j=0; j<(last+1); j++)     //scorre fino all'ultimo elemento della riga
				(*m)[*row][ map[i][j] ] = (*m)[riga][j];  //shift nella posizione corretta indicata dalla mappa
			*row = *row + 1;			//aumento del conteggio delle righe
		}
	}
}


int grado_monomio(int posizione, int *degree_position){
	
	int i,grado;
	grado = max_degree;
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


void moltiplica_matrice(long long ***m, int *row, int col, int **map, int *degree, int *degree_position){
	
	int n,i;
	n = *row;    //n conta il numero di righe della matrice di partenza che devo moltiplicare
	for(i=0; i<n; i++)
		moltiplica_riga(m,row,col,i,map,degree,degree_position);
}

void matrix_degree(long long **m, int row, int col, int *m_deg, int *degree_position){
	
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
	for(i=0; i<max_degree+1; i++)
		if( m_deg[i] != 0 )	printf(" %d ",i);		
	printf("}\n");
}



int null_rows(long long **m, int row, int col){
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

	int null_row = null_rows(*m,*row,col);
	*m = realloc( *m , (*row - null_row ) * sizeof (long long *));
	*row = *row - null_row;
}


void execute(long long ***m, int *d_row, int col, int **map, int *degree, int *degree_position){

	int *m_deg = calloc(max_degree+1, sizeof(int));
	printf("Inizio procedura\n");
	matrix_degree(*m,*d_row,col,m_deg,degree_position);

	int flag,old,new;
	flag = old = new = 0;
	old = *d_row;
	while( flag != 1 ){
		printf("\n -Eseguo moltiplicazione, ");
		moltiplica_matrice(m,d_row,col,map,degree,degree_position);  //moltiplico la matrice per tutti i monomi possibili
		printf("numero righe: %d\n", *d_row);
		printf(" -Eseguo Gauss, ");	
		gauss(*m, *d_row, col);                                     //applico la riduzione di Gauss
		eliminate_null_rows(m,d_row,col);							//elimino le righe nulle della matrice
  		matrix_degree(*m,*d_row,col,m_deg,degree_position);
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

	scanf("%d",&module); //leggo il modulo
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
	for (int i=0; i<row; i++)      
		free((*m)[i]);
	free(*m);	
}

void matrix_alloc_int(int ***m, int row, int col){
	*m = malloc(row * sizeof (int *) );            // allocazione della matrice dei coefficienti
	if( *m != NULL )
		for (int i=0; i<row; i++)
			(*m)[i] = calloc(col , sizeof (int) );	
}

void matrix_free_int(int ***m, int row, int col){
	for (int i=0; i<row; i++)      
		free((*m)[i]);
	free(*m);	
}

void parse(int num_var, char *vet, long long **m, int **vet_grd, int len){
	
	int pos_pol = 0,i,cof,col;
	char c,* mon;
	cof = 0;
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
			parse_mon(mon,i,&cof,num_var,vet,grade,pos_pol);
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

void parse_mon(char * mon, int len,int * val, int num_var, char *vet, int *grade, int pos_pol){

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
		*val = atoi(cof);
	}else{
		*val = 1;
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
						}
					}else{
						grade[pos_var] = 1;	
					}
				}else{
					grade[pos_var] = 1;
				}								
			}
			i++;
		}		
	}
	free(exp);
	free(cof);
}