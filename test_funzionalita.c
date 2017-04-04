#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

int max_degree = 0;

int combination(int n, int k){
	
	double a,b,c;
	a = gsl_sf_fact(n+k-1);
	b = gsl_sf_fact(k);
	c = gsl_sf_fact((n+k-1)-k);
	return (int) a/(c*b);
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


void print_matrix (gsl_matrix * m, int row, int col){
	
	int i,j;	
	
	printf("Matrice di Gauss\n");
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++)
		{
			printf("%g ",gsl_matrix_get(m,i,j));						
	
		}
		printf("\n");
	}
	printf("\n");

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


void moltiplica_riga(gsl_matrix * m, int row, int col, int riga, int * degree, int * degree_position){
//moltiplica la riga indicata per ogni monomio in modo tale che il prodotto abbia grado <= del grado massimo

	int grado_massimo_riga, grado_massimo_monomio;
	int i,j,k,last;
	double v;	
	
    v = 0.0;
	last = -1;
	//cerco la posizione dell'ultimo coefficiente non nullo del polinomio rappresentato nella riga, grado più alto
	for(i=col-1; i>0; i--){
		v = gsl_matrix_get(m,riga,i);		
		if( v != 0.0 ){
			last = i;
		}
	}
	
	//risalgo al grado del monomio appena trovato
	//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado
	if( last != -1 ){
		grado_massimo_riga = grado_monomio(last,degree_position);
		
		//calcolo il grado massimo che deve avere il monomio per cui moltiplicare		
		grado_massimo_monomio = max_degree - grado_massimo_riga;		
		
		printf("grado riga %d grado monomio %d\n", grado_massimo_riga,grado_massimo_monomio);
		//moltiplico la riga per ogni monomio possibile
		for(i=1; i<grado_massimo_monomio; i++){
			for(j=0; j<degree[i]; j++){
				for(k=0; k<col; k++){
					
					
				}
			
			}
		}

	}


}

int main(void){
	
	max_degree = 7;
	
	int i,*row_p,col,row,num_var,row_max;

	int degree[max_degree+1],degree_position[max_degree+1];

	row_max = 100;
	row = 4;
	row_p = &row;   //numero di righe variabile
	col = 10;        //numero di colonne fisse
	num_var = 3;    //numero di variabili 3 -> x,y,z
	

	gsl_matrix * m = gsl_matrix_alloc (row_max, col);   //allocazione matrice

	init_degree_vector(degree,degree_position,num_var); //inizializza i vettori che vanno utilizzati per la moltiplicazione della matrice

	gsl_matrix_set(m,0,7,5.0);

	print_matrix(m,*row_p,col);		


	moltiplica_riga(m,*row_p,col,0,degree,degree_position);

	return 0;
}
