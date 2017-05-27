#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "scan.h"
#include "linalg.h"

void allocation(long long ***m, int *row, int *col, int *num_var, char **v, int *n, long long *module, int *max_degree){
/*
Legge da input le seguenti informazioni:
	- modulo dei coefficienti
	- grado massimo
	- numero dei polinomi di partenza
	- tipo di ordinamento
	- variabili utilizzate nei polinomi


con queste informazioni alloca la matrice principale (matrice che conterrà i polinomi) e stabilisce il numero di variabili utilizzate.
*/
	scanf("%lli",module); //leggo il modulo
	getchar();
	scanf("%d",max_degree); //leggo il grado massimo
	getchar();
	scanf("%d",row);  //leggo numero dei polinomi di partenza
	getchar();
	scanf("%d",n);  //leggo tipo di ordinamento
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

	*col = 1+monomial_combinations(*num_var, *max_degree);

	*m = malloc((*row) * sizeof (long long *) );            // allocazione della matrice dei coefficienti
	if( *m != NULL )
		for (int i=0; i<(*row); i++)
			(*m)[i] = calloc((*col) , sizeof (long long) );	


}



int parse(int num_var, char *vet, long long **m, int **vet_grd, int len, long long module, int (*ord) (const void *, const void *, void*) ){
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
			if( parse_mon(mon,i,&cof,num_var,vet,grade,pos_pol,module) == -1 ){
				return -1;
			}
			//inserire monomio in posizione corretta
			col = (int **)(bsearch_r((void *) &grade, (void *) vet_grd, len, (sizeof(int*)), ord, &num_var)) - vet_grd;
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

int parse_mon(char * mon, int len,long long * val, int num_var, char *vet, int *grade, int pos_pol, long long module){
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
	if( isdigit(mon[i]) != 0 ){  // se il primo carattere letto è un numero
		i = 0;
		while( isdigit(mon[i]) && i<len){ //procedo alla lettura del coefficiente
			cof = realloc(cof, (i+1) * sizeof(char));
			cof[i] = mon[i];
			i++;                   
		}
		*val = mod(atoll(cof),module);  //a questo punto ho letto il coefficiente e lo riduco in modulo se necessario

		
	}else{ // altrimenti se il primo carattere letto non è un numero ma una lettera
		if( isalpha(mon[i]) != 0 ){
			*val = 1;   // significa che il coefficiente del monomio è 1
			
		}else{ 
			// se il primo carattere non è ne un numero ne una lettera allora è un carattere non valido -> ERRORE
			return -1;  //terminazione parse con codice di errore -1
		}	
	}
	if( i < len ){   
		while( i < len ){
			if( mon[i] == '*' ){ 
				i++;    
			}
			if( i<len && (isalpha(mon[i]) != 0) ){  //cerco le variabili presenti
				for(k=0; k<num_var; k++){
					if( vet[k] == mon[i] ){
						pos_var = k;
						break;			//quando trovo una variabile salvo la sua posizione e interrompo
					}
				}
				i++;
				if( i<len ){  //controllo se la variabile ha un grado
					if( mon[i] == '^' ){ //se trovo il carattere di elevamento a potenza allora c'è un grado
						i++;
						if( isdigit(mon[i]) != 0 ){
							k = 0;
							while( isdigit(mon[i]) != 0 && i < len ){
								exp = realloc(exp, (k+1) * sizeof(char));
								exp[k] = mon[i];
								i++;
								k++;
							}
						grade[pos_var] = atoi(exp); //leggo il grado della variabile e lo salvo nel 
													//vettore che rappresenta tutte le possibili variabili nella posizione adeguata
						}else{
							// se il carattere successivo al carattere di elevamento a potenza non è un numero allora errore di input
							return -1;									
						}
					}else{  // se non esiste il simbolo di elevamento a potenza allora significa che il grado della variabile è 1
						grade[pos_var] = 1;	
					}
				}else{ // se termino il monomio con una varibile senza grado allora il suo grado è 1
					grade[pos_var] = 1;
				}								
			}else{
				// se dopo il simbolo * non trovo un carattere alfabetico che indica una variabile errore di input
				return -1;				
			}
			i++;
		}		
	}
	free(exp);
	free(cof);
}


int order(int (**ord) (const void *, const void *, void*), int n){
//inizializza il puntatore ord alla funzione di ordinamento adeguata. Il numero n indica quale funzione scegliere.

	switch(n){

		case 0:
			*ord = grevlex_comparison;
			return 0;
			break;

		default:
			return -1;
			break;	

	}
}


