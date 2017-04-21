#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


void parse_mon(char * mon, int len,int * val, int num_var, char *vet, int *grade, int pos_pol){

	int i,k,pos_var;
	char c,* cof,*var;
	cof = malloc( sizeof(char) );
	var = malloc( sizeof(char) );
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
				//printf("* in posizione %d\n", i);         
				i++;    
			}
			if( i<len && (isalpha(mon[i]) != 0) ){
				for(k=0; k<num_var; k++){
					if( vet[k] == mon[i] ){
						pos_var = k;
						//printf("posizione varaibile %d\n",pos_var);
						break;
					}
				}
				i++;
				if( i<len ){
					if( mon[i] == '^' ){ //ho trovato il grado della variabile
						i++;
						if( isdigit(mon[i]) != 0 ){
							var[0] = mon[i];
							grade[pos_var] = atoi(var);
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
	free(var);
	free(cof);
}


void parse(int num_var, char *vet){
	
	int pos_pol = 0,i,cof;
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
			printf("%s   cof:%d   grade: %d %d %d \n", mon,cof,grade[0],grade[1],grade[2]);

			//inserire monomio in posizione corretta


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





int main (void){

	int i,j,k,num_var = 0,pos_pol;
	char c,*v;
	i=0;
	pos_pol = 0;
	v = malloc(sizeof(char));
	c = getchar();
	while( c != '\n' ){
		v[i] = c;
		i++;
		num_var++;
		v = realloc(v, (i+1)*sizeof(char) );
		c = getchar();
	}
	parse(num_var,v);
	
	free(v);

	return 0;
}