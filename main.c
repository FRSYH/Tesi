#include <stdio.h>
#include <gsl/gsl_linalg.h>

int max_degree = 0;
int module = 0;

void swap_rows(gsl_matrix * m, int row, int col, int j, int i){
	
	int k;
	double tmp;

	for(k=0;k<row;k++){
		tmp = gsl_matrix_get(m,i,k);
		gsl_matrix_set(m,i,k, gsl_matrix_get(m,j,k));
		gsl_matrix_set(m,j,k, tmp);
	}

}

void sort_matrix (gsl_matrix * m, int row, int col){

	int i,j,*v,k;	
	double pivot;	

 	v = (int *) malloc(sizeof(int) * row);
	
	for(i=0;i<row;i++){
		j=0;
		do{
			pivot = gsl_matrix_get(m,i,j);
			j++;
		}while(pivot == 0.0);
		j--;
		v[i] = j;
		
	}
	
	for(i=0;i<row-1;i++){
		for(j=i+1;j<row;j++){
			if(v[j]<v[i]){
				k=v[i];
				v[i]=v[j];
				v[j]=k;
				swap_rows(m,row,col,j,i);
			}
		}	
	}
	free(v);
}


void print_matrix (gsl_matrix * m, int row, int col){
	
	int i,j;	
	
	printf("Matrice di Gauss\n");
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++)
		{
			if(j==col-1)
			{
				printf("| %g ",gsl_matrix_get(m,i,j));						
			}else{
				printf("%g ",gsl_matrix_get(m,i,j));						
			}	
		}
		printf("\n");
	}
	printf("\n");

}


void gauss_riduzione(gsl_matrix * m, int row, int col){

	int i,j,k,n,o;
	double x,y,q,pivot,tmp;
	
	sort_matrix(m,row,col);

	if( col < row ){
		n = col-1;
		o = row;
	}else{
		n = row;
		o = col-1;
	}

	for(j=0; j<n;j++) 
		{
			for(i=0; i<o; i++)
			{
				if(i>j && i<row)
				{         
					q = gsl_matrix_get(m,i,j) / gsl_matrix_get(m,j,j); //q=M[i][j]/M[j][j];
					for(k=0; k<col; k++)
					{
						x = gsl_matrix_get(m,i,k);
						y = gsl_matrix_get(m,j,k);
						tmp = q*y;
						gsl_matrix_set (m,i,k, x-tmp);      //M[i][k]=M[i][k]-q*M[j][k];
		
					}
				}
			}
		}


	//riduce i pivot a 1
	for(i=0;i<row;i++){
		j=0;
		do{
			pivot = gsl_matrix_get(m,i,j);
			j++;
		}while(pivot == 0.0);
		if(pivot != 1.0){
			for(k=j-1;k<col;k++){
				gsl_matrix_set (m,i,k, gsl_matrix_get(m,i,k)/pivot);			
			}
		}
	}


}

void riduzione(gsl_matrix * m, int row, int col, int riga_pivot, int j){
	
	int r,k;
	double s,tmp,a;
	printf("entro riduzione\n");
	printf("riga_pivot %d\n", riga_pivot);
	for( r=riga_pivot+1; r<row; r++ ){
		if( gsl_matrix_get(m,r,j) != 0.0 ){
			s = (1/gsl_matrix_get(m,riga_pivot,j))*gsl_matrix_get(m,r,j);
			for( k=0; k<col; k++ ){
				tmp = s*gsl_matrix_get(m,riga_pivot,k); 
				printf("tmp %g\n", tmp);
				a = gsl_matrix_get(m,r,k)-tmp;
				gsl_matrix_set(m,r,k, a);
			}
		}
	}
	
}


void gauss(gsl_matrix * m, int row, int col){
	
	int pivot_riga, pivot_colonna,righe_trovate,j;
	
	righe_trovate = -1;
	printf("entro gauss\n");
	for( j=0; j<col; j++){
		pivot_colonna = j;
		for( pivot_riga=(righe_trovate+1); pivot_riga<row; pivot_riga++ ){
			if( gsl_matrix_get(m,pivot_riga,pivot_colonna) != 0.0 ){
				riduzione(m,row,col,pivot_riga,pivot_colonna);
				righe_trovate++;
				swap_rows(m,row,col,righe_trovate,pivot_riga);
				break;
			}
		}
	}

}





int main (void)
{
	double a_data[] = { 1.0, 2.0, -1.0, 3.0, 7.0, 6.0, 
                        4.0, 3.0, 1.0, 2.0, 1.0, 1.0,
                        2.0, 2.0, 3.0, 8.0, 2.0, 4.0,
						3.0, 5.0, 1.0, 1.0, 2.0, 1.0,
						5.0, 3.0, 3.0, 4.0, 1.0, 2.0    };
					  	
	int row,col;

	row = 5;
	col = 6;	

	gsl_matrix_view m = gsl_matrix_view_array (a_data, row, col);

	gauss(&m.matrix, row, col);
	
	print_matrix(&m.matrix, row, col);

	return 0;
}
