#include <stdio.h>
#include <gsl/gsl_linalg.h>

int main (void)
{
	double a_data[] = { 1.0, 2.0, -1.0, 2.0, 
                        4.0, 3.0, 1.0,  3.0,
                        2.0, 2.0, 3.0, 5.0};
					  
	int s,r,c,i,j,row,col,n,k;
	double x,y,q;

	row = 3;
	col = 4;
	n=row;	
	gsl_matrix_view m = gsl_matrix_view_array (a_data, row, col);

	for(j=0; j<n; j++) 
		{
			for(i=0; i<n; i++)
			{
				if(i>j)
				{         
					q = gsl_matrix_get(&m.matrix,i,j) / gsl_matrix_get(&m.matrix,j,j); //q=M[i][j]/M[j][j];
					for(k=0; k<n+1; k++)
					{
						x = gsl_matrix_get(&m.matrix,i,k);
						y = gsl_matrix_get(&m.matrix,j,k);
						gsl_matrix_set (&m.matrix,i,k, x-(q*y));      //M[i][k]=M[i][k]-q*M[j][k];
					}
				}
			}
		}


	printf("Matrice di Gauss\n");
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++)
		{
			if(j==col-1)
			{
				printf("| %f ",gsl_matrix_get(&m.matrix,i,j));						
			}else{
				printf("%f ",gsl_matrix_get(&m.matrix,i,j));						
			}	
		}
		printf("\n");
	}
	printf("\n");
 
	return 0;
}
