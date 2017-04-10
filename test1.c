#define _GNU_SOURCE 
//define needed to use qsort_R without warnings from the compiler (no c99 ansi)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>

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


int main(void) {


	int n = 3, m = 7, **vet, len, **map;
	
	
	len = 1+monomial_combinations(n, m);
	
	
	vet = malloc(sizeof(int*) * len);
	for (int i = 0; i < len; i++)
		vet[i] = malloc(n * sizeof(int));
	
	int *mon = malloc(n*sizeof(int));
	
	monomial_computation(n, m, vet, 0, mon);

	free(mon);

	qsort_r(vet, len, sizeof(int*), grevlex_comparison, &n);
	
	
	map = malloc(len * sizeof(int *));
	for (int i = 0; i < len; i++)
		map[i] = malloc(len * sizeof(int));

	
	
	
	/*
	for (int i = 0; i < len; i++) {
		for (int k = 0; k < n; k++)
			printf("%d ",vet[i][k] );
		printf("\n");
	}*/
	
	setup_map(map, vet, len, n, m);
	
	for (int row = 0; row < len; row++){
		for (int col = 0; col < len; col++)
			printf("%d ", map[row][col]);
		printf("\n");
	}
	
	

	for (int i = 0; i < len; i++) {
		free(vet[i]);
	}
	free(vet);

	return 0;
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
