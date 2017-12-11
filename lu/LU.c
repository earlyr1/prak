#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>


long double func(int i, int j, int n) {
	return (long double)(1.0 / (2 * n - i - j - 1.0)); 
}

long double L(long double ** A, int i, int j) {
    if (i > j) return A[i][j];
    if (i == j) return 1;
    return 0;
}

long double U(long double ** A, int i, int j) {
	if (i <= j) return A[i][j];
	return 0;
}

void forth(int n, long double **L, long double *y, long double *b) {
	long double k;
	y[0]=b[0] / L[0][0];
	for(int i = 1; i < n; i++) {
		k = 0;
		for(int j = 0; j <= i - 1; j++) {
			k += y[j] * L[i][j];
		}
		y[i] = (b[i] - k) / L[i][i];
	}
}

void back(int n, long double **U, long double *x, long double *y) {
	long double k;
	x[n - 1] = y[n - 1];
	for(int i = n - 2; i >= 0; i--) {
	    k = 0;
	    for(int j = n - 1; j > i; j--) {
	    	k += x[j] * U[i][j];
	    }
		x[i] = y[i] - k;
	}
	for(int j = 0; j < n; j++) {
		printf("%Lf\n", x[j]);
	}
}


void input_matrix(int * n, long double *** A, long double ** b) {
	printf("From file? Print y/n:\n");
	char c; 
	scanf("%c", &c);
	FILE * f;
	if (c == 'y') {
		printf("Enter filename:\n");
		char filename[100];
		scanf("%s", filename);
		f = fopen(filename, "r");
		fscanf(f, "%d", n);
	} else {
		printf("Enter n:\n");
		scanf("%d", n);
	}
	*A = malloc((*n) * sizeof(long double*));
	*b = malloc((*n) * sizeof(long double));
	for (int i = 0; i < (*n); i++) {
		(*A)[i] = malloc((*n) * sizeof(long double)); 
	}
	if (c == 'y') {
		for(int i = 0; i < (*n); ++i) {
			for(int j = 0; j < *n; ++j) { 
				fscanf(f, "%Lf", &(*A)[i][j]);
			}
			fscanf(f, "%Lf", &(*b)[i]);
		}
	} else {
		for(int i = 0; i < (*n); ++i) {
			long double sum = 0;
			for(int j = 0; j < *n; ++j) { 
				(*A)[i][j] = func(i, j, *n);
				sum += func(i, j, *n);
			}
			(*b)[i] = sum;
		}
	}
	for(int i = 0; i < (*n); ++i) {
			for(int j = 0; j < *n; ++j) { 
				printf("%Lf ", (*A)[i][j]);
			}
			printf("%Lf\n", (*b)[i]);
		}
	printf("\n");
}


void LU_decomp(int n, long double **A, long double *b, long double * b1, long double *x) {
	long double * y = malloc(n * sizeof(long double));
  	for (int k = 0; k < n; k++) {
      	for (int i = 0; i < n; i++) {
      		//printf(" %Lf ", A[i][j]);
          	if (i >= k) {
          		long double sum = 0;
              	for (int j = 0; j < k; j++) {
                  	sum += A [i][j] * A [j][k]; 
              	} 
              	A[i][k] -= sum;
          	} 
       	}
       	int idmax = k;
		for(int j = k; j < n; j++) {
			if (fabs(A[j][k]) > fabs(A[idmax][k])) {
				idmax = j;
			}
		} 
		/*
		for(int j = i; j < n; j++) {
			long double tmp = A[i][j];
			A[i][j] = A[idmax][j];	
			A[idmax][j] = tmp;
		}*/
        long double *tmp = A[k];
		A[k] = A[idmax];	
		A[idmax] = tmp;

		long double tmp1 = b[k];
		b[k] = b[idmax];
		b[idmax] = tmp1;

		long double tmp2 = b1[k];
		b1[k] = b1[idmax];
		b1[idmax] = tmp1;
      	for(int i = 0; i < n; i++) {
      		if (k < i) {
   				long double sum = 0;
            	for (int j = 0; j < k; j++) {
                	sum += A [k][j] * A [j][i];
            	}
            	A[k][i] -= sum;
				if (A[k][k] != 0) {
					A [k][i] /= A [k][k];
				} else {
					printf("Wrong matrix");
					exit(0);
				}
			}
   		}      	
   	}


   	printf("\n\n");
   	for (int i = 0; i < n; i++) {
     	for (int j = 0; j < n; j++)
    		printf("  %Lf", A[i][j]);
     	printf("\n\n");
   	}
}


int main() {
	int n;
	long double ** A;
	long double * b;
	input_matrix(&n, &A, &b);
	long double * b1 = malloc(n * sizeof(long double));
	for(int i = 0; i < n; i++) {
		b1[i] = b[i];
	}
	long double * x = malloc(n * sizeof(long double));
	long double * y = malloc(n * sizeof(long double));
	
	LU_decomp(n, A, b, b1, x);
	forth(n, A, y, b);
	back(n, A, x, y);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(i == j) { 
				y[i] +=x[j];
			}
			else if (j > i) {
				y[i]+=A[i][j]*x[j];
			}
		}
	}


	for(int i = 0; i < n; i++) {
		x[i] = 0;
		for(int j = 0; j < n; j++) {
			if (i >= j) {
				x[i] += A[i][j] * y[j];
			}
		}
	}
	
	
	long double s = 0;
	for(int i = 0; i < n; i++)
	{

		x[i] -= 2 * b1[i];
		s+= x[i] * x[i];
	}
	printf("Невязка %Lf", (s));
// ДОДЕЛАТЬ норм невязку и чтоб не расходилось так сильно
}



