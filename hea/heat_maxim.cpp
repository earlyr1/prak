#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

//using namespace std;

double f( double x, double t, double a ){
	return 2*(t*(x*x-x) - a*a*(t*t+1)) ;
}

double U( double x, double t ){
	return (t*t+1) * (x*x-x);
}

double bordU( int x, double t ){
	if ( x == 0 )
		return U(0,t);
	if ( x == 1 )
		return U(1,t);
	else
		printf("Неправильная граница\n");
	exit(1);
}

double U0 ( double x ){
	if ( x > 1 || x < 0 )
		printf("Выход за границы [0,1]\n");
	return U(x,0);
}

int main ( int argc, char ** argv ){
	double h, a, Thau, M, N;
	double * prevLayer, * currLayer, * endLayer;
	FILE * err, * dat;
	int filter;
	int mode;
	if ( !argv[1] ){
		printf("Не указан файл для вывода оценки погрешности\n");
		exit(1);
	}
	if ( !argv[2] ){
		printf("Не указан файл для записи последовательных приближений ф-ии\n");
		exit(1);
	}
	err = fopen(argv[1], "w" );
	dat = fopen(argv[2], "w" );
	printf("Введите M\n");
	scanf("%lf",&M);
	h = 1/M;
	prevLayer = (double *) calloc(M+1, sizeof(double));
	currLayer = (double *) calloc(M+1, sizeof(double));
	endLayer = (double *) calloc(M+1, sizeof(double));

	printf("Введите N\n");
	scanf("%lf",&N);
	Thau = 1/N;
	filter = N/100;
	printf("Введите a\n");
	scanf("%lf",&a);
	printf("f(t,x) = 2*t*(x*x-x)) - 2*%lf*%lf(t*t+1)\n", a, a);
	if ( 2*M*M*a*a < N ){
		mode = 1;
		printf("Автоматически выбрана явная схема\n");
	}
	else{
		printf("Автоматически выбрана неявная схема\n");
		mode = 2;
	} 
	for ( int m = 0; m <= M; m++ ){
		currLayer[m] = U0(m/M);
		//printf("%lf ", currLayer[m]);
		endLayer[m] = U(m/M,1);
	}
	//printf("\n");
	if ( mode == 1 ){//явная схема
		for ( int n = 0; n < N; n++ ){
			for ( int m = 0; m <= M; m++ ){
				prevLayer[m] = currLayer[m];
			}
			currLayer[0] = bordU(0,n);
			currLayer[(int)M] = bordU(1,n);
			for ( int m = 1; m < M; m++ ){
				currLayer[m] = prevLayer[m]+ Thau*(a*a*(prevLayer[m+1] - 2*prevLayer[m] + prevLayer[m-1])/(h*h) + f(m*h,n*Thau,a));		
			}
			if ( (n % filter) == 0 ){
				for ( int m = 0; m < M + 1; m++ ){
					fprintf(dat, "%lf %lf %lf\n", m/M, currLayer[m], endLayer[m]);
				}
				fprintf(dat, "\n\n");
			}
		}
	}
	if ( mode == 2 ){//неявная схема
		double * alpha, * betta;
		alpha = (double *) calloc(M+1, sizeof(double));
		betta = (double *) calloc(M+1, sizeof(double));

		for ( int n = 0; n < N; n++ ){
			int m;
			for ( int m = 0; m <= M; m++ ) {
				prevLayer[m] = currLayer[m];
			}
			alpha[0] = 0;
			betta[0] = 0;
			int xi1 = bordU(0,n)/prevLayer[1];
			int mu1 = bordU(0,n) - prevLayer[1]*xi1;
			alpha[1] = xi1;
			betta[1] = mu1;
			for ( m = 1; m < M; m++ ){
				alpha[m+1] = Thau*a*a/(h*h) / ( 1 + 2*Thau*a*a/(h*h) - alpha[m]*Thau*a*a/(h*h) );
				betta[m+1] = ( Thau*f(m*h,n*Thau,a) + prevLayer[m] + betta[m]*(Thau*a*a)/(h*h) ) / ( 1 + 2*Thau*a*a/(h*h) - alpha[m]*Thau*a*a/(h*h) );
			}
			int xi2 = bordU(1,n)/prevLayer[m-1];
			int mu2 = bordU(1,n) - prevLayer[m-1]*xi2;
			currLayer[(int)M] = (xi2*betta[m]+mu2)/(1-xi2*alpha[m]);
			for ( int m = M; m > 0/*1*/; m-- ){
				currLayer[m-1] = alpha[m]*currLayer[m]+betta[m];
			}
			//currLayer[0] = bordU(0,n);
			if ( (n % filter) == 0 ){
				for ( int m = 0; m < M + 1; m++ ) {
					fprintf(dat, "%lf %lf %lf\n", m/M, currLayer[m], endLayer[m]);
				}
				fprintf(dat, "\n\n");
			}
		}
	}

	double max = 0, maxDiff = 0, sum = 0, sumDiff = 0;
	for ( int m = 1; m < M; m++ ){
		double tmp;
		if ( (tmp = fabs(currLayer[m] - U(m*h,1)) ) > maxDiff )
			maxDiff = tmp;
		sumDiff += tmp*tmp;
		if ( ( tmp = fabs(U(m*h,1)) ) > max )
			max = tmp;
		sum += tmp*tmp;
	}
	fprintf(err, "%lf ", maxDiff);
	fprintf(err, "%lf ", maxDiff/max);
	fprintf(err, "%lf ", sqrt(sumDiff));
	fprintf(err, "%lf\n", sqrt(sumDiff)/sqrt(h*sum));
	/*
	for ( int m = 0; m < M+1; m++ )
		printf("%lf ", currLayer[m]);
	printf("\n");
	*/
	free(prevLayer);
	free(currLayer);
	fclose(err);
	fclose(dat);
	return 0;
}