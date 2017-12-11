#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int M, P, T;
double dx, dt;
double pi = 3.1415926;
double a = 1.0;
int edge;
double * prevLayer;
double * currLayer;

double f(double x, double t) {
	return 2 * t * cos(pi * x) + a * a * pi * pi * (t * t + 1) * cos(pi * x);
}
double U(double x, double t) {
	return (t * t + 1) * cos(pi * x);
}

double dUdx(double x, double t) {
	return -(t * t + 1) * pi * sin(pi * x);
}

double phi(double x) {
	return U(x, 0);
}

double a1(double t) {
	return U(0, t);
}

double a2(double t) {
	return U(1, t);
}



int main(int argc, char ** argv) {
	ofstream ofs("data.txt");
	int mode = 2;
	M = atoi(argv[1]); // N per x 
	P = atoi(argv[2]); // N per time	
	T = atof(argv[3]); // complete time
	char * edge1 = argv[4]; //0 if values; 1 if derivative
	edge = argv[4][0] == 'd'? 1:0;
	dx = 1.0 / M;
	dt = (T + 0.0) / P;
	prevLayer = new double[M + 1](); 
	currLayer = new double[M + 1]();
	vector<double> ans(M + 1, 0.0);
	for(int i = 0; i <= M; i++) {
		ans[i] = U(dx * i, T);
	}
	if (mode == 1) {
		
		double * tmp;
		for(int i = 0; i <= M; i++) prevLayer[i] = phi(dx * i);		
		for(int t = 1; t <= P; t++) {
			for(int x = 0; x <= M; x++) {
				if (x == 0) currLayer[x] = a1(t * dt);
				else if (x == M) currLayer[x] = a2(t * dt);
				else currLayer[x] = (f(dx * x, dt * (t - 1)) + a * a * ((prevLayer[x - 1] - 2 * prevLayer[ x] + prevLayer[ x + 1]) / (dx * dx))) * dt + prevLayer[x];
			}
			tmp = prevLayer;
			prevLayer = currLayer;
			currLayer = tmp;
		}
		for(int i = 0; i < M + 1; i++) {
			ofs << (i + 0.) / M << " " << prevLayer[i] << " " <<  ans[i] << endl;
		}
	} else {
		for(int i = 0; i < M + 1; i++) {
			prevLayer[i] = currLayer[i];
		}          
		double * alpha, * beta, *tmp;
		alpha = new double[M + 1]();
		beta = new double[M + 1]();
		for(int i = 0; i <= M; i++) currLayer[i] = phi(dx * i);
		for ( int t = 0; t < P; t++ ){
			tmp = prevLayer;
			prevLayer = currLayer;
			currLayer = tmp;
			int x;
			alpha[0] = 0;
			beta[0] = 0;
			int xi1 = a1(t * dt) / prevLayer[1];
			int mu1 = a1(t * dt) - prevLayer[1]*xi1;
			alpha[1] = xi1;
			beta[1] = mu1;
			for ( x = 1; x < M; x++ ){
				alpha[x+1] = dt*a*a/(dx*dx) / ( 1 + 2*dt*a*a/(dx*dx) - alpha[x]*dt*a*a/(dx*dx) );
				beta[x+1] = ( dt*f(x*dx,t*dt) + prevLayer[x] + beta[x]*(dt*a*a)/(dx*dx) ) / ( 1 + 2*dt*a*a/(dx*dx) - alpha[x]*dt*a*a/(dx*dx) );
			}
			int xi2 = a2(t * dt) / prevLayer[x-1];
			int mu2 = a2(t * dt) - prevLayer[x-1]*xi2;
			currLayer[M] = (xi2*beta[x]+mu2)/(1-xi2*alpha[x]);
			for ( int x = M; x > 0/*1*/; x-- ){
				currLayer[x-1] = alpha[x]*currLayer[x]+beta[x];
			}
			
		}
		for(int i = 0; i < M + 1; i++) {
			ofs << (i + 0.) / M << " " << currLayer[i] << " " <<  ans[i] << endl;
		}
	}
}
