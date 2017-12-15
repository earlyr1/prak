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
	int mode = atoi(argv[4]); //1 if явная, 0 else
	M = atoi(argv[1]); // N per x 
	P = atoi(argv[2]); // N per time	
	T = 1.0; // complete time
	char * edge1 = argv[3]; //0 if values; 1 if derivative
	edge = argv[4][0] == 'd'? 1:0;
	ofstream ofs("data" + to_string(mode) + to_string(edge) + ".txt");
	dx = 1.0 / M;
	dt = (T + 0.0) / P;
	prevLayer = new double[M + 1](); 
	currLayer = new double[M + 1]();
	vector<double> ans(M + 1, 0.0);
	for(int i = 0; i <= M; i++) {
		ans[i] = U(dx * i, T);
	}
	if (mode) {
		
		double * tmp;
		for(int i = 0; i <= M; i++) prevLayer[i] = phi(dx * i);	
		if (!edge) 	{
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
		} else {
			for(int t = 1; t <= P; t++) {
				for(int x = 1; x < M; x++) {
					currLayer[x] = (f(dx * x, dt * (t - 1)) + a * a * ((prevLayer[x - 1] - 2 * prevLayer[ x] + prevLayer[ x + 1]) / (dx * dx))) * dt + prevLayer[x];
				}
				currLayer[0] = (4 * currLayer[1] - currLayer[2]) / 3;
				currLayer[M] = (4 * currLayer[M - 1] - currLayer[M - 2]) / 3;

				tmp = prevLayer;
				prevLayer = currLayer;
				currLayer = tmp;
			}
		}
		double error = 0.;
		for(int i = 0; i < M + 1; i++) {
			ofs << (i + 0.) / M << " " << prevLayer[i] << " " <<  ans[i] << endl;
			error += (ans[i] - prevLayer[i]) * (ans[i] - prevLayer[i]);
		}
		cout << mode << " " << edge << " " << error << endl;
	} else {
		for(int i = 0; i < M + 1; i++) {
			prevLayer[i] = currLayer[i];
		}          
		double * alpha, * beta, *tmp;
		alpha = new double[M + 1]();
		beta = new double[M + 1]();
		for(int i = 0; i <= M; i++) currLayer[i] = phi(dx * i);
		if (!edge) {
			for ( int t = 0; t < P; t++ ){
				tmp = prevLayer;
				prevLayer = currLayer;
				currLayer = tmp;
				int x;
				alpha[0] = 0;
				beta[0] = 0;
				int X1 = a1(t * dt) / prevLayer[1];
				int mew1 = a1(t * dt) - prevLayer[1]*X1;
				alpha[1] = X1;
				beta[1] = mew1;
				for ( x = 1; x < M; x++ ){
					alpha[x+1] = dt*a*a/(dx*dx) / ( 1 + 2*dt*a*a/(dx*dx) - alpha[x]*dt*a*a/(dx*dx) );
					beta[x+1] = ( dt*f(x*dx,t*dt) + prevLayer[x] + beta[x]*(dt*a*a)/(dx*dx) ) / ( 1 + 2*dt*a*a/(dx*dx) - alpha[x]*dt*a*a/(dx*dx) );
				}
				int X2 = a2(t * dt) / prevLayer[x-1];
				int mew2 = a2(t * dt) - prevLayer[x-1]*X2;
				currLayer[M] = (X2*beta[x]+mew2)/(1-X2*alpha[x]);
				for ( int x = M; x > 0/*1*/; x-- ){
					currLayer[x-1] = alpha[x]*currLayer[x]+beta[x];
				}		
			}
		} else {
			for ( int t = 0; t < P; t++ ){
				int X1, mew1, X2, mew2;
				tmp = prevLayer;
				prevLayer = currLayer;
				currLayer = tmp;
				int x;
				alpha[0] = 0;
				beta[0] = 0;
				if (!edge) {
					X1 = 0;
					mew1 = a1(t * dt);
				} else {
					X1 = 1;
					mew1 = 0;
				}
				alpha[1] = X1;
				beta[1] = mew1;
				for ( x = 1; x < M; x++ ){
					alpha[x+1] = dt*a*a/(dx*dx) / ( 1 + 2*dt*a*a/(dx*dx) - alpha[x]*dt*a*a/(dx*dx) );
					beta[x+1] = ( dt*f(x*dx,t*dt) + prevLayer[x] + beta[x]*(dt*a*a)/(dx*dx) ) / ( 1 + 2*dt*a*a/(dx*dx) - alpha[x]*dt*a*a/(dx*dx) );
				}
				if (!edge) {
					X2 = 0;
					mew2 = a2(t * dt);
				} else {
					X2 = 1;
					mew2 = 0;
				}
				currLayer[M] = (X2*beta[x]+mew2)/(1-X2*alpha[x]);
				for ( int x = M; x > 1/*1*/; x-- ){
					currLayer[x-1] = alpha[x]*currLayer[x]+beta[x];
				}
				currLayer[0] = (4 * currLayer[1] - currLayer[2]) / 3;
				currLayer[M] = (4 * currLayer[M - 1] - currLayer[M - 2]) / 3;
			}
		}
		double error = 0.;
		for(int i = 0; i < M + 1; i++) {
			ofs << (i + 0.) / M << " " << currLayer[i] << " " <<  ans[i] << endl;
			error += (ans[i] - currLayer[i]) * (ans[i] - currLayer[i]);
		}
		cout << mode << " " << edge << " " << error << endl;
	}
}
