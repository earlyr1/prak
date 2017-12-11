#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <ctime>

int n;
int cnt;
int cnt1;
char t;

double f(int i, int j) {
	if (i == j && i == n - 1) return n;
	if (i == j) return 1;
	if (i == n - 1) return j + 1;
	if (j == n - 1) return i + 1;
	return 0;
}

void way(char &t) {
	std::cout << "How to choose maximum? Print u for maximal absolute value, m for maximal row and c for cycle" << std::endl;
	std::cin >> t;
}

double iteration(std::vector< std::vector<double> > &A) {
	
	int maxi, maxj;
	double maxk, m;
	//usual maximum
	

	if (t == 'u') {
		maxi = 1;
		maxj = 0;
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				if (fabs(A[i][j]) > fabs(A[maxi][maxj]) && i > j) {
					maxi = i;
					maxj = j;
				}
			}
		}
		maxk = fabs(A[maxi][maxj]);
	}
	
	//max-row
	
	if (t == 'm') {
		maxi = 0;
		double smax = 0;
		for(int i = 0; i < n; i++) {
			double s = 0;
			for(int j = 0; j < n; j++) {
				if (i == j) continue;
				s += A[i][j] * A[i][j];
			}
			if (s > smax) {
				smax = s;
				maxi = i;
			}
		}
		maxj = maxi == 0? 1: 0;
		for(int j = 0; j < n; j++) {
			if (fabs(A[maxi][j]) > fabs(A[maxi][maxj]) && j != maxi) {
				maxj = j;
			}
		}
		if (maxi < maxj) {
			int tmp = maxi;
			maxi = maxj;
			maxj = tmp;
		}
		maxk = fabs(A[maxi][maxj]);
	}
	
	//cycle
	if (t == 'c') {

		m = 0;
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				if (fabs(A[i][j]) > m && i > j) {
					m = fabs(A[i][j]);
				}
			}
		}
		int count = cnt1 % (n * n);
		maxi = count / n;
		maxj = count % n;
		while (maxi <= maxj) {
			cnt1++;
			count = cnt1 % (n * n);
			maxi = count / n;
			maxj = count % n;
		}
		maxk = m; 
	}
	//double maxk = fabs(A[maxi][maxj]);
	std::cout << maxi << " " << maxj << " " << maxk << " " << std::endl;
	double i = maxi;
	double j = maxj;
	double phi;
	if (fabs(A[i][j]) < 1e-8) phi = 0;
	else if (A[i][i] - A[j][j] < 1e-8) phi = atan(2 * A[i][j] / (A[i][i] - A[j][j])) / 2;
	else phi = 3.1415926 / 4;
	double c = cos(phi);
	double s = sin(phi);
	
	std::vector<double> v1(n, 0.0);
    std::vector<double> v2(n, 0.0);
    for (int m = 0; m < n; m++)
        v1[m] = A[i][m] * c + A[j][m] * s;
    for (int m = 0; m < n; m++)
        v2[m] = A[i][m] * (-s) + A[j][m] * c;
    for (int m = 0; m < n; m++) {
        A[i][m] = v1[m];
        A[j][m] = v2[m];
    }

    v1 = std::vector<double> (n, 0.0);
    v2 = std::vector<double> (n, 0.0);
    for (int m = 0; m < n; m++)
        v1[m] = A[m][i] * c + A[m][j] * s;
    for (int m = 0; m < n; m++)
        v2[m] = A[m][i] * (-s) + A[m][j] * c;
    for (int m = 0; m < n; m++) {
        A[m][i] = v1[m];
        A[m][j] = v2[m];
    }

    cnt++; cnt1++;
	return maxk;
}

int main() {
	std::vector< std::vector<double> > A; 
	std::cout << "From file? Print y/n\n";
	char c; std::cin >> c;
	if (c == 'y') {
		std::cout << "Print filename\n";
		std::string fname;
		std::cin >> fname;
		std::ifstream ifs;
		ifs.open(fname);
		way(t);
		ifs >> n;
		A.resize(n);
		for(auto &i: A) {
			i.resize(n);
			for(auto &j: i) {
				ifs >> j;
			}
		}
	} else {
		std::cout << "Enter N\n";
		std::cin >> n;
		way(t);
		A.resize(n);
		for(int i = 0; i < n; i++) {
			A[i].resize(n);
			for(int j = 0; j < n; j++) {
				A[i][j] = f(i, j);
			}
		}
	}
	double maxi = 1e10;
	double tm = clock();
	while (maxi > 1e-8) {
		maxi = iteration(A);
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				//std::cout << A[i][j] << " ";
			}
			//usleep(100000);
			//std::cout << std::endl;

		}
		//std::cout << std::endl;
	}
	tm = (clock() - tm) / CLOCKS_PER_SEC;
	std::cout << std::endl;
	std::cout << "Number of iterations: " << cnt << std::endl;
	std::cout << "Time: " << tm << std::endl;
	for(int i = 0; i < n; i++) std::cout << A[i][i] << std::endl;
	std::cout << std::endl;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			//std::cout << A[i][j] << " ";
		}
		//std::cout << std::endl;
	}
}













