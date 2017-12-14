all:
	g++ -std=c++11 heat.cpp -o one
	./one 50 10000 1 d 1
	gnuplot instruction11
	./one 500 1000 1 d 2
	gnuplot instruction12
	./one 50 10000 1 v 1
	gnuplot instruction21
	./one 500 1000 1 v 2
	gnuplot instruction22
test:
	g++ -std=c++11 heat.cpp -o one
	./one 50 10000 1 v 2
	./one 100 40000 1 v 2