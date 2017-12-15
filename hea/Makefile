all:
	g++ -std=c++11 heat.cpp -o one
	./one 500 1000 v 0
	gnuplot instruction11
	./one 50 10000 v 1
	gnuplot instruction12
	./one 500 1000 d 0
	gnuplot instruction21
	./one 50 10000 d 1
	gnuplot instruction22
test:
	g++ -std=c++11 heat.cpp -o one
	./one 50 10000 v 2
	./one 100 40000 v 2 