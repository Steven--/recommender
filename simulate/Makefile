all: simulate simulate-bounded-100
	

simulate: simulate.cpp
	g++ -O3 -march=native -mtune=native -Wall simulate.cpp -lm -lgsl -lgslcblas -lpthread -osimulate

simulate-bounded-100: simulate.cpp
	g++ -DBOUNDED_HISTORY -DHISTORY_BOUND=100 -O3 -march=native -mtune=native -Wall simulate.cpp -lm -lgsl -lgslcblas -lpthread -o$@

