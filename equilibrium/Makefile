all: equilibrium equilibrium-hub
	

equilibrium: equilibrium.cpp
	g++ -DDIRECTED equilibrium.cpp rtnorm.cpp --std=c++11 -I . -lgsl -lgslcblas -O3 -o$@
	
equilibrium-hub: equilibrium.cpp
	g++ -DDIRECTED -DADD_HUB equilibrium.cpp rtnorm.cpp --std=c++11 -I . -lgsl -lgslcblas -O3 -o$@
