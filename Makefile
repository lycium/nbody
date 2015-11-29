all:nbody

nbody: nbody.cpp
	g++ -O3 -lGL -lglut -lGLU -lpthread -std=c++11 nbody.cpp -o nbody
