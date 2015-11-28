all:nbody

nbody: nbody.cpp
	g++ -lGL -lglut -lGLU -lpthread -std=c++11 nbody.cpp -o nbody
