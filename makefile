
FLAG = -O3 -m64 -Wall -pipe -mmmx -msse -msse2 -msse3 -mcmodel=large --std=c++11 -W -Wno-sign-compare

all: EBR

EBR.o: EBR.cpp ./christofides-algorithm/MST.h ./christofides-algorithm/Matching/Matching.h ./christofides-algorithm/Matching/Graph.h ./christofides-algorithm/Christofides.h
	g++ $(FLAG) -c EBR.cpp -o EBR.o

EBR: ./christofides-algorithm/Matching.o ./christofides-algorithm/BinaryHeap.o ./christofides-algorithm/Graph.o EBR.o
	g++ $(FLAG) ./christofides-algorithm/Matching.o ./christofides-algorithm/BinaryHeap.o ./christofides-algorithm/Graph.o EBR.o -o EBR