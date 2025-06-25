all: ctyper2

ctyper2: CramReader.cpp FastaReader.cpp FastqReader.cpp Genotyper.cpp KmerCounter.cpp KmerHash.cpp KmerMatrix.cpp KmerWindow.cpp KtableReader.cpp main.cpp PriorData.cpp Processor.cpp Regression.cpp TreeRound.cpp
	g++ -std=c++17 -O3 -flto -march=native -funroll-loops  -fno-rtti \
 -Wall -Wextra -Wshadow -fdiagnostics-color=always -I ${CONDA_PREFIX}/include -I ${CONDA_PREFIX}/include/eigen3 -L${CONDA_PREFIX}/lib $^ -lz -pthread -lhts -o ctyper2
