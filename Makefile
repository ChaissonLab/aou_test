all: ctyper2
objs=CramReader.o FastaReader.o FastqReader.o Genotyper.o KmerCounter.o KmerHash.o KmerMatrix.o KmerWindow.o KtableReader.o main.o PriorData.o Processor.o Regression.o TreeRound.o
%.o: %.cpp
	g++ -std=c++17 -O3 -flto -march=native -funroll-loops  -fno-rtti \
 -Wall -Wextra -Wshadow -fdiagnostics-color=always -I ${CONDA_PREFIX}/include -I ${CONDA_PREFIX}/include/eigen3 -L${CONDA_PREFIX}/lib $^ -lz -pthread -lhts -c -o $@

ctyper2: ${objs}
	g++ -std=c++17 -O3 -flto -march=native -funroll-loops  -fno-rtti \
 -Wall -Wextra -Wshadow -fdiagnostics-color=always -I ${CONDA_PREFIX}/include -I ${CONDA_PREFIX}/include/eigen3 -L${CONDA_PREFIX}/lib $^ -lz -pthread -lhts -o ctyper2
