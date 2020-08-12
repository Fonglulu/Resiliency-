EXECS=RandomWalker_MPI
MPICXX?=mpicxx

all: ${EXECS}

RandomWalker_MPI: RandomWalker_MPI.cpp
	${MPICXX} -o RandomWalker_MPI RandomWalker_MPI.cpp

clean:
	rm -f ${EXECS}
