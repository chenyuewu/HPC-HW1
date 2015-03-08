EXECUTABLE = int_ring jacobi-mpi
COMPILER = mpicc
FLAG = -lrt -O3 -Wall

all: $(EXECUTABLE)

int_ring: int_ring.c
	$(COMPILER) $< -o $@ $(FLAG)

jacobi-mpi: jacobi-mpi.c
	$(COMPILER) $< -o $@ $(FLAG)

clean:
	rm -rf $(EXECUTABLE)
