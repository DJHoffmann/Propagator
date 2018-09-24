# GCC Makefile
# =================

TARGET = Propagator
DIR = $(HOME)/Code/Propagator/
INCLUDE = $(DIR)

it: $(DIR)$(TARGET).cpp
	cp $(DIR)$(TARGET).cpp ./ ; \
	g++ -O2 $(TARGET).cpp -o X$(TARGET).x -I $(INCLUDE) -lfftw3 -lm

mpi: $(DIR)$(TARGET).cpp
	cp $(DIR)$(TARGET).cpp ./ ; \
	icc $(MPI_LIBS) -O2 -axT -xT -i-static $(TARGET).cpp -o $(TARGET).x -I$(INCLUDE) -I/apps/fftw/3.2.2-double/include/ -L/usr/lib64/ -L/apps/fftw/3.2.2-double/lib/ -I$(MPI_INCLUDE) -lfftw3 -lm -DuseMPI 

mpiRed: $(DIR)$(TARGET).cpp
	cp $(DIR)$(TARGET).cpp ./ ; \
	icc $(MPI_LIBS) -O2 -axT -xT -i-static $(TARGET).cpp -o $(TARGET).x -I$(INCLUDE) -I/apps/fftw/3.2.2-double/include/ -L/usr/lib64/ -L/apps/fftw/3.2.2-double/lib/ -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/  -L/usr/lib/gcc/x86_64-redhat-linux6E/4.4.4/ -I$(MPI_INCLUDE) -lfftw3 -lm -DuseMPI 

go: X$(TARGET).x
	./X$(TARGET).x

crun: $(DIR)$(TARGET).cpp
	make it ; \
	make go
