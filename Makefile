CC = mpicc
FFTW2 ?= /usr
CFLAGS = -O2 -I$(FFTW2)/include
LDFLAGS = -L$(FFTW2)/lib
LDLIBS = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lm

MPI_PROGRAMS = genus genus2 genus3 genus_plate genus_rod genus_rod2 genus_rod3 genus_subcubic sigma12
SERIAL_PROGRAMS = statsubg statsubg2 checkplatevol

all: $(MPI_PROGRAMS) $(SERIAL_PROGRAMS)

genus3: genus3.c mpiqsort.c mpiqsort.h
	$(CC) $(CFLAGS) -o $@ genus3.c mpiqsort.c $(LDFLAGS) $(LDLIBS)

$(filter-out genus3,$(MPI_PROGRAMS)): %: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

$(SERIAL_PROGRAMS): %: %.c
	$(CC) -O2 -o $@ $< -lm

clean:
	rm -f $(MPI_PROGRAMS) $(SERIAL_PROGRAMS) *.o

.PHONY: all clean
