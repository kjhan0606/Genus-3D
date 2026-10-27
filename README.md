# Genus-3D

MPI C programs that measure the genus topology of three-dimensional density
fields, with smoothing through FFTW 2.x MPI transforms.

## Programs

Built by `make`. All take command-line arguments; nothing is read from stdin.

| Program | Command line |
| --- | --- |
| genus3 | `./genus3 indenfile boxsize outnum p/i nsub` |
| genus2 | `./genus2 indenfile boxsize outnum p/i nsub` |
| genus | `./genus indenfile boxsize outnum p/i` |
| genus_plate | `./genus_plate indenfile boxsize outnum p/i nsub` |
| genus_rod | `./genus_rod indenfile boxsize outnum p/i nsub` |
| genus_rod2 | `./genus_rod2 indenfile boxsize outnum p/i nsub` |
| genus_rod3 | `./genus_rod3 indenfile boxsize outnum p/i nsub` |
| genus_subcubic | `./genus_subcubic indenfile boxsize outnum p/i nsub` |
| sigma12 | `./sigma12 indenfile boxsize outnum p/i nsub` |
| statsubg | `./statsubg ithre` |
| statsubg2 | `./statsubg2 ithre` |
| checkplatevol | `./checkplatevol n nsub` |

`genus3` finds the density threshold exactly with a parallel sort
(`mpiqsort.c`); `genus2` estimates it by sampling. `statsubg` and
`statsubg2` combine `subgenus*` outputs into `statsubgenus*` summaries.
See the sources for argument details.

## Requirements

- An MPI C compiler (`mpicc`, or `mpiicx` for Intel MPI)
- FFTW 2.x with MPI (`srfftw_mpi.h`, `libsrfftw_mpi`, `libsfftw_mpi`,
  `libsrfftw`, `libsfftw`)

## Build

```sh
make
make FFTW2=/path/to/fftw2
make CC=mpiicx FFTW2=/path/to/fftw2
make clean
```

## Run

```sh
mpirun -np 4 ./genus3 indenfile boxsize outnum p/i nsub
```

`indenfile` is a binary file starting with three `int` grid dimensions
(nx, ny, nz) followed by the float density array. The subcube programs
also read a text `contour.inp` in the working directory.
