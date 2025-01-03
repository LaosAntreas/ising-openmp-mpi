# ising-openmp-mpi


Parallelising the 2D Ising model with OpenMP and OpenMPI


---
## The Ising Model

The energy of the System:

$$H(s) = -\sum_{<ij>}{s_i s_j}$$

Where $s_i$ are the spins that can have a value of either $1$ or $-1$ and $<ij>$ is the sum over the nearest neighbours of $s$.

The goal of the model is to minimise the energy using the Metropolis-Hastings algorithm.

The magnetisation of the system: 

$$M = {1\over N}\sum_{i=0}^N{s_i}$$

More information [here](https://en.wikipedia.org/wiki/Ising_model)

---

## ising_omp.c

#### Compiling

```
gcc ising_omp.c -o omp.out -lm -fopenmp
```

#### Running

```
./omp.out <N> <beta>
```

  

In this implementation, the 2D lattice is represented in a single 1D array so that the problem can be distributed between the threads efficiently.

#### Scaling 

<img width="357" alt="Screenshot 2024-12-19 at 12 19 55 PM" src="https://github.com/user-attachments/assets/2f8751c4-651f-4d81-b0b9-b1822c4df370" />


---

## ising_mpi.c

#### Compiling

```
mpicc ising_mpi.c -o mpi.out -lm
```

#### Running

```
mpirun -np <processes> ./mpi.out <N> <beta>
```

In this implementation, the 2D lattice is split into sublattices(split by rows) that share halo regions with neighbouring processes. At each step, each process updates half of the lattice and shares the halo region with its neighbour, then updates the other half and shares the other halo region with its other neighbour.

#### Scaling

<img width="355" alt="Screenshot 2024-12-19 at 12 20 46 PM" src="https://github.com/user-attachments/assets/38bb38ed-3478-4b7f-84e0-a6889cabe2e6" />
