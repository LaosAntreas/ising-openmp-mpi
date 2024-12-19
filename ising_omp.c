#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NO_SWEEPS 1000
#define BURN_IN 100

int *lattice;


int *init_lattice(int N, unsigned int seed){
    int *lattice = (int *)malloc(N*N*sizeof(int));
    #pragma omp parallel private(seed)
    {
    seed += omp_get_thread_num();
    #pragma omp for
    for(int i = 0; i < N * N; i++){
        lattice[i] = 2*(rand_r(&seed) & 1) - 1;
    }
    }
    return lattice;
}

double magn(int N){
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < N * N; i++){
            sum += lattice[i];
    }
    sum /= (N*N);
    return sum;
}

int energy(int N){
    int sum = 0;
    int spin_sum = 0;
    int right, left, up, down;
    #pragma omp parallel for private(right, left, up, down, spin_sum) reduction(+:sum)
    for(int i = 0; i < N*N; i++){
        right = (i/N)*N + (i+1)%N;
        left = (i/N)*N + (i-1+N)%N;
        up = ((i/N) - 1 + N) % N * N + i % N;
        down = ((i/N) + 1) % N * N + i % N;

        spin_sum = lattice[right] + lattice[left] + lattice[up] + lattice[down];
        sum += -1 * lattice[i] * spin_sum;
    }
}

int energy_diff(int i, int N){
    int right = (i/N)*N + (i+1)%N;
    int left = (i/N)*N + (i-1+N)%N;
    int up = ((i/N) - 1 + N) % N * N + i % N;
    int down = ((i/N) + 1) % N * N + i % N;

    int spin_sum = lattice[right] + lattice[left] + lattice[up] + lattice[down];

    int diff = -1 * ((-1 * lattice[i]) * spin_sum - lattice[i] * spin_sum ) ;
    return diff;
}

int sweep(int i_0, int N, double beta, unsigned int seed){
    int diff;
    double p;
    double rnd;
    int right, left, up, down;
    int spin_sum;

    #pragma omp parallel private(diff, p, rnd, left, right, up, down, spin_sum, seed) 
    {
    seed += omp_get_thread_num();
    #pragma omp for
    for(int i = i_0; i < N*N; i+=2){
        
        diff = energy_diff(i, N);
        if (diff <= 0){
            lattice[i] = - lattice[i];
        } else {
            p = exp( - beta * diff);
            rnd = rand_r(&seed) / ((double) (RAND_MAX)); 
            if ( rnd < p){
                lattice[i] = - lattice[i];
            }
        }
    }
    }
}

int main(int argc, char *argv[]){
    if(argc != 3){
        printf("Usage: %s <N> <beta>\n", argv[0]);
        exit(1);
    }
    int N = atoi(argv[1]);
    double beta = atof(argv[2]);

    unsigned int seed = 42;
    int i, j, r;

    lattice = init_lattice(N, seed);

    int energy_arr[NO_SWEEPS];
    double magn_arr[NO_SWEEPS];

    for(int i = 0; i < BURN_IN; i++){
        sweep(0, N, beta, seed);
        sweep(1, N, beta, seed);

    }

    for(int i = 0; i < NO_SWEEPS; i++){
        sweep(0, N, beta, seed);
        sweep(1, N, beta, seed);

        int e = energy(N);
        double m = magn(N);

        printf("Total Energy: %d | Total Magnetization: %.2f\n", e, m);

        energy_arr[i] = e;
        magn_arr[i] = m;
    }

    double sum_m = 0.0;
    double sum_e = 0.0;

    for(int i = 0; i < NO_SWEEPS; i++){
        sum_m += magn_arr[i];
        sum_e += energy_arr[i];
    }

    sum_m /= NO_SWEEPS;
    sum_e /= NO_SWEEPS;

    double tmp_m = 0.0;
    double tmp_e = 0.0;

    for(int i = 0; i < NO_SWEEPS; i++){
        tmp_m += (sum_m - magn_arr[i]) * (sum_m - magn_arr[i]);
        tmp_e += (sum_e - energy_arr[i]) * (sum_e - energy_arr[i]);
    }

    tmp_m /= NO_SWEEPS;
    tmp_e /= NO_SWEEPS;

    tmp_m = sqrt(tmp_m);
    tmp_e = sqrt(tmp_e);

    printf("Average Energy %f +- %f, Average Magn %f +- %f\n", sum_e, tmp_e, sum_m, tmp_m);

    return EXIT_SUCCESS;
}
