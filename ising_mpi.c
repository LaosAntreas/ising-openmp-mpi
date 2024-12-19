#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>

#define NO_SWEEPS 1000
#define BURN_IN 100


int** init_sublattice(int columns, int rows, int seed){
    int **sublattice = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++){
        sublattice[i] = (int *)malloc(columns * sizeof(int));
        for (int j = 0; j < columns; j++){
            sublattice[i][j] = 2*(rand_r(&seed) & 1) - 1;
            //sublattice[i][j] = 1;
        }
    }

    return sublattice;
}

void print_sublattice(int **sublattice, int columns, int rows){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            if (sublattice[i][j] == 1){
                printf("_|");
            }else{
                printf("X|");
            }
        }
        printf("\n");

    }
}

double magn(int **sublattice, int columns, int rows){
    double sum = 0.0;
    for(int i = 1; i < rows-1; i++){
        for(int j = 0; j < columns; j++){
            sum += sublattice[i][j];
        }
    }
    sum /= (columns*(rows-2));
    return sum;
}

int subenergy(int **sublattice, int columns, int rows){
    int sum = 0;

    for(int i = 1; i < rows - 1; i++){
        for(int j = 0; j < columns; j++){
            int up = sublattice[i-1][j]; //periodic boundary conditions using the halo regions
            int down = sublattice[i+1][j]; //periodic boundary conditions using the halo regions
            int left = sublattice[i][(j-1+ columns)%columns];
            int right = sublattice[i][(j+1)%columns];

            int spin_sum = up + down + left + right;
            sum += -1 * sublattice[i][j] * spin_sum;

        }
    }
    return sum;

}

int energy_diff(int **sublattice, int i, int j, int columns, int rows){
    int up = sublattice[i-1][j]; //periodic boundary conditions using the halo regions
    int down = sublattice[i+1][j]; //periodic boundary conditions using the halo regions
    int left = sublattice[i][(j-1+ columns)%columns];
    int right = sublattice[i][(j+1)%columns];

    int spin_sum = up + down + left + right;
    int diff = -1 * ((-1 * sublattice[i][j]) * spin_sum - sublattice[i][j] * spin_sum);
    return diff;
}

void sweep(int **sublattice, int columns, int rows, double beta, int seed, int up, int down){
    
    int *send_up = (int *)malloc(columns * sizeof(int));
    int *recv_up = (int *)malloc(columns * sizeof(int));
    int *send_down = (int *)malloc(columns * sizeof(int));
    int *recv_down = (int *)malloc(columns * sizeof(int));


    send_up = sublattice[1];
    
    MPI_Request requests[4];
    MPI_Isend(send_up, columns, MPI_INT, up, 0, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(recv_down, columns, MPI_INT, down, 0, MPI_COMM_WORLD, &requests[1]);
   
    MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
    MPI_Wait(&requests[1], MPI_STATUS_IGNORE);

    sublattice[rows-1] = recv_down;

    
    for(int i = rows / 2; i < rows - 1 i++){
        for(int j = 0; j < columns; j++){
            int diff = energy_diff(sublattice, i, j, columns, rows);
            
            if(diff < 0){
                if (sublattice[i][j] == 1){
                    sublattice[i][j] = -1;
                }else if (sublattice[i][j] == -1){
                    sublattice[i][j] = 1;
                }else{
                    printf("Error at [%d,%d]\n", i, j);
                }
                
            }else{
                double p = exp(-beta * diff);
                double rnd = (double)rand_r(&seed)/RAND_MAX;
                if (rnd < p){
                    if (sublattice[i][j] == 1){
                        sublattice[i][j] = -1;
                    }else if (sublattice[i][j] == -1){
                        sublattice[i][j] = 1;
                    }else{
                        printf("Error at [%d,%d]\n", i, j);
                    }
                }
            }

        }
    }
    
    send_down = sublattice[rows - 2];

    MPI_Isend(send_down, columns, MPI_INT, down, 1, MPI_COMM_WORLD, &requests[2]);
    MPI_Irecv(recv_up, columns, MPI_INT, up, 1, MPI_COMM_WORLD, &requests[3]);

    MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
    MPI_Wait(&requests[3], MPI_STATUS_IGNORE);

    sublattice[0] = recv_up;

    //print received row

    // printf("\nReceived row\n");
    // for(int i = 0; i < columns; i++){
    //     printf("%d| ", recv_up[i]);
    // }


    for(int i = 1; i < rows / 2; i++){
        for(int j = 0; j < columns; j++){
            int diff = energy_diff(sublattice, i, j, columns, rows);
            double p = exp(-beta * diff);
            double rnd = (double)rand_r(&seed)/RAND_MAX;
            if (rnd < p){
                sublattice[i][j] *= -1;
            }
        }
    }

    

}

int main(int argc, char * argv[]){

    MPI_Init(&argc, &argv);
    

    int n = 0 ;
    MPI_Comm_size(MPI_COMM_WORLD, &n);

    int id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if(argc != 3){
        if(id == 0){
            printf("Usage: %s <N> <beta>\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    
    double beta = atof(argv[2]);
    int N = atoi(argv[1]);

    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if ((N % n) % 2 != 0){
        if(id == 0){
            printf("Incorrect number of processes\n");

        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }


    unsigned int seed = 53 + id;

    int rows = N/n + 2; //sublattice plus halo region
    int columns = N;

    int up = (id + 1)%n; //neighbouring process above
    int down = (id - 1 + n)%n; //neighbouring process below


    int **sublattice = init_sublattice(columns, rows, seed);
    //print_sublattice(sublattice, columns, rows);
    
    //burn in
    for(int i = 0; i < BURN_IN; i++){
        sweep(sublattice, columns, rows, 1.0, seed, up, down);
    }

    int energy_arr[NO_SWEEPS];
    double magn_arr[NO_SWEEPS];

    for(int i = 0; i < NO_SWEEPS; i++){
        sweep(sublattice, columns, rows, beta, seed, up, down);


        double m = magn(sublattice, columns, rows);
        int e = subenergy(sublattice, columns, rows);
        //printf("\nMagn %.2f id %d\n", m, id);

        double total_magn = 0.0;
        int total_energy = 0;
        

        MPI_Reduce(&m, &total_magn, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&e, &total_energy, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if(id == 0){
            //printf("Magnetization: %.2f\n", total_magn/n);
            //printf("Energy: %d\n", total_energy/n);
            energy_arr[i] = total_energy/n;
            magn_arr[i] = total_magn/n;

        }
        
    }

    if(id == 0){
        double avg_m = 0.0;
        double avg_e = 0.0;

        for(int i = 0; i < NO_SWEEPS; i++){
            avg_m += magn_arr[i];
            avg_e += energy_arr[i];
        }

        avg_m /= NO_SWEEPS;
        avg_e /= NO_SWEEPS;

        double std_m = 0.0;
        double std_e = 0.0;
        for(int i = 0; i < NO_SWEEPS; i++){
            std_m += (avg_m - magn_arr[i]) * (avg_m - magn_arr[i]);
            std_e += (avg_e - energy_arr[i]) * (avg_e - energy_arr[i]);
        }

        std_m /= NO_SWEEPS;
        std_e /= NO_SWEEPS;

        printf("Average Energy %f +- %f, Average Magn %f +- %f\n", avg_e, std_e, avg_m, std_m);
    }






    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}