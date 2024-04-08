/*
This file contains methods, logic, and code in order to calculate the primality of all numbers from 1 to some max using the Fermat Primality Test,
the Miller-Rabin Primality Test, and the Solovay-Strassen Primality test. MPI is used to create three processes, with each process performing one
of the previously mentioned primality tests. OpenMP is then used on each process to partition the range of numbers based on the number of threads 
given in the command line and calculate the primality of its given set of the range of numbers. The total number of primes calculated by each thread
is then collected and combined in order to get the total number of primes each process calculated, then each total is printed along with how many
each process/test should have calculated.

NOTE: This program requires exactly 3 MPI processes to be generated and will not run if asked to create more or less then 3.

Compile Command : mpicc -fopenmp PT.c -o PT
Run Command     : mpiexec -n 3 ./PT <Number of Threads>

Class : AUCSC 450
Name  : Zachary Kelly
ID    : 1236421
Date  : April 18th, 2023
*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

/* Shared Globals */

//number of iterations each primality test will perform on each number being checked for primality.
//Increase for more certainty of primality.
int numberOfIterations = 50; 

//Checks for the number of probable primes each algorithm should find when testing to some power of ten.
int numberOfPrimesBelow100000     = 9592;
int numberOfPrimesBelow1000000    = 78498;
int numberOfPrimesBelow10000000   = 664579;
int numberOfPrimesBelow100000000  = 5761455;
int numberOfPrimesBelow1000000000 = 50847534;

int threadCount;

//The range of numbers to be tested for primality.
//The program will test from 1 to max each number for primality.
//Only one 'max' should be uncommented at run time.
int max = 100000;
//int max = 1000000;
//int max = 10000000;
//int max = 100000000;
//int max = 1000000000;

//First 51 Mersenne Primes for further testing of the algorithms/
long long mersennePrimes[51] = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 
                                9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 
                                1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 
                                32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933};
//###############################################################################
/* Fermat Primality Test Globals */
//###############################################################################

int fermatInitializer = 0;
double fermatTimeInitializer = 0.0;

int* fermatTotal = &fermatInitializer;
double* fermatTime = &fermatTimeInitializer;

//###############################################################################
/* Miller-Rabin Primality Test Globals */
//###############################################################################

int millerInitializer = 0;
double millerTimeInitializer = 0.0;

int *millerTotal = &millerInitializer;
double* millerTime = &millerTimeInitializer;

//###############################################################################
/* Solovay-Strassen Primality Test Globals*/
//###############################################################################

int solovayInitializer = 0;
double solovayTimeInitializer = 0.0;

int *solovayTotal = &solovayInitializer;
double* solovayTime = &solovayTimeInitializer;

//###############################################################################
/* Shared Methods*/
//###############################################################################

/*
modExp:
Performs modular exponentiation (c = (base ^ exponent) mod modular).
https://www.sanfoundry.com/c-program-perform-fermat-primality-test/

base     : The base number that will be raised to some power (the exponent)
exponent : The number to which we will be raising the base
modular  : The number to which we will mod the base^exponent result
*/

long long modExp(long long base, long long exponent, long long modular){

    long long x = 1;

    while(exponent > 0){

        if((exponent % 2) == 1){

            x = (x * base) % modular;

        }

        base = (base * base) % modular;

        exponent = exponent / 2;

    }

    return x % modular;

}

/*
multiplyMod:
Calculates (a * b) % c taking into account that a * b might overflow.
https://www.sanfoundry.com/c-program-perform-fermat-primality-test/

a       : First number to be multiplied
b       : Second number to be multiplied
modular : Number that will be used to mod the result of (a * b)
*/

long long multiplyMod(long long a, long long b, long long modular){

    long long x = 0;
    long long y = a % modular;

    while(b > 0){

        if((b % 2) == 1){

            x = (x + y) % modular;

        }

        y = (y * 2) % modular;

        b = b / 2;

    }

    return x % modular;

}

//###############################################################################
/* Fermat Primality Test Methods */
//###############################################################################

/*
fermatPrimalityTest:
Performs the Fermat Primality Test
https://www.sanfoundry.com/c-program-perform-fermat-primality-test/

possiblePrime : The number to be tested for primality
iterations    : The number of iterations that possiblePrime should be tested for primality
*/

bool fermatPrimalityTest(long long possiblePrime, int iterations){

    if(possiblePrime == 1){ //0 and 1 are not a prime numbers

        return false;

    }

    if((possiblePrime % 2) == 0 && possiblePrime != 2){ //If not 2, and possiblePrime % 2 = 0 is must be an even number

        return false;

    }

    if(possiblePrime == 3){ // If less then three must be prime (1, 2, and 3 all prime)

        return true;

    }

    for(int i = 0; i < iterations; i++){ //Perform the Fermat Primality test

        long long a = rand() % (possiblePrime - 1) + 1;

        if(modExp(a, (possiblePrime - 1), possiblePrime) != 1){

            return false;

        }
    }

    return true;

}

//###############################################################################
/* Miller-Rabin Primality Test Methods */
//###############################################################################

/*
millerRabinPrimalityTest:
Performs the Miller-Rabin Primality Test
https://www.sanfoundry.com/c-program-implement-rabin-miller-primality-test-check-number-prime/

possiblePrime : The number to be tested for primality
iterations    : The number of iterations that possiblePrime should be tested for primality
*/

bool millerRabinPrimalityTest(long long possiblePrime, int iterations){

    long long s = possiblePrime - 1;

    if(possiblePrime < 2){ //0 and 1 are not a prime numbers

        return false;

    }

    if(possiblePrime != 2 && (possiblePrime % 2) == 0){ //If not 2, and possiblePrime % 2 = 0 is must be an even number

        return false;

    }

    while(s % 2 == 0){ 

        s = s / 2;

    }

    for(int i = 0; i < iterations; i++){ //Perform the Miller-Rabin primality test

        long long a = rand() % (possiblePrime - 1) + 1;
        long long temp = s;
        long long modular = modExp(a, temp, possiblePrime);

        while(temp != (possiblePrime - 1) && modular != 1 && modular != (possiblePrime - 1)){

            modular = multiplyMod(modular, modular, possiblePrime);
            temp = temp * 2;

        }

        if(modular != (possiblePrime - 1) && (temp % 2) == 0){

            return false;

        }
    }

    return true;

}

//###############################################################################
/* Solovay-Strassen Primality Test Methods */
//###############################################################################

/*
calcJacobian:
Calculates the Jacobi Symbol (a/n) for the Solovay-Strassen Primality Test
https://www.sanfoundry.com/cpp-program-implement-solovay-strassen-primality-test/

a : A random number less then b
b : The number being tested for primality. Inherently must be an odd number
*/

long long calcJacobian(long long a, long long b){

    long long answer;
    long long temp;

    if(a == 0){

        return 0;

    }//if(a == 0)

    answer = 1;

    if(a < 0){

        a = -a;

        if(b % 4 == 3){

            answer = -answer;

        }//if(b % 4 == 3)

    }//if(a < 0)

    if(a == 1){

        return answer;

    }//if(a == 1)

    while(a != 0){

        if (a < 0){

            a = -a;

            if(b % 4 == 3){

                answer = -answer;

            }//if(b % 4 == 3)

        }//if (a < 0)

        while(a % 2 == 0){

            a = a / 2;

            if(b % 8 == 3 || b % 8 == 5){

                answer = -answer;

            }//if(b % 8 == 3 || b % 8 == 5)

        }//while(a % 2 == 0)

        temp = a;
        a = b;
        b = temp;

        if(a % 4 == 3 && b % 4 == 3){

            answer = -answer;

        }//if(a % 4 == 3 && b % 4 == 3)

        a = a % b;

        if(a > (b / 2)){

            a = a - b;

        }//if(a > (b / 2))

    }//while(a != 0)

    if(b == 1){

        return answer;

    }

    return 0;

}

/*
solovayStrassenPrimalityTest:
Performs the Solovay-Strassen Primality Test
https://www.sanfoundry.com/cpp-program-implement-solovay-strassen-primality-test/

possiblePrime : The number to be tested for primality
iterations    : The number of iterations that possiblePrime should be tested for primality
*/

bool solovayStrassenPrimalityTest(long long possiblePrime, int iterations){

    if(possiblePrime < 2){ //0 and 1 are not a prime numbers

        return false;

    }

    if(possiblePrime != 2 && (possiblePrime % 2) == 0){ //If not 2, and possiblePrime % 2 = 0 is must be an even number

        return false;

    }

    for(int i = 0; i < iterations; i++){ //Perform the solovay-strassen test.

        long long a = rand() % (possiblePrime - 1) + 1;
        long long jacobian = (possiblePrime + calcJacobian(a, possiblePrime)) % possiblePrime;
        long long modular = modExp(a, ((possiblePrime - 1) / 2), possiblePrime);

        if(jacobian == 0 || modular != jacobian){

            return false;

        }
    }

    return true;

}

//###############################################################################
/* Main */
//###############################################################################

int main(int argc, char *argv[]){

    int my_rank; //the rank of the processor
    int comm_sz; //the number of processors

    threadCount = atoi(argv[1]); //The number of threads to be used

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    srand(time(NULL)); //Initializing seed for random number generation

    if(comm_sz != 3){ //Make sure that the appropriate number of processes has been initialized

        if(my_rank == 0){

            printf("Sorry this program requires exactly 3 processes\n");
            MPI_Finalize();
            exit(0);
        }
        else{

            exit(0);

        }
    }

    if(my_rank == 0){ //Fermat

        MPI_Barrier(MPI_COMM_WORLD);

        printf("Process %d started\n", my_rank);

        int i;

        //spool up threads
        #pragma omp parallel num_threads(threadCount) default(none) shared(fermatTotal, fermatTime, max, numberOfIterations, mersennePrimes) private(i)
        {

            int rank = omp_get_thread_num();
            int partialSum = 0;

            //Barrier threads so none start early
            #pragma omp barrier

            double fTime = omp_get_wtime(); //start timing

            //Spilt up the range of numbers and have each thread calculate its own portion
            #pragma omp for

                for(i = 0; i < max; i++){

                    if(fermatPrimalityTest(i, numberOfIterations) == true){

                        printf("Rank %d checking %d\n", rank, i);
                        partialSum = partialSum + 1;

                    }//if
                }//for
            
            //Barrier threads so that all end at the same time. Perserves total time, not just single thread time
            #pragma omp barrier

            fTime = omp_get_wtime() - fTime; //Get total time elapsed

            //algomate the results from all threads
            #pragma omp critical 
            {

                *fermatTotal = *fermatTotal + partialSum;

            }//omp critical

            if(rank == 0){ //Assume minimal difference in start/stop times between threads and get rank 0's time

                *fermatTime = fTime;

            }
        }//omp parallel

        MPI_Barrier(MPI_COMM_WORLD);

        printf("Process %d done with Fermat Primality Test\n", my_rank);

    }

    if(my_rank == 1){ //Miller

        MPI_Barrier(MPI_COMM_WORLD);

        printf("Process %d started\n", my_rank);

        int i;

        //spool up threads
        #pragma omp parallel num_threads(threadCount) default(none) shared(millerTotal, millerTime, max, numberOfIterations, mersennePrimes) private(i)
        {

            int rank = omp_get_thread_num();
            int partialSum = 0;

            //Barrier threads so none start early
            #pragma omp barrier

            double mTime = omp_get_wtime(); //start timing

            //Spilt up the range of numbers and have each thread calculate its own portion
            #pragma omp for

                for(i = 0; i < max; i++){

                    if(millerRabinPrimalityTest(i, numberOfIterations) == true){

                        partialSum = partialSum + 1;

                    }//if
                }//for

            //Barrier threads so that all end at the same time. Perserves total time, not just single thread time
            #pragma omp barrier

            mTime = omp_get_wtime() - mTime; //Get total time elapsed

            //algomate the results from all threads
            #pragma omp critical
            {

                *millerTotal = *millerTotal + partialSum;

            }//omp critical

            if(rank == 0){ //Assume minimal difference in start/stop times between threads and get rank 0's time

                *millerTime = mTime;

            }
        }//omp parallel

        MPI_Barrier(MPI_COMM_WORLD);

        printf("Process %d done with Miller-Rabin Primality Test\n", my_rank);

    }

    if(my_rank == 2){ //Solovay

        MPI_Barrier(MPI_COMM_WORLD);

        printf("Process %d started\n", my_rank);

        int i;

        //spool up threads
        #pragma omp parallel num_threads(threadCount) default(none) shared(solovayTotal, solovayTime, max, numberOfIterations, mersennePrimes) private(i)
        {

            int rank = omp_get_thread_num();
            int partialSum = 0;

            //Barrier threads so none start early
            #pragma omp barrier

            double sTime = omp_get_wtime(); //start timing

            //Spilt up the range of numbers and have each thread calculate its own portion
            #pragma omp for

                for(int i = 0; i < max; i++){

                    if(solovayStrassenPrimalityTest(i, numberOfIterations) == true){

                        partialSum = partialSum + 1;

                    }
                }
            
            //Barrier threads so that all end at the same time. Perserves total time, not just single thread time
            #pragma omp barrier

            sTime = omp_get_wtime() - sTime; //Get total time elapsed

            //algomate the results from all threads
            #pragma omp critical
            {

                *solovayTotal = *solovayTotal + partialSum;

            }//omp critical

            if(rank == 0){ //Assume minimal difference in start/stop times between threads and get rank 0's time

                *solovayTime = sTime;

            }
        }//omp parallel

        MPI_Barrier(MPI_COMM_WORLD);

        printf("Process %d done with Solovay-Strassen Primality Test\n", my_rank);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Sending totals back to process 0 */
    MPI_Bcast(millerTotal, 1, MPI_INT, 1, MPI_COMM_WORLD);
    MPI_Bcast(millerTime, 1, MPI_DOUBLE, 1, MPI_COMM_WORLD);

    MPI_Bcast(solovayTotal, 1, MPI_INT, 2, MPI_COMM_WORLD);
    MPI_Bcast(solovayTime, 1, MPI_DOUBLE, 2, MPI_COMM_WORLD);

    if(my_rank == 0){ //Print out results for each Primality Test

        printf("\n-------------------------------------------------------------------------\n");
        printf("Number of threads = %d\n",threadCount);
        printf("-------------------------------------------------------------------------\n");

        printf("\n-------------------------------------------------------------------------\n");

        if(max == 100000){

            printf("Number of primes below 100,000:                                  %d\n\n", numberOfPrimesBelow100000);

        }
        else if(max == 1000000){

            printf("Number of primes below 1,000,000:                                %d\n\n", numberOfPrimesBelow1000000);

        }
        else if(max == 10000000){

            printf("Number of primes below 10,000,000:                               %d\n\n", numberOfPrimesBelow10000000);

        }
        else if(max == 100000000){

            printf("Number of primes below 100,000,000:                              %d\n\n", numberOfPrimesBelow100000000);

        }
        else if(max == 1000000000){

            printf("Number of primes below 1,000,000,000:                            %d\n\n", numberOfPrimesBelow1000000000);

        }
    
        printf("Number of primes found with the Fermat Primality Test:           %d\n", *fermatTotal);
        printf("Number of primes found with the Miller-Rabin Primality Test:     %d\n", *millerTotal);
        printf("Number of primes found with the Solovay-Strassen Primality Test: %d\n", *solovayTotal);
        printf("-------------------------------------------------------------------------\n");

        printf("\n-------------------------------------------------------------------------\n");
        printf("Time taken by the Fermat Primality Test:           %.2f\n", *fermatTime);
        printf("Time taken by the Miller-Rabin Primality Test:     %.2f\n", *millerTime);
        printf("Time taken by the Solovay-Strassen Primality Test: %.2f\n", *solovayTime);
        printf("-------------------------------------------------------------------------\n");
    
    }

    MPI_Finalize();

    return 0;

}