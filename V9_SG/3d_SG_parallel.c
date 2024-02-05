//Spin Glass, Parallel tempering
//L.H.Miranda-Filho., ITP-CAC, lucmiranda@gmail.com, 2023
//**************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <omp.h>

#define SCALE 2.328306436e-10


// static int t_PT = 400;
static int t_PT = 400;
int num_neighbours = 6;


//######################################################################################
//Structure definition

#include "global.h"


//######################################################################################
//random number generator
double xor64(int ii, Seed* seeds ){
    unsigned int t = (seeds[ii].ux^(seeds[ii].ux<<8));
    seeds[ii].ux = seeds[ii].uy;
    seeds[ii].uy = (seeds[ii].uy^(seeds[ii].uy>>22))^(t^(t>>9));
    return ((double) seeds[ii].uy) * SCALE;
}

//######################################################################################

//Paralellized functions

#include "mc.c"

void linspace(double start, double end, int numPoints, Sample* samples) {
    double step = (end - start) / (numPoints - 1);

    double bac;

    for (int ii = 0; ii < numPoints; ii++) {
        bac = start + ii * step;

        samples[ii].beta = bac;
    }
}
//######################################################################################

void printScrollingValue(double value, double targetValue) {
    printf("\rCurrent/Target: %.4g/%g ", value, targetValue);
    fflush(stdout);
}

int main(int argc, char** argv) {

    int ii,jj,kk, label;
    char file_name[50];

    sscanf (argv[1],"%d",&label);
//###########################################################################
//variables declaration

    L = 12;
    L2 = L*L;
    N = L*L*L;


    FILE *fin0;
    sprintf(file_name,"par.dat");
    fin0 = fopen(file_name, "r");

    fscanf(fin0,"%d\n", &tobs);
    fscanf(fin0,"%d\n", &t_check);
    fscanf(fin0,"%d\n", &t_rlx);
    fscanf(fin0,"%d\n", &num_threads);
    fscanf(fin0,"%d\n", &num_cores);

    fscanf(fin0,"%lf\n", &beta0);
    fscanf(fin0,"%lf\n", &betaf);

    fclose(fin0);

    fac = 0.5*num_threads*(num_threads-1);

//######################################################################################
//initial conditions definitions

    Sample* samples = (Sample*)malloc(num_threads * sizeof(Sample));
    Seed* seeds = (Seed*)malloc((num_threads+1) * sizeof(Seed));


    for (ii = 0; ii < num_threads; ii++) {
        samples[ii].particles = (Particle*)malloc(N * sizeof(Particle));
        samples[ii].time = (Time*)malloc(tobs*t_check * sizeof(Time));
    }

    double* q = (double*)malloc((fac) * sizeof(double));
    double* chi_plot = (double*)malloc((tobs) * sizeof(double));



    //random number initialization
    srand( (int)time(NULL) );
//     srand( 9090911 );
    for(jj=0; jj< num_threads+1; ++jj){

        for (ii=0; ii<100; ++ii) seeds[jj].ux = rand();
        for (ii=0; ii<100; ++ii) seeds[jj].uy = rand();
    }


    for(ii=0; ii < N; ++ii)
        samples[0].particles[ii].x0 = (xor64(num_threads, seeds) >= 0.5) ? 1 : -1;


    for(ii=0; ii < num_threads; ++ii)
        for(jj=0; jj < N; ++jj)
            samples[ii].particles[jj].x = samples[0].particles[jj].x0;


    //seting index_list
    for(ii=0; ii< num_threads; ++ii)
        samples[ii].index_list = ii;

    for (ii = 0; ii < num_threads; ii++)
        for (jj = 0; jj < N; ++jj)
            samples[ii].particles[jj].J =  (int*)malloc(num_neighbours * sizeof(int));

    srand( 73857 );
    for (ii=0; ii<100; ++ii) seeds[num_threads].ux = rand();
    for (ii=0; ii<100; ++ii) seeds[num_threads].uy = rand();

    for(ii=0; ii<N; ++ii){

        samples[0].particles[ii].J[0] = (xor64(num_threads, seeds) >= 0.5) ? 1 : -1; //Right
        samples[0].particles[ nl(ii,0) ].J[1] = samples[0].particles[ii].J[0];

        samples[0].particles[ii].J[2] = (xor64(num_threads, seeds) >= 0.5) ? 1 : -1; //Up
        samples[0].particles[ nl(ii,2) ].J[3] = samples[0].particles[ii].J[2];

        samples[0].particles[ii].J[4] = (xor64(num_threads, seeds) >= 0.5) ? 1 : -1; //Up2
        samples[0].particles[ nl(ii,4) ].J[5] = samples[0].particles[ii].J[4];

    }

//     for(ii=0; ii<N; ++ii)
//         for(jj=0; jj<num_neighbours; ++jj)
//             samples[0].particles[ii].J[jj] = 1;

//##################################
//#initial energy

    for(ii=0; ii < num_threads; ++ii)
        energy(ii, samples);

//##################################
// //#beta input from par.dat

//     if(beta0 == betaf)
//         for(ii=0; ii < num_threads; ++ii)
//             samples[ii].beta = beta0;
//
//     else
//         linspace(beta0, betaf, num_threads, samples);
//##################################

//##################################
// //#beta input from beta_in.dat
    FILE *fin1;
    sprintf(file_name,"beta_in.dat");
    fin1 = fopen(file_name, "r");

    for(ii=0; ii < num_threads; ++ii)
        fscanf(fin1,"%lf\n", &samples[ii].beta);
    fclose(fin1);
//##################################


    //seting counter
    for(ii=0; ii< num_threads; ++ii)
        samples[ii].accept_num_PT_bac = samples[ii].accept_num_PT = 0;

//##################################
//#Dynamics

    printf("\n#############################\n");
    printf("label = %d, N(L3) = %d\n\n", label, N);

    printf("num_threads = %d\n\n", num_threads);

    printf("tobs = %.1e, t_check = %d, t_rlx = %d, t_PT = %d\n", (double)tobs, t_check, t_rlx, t_PT);

    printf("\nx[0]=\n");
    printf("%d %d %d %d", samples[0].particles[0].x, samples[0].particles[1].x, samples[0].particles[2].x, samples[0].particles[3].x );
    printf("...\n#############################\n");


    printf("\nJ=\n");
    for(ii=0; ii< 3; ++ii){
        printf("%d %d %d ", samples[0].particles[ii].J[0], samples[0].particles[ii].J[1], samples[0].particles[ii].J[2] );
        printf("%d %d %d\n", samples[0].particles[ii].J[3], samples[0].particles[ii].J[4], samples[0].particles[ii].J[5] );
    }

    printf("...\n#############################\n\n");


//relaxation

    omp_set_num_threads(num_cores);

//     jj = 0;
    for (ii = 0; ii < t_rlx; ii++){

        #pragma omp parallel
        {

            int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
            int thread_id = omp_get_thread_num();

            for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++)
                MC(ll, samples, seeds);
        }

        printScrollingValue((double)ii, (double)t_rlx);


//         if(ii%10==0){
//             MC_PT(num_threads, samples, jj, seeds);
//
//             if((jj+1)%t_PT==0)
//                 optm_T(samples);
//             jj+=1;
//         }

    }

    printf("\n\n");
    printf("beta\t\te\t\tindex\t\taccept\n");
    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        printf("%.2lf\t\t%.2e\t\t%d\t\t%.2g\n", samples[ii].beta, (double)samples[jj].e, samples[ii].index_list, samples[ii].accept_num_PT_bac);
    }


    printf("\n\n");


    for(ii=0; ii < num_threads; ++ii)
        for(jj=0; jj < N; ++jj)
            samples[ii].particles[jj].x0 = samples[ii].particles[jj].x; //redefinition x0 (after wait time t_rlx)


    jj = 0;
    for (ii = 0; ii < tobs; ii++){

        #pragma omp parallel
        {

            int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
            int thread_id = omp_get_thread_num();

            for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++){
                MC(ll, samples, seeds);
                samples[ ll ].time[ii].q_plot = q_tcf(ll, samples);
            }
        }

//         printScrollingValue((double)ii, (double)tobs);
//         chi_plot[ii] = chi(q, samples);


        if(ii%10==0){
            MC_PT(num_threads, samples, jj, seeds);

            if((jj+1)%t_PT==0){
                optm_T_mod(samples);
                printScrollingValue((double)ii, (double)tobs);
            }
            jj+=1;
        }
    }

    #pragma omp parallel
    {

        int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
        int thread_id = omp_get_thread_num();

        for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++)
            MC(ll, samples, seeds);

    }

    printf("\n\n");


// ##################################
// #plot

    double av;

//     FILE *fout2;
//     sprintf(file_name,"./output/chi_%.2d.dat",label);
//     fout2 = fopen(file_name, "w");
//
//     for (ii = 1; ii < tobs; ii++)
//         fprintf(fout2,"%d %g\n", ii*t_check, chi_plot[ii]);
//
//     fclose(fout2);

// ##################################
// #time autocorrelation function

//     FILE *fout3;
//     sprintf(file_name,"./output/q_tcf_%.2d.dat",label);
//     fout3 = fopen(file_name, "w");
//
//
//
// //     samples with same beta value (no replica exchanges)
//     for (ii = 0; ii < tobs; ii++){
//
//         av = 0;
//
//         for (kk = 0; kk < num_threads; kk++)
//             av += samples[ kk ].time[ii].q_plot;
//
//         fprintf(fout3,"%d %g\n", ii, av/(double)num_threads);
//
//     }
//
//     fclose(fout3);

    FILE *fout3;
    sprintf(file_name,"./output/q_tcf_%.2d.dat",label);
    fout3 = fopen(file_name, "w");

    for (ii = 0; ii < num_threads; ii++){


        for (kk = 0; kk < tobs; kk++)
            fprintf(fout3,"%d %g\n", kk*t_check, samples[ ii ].time[kk].q_plot);

//         fprintf(fout3,"\n");

    }

    fclose(fout3);

// ##################################

    FILE *fout4;
    sprintf(file_name,"./output/report_%.2d.dat",label);
    fout4 = fopen(file_name, "w");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        fprintf(fout4,"%.2lf %d %d %.2g %d %d\n", samples[ii].beta, samples[jj].e, samples[ii].index_list, samples[ii].accept_num_PT_bac, tobs, t_rlx);
    }

    fclose(fout4);

//     FILE *fout5;
//     sprintf(file_name,"./output/EA_dist_%.2d.dat",label);
//     fout5 = fopen(file_name, "w");
//     for (ii = 0; ii < num_threads; ii++)
//         for (jj = ii+1; jj < num_threads; jj++)
//             fprintf(fout5,"%g\n", q_tcf_mod(ii,jj,samples));
//     fclose(fout5);

//     kk = 28;
//     printf("sample target, beta_index:%d   samples_index:%d,   beta:%g\n", kk, samples[kk].index_list, samples[kk].beta);
//     FILE *fout5;
//     sprintf(file_name,"./output/conf_%.2d.dat",label);
//     fout5 = fopen(file_name, "w");
//
//     jj = samples[kk].index_list;
//     for (ii = 0; ii < N; ii++)
//         fprintf(fout5,"%d\n", samples[jj].particles[ ii ].x);
//
//     fclose(fout5);

    printf("\n\n");
    printf("beta\t\te\t\tindex\t\taccept\n");
    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        printf("%.2lf\t\t%.2e\t\t%d\t\t%.2g\n", samples[ii].beta, (double)samples[jj].e, samples[ii].index_list, samples[ii].accept_num_PT_bac);
    }

    printf("\n\n");

    return 0;
}

