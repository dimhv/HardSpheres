//Hard Sphere System, Monte Carlo
//L.H.Miranda-Filho., ITP-CAS, lucmiranda@gmail.com, 2023
//**************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <omp.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#pragma GCC diagnostic ignored "-Wunused-result"

//######################################################################################
//Variables definitio//#################################################################

static int L = 10;
static int L2 = 100;

static double beta = 1.;

double mesh_fac = 2.; //open for adjustment
int t_PT = 400;

int nhis = 8000;
// double bins[8000];


#include "global.h"

//######################################################################################
//random number generator//#############################################################

#define SCALE 2.328306436e-10
double xor64(Seed* s){
    unsigned int t = (s->ux^(s->ux<<8));
    s->ux = s->uy;
    s->uy = (s->uy^(s->uy>>22))^(t^(t>>9));
    return ((double) s->uy) * SCALE;
}

//######################################################################################
//functions//###########################################################################
#include "mc.c"

//######################################################################################
//Parallel Tempering//##################################################################
#include "optm_T.c"
#include "MC_PT.c"

//######################################################################################
//Utilities
#include "testz.c"
#include "ovito_plot.c"
#include "utilities.c"

//###########################################################################

void printScrollingValue(double value, double targetValue) {
    printf("\rCurrent/Target: %.4g/%g ", value, targetValue);
    fflush(stdout);
}

void p_evaluation(int thread_id, Sample* samples, Seed* seeds){

    int ii, jj, kk;
    kk = 500;

    for(jj=0; jj < nhis; ++jj)
        samples[ samples[thread_id].index_list ].bins[jj] = 0.;

    for(ii=0; ii< kk; ++ii){

//         printf("%d\n", thread_id);

        (void)legal_configuration_dyn(thread_id, samples, seeds);
        g_r_dist_count(thread_id, samples);
    }

    samples[ samples[thread_id].index_list ].p_measured = g_r_dist_return(thread_id, samples, (double)kk);

    return;
}

int main(int argc, char** argv) {

    int label;
    sscanf (argv[1],"%d",&label);

    char file_name[50];

    FILE *fin0;
    sprintf(file_name,"par.dat");
    fin0 = fopen(file_name, "r");

    fscanf(fin0,"%d\n", &N);
    fscanf(fin0,"%d\n", &num_threads);
    fscanf(fin0,"%d\n", &num_cores);
    fscanf(fin0,"%d\n", &t_check);
    fscanf(fin0,"%d\n", &tobs);
    int tw;
    fscanf(fin0,"%d\n", &tw);
    fscanf(fin0,"%lf\n", &acept_Vbac);

    acept_Xbac = 0.;

    fclose(fin0);

    int input_index;//initial condition (independent)
    sscanf (argv[2],"%d",&input_index);

//     int input_index2;//initial condition (replica)
//     sscanf (argv[3],"%d",&input_index2);
//label, initial condition (replica)

    L_mesh = L/mesh_fac;
    L2_mesh = L_mesh*L_mesh;
    cell_num = L2_mesh*L_mesh;


//######################################################################################
//loop counters//#######################################################################

    int ii,jj,kk;

//######################################################################################
//structures initialization//###########################################################

    Sample* samples = (Sample*)malloc(num_threads * sizeof(Sample));
    Seed* seeds = (Seed*)malloc( (num_threads+1) * sizeof(Seed));

    FILE *fin[num_threads];
    FILE *fin2[num_threads];

    double A0, B0;
    double ampF_0;
    for (ii = 0; ii < num_threads; ii++) {

//         printf("\nCHECK %d\n", ii);
        sprintf(file_name,"./input/output_%.2d/conf/par_%.3d_%.3d.dat", input_index, ii, label);


        fin[ii] = fopen(file_name, "r");
        fscanf(fin[ii],"%lf\n", &ampF_0);
        fscanf(fin[ii],"%lf\n", &A0);
        fscanf(fin[ii],"%lf\n", &B0);
        fscanf(fin[ii],"%lf\n", &samples[ii].P);
        fclose(fin[ii]);


        samples[ii].particles = (Particle*)malloc(N * sizeof(Particle));
        samples[ii].time = (Time*)malloc(tobs * sizeof(Time));

        samples[ii].A = 0.08;
        samples[ii].B = 0.005;

        samples[ii].amp_fac = ampF_0;

        samples[ii].V = L2*L/ampF_0/ampF_0/ampF_0;
//         samples[ii].V = L2/ampF_0/ampF_0;

        samples[ii].cell_ind = (int*)malloc(cell_num * sizeof(int));
        samples[ii].cell = (int**)malloc(cell_num * sizeof(int*));
        for (int jj = 0; jj < cell_num; jj++) {
            samples[ii].cell[jj] = (int*)malloc(N * sizeof(int));
        }

        samples[ii].bins = (int*)malloc(nhis * sizeof(int));

        sprintf(file_name,"./input/output_%.2d/conf/c_%.3d_%.3d.dat", input_index, ii, label);

        fin2[ii] = fopen(file_name, "r");
        for (jj = 0; jj < N; jj++)
            fscanf(fin2[ii],"%lf %lf %lf %lf %d %d %d\n", &samples[ii].particles[jj].x0, &samples[ii].particles[jj].y0, &samples[ii].particles[jj].z0, &samples[ii].particles[jj]._2r, &samples[ii].particles[jj].x0_pbc, &samples[ii].particles[jj].y0_pbc, &samples[ii].particles[jj].z0_pbc);
//             fscanf(fin2[ii],"%lf %lf %lf %lf\n", &samples[ii].particles[jj].x0, &samples[ii].particles[jj].y0, &samples[ii].particles[jj].z0, &samples[ii].particles[jj]._2r);

            fclose(fin2[ii]);
    }


//######################################################################################
//seeds random numbers settings//#######################################################

// long int bac_seed = 27278911;
long int bac_seed =  (int)time(NULL);

//         srand( (int)time(NULL) );
    srand( bac_seed );

    for (jj = 0; jj < num_threads+1; jj++) {

        for (ii=0; ii<100; ii++) seeds[jj].ux = rand();
        for (ii=0; ii<100; ii++) seeds[jj].uy = rand();


    }

//######################################################################################
//initial configuration setting//#######################################################

    for (ii = 0; ii < num_threads; ii++)
        for (jj = 0; jj < N; jj++){
            samples[ii].particles[jj].x = samples[ii].particles[jj].x0;
            samples[ii].particles[jj].y = samples[ii].particles[jj].y0;
            samples[ii].particles[jj].z = samples[ii].particles[jj].z0;

            samples[ii].particles[jj].x_pbc = samples[ii].particles[jj].x0_pbc;
            samples[ii].particles[jj].y_pbc = samples[ii].particles[jj].y0_pbc;
            samples[ii].particles[jj].z_pbc = samples[ii].particles[jj].z0_pbc;

        }

    for (ii = 0; ii < num_threads; ii++)
        for (jj = 0; jj < cell_num; jj++)
            samples[ii].cell_ind[jj] = 0;

    for (kk = 0; kk < num_threads; kk++)
        for(ii=0; ii < N; ii++){

            jj = (int)(samples[kk].particles[ii].x/mesh_fac)+L_mesh*(int)(samples[kk].particles[ii].y/mesh_fac)+L2_mesh*(int)(samples[kk].particles[ii].z/mesh_fac);

            samples[kk].cell[ jj ][ samples[kk].cell_ind[jj] ] = ii;
            samples[kk].cell_ind[jj] += 1;
        }

//######################################################################################
//seting pressure distribution//########################################################

//     linspace(1./samples[0].P, 1./samples[num_threads-1].P, num_threads, samples);
//     logspace(P0, Pf, num_threads, samples);

    for(ii=0; ii< num_threads; ++ii){

        samples[ii].P = (1./samples[ii].P)*(double)N/samples[ii].V;

    }

//######################################################################################
//reseting index list//#################################################################

    for(ii=0; ii< num_threads; ++ii)
        samples[ii].index_list = ii;

//######################################################################################
//reseting counters//###################################################################

    for(ii=0; ii< num_threads; ++ii){
        samples[ii].accept_num_PT_bac = samples[ii].accept_num_PT = 0;
        samples[ii].accept_mov = samples[ii].accept_V = 0;
    }
//######################################################################################
//#Dynamics//###########################################################################


    printf("\n#############################\n");
    printf("label = %d, N = %d, L = %d\n\n", label, N, L);
    printf("V[0] = %g, mesh_fac = %g, cell_num = %d, P[0] = %g\n", samples[0].V, mesh_fac,cell_num, samples[0].P);

    printf("num_threads = %d, seed = %ld\n\n", num_threads, bac_seed);

    printf("tobs = %.1e, t_check = %d, t_PT = %d, acept_Vbac = %g, acept_Xbac = %g, tw = %d\n", (double)tobs, t_check, t_PT, acept_Vbac, acept_Xbac, tw);
    printf("A0 = %g, B0 = %g\n", samples[0].A, samples[0].B);
    printf("#############################\n\n");


    omp_set_num_threads(num_cores);

    #pragma omp parallel
    {
        int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
        int thread_id = omp_get_thread_num();

        for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++)
            testz(ll, samples);
    }

    printf("\n");


    #pragma omp parallel

    {
        int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
        int thread_id = omp_get_thread_num();

        for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++){
            V_hard_sphere(ll, samples);
        }
    }

    printf("\n\ntw\n");

    for (ii = 0; ii < tw; ii++){

        #pragma omp parallel
        {
            int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
            int thread_id = omp_get_thread_num();

            for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++)
                MC(ll, samples, seeds);

        }

        printScrollingValue((double)ii, (double)tw);
    }


    for(ii=0; ii< num_threads; ++ii){
        samples[ii].accept_num_PT_bac = samples[ii].accept_num_PT = 0;
        samples[ii].accept_mov = samples[ii].accept_V = 0;

        samples[ii].accept_av_PT = 0;


        com(ii, samples);
    }


    for (ii = 0; ii < num_threads; ii++){
        for (jj = 0; jj < N; jj++){
            samples[ii].particles[jj].x0 = samples[ii].particles[jj].x;
            samples[ii].particles[jj].y0 = samples[ii].particles[jj].y;
            samples[ii].particles[jj].z0 = samples[ii].particles[jj].z;
        }
    }

    for (ii = 0; ii < num_threads; ii++)
        com0(ii, samples);

    printf("\n\ntobs\n");

    int frame = 1000;

    for (ii = 0; ii < tobs; ii++){

        #pragma omp parallel
        {
            int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
            int thread_id = omp_get_thread_num();

            for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++)
                MC(ll, samples, seeds);


        }

        if(ii%frame==0){
            #pragma omp parallel
            {
                int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
                int thread_id = omp_get_thread_num();

                for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++){

                    p_evaluation(ll, samples, seeds);
                    samples[samples[ll].index_list].time[ii/frame].p_plot = samples[ samples[ll].index_list ].p_measured;
                    samples[samples[ll].index_list].time[ii/frame].pf_plot = samples[ samples[ll].index_list ].V_hs/samples[ samples[ll].index_list ].V;
                    samples[ll].time[ii/frame].volume_plot = samples[samples[ll].index_list].V;


                }
            }
        }

        printScrollingValue((double)ii, (double)tobs);

        MC_PT(samples, ii, seeds);
        if((ii+1)%(t_PT)==0)
            optm_T( samples );
    }

    printf("\n\n");


    #pragma omp parallel
    {
        int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
        int thread_id = omp_get_thread_num();

        for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++)
            testz(ll, samples);
    }

    printf("\n\n");

    #pragma omp parallel

    {
        int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
        int thread_id = omp_get_thread_num();

        for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++){
            V_hard_sphere(ll, samples);
            p_evaluation(ll, samples, seeds);

        }
    }
//######################################################################################
//plot//################################################################################

    printf("\n");


    FILE *fout;
    sprintf(file_name,"./output/report_%.2d.dat",label);
    fout = fopen(file_name, "w");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        fprintf(fout,"%g %g %g %g %g %d %g %g\n", samples[jj].amp_fac, samples[ii].P, samples[jj].V, samples[jj].A, samples[jj].B, samples[ii].index_list, samples[ii].accept_num_PT_bac, samples[jj].p_measured);
    }

    fclose(fout);


    FILE *fout2;
    sprintf(file_name,"./output/volume_%.2d.dat",label);
    fout2 = fopen(file_name, "w");

    double bac_plot;


    for (jj = 0; jj < num_threads; jj++){

        for (ii = 0; ii < tobs/frame; ii++)
            fprintf(fout2,"%d %g\n", ii*t_check*frame, samples[jj].time[ii].volume_plot);

        fprintf(fout2,"\n");

    }


    fclose(fout2);

    FILE *fout3;
    sprintf(file_name,"./output/prss_%.2d.dat",label);
    fout3 = fopen(file_name, "w");


    for (jj = 0; jj < num_threads; jj++){

        for (ii = 0; ii < tobs/frame; ii++)
            fprintf(fout3,"%d %g\n", ii*t_check*frame, samples[jj].time[ii].p_plot);

        fprintf(fout3,"\n");

    }


    fclose(fout3);

    FILE *fout4;
    sprintf(file_name,"./output/accept_av_plot_%.2d.dat",label);
    fout4 = fopen(file_name, "w");


    for (jj = 0; jj < num_threads; jj++)
        fprintf(fout4,"%d %g\n", jj, samples[jj].accept_av_PT/( tobs/t_PT) );


    fclose(fout4);

    FILE *fout5;
    sprintf(file_name,"./output/pf_%.2d.dat",label);
    fout5 = fopen(file_name, "w");


    for (jj = 0; jj < num_threads; jj++){

        for (ii = 0; ii < tobs/frame; ii++)
            fprintf(fout5,"%d %g\n", ii*t_check*frame, samples[jj].time[ii].pf_plot);

        fprintf(fout5,"\n");

    }

    fclose(fout5);

    printf("amp\t\tP\t\tV\t\tA\t\tacX\t\tB\t\tacV\t\tid\taccept\n");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        printf("%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%d\t%g\n", samples[jj].amp_fac, samples[ii].P, samples[jj].V, samples[jj].A, samples[jj].accept_mov, samples[jj].B, samples[jj].accept_V, samples[ii].index_list, samples[ii].accept_num_PT_bac);
    }
    printf("\n");

    for (ii = 0; ii < num_threads; ii++)
        ovito_plot(ii, samples, label, 0);

    return 0;
}

