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

//######################################################################################
//Variables definitio//#################################################################

static int L = 10;
static int L2 = 100;

static double beta = 1.;
static double k = 5.25;

double mesh_fac = 2.; //open for adjustment
int t_PT = 400;

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

int main(int argc, char** argv) {
    #pragma GCC diagnostic ignored "-Wunused-result"

    int label;
    sscanf (argv[1],"%d",&label);

    char file_name[50];

    FILE *fin0;
    sprintf(file_name,"par.dat");
    fin0 = fopen(file_name, "r");
//     fscanf(fin,"%d %lf %lf %lf %lf %lf\n", &N, &amp_fac0, &A0, &B0, &P0, &Pf);
    fscanf(fin0,"%d\n", &N);
    fscanf(fin0,"%lf\n", &P0);
    fscanf(fin0,"%lf\n", &Pf);
    fscanf(fin0,"%d\n", &num_threads);
    fscanf(fin0,"%d\n", &num_cores);
    fscanf(fin0,"%d\n", &t_check);
    fscanf(fin0,"%d\n", &tobs);
    fscanf(fin0,"%lf\n", &acept_Vbac);

    int input_index;
    fscanf(fin0,"%d\n", &input_index);
//     sscanf (argv[2],"%d",&input_index);

    fclose(fin0);

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


    double A0, B0, ampF_0;

    for (ii = 0; ii < num_threads; ii++) {

//         sprintf(file_name,"par_xyz.dat");
        sprintf(file_name,"./input/par_%.5d_000_000.dat",input_index);

        fin[ii] = fopen(file_name, "r");
        fscanf(fin[ii],"%lf\n", &ampF_0);
        fscanf(fin[ii],"%lf\n", &A0);
        fscanf(fin[ii],"%lf\n", &B0);
        fclose(fin[ii]);

        samples[ii].particles = (Particle*)malloc(N * sizeof(Particle));
        samples[ii].time = (Time*)malloc(tobs * sizeof(Time));

        samples[ii].A = 1.;
        samples[ii].B = B0;

        samples[ii].amp_fac = ampF_0;

        samples[ii].V = L2*L/ampF_0/ampF_0/ampF_0;
//         samples[ii].V = L2/ampF_0/ampF_0;

        samples[ii].cell_ind = (int*)malloc(cell_num * sizeof(int));
        samples[ii].cell = (int**)malloc(cell_num * sizeof(int*));
        for (int jj = 0; jj < cell_num; jj++) {
            samples[ii].cell[jj] = (int*)malloc(N * sizeof(int));
        }

//         sprintf(file_name,"xyz.dat");
        sprintf(file_name,"./input/c_%.5d_000_000.dat",input_index);

        fin2[ii] = fopen(file_name, "r");
        for (jj = 0; jj < N; jj++)
            fscanf(fin2[ii],"%lf %lf %lf %lf\n", &samples[ii].particles[jj].x0, &samples[ii].particles[jj].y0, &samples[ii].particles[jj].z0, &samples[ii].particles[jj]._2r);
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
        }


    for (kk = 0; kk < num_threads; kk++)
        for(ii=0; ii < N; ii++){

            jj = (int)(samples[kk].particles[ii].x/mesh_fac)+L_mesh*(int)(samples[kk].particles[ii].y/mesh_fac)+L2_mesh*(int)(samples[kk].particles[ii].z/mesh_fac);

            samples[kk].cell[ jj ][ samples[kk].cell_ind[jj] ] = ii;
            samples[kk].cell_ind[jj] += 1;
        }

//######################################################################################
//seting pressure distribution//########################################################

//     linspace(P0, Pf, num_threads, samples);
//     logspace(P0, Pf, num_threads, samples);

    for(ii=0; ii< num_threads; ++ii){
        samples[ii].P = 0;
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

        for (jj = 0; jj < N; jj++){
            samples[ii].particles[jj].x_pbc = 0;
            samples[ii].particles[jj].y_pbc = 0;
            samples[ii].particles[jj].z_pbc = 0;
        }
    }
//######################################################################################
//#Dynamics//###########################################################################


    printf("\n#############################\n");
    printf("label = %d (input_file: %d), N = %d, L = %d, P0,Pf = %g,%g\n\n", label, input_index, N, L, P0,Pf);
    printf("V[0] = %g, mesh_fac = %g, cell_num = %d\n", samples[0].V, mesh_fac,cell_num);

    printf("num_threads = %d, seed = %ld\n\n", num_threads, bac_seed);

    printf("tobs = %.1e, t_check = %d, t_PT = %d, acept_Vbac = %g\n", (double)tobs, t_check, t_PT, acept_Vbac);
    printf("A0 = %g, B0 = %g\n", A0, B0);
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

    int trlx;
    trlx = 1000/t_check;

    printf("\n\ntrlx\n");

    for (ii = 0; ii < trlx; ii++){

        #pragma omp parallel
        {
            int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
            int thread_id = omp_get_thread_num();

            for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++)
                MC(ll, samples, seeds);

        }

        printScrollingValue((double)ii, (double)trlx);
    }


    for(ii=0; ii< num_threads; ++ii){
        samples[ii].accept_num_PT_bac = samples[ii].accept_num_PT = 0;
        samples[ii].accept_mov = samples[ii].accept_V = 0;

        for (jj = 0; jj < N; jj++){
            samples[ii].particles[jj].x_pbc = 0;
            samples[ii].particles[jj].y_pbc = 0;
            samples[ii].particles[jj].z_pbc = 0;
        }

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

    for (ii = 0; ii < num_threads; ii++)
        ovito_plot(ii, samples, label, 0);

    printf("\n\ntobs\n");

    int frame = 200;
    int t_plot = 0.9*tobs;

    for (ii = 0; ii < tobs; ii++){

        #pragma omp parallel
        {
            int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
            int thread_id = omp_get_thread_num();

            for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++){

                MC(ll, samples, seeds);
                samples[ll].time[ii].volume_plot = samples[samples[ll].index_list].V;

                samples[ll].time[ii].msd_plot = msd(ll, samples);
                samples[ll].time[ii].com_plot = samples[ll].COM;

            }
        }

        printScrollingValue((double)ii, (double)tobs);

//         MC_PT(samples, ii, seeds);
//         if((ii+1)%t_PT==0)
//             optm_T( samples );
    }

    printf("\n");


    #pragma omp parallel
    {
        int num_iterations_per_thread = (num_threads + (num_cores-1)) / num_cores; // Round up division
        int thread_id = omp_get_thread_num();

        for (int ll = thread_id * num_iterations_per_thread; ll < (thread_id + 1) * num_iterations_per_thread && ll < num_threads; ll++)
            testz(ll, samples);
    }

//######################################################################################
//plot//################################################################################

    printf("\n");

    FILE *fout0;
    sprintf(file_name,"./output/p_pacF_%.2d.dat",label);
    fout0 = fopen(file_name, "w");

    for (ii = 0; ii < num_threads; ii++){

        jj = samples[ii].index_list;
        fprintf(fout0,"%g\n", samples[jj].V_hs/samples[jj].V);
//         fprintf( fout0,"%g %g\n", samples[jj].V_hs/(L2*L/pow(samples[jj].amp_fac,3)), samples[ii].P*(L2*L/pow(samples[jj].amp_fac,3))/(double)N);
    }
    fclose(fout0);

    FILE *fout1;
    sprintf(file_name,"./output/r_%.2d.dat",label);
    fout1 = fopen(file_name, "w");

    for (ii = 0; ii < N; ii++)
        fprintf(fout1,"%g\n", samples[0].particles[ii]._2r);

    fclose(fout1);

    FILE *fout2;
    sprintf(file_name,"./output/report_%.2d.dat",label);
    fout2 = fopen(file_name, "w");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        fprintf(fout2,"%g %g %g %g %g %d %g\n", samples[jj].amp_fac, samples[ii].P, samples[jj].V, samples[jj].A, samples[jj].B, samples[ii].index_list, samples[ii].accept_num_PT_bac);
    }

    fclose(fout2);

    FILE *fout3;
    sprintf(file_name,"./output/volume_%.2d.dat",label);
    fout3 = fopen(file_name, "w");

    for (ii = 0; ii < num_threads; ii++){

        kk = samples[ii].index_list;

        for (jj = 0; jj < tobs; jj++)
            fprintf(fout3,"%g\n", samples[kk].time[jj].volume_plot);

        fprintf(fout3,"\n");
    }
    fclose(fout3);


    FILE *fout4;
    sprintf(file_name,"./output/msd_av_%.2d.dat",label);
    fout4 = fopen(file_name, "w");

    double bac_plot;

    bac_plot=0;

    for (jj = 0; jj < num_threads; jj++)
        bac_plot += samples[jj].time[0].msd_plot;

    fprintf(fout4,"%d %g\n", 1, bac_plot/(double)num_threads);

    for (ii = 1; ii < tobs; ii++){

        bac_plot=0;

        for (jj = 0; jj < num_threads; jj++)
            bac_plot += samples[jj].time[ii].msd_plot;

        fprintf(fout4,"%d %g\n", ii*t_check, bac_plot/(double)num_threads);

    }
    fclose(fout4);

    FILE *fout5;
    sprintf(file_name,"./output/msd_%.2d.dat",label);
    fout5 = fopen(file_name, "w");


    for (jj = 0; jj < num_threads; jj++){

        for (ii = 0; ii < tobs; ii++)
            fprintf(fout5,"%d %g\n", (ii+1)*t_check, samples[jj].time[ii].msd_plot);

    fprintf(fout5,"\n");

    }
    fclose(fout5);

    FILE *fout6;
    sprintf(file_name,"./output/com%.2d.dat",label);
    fout6 = fopen(file_name, "w");


    for (jj = 0; jj < num_threads; jj++){

        for (ii = 0; ii < tobs; ii++)
            fprintf(fout6,"%d %g\n", (ii+1)*t_check, samples[jj].time[ii].com_plot);

    fprintf(fout6,"\n");

    }
    fclose(fout6);

    printf("amp\tP\tV\tA\tacX\tB\tacV\tid\taccept\n");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        printf("%.4g\t%.2g\t%.3g\t%.3g\t%.2g\t%.2g\t%.2g\t%d\t%g\n", samples[jj].amp_fac, samples[ii].P, samples[jj].V, samples[jj].A, samples[jj].accept_mov, samples[jj].B, samples[jj].accept_V, samples[ii].index_list, samples[ii].accept_num_PT_bac);
    }
    printf("\n");

    for (ii = 0; ii < num_threads; ii++)
        ovito_plot(ii, samples, label, 1);

    return 0;
}

