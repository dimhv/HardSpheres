//Hard Sphere System, Monte Carlo
//L.H.Miranda-Filho., ITP-CAS, lucmiranda@gmail.com, 2023
//**************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


//######################################################################################
//Variables definitio//#################################################################

static int L = 10;
static int L2 = 100;

static double beta = 1.;
static double k = 5.25;

static int num_threads = 1;

static int t_max = (int)1e8;

double mesh_fac = 2.5; //open for adjustment

#include "global.h"

//######################################################################################
//random number generator//#############################################################

#define SCALE 2.328306436e-10

double xor64(int ii, Seed* seeds ){
    unsigned int t = (seeds[ii].ux^(seeds[ii].ux<<8));
    seeds[ii].ux = seeds[ii].uy;
    seeds[ii].uy = (seeds[ii].uy^(seeds[ii].uy>>22))^(t^(t>>9));
    return ((double) seeds[ii].uy) * SCALE;
}
//######################################################################################
//functions//###########################################################################
#include "mc.c"

//######################################################################################
//Parallel Tempering//##################################################################

//######################################################################################
//Utilities
#include "testz.c"
#include "ovito_plot.c"
#include "utilities.c"

//###########################################################################

void printScrollingValue(double value0, double targetValue0, double value1, double targetValue1, double value2, double targetValue2, double value3, double value4, double value5) {
    printf("\rp: %.3g / %g ## pc_fr: %.3g / %g ## steps: %.2g / %.2g ## B = %.2g ## accept_V = %.2g ## P = %.2g  ", value0, targetValue0, value1, targetValue1, value2, targetValue2, value3, value4, value5);
    fflush(stdout);
}

int main(int argc, char** argv) {

    int label;
    sscanf (argv[1],"%d",&label);

    char file_name[50];
    FILE *fin0;
    sprintf(file_name,"par.dat");
    fin0 = fopen(file_name, "r");

    double rate_bac, P0_label;

    fscanf(fin0,"%d\n", &N);
    fscanf(fin0,"%lf\n", &rate_bac);
    fscanf(fin0,"%d\n", &t_check);
    fscanf(fin0,"%lf\n", &p_max);
    fscanf(fin0,"%lf\n", &pc_fr_max);
    fscanf(fin0,"%lf\n", &B_accept);
    fscanf(fin0,"%lf\n", &P0_label);

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


        sprintf(file_name,"par_xyz.dat");
        fin[ii] = fopen(file_name, "r");
        fscanf(fin[ii],"%lf\n", &ampF_0);
        fscanf(fin[ii],"%lf\n", &A0);
        fscanf(fin[ii],"%lf\n", &B0);
        fclose(fin[ii]);


        samples[ii].particles = (Particle*)malloc(N * sizeof(Particle));

        samples[ii].A = A0;
        samples[ii].B = B0;

        samples[ii].amp_fac = ampF_0;

        samples[ii].V = L2*L/ampF_0/ampF_0/ampF_0;
//         samples[ii].V = L2/ampF_0/ampF_0;

        samples[ii].cell_ind = (int*)malloc(cell_num * sizeof(int));
        samples[ii].cell = (int**)malloc(cell_num * sizeof(int*));
        for (int jj = 0; jj < cell_num; jj++) {
            samples[ii].cell[jj] = (int*)malloc(N * sizeof(int));
        }

        sprintf(file_name,"xyz.dat");
        fin2[ii] = fopen(file_name, "r");
        for (jj = 0; jj < N; jj++)
            fscanf(fin2[ii],"%lf %lf %lf %lf\n", &samples[ii].particles[jj].x0, &samples[ii].particles[jj].y0, &samples[ii].particles[jj].z0, &samples[ii].particles[jj]._2r);
        fclose(fin2[ii]);



        samples[ii].rate = rate_bac;

    }

//######################################################################################
//seeds random numbers settings//#######################################################
int seed_bac = (int)time(NULL);
srand( seed_bac );
//     srand( 27278911 );

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
//reseting index list//#################################################################

    for(ii=0; ii< num_threads; ++ii)
        samples[ii].index_list = ii;

//######################################################################################
//reseting counters//###################################################################

    for(ii=0; ii< num_threads; ++ii)
       samples[ii].accept_num_PT_bac = samples[ii].accept_num_PT = 0;

//######################################################################################
//seting pressure distribution//########################################################

    double bac_double = 0.;
    for(ii=0; ii< 100; ++ii){
        (void)legal_configuration_dyn(0, samples, seeds);
        bac_double += g_r_dist(0, samples);
    }

    bac_double /= 100.;
    P0 = bac_double*samples[0].V/(double)N;

    if(P0_label==0)
        samples[0].P = P0;

    else
        samples[0].P = 1000;
//######################################################################################
//#Dynamics//###########################################################################


    printf("\n#############################\n");
    printf("label = %d, N = %d, L = %d, p0(sample) = %g\n\n", label, N, L, bac_double);
    printf("V[0] = %g, mesh_fac = %g, cell_num = %d\n", samples[0].V, mesh_fac,cell_num);

    printf("num_threads = %d, seed = %d\n\n", num_threads, seed_bac);

    printf("t_check = %1.1g, rate = %g, B_accept = %.2g\n", (double)t_check, samples[0].rate, B_accept);
    printf("t_max = %1.1g, p_max = %g, pc_fr_max = %g\n", (double)t_max, p_max, pc_fr_max);
    printf("A = %1.1g, B = %g\n", samples[0].A, samples[0].B);

    printf("#############################\n\n");

    clock_t start_time, end_time;
    double cpu_time_used;

    start_time = clock(); // Record the start time


    printf("\n");
    testz(0,samples);
    printf("\n\n");

    V_hard_sphere(0, samples);
    printf("\n\n");


    FILE *fout0;
    sprintf(file_name,"./output/p_pacF_%.2d.dat",label);
    fout0 = fopen(file_name, "w");

//     fprintf(fout0,"%g %g %g\n", samples[0].V_hs/samples[0].V, bac_double, samples[0].V);
    fprintf(fout0,"%g %g %g\n", samples[0].V_hs/samples[0].V, 1./bac_double, samples[0].V);

    ovito_plot(0, samples, label, ii);


    samples[0].p_measured = 0.;
    ii=1;
    while(samples[0].p_measured < p_max && samples[0].V_hs/samples[0].V < pc_fr_max && ii < t_max){

        MC(0, samples, seeds);
//         fprintf(fout0,"%g %g %g\n", samples[0].V_hs/samples[0].V_bac, samples[0].p_measured, samples[0].V_bac);
        fprintf(fout0,"%g %g %g\n", samples[0].V_hs/samples[0].V_bac, 1./samples[0].p_measured, samples[0].V_bac);

        printScrollingValue(samples[0].p_measured, p_max, samples[0].V_hs/samples[0].V, pc_fr_max, (double)ii, (double)t_max, samples[0].B, samples[0].accept_V, samples[0].P);

        ii+=1;

        ovito_plot(0, samples, label, ii);
    }
    fclose(fout0);

    end_time = clock(); // Record the end time
    // Calculate the CPU time used in seconds
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\n\n");
    testz(0,samples);
    printf("\n\n");

    V_hard_sphere(0, samples);
    printf("\n\n");



//######################################################################################
//plot//################################################################################

    printf("\n");

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
    sprintf(file_name,"./output/par_%.2d.dat",label);
    fout3 = fopen(file_name, "w");

    fprintf(fout3,"%d\n", N);
    fprintf(fout3,"%g\n", rate_bac);
    fprintf(fout3,"%d\n", t_check);
    fprintf(fout3,"%g\n", p_max);
    fprintf(fout3,"%g\n", pc_fr_max);
    fprintf(fout3,"%g\n", B_accept);
    fprintf(fout3,"%g\n", P0_label);
    fprintf(fout3,"%g\n", cpu_time_used);

    fclose(fout3);


    printf("amp\t\tP\t\tV\t\tA\t\tacX\t\tB\t\tacV\t\tid\t\taccept\n");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        printf("%.4g\t\t%.2g\t\t%.4g\t\t%.2g\t\t%.2g\t\t%.2g\t\t%.2g\t\t%d\t\t%g\n", samples[jj].amp_fac, samples[ii].P, samples[jj].V, samples[jj].A, samples[jj].accept_mov, samples[jj].B, samples[jj].accept_V, samples[ii].index_list, samples[ii].accept_num_PT_bac);
    }
    printf("\n");

    for (ii = 0; ii < num_threads; ii++)
        ovito_plot(ii, samples, label, 0);

    return 0;
}

