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

double mesh_fac = 2.5; //open for adjustment
int t_PT = 400;

static int t_max = (int)1e8;
static int num_threads = 1;


int len_pf_target = 64;
double pf_target[64];

int nhis = 8000;
double bins[8000];

double p_max = 50000;

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

void printScrollingValue_2(double value1, double targetValue1, double value2, double targetValue2, double value4, double value5, double value6) {
    printf("\rpc_fr:%.5g/%g##steps:%.2g/%.2g##accept_V:%.1e/%.1e##P:%.2g  ", value1, targetValue1, value2, targetValue2, value4, value6, value5);
    fflush(stdout);
}

double logpace(double start, double end, int numPoints, int Point) {

    double base = exp(log(end / start) / (numPoints - 1));

    if(start == end)
        return end;

    return start * pow(base, Point);
}

void p_evaluation(int thread_id, Sample* samples, Seed* seeds){

    int ii, jj, kk;
    kk = 500;


    for(jj=0; jj < nhis; ++jj)
        bins[jj] = 0.;

    for(ii=0; ii< kk; ++ii){

        (void)legal_configuration_dyn(thread_id, samples, seeds);
        g_r_dist_count(thread_id, samples);
    }

    samples[thread_id].p_measured = g_r_dist_return(thread_id, samples, (double)kk);

    return;
}

int main(int argc, char** argv) {

    int label;
    sscanf (argv[1],"%d",&label);

    char file_name[50];
    FILE *fin0;
    sprintf(file_name,"par.dat");
    fin0 = fopen(file_name, "r");

    double rate_bac;

    fscanf(fin0,"%d\n", &N);
    fscanf(fin0,"%lf\n", &rate_bac);
    fscanf(fin0,"%d\n", &t_check);

    double pc_fr_min, pc_fr_max;
    fscanf(fin0,"%lf\n", &pc_fr_min);
    fscanf(fin0,"%lf\n", &pc_fr_max);
    fscanf(fin0,"%lf\n", &acept_Vbac);
//     fscanf(fin0,"%d\n", &input_index);

    sscanf (argv[2],"%d",&input_index);


    fclose(fin0);

    L_mesh = L/mesh_fac;
    L2_mesh = L_mesh*L_mesh;
    cell_num = L2_mesh*L_mesh;

//######################################################################################
//loop counters, pressure target inverval//#############################################
    int ii,jj,kk;
    printf("\n\ntrlx\n");


    for (ii = 0; ii < len_pf_target; ii++){
        pf_target[ii] = logpace(pc_fr_min, pc_fr_max, len_pf_target, ii);
        printf("pf_target, %d: %g\n", ii, pf_target[ii]);
    }
//######################################################################################
//structures initialization//###########################################################

    Sample* samples = (Sample*)malloc(num_threads * sizeof(Sample));
    Seed* seeds = (Seed*)malloc( (num_threads+1) * sizeof(Seed));

    FILE *fin[num_threads];
    FILE *fin2[num_threads];
    double A0, B0, ampF_0;

//     kk = 100;
    kk = 0;

    for (ii = 0; ii < num_threads; ii++) {

        sprintf(file_name,"./input/par_%.4d_%.4d.dat", input_index,kk);

        fin[ii] = fopen(file_name, "r");
        fscanf(fin[ii],"%lf\n", &ampF_0);
        fscanf(fin[ii],"%lf\n", &A0);
        fscanf(fin[ii],"%lf\n", &B0);
        fclose(fin[ii]);

        samples[ii].particles = (Particle*)malloc(N * sizeof(Particle));
        samples[ii].time = (Time*)malloc(tobs * sizeof(Time));

        samples[ii].A = A0;
        samples[ii].B = B0;

        samples[ii].amp_fac = ampF_0;

        samples[ii].V = L2*L/ampF_0/ampF_0/ampF_0;
//         samples[ii].V = L2/ampF_0/ampF_0;

        samples[ii].rate = rate_bac;

        samples[ii].cell_ind = (int*)malloc(cell_num * sizeof(int));
        samples[ii].cell = (int**)malloc(cell_num * sizeof(int*));
        for (int jj = 0; jj < cell_num; jj++) {
            samples[ii].cell[jj] = (int*)malloc(N * sizeof(int));
        }

        sprintf(file_name,"./input/c_%.4d_%.4d.dat", input_index,kk);

        fin2[ii] = fopen(file_name, "r");
        for (jj = 0; jj < N; jj++)
            fscanf(fin2[ii],"%lf %lf %lf %lf\n", &samples[ii].particles[jj].x0, &samples[ii].particles[jj].y0, &samples[ii].particles[jj].z0, &samples[ii].particles[jj]._2r);
        fclose(fin2[ii]);

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

            samples[ii].particles[jj].x0 = samples[ii].particles[jj].x;
            samples[ii].particles[jj].y0 = samples[ii].particles[jj].y;
            samples[ii].particles[jj].z0 = samples[ii].particles[jj].z;
        }
    }
//seting pressure distribution//########################################################


    p_evaluation(0, samples, seeds);
    samples[0].P = samples[0].p_measured*(double)N/samples[0].V;
//     samples[0].P = 40;


//######################################################################################
//#Dynamics//###########################################################################


    printf("\n#############################\n");
    printf("label = %d, N = %d, L = %d, p0(sample) = %g\n\n", label, N, L, samples[0].p_measured);
    printf("V[0] = %g, mesh_fac = %g, cell_num = %d\n", samples[0].V, mesh_fac,cell_num);

    printf("num_threads = %d, seed = %d\n\n", num_threads, seed_bac);

    printf("t_check = %1.1g, rate = %g\n", (double)t_check, samples[0].rate);
    printf("t_max = %1.1g, p_max = %g, pc_fr_max = %g\n", (double)t_max, p_max, pc_fr_max);
    printf("A = %.2g, B = %.2g, acept_Vbac = %.2g\n", samples[0].A, samples[0].B, acept_Vbac);
    printf("#############################\n\n");


    for (ii = 0; ii < num_threads; ii++)
        V_hard_sphere(ii, samples);

    FILE *fout0;
    sprintf(file_name,"./output/p_pacF_%.2d.dat", label);

    fout0 = fopen(file_name, "w");


    printf("\n\n");
    testz(0,samples);

    printf("\n\ntobs\n");


    ii=0;
    jj=0;
    while(samples[0].p_measured < p_max && samples[0].V_hs/samples[0].V < pc_fr_max && ii < t_max){

        if(ii%10==0)
            printScrollingValue_2(samples[0].V_hs/samples[0].V, pc_fr_max, (double)ii, (double)t_max, samples[0].accept_V, samples[0].P, samples[0].B);

        MC(0, samples, seeds);


        ii+=1;

        if( (samples[0].V_hs/samples[0].V - pf_target[jj]) > 0){

            p_evaluation(0, samples, seeds);

            printf("\npf_plt = %g, prss = %g\n", samples[0].V_hs/samples[0].V, samples[0].p_measured);
            fprintf(fout0,"%g %g %g %g %g\n", samples[0].V_hs/samples[0].V, 1./samples[0].p_measured, samples[0].V, samples[0].p_measured*(double)N/samples[0].V, samples[0].P);
//             fprintf(fout0,"%d %g %g %g %g %g\n", ii, samples[0].V_hs/samples[0].V, samples[0].V, samples[0].p_measured, samples[0].p_measured*(double)N/samples[0].V, samples[0].P);

            ovito_plot(0, samples, label, jj);
            jj+=1;
        }

        if(jj==len_pf_target)
            break;
    }
    fclose(fout0);

    printf("\n\n");
    testz(0,samples);
    printf("\n\n");

    V_hard_sphere(0, samples);
    printf("\n\n");
//######################################################################################
//plot//################################################################################

    printf("\n");

    FILE *fout1;
    sprintf(file_name,"./output/output_%.2d/r_%.2d.dat",input_index, label);
//     sprintf(file_name,"./output/r_%.2d.dat", label);

    fout1 = fopen(file_name, "w");

    for (ii = 0; ii < N; ii++)
        fprintf(fout1,"%g\n", samples[0].particles[ii]._2r);

    fclose(fout1);

    FILE *fout2;
    sprintf(file_name,"./output/output_%.2d/report_%.2d.dat",input_index,label);
//     sprintf(file_name,"./output/report_%.2d.dat", label);
    fout2 = fopen(file_name, "w");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        fprintf(fout2,"%g %g %g %g %g %d %g\n", samples[jj].amp_fac, samples[ii].P, samples[jj].V, samples[jj].A, samples[jj].B, samples[ii].index_list, samples[ii].accept_num_PT_bac);
    }

    fclose(fout2);

    printf("amp\tP\tV\tA\tacX\tB\tacV\tid\taccept\n");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        printf("%.4g\t%.2g\t%.4g\t%.3g\t%.2g\t%.2g\t%.2g\t%d\t%g\n", samples[jj].amp_fac, samples[ii].P, samples[jj].V, samples[jj].A, samples[jj].accept_mov, samples[jj].B, samples[jj].accept_V, samples[ii].index_list, samples[ii].accept_num_PT_bac);
    }
    printf("\n");

//     for (ii = 0; ii < num_threads; ii++)
//         ovito_plot(ii, samples, label, 1);

    return 0;
}

