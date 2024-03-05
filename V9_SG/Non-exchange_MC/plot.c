void plot_conf(int time, int thread_id, Sample* samples){


    int II = samples[thread_id].index_list;
    int ii;

    FILE *fout;
    sprintf(file_name,"./output/conf_%.3d_%.4d.dat",label, time);
    fout = fopen(file_name, "w");

    for(ii=0; ii< N; ++ii)
        fprintf(fout,"%d\n", samples[II].particles[ii].x);

    fclose(fout);

    return;
}


void plot_NO_EXCHANGE_chi(Sample* samples, double* chi_plot){

    FILE *fout;
    sprintf(file_name,"./output/chi_%.2d.dat",label);
    fout = fopen(file_name, "w");

    int ii;

    for (ii = 1; ii < tobs; ii++)
        fprintf(fout,"%d %g\n", ii*t_check, chi_plot[ii]);

    fclose(fout);

    return;
}


void plot_NO_EXCHANGE_av_q_tcf(Sample* samples){

    int ii, kk;
    double av;

    FILE *fout;
    sprintf(file_name,"./output/q_tcf_%.2d.dat",label);
    fout = fopen(file_name, "w");

    for (ii = 0; ii < tobs; ii++){
        av = 0;
        for (kk = 0; kk < num_threads; kk++)
            av += samples[ kk ].time[ii].q_plot;
        fprintf(fout,"%d %g\n", ii*t_check, av/(double)num_threads);
    }
    fclose(fout);

    return;
}


void plot_PT_q_tcf(Sample* samples){

    int ii, kk;

    FILE *fout;
    sprintf(file_name,"./output/q_tcf_%.2d.dat",label);
    fout = fopen(file_name, "w");
    for (ii = 0; ii < num_threads; ii++){

        for (kk = 0; kk < tobs; kk++)
            fprintf(fout,"%d %g\n", kk*t_check, samples[ ii ].time[kk].q_plot);

//         fprintf(fout,"\n");
    }
    fclose(fout);

    return;
}


void plot_report(Sample* samples){

    int ii, jj;

    FILE *fout;
    sprintf(file_name,"./output/report_%.2d.dat",label);
    fout = fopen(file_name, "w");

    for (ii = 0; ii < num_threads; ii++){
        jj = samples[ii].index_list;
        fprintf(fout,"%.2lf %d %d %.2g %d %d\n", samples[ii].beta, samples[jj].e, samples[ii].index_list, samples[ii].accept_num_PT_bac, tobs, t_rlx);
    }

    fclose(fout);

    return;
}

void plot_NO_EXCHANGE_EA(Sample* samples){

    int ii, jj;

    FILE *fout;
    sprintf(file_name,"./output/EA_dist_%.2d.dat",label);
    fout = fopen(file_name, "w");
    for (ii = 0; ii < num_threads; ii++)
        for (jj = ii+1; jj < num_threads; jj++)
            fprintf(fout,"%g\n", q_tcf_mod(ii,jj,samples));
    fclose(fout);

    return;
}


