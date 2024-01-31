void optm_T(Sample* samples){

    int jj;
    double bac_double, delta_beta;

    for (jj = 0; jj < num_threads-1; ++jj){
        samples[jj].accept_num_PT_bac = samples[jj].accept_num_PT/( (double)(t_PT*0.5) );
        samples[jj].accept_av_PT += samples[jj].accept_num_PT/( (double)(t_PT*0.5) );
    }

//     bac_double = 0.0;
//     for (jj = 0; jj < num_threads-1; ++jj)
//         bac_double += samples[jj].accept_num_PT_bac;
//
//     bac_double /= (double)(num_threads-1);
//
//     for (jj = 0; jj < num_threads; ++jj)
//         samples[jj].P_prime = samples[jj].P;
//
//
//     for (jj = 1; jj < num_threads; ++jj){
//
//         delta_beta = samples[jj].P_prime - samples[jj-1].P_prime;
//         samples[jj].P = samples[jj-1].P + delta_beta*(samples[jj-1].accept_num_PT_bac/bac_double);
//     }

    for (jj = 0; jj < num_threads; ++jj)
        samples[jj].accept_num_PT = 0.;

    return;
}
