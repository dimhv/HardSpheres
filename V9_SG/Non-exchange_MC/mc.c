//######################################################################################
//PT

void optm_T(Sample* samples){

    int jj;
    double bac_double, delta_beta;

    for (jj = 0; jj < num_threads-1; ++jj)
        samples[jj].accept_num_PT_bac = samples[jj].accept_num_PT/( (double)(t_PT*0.5) );

    bac_double = 0.0;
    for (jj = 0; jj < num_threads-1; ++jj)
        bac_double += samples[jj].accept_num_PT_bac;

    bac_double /= (double)(num_threads-1);

    for (jj = 0; jj < num_threads; ++jj)
        samples[jj].beta_prime = samples[jj].beta;


    for (jj = 1; jj < num_threads; ++jj){

        delta_beta = samples[jj].beta_prime - samples[jj-1].beta_prime;
        samples[jj].beta = samples[jj-1].beta + delta_beta*(samples[jj-1].accept_num_PT_bac/bac_double);
    }

    for (jj = 0; jj < num_threads; ++jj)
        samples[jj].accept_num_PT = 0.;

    return;
}

void optm_T_mod(Sample* samples){

    int jj;

    for (jj = 0; jj < num_threads-1; ++jj)
        samples[jj].accept_num_PT_bac = samples[jj].accept_num_PT/( (double)(t_PT*0.5) );

//     bac_double = 0.0;
//     for (jj = 0; jj < num_threads-1; ++jj)
//         bac_double += samples[jj].accept_num_PT_bac;
//
//     bac_double /= (double)(num_threads-1);
//
//     for (jj = 0; jj < num_threads; ++jj)
//         samples[jj].beta_prime = samples[jj].beta;
//
//
//     for (jj = 1; jj < num_threads; ++jj){
//
//         delta_beta = samples[jj].beta_prime - samples[jj-1].beta_prime;
//         samples[jj].beta = samples[jj-1].beta + delta_beta*(samples[jj-1].accept_num_PT_bac/bac_double);
//     }

    for (jj = 0; jj < num_threads; ++jj)
        samples[jj].accept_num_PT = 0.;

    return;
}

void MC_PT(int num_threads, Sample* samples, int time_bac, Seed* seeds){

    int ii, II, II0;
    double delta;

    if(time_bac%2==0)
        for(ii=1; ii < num_threads; ii+=2){

            II = samples[ii].index_list;
            II0 = samples[ii-1].index_list;

            delta = (samples[ii].beta - samples[ii-1].beta)*(double)(samples[II0].e - samples[II].e);
            if( delta <= 0. || xor64(num_threads, seeds) < exp(-delta) ){

                samples[ii].index_list = II0;
                samples[ii-1].index_list = II;

                samples[ii-1].accept_num_PT += 1.;
            }
        }

    else
        for(ii=2; ii < num_threads; ii+=2){

            II = samples[ii].index_list;
            II0 = samples[ii-1].index_list;

            delta = (samples[ii].beta - samples[ii-1].beta)*(double)(samples[II0].e - samples[II].e);
            if( delta <= 0. || xor64(num_threads, seeds) < exp(-delta) ){

                samples[ii].index_list = II0;
                samples[ii-1].index_list = II;

                samples[ii-1].accept_num_PT += 1.;
            }
        }

    return;
}

//######################################################################################
//neighbour list
int nl(int ii, int jj){

    int bac;

    switch (jj){

        case 0:
            //right
            bac = ii + 1;
            if( bac%L == 0 )
                bac -= L;
            return bac;

        case 1:
            //left
            bac = ii - 1;
            if( ii%L == 0 )
                bac += L;
            return bac;

        case 2:
            //up
            bac = ii + L;
            if( ii%L2 >= L*(L-1) )
                bac -= L2;
            return bac;

        case 3:
            //down
            bac = ii - L;
            if( ii%L2 < L )
                bac += L2;
            return bac;

        case 4:
            //up2
            bac = ii + L2;
            if(ii%N >= L2*(L-1))
                bac -= N;
            return bac;

        case 5:
            //down2
            bac = ii - L2;
            if(ii%N < L2)
                bac += N;
            return bac;
    }
}




//######################################################################################
// mc


void MC(int thread_id, Sample* samples, Seed* seeds){

    int ii, jj, ll, kk,delta_e;
    int II = samples[thread_id].index_list;

    for(ii=0; ii < t_check; ii++){

        for(ll=0; ll < N; ll++){

            jj = (int)(N*xor64(thread_id, seeds));
            samples[II].particles[jj].x *= -1;

            delta_e = 0;
            for(kk = 0; kk < num_neighbours; ++kk)
                delta_e -= 2*samples[0].particles[jj].J[kk]*samples[II].particles[ nl(jj,kk) ].x;

            delta_e *= samples[II].particles[ jj ].x;

            if( delta_e <= 0. || xor64(thread_id, seeds) < exp( -samples[thread_id].beta*(double)delta_e) )
                samples[II].e += delta_e;

            else
                samples[II].particles[jj].x *= -1;
        }
    }
    return;

}

void MC_rlx(int thread_id, Sample* samples, Seed* seeds){

    int ii, jj, ll, kk,delta_e;
    int II = samples[thread_id].index_list;

    for(ii=0; ii < t_rlx; ii++){

        for(ll=0; ll < N; ll++){

            jj = (int)(N*xor64(thread_id, seeds));
            samples[II].particles[jj].x *= -1;

            delta_e = 0;
            for(kk = 0; kk < num_neighbours; ++kk)
                delta_e -= 2*samples[0].particles[jj].J[kk]*samples[II].particles[ nl(jj,kk) ].x;

            delta_e *= samples[II].particles[ jj ].x;

            if( delta_e <= 0. || xor64(thread_id, seeds) < exp( -samples[thread_id].beta*(double)delta_e) )
                samples[II].e += delta_e;

            else
                samples[II].particles[jj].x *= -1;
        }
    }
    return;

}


//######################################################################################
// measurements

double mag(int thread_id, Sample* samples){

    int ii;
    int II = samples[thread_id].index_list;


    double bac_double = 0;
    for(ii=0; ii < N; ++ii)
        bac_double += samples[II].particles[ii].x;

    return fabs(bac_double/(double)N);
}

double q_tcf(int thread_id, Sample* samples){

    int ii;

    int II = samples[thread_id].index_list;


    double bac_double = 0;
    for(ii=0; ii < N; ++ii)
        bac_double += samples[II].particles[ii].x0*samples[II].particles[ii].x;

    return bac_double/(double)N;
}

void energy(int thread_id, Sample* samples){

    int ii,jj, bac;

    int II = samples[thread_id].index_list;
    samples[II].e = 0;

    for (ii = 0; ii < N; ii++){

        bac = 0;

        for (jj = 0; jj < num_neighbours; jj+=2)
            bac -= samples[0].particles[ii].J[jj]*samples[II].particles[ nl(ii,jj) ].x;

        bac *= samples[II].particles[ii].x;
        samples[II].e += bac;
    }

    return;
}

double q_tcf_mod(int index_A, int index_B, Sample* samples){

    int ii;

    double bac_double = 0;
    for(ii=0; ii < N; ++ii)
        bac_double += samples[index_A].particles[ii].x*samples[index_B].particles[ii].x;

    return bac_double/(double)N;
}

double chi(double* q, Sample* samples){

    int ii, jj, kk;

    kk = 0;
    for(ii=0; ii < num_threads; ii++)
        for(jj=ii+1; jj < num_threads; jj++){

            q[kk] = q_tcf_mod(ii, jj, samples);
            kk+=1;

        }


    double msd = 0;

    for(ii=0; ii < fac; ii++)
        msd += q[ii]*q[ii];

    return (double)N*msd/(double)fac;

}

