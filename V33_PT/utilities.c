void logspace(double start, double end, int numPoints, Sample* samples) {
    double base = exp(log(end / start) / (numPoints - 1));

    for (int ii = 0; ii < numPoints; ii++)
        samples[ii].P = start * pow(base, ii);
}


void linspace(double start, double end, int numPoints, Sample* samples) {
    double step = (end - start) / (numPoints - 1);

    for (int ii = 0; ii < numPoints; ii++) {
        samples[ii].P = start + ii * step;
    }
}

void V_hard_sphere(int thread_id, Sample* samples){

    int ii;
    int II = samples[thread_id].index_list;

    samples[II].V_hs = 0.;
    for(ii=0; ii < N; ++ii)
        samples[II].V_hs += pow(samples[II].particles[ii]._2r,3);
    samples[II].V_hs *= (4./3.)*M_PI/8.;

//     for(ii=0; ii < N; ++ii)
//         samples[II].V_hs += pow(samples[II].particles[ii]._2r,2);
//     samples[II].V_hs *= M_PI/4.;

    printf("sample: %.2d,\tPF = %.3g,\tamp_fac = %.3g \t\tP = %.3g\n", thread_id, samples[II].V_hs/(L2*L/pow(samples[II].amp_fac,3)), samples[II].amp_fac, samples[II].P);
//     printf("sample: %d,\t\tPF = %g,\t\tamp_fac = %g\n", thread_id, samples[II].V_hs/(L2/pow(samples[II].amp_fac,2)), samples[II].amp_fac);

}

void com(int kk, Sample* samples){

    int jj;


    samples[kk].x_com = samples[kk].y_com = samples[kk].z_com = 0.;
    for (jj = 0; jj < N; jj++){
        samples[kk].x_com += samples[kk].particles[jj].x + (double)(samples[kk].particles[jj].x_pbc*L);
        samples[kk].y_com += samples[kk].particles[jj].y + (double)(samples[kk].particles[jj].y_pbc*L);
        samples[kk].z_com += samples[kk].particles[jj].z + (double)(samples[kk].particles[jj].z_pbc*L);
    }

    samples[kk].x_com /= (double)N;
    samples[kk].y_com /= (double)N;
    samples[kk].z_com /= (double)N;

    samples[kk].COM = samples[kk].x_com*samples[kk].x_com + samples[kk].y_com*samples[kk].y_com + samples[kk].z_com*samples[kk].z_com;
    samples[kk].COM = sqrt(samples[kk].COM);

    return;
}

void com0(int kk, Sample* samples){

    int jj;

    samples[kk].x0_com = samples[kk].y0_com = samples[kk].z0_com = 0.;
//     for (jj = 0; jj < N; jj++){
//         samples[kk].x0_com += samples[kk].particles[jj].x0;
//         samples[kk].y0_com += samples[kk].particles[jj].y0;
//         samples[kk].z0_com += samples[kk].particles[jj].z0;
//     }

    for (jj = 0; jj < N; jj++){
        samples[kk].x0_com += samples[kk].particles[jj].x0 + (double)(samples[kk].particles[jj].x0_pbc*L);
        samples[kk].y0_com += samples[kk].particles[jj].y0 + (double)(samples[kk].particles[jj].y0_pbc*L);
        samples[kk].z0_com += samples[kk].particles[jj].z0 + (double)(samples[kk].particles[jj].z0_pbc*L);
    }

    samples[kk].x0_com /= (double)N;
    samples[kk].y0_com /= (double)N;
    samples[kk].z0_com /= (double)N;

    samples[kk].COM0 = samples[kk].x0_com*samples[kk].x0_com + samples[kk].y0_com*samples[kk].y0_com + samples[kk].z0_com*samples[kk].z0_com;
    samples[kk].COM0 = sqrt(samples[kk].COM0);

    return;
}

double msd(int thread_id, Sample* samples){

    double dx, dy, dz;
    int jj, kk;

    kk = samples[thread_id].index_list;
    double msd_bac = 0.;

    com(kk, samples);

    for (jj = 0; jj < N; jj++){

        dx = samples[kk].particles[jj].x + (double)(samples[kk].particles[jj].x_pbc*L) - ( samples[kk].particles[jj].x0 + (double)(samples[kk].particles[jj].x0_pbc*L) );
        dy = samples[kk].particles[jj].y + (double)(samples[kk].particles[jj].y_pbc*L) - ( samples[kk].particles[jj].y0 + (double)(samples[kk].particles[jj].y0_pbc*L) );
        dz = samples[kk].particles[jj].z + (double)(samples[kk].particles[jj].z_pbc*L) - ( samples[kk].particles[jj].z0 + (double)(samples[kk].particles[jj].z0_pbc*L) );

        dx -= samples[kk].x_com - samples[kk].x0_com;
        dy -= samples[kk].y_com - samples[kk].y0_com;
        dz -= samples[kk].z_com - samples[kk].z0_com;

        msd_bac += dx*dx + dy*dy + dz*dz;
    }

    msd_bac /= ( (double)N*samples[kk].amp_fac*samples[kk].amp_fac );
    return msd_bac;
}

double msd_AB(int A_ind, int B_ind, Sample* samples){

    double dx, dy, dz;
    int jj;

    int kk = samples[A_ind].index_list;
    int ll = samples[B_ind].index_list;

    double msd_bac = 0.;

    com(kk, samples);
    com(ll, samples);

    for (jj = 0; jj < N; jj++){

        dx = samples[kk].particles[jj].x + (double)(samples[kk].particles[jj].x_pbc*L) - ( samples[ll].particles[jj].x + (double)(samples[ll].particles[jj].x_pbc*L) );
        dy = samples[kk].particles[jj].y + (double)(samples[kk].particles[jj].y_pbc*L) - ( samples[ll].particles[jj].y + (double)(samples[ll].particles[jj].y_pbc*L) );
        dz = samples[kk].particles[jj].z + (double)(samples[kk].particles[jj].z_pbc*L) - ( samples[ll].particles[jj].z + (double)(samples[ll].particles[jj].z_pbc*L) );

        dx -= samples[kk].x_com - samples[ll].x_com;
        dy -= samples[kk].y_com - samples[ll].y_com;
        dz -= samples[kk].z_com - samples[ll].z_com;

        msd_bac += dx*dx + dy*dy + dz*dz;
    }

    msd_bac /= ( (double)N*samples[kk].amp_fac*samples[kk].amp_fac );
    return msd_bac;
}

void g_r_dist_count(int thread_id, Sample* samples){

    int ii, jj, ig;
    double delta_x, delta_y, delta_z, d, D;

    double delg = (double)L/(2.*nhis);

    int II = samples[thread_id].index_list;

    for(ii=0; ii < N; ++ii){


        for(jj=ii+1; jj < N; ++jj){

            delta_x = fabs(samples[II].particles[ii].x - samples[II].particles[jj].x);
            if( delta_x > 0.5*(double)L )
                delta_x = (double)L - delta_x;


            delta_y = fabs(samples[II].particles[ii].y - samples[II].particles[jj].y);
            if( delta_y > 0.5*(double)L )
                delta_y = (double)L - delta_y;

            delta_z = fabs(samples[II].particles[ii].z - samples[II].particles[jj].z);
            if( delta_z > 0.5*(double)L )
                delta_z = (double)L - delta_z;

            d = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;
            d = sqrt(d);

            D = samples[II].amp_fac*0.5*(samples[II].particles[ii]._2r + samples[II].particles[jj]._2r);
            d /= D;

            ig = (int)( d/delg );

            if(ig<nhis)
                samples[II].bins[ig] += 2.;
        }

    }

    return;
}

double g_r_dist_return(int thread_id, Sample* samples, double norm_fact){

    int ii;
    double d, nid, vb;
    double delg = (double)L/(2.*nhis);

    double bac_bins;

    int II = samples[thread_id].index_list;
    double rho = (double)N/( samples[II].V );


    for(ii=0; ii < nhis; ++ii){

        d = delg * (ii + 0.5);
        vb = ( (ii + 1) * (ii + 1) * (ii + 1) - ii * ii * ii ) *delg * delg * delg;
        nid =  (4./3.)*M_PI*vb*rho;

        samples[II].bins[ii] /= ((double)N * nid * norm_fact);

        if(d > 1){
            bac_bins = 1. + 0.5*(4./3.)*M_PI*rho*samples[II].bins[ii];
            return bac_bins;
        }
    }

    return 0.0;
}
