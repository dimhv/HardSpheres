void testz(int thread_id, Sample* samples){

    double d, delta_x, delta_y, delta_z, D;
    int ii,jj;

    int II = samples[thread_id].index_list;


    for(ii=0; ii < N; ++ii)
        for(jj=ii+1; jj < N; ++jj){

            delta_x = samples[II].particles[ii].x - samples[II].particles[jj].x;
            delta_y = samples[II].particles[ii].y - samples[II].particles[jj].y;
            delta_z = samples[II].particles[ii].z - samples[II].particles[jj].z;

            delta_x = fabs(delta_x);
            delta_y = fabs(delta_y);
            delta_z = fabs(delta_z);

            delta_x = fmin( (double)L-delta_x, delta_x );
            delta_y = fmin( (double)L-delta_y, delta_y );
            delta_z = fmin( (double)L-delta_z, delta_z );

            d = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;
            d = sqrt(d);

            D = 0.5*samples[II].amp_fac*(samples[II].particles[ii]._2r + samples[II].particles[jj]._2r);

            if(d < D){
                printf("Overlap Warning: d=%g, D = %g, ii = %d jj = %d sample = %d\n", d, D, ii, jj, II);
                return;

            }
        }

    printf("No Overlap detected, sample: %d\n", II);
    return;
}
