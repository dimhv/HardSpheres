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

    printf("sample: %.2d,\tPF = %.3g,\tamp_fac = %.3g, \tP(MC) = %.3g\n", thread_id, samples[II].V_hs/(L2*L/pow(samples[II].amp_fac,3)), samples[II].amp_fac, samples[II].P);
//     printf("sample: %d,\t\tPF = %g,\t\tamp_fac = %g\n", thread_id, samples[II].V_hs/(L2/pow(samples[II].amp_fac,2)), samples[II].amp_fac);

}
