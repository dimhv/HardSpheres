void ovito_plot(int thread_id, Sample* samples, int label, int tt){

    int ii;
    int II = samples[thread_id].index_list;

    char file_name[50];

//     FILE *fout_conf;
//     sprintf(file_name,"./output/conf/c_%.3d_%.3d.xyz", tt,label);
//     fout_conf = fopen(file_name, "w");
//
//
//     fprintf( fout_conf,"ITEM: TIMESTEP\n");
//     fprintf( fout_conf,"%d\n", tt);
//
//     fprintf( fout_conf,"ITEM: NUMBER OF ATOMS\n");
//     fprintf( fout_conf,"%d\n", N);
//
//     fprintf( fout_conf,"ITEM: BOX BOUNDS pp pp pp\n");
//     fprintf( fout_conf,"%d %g\n", 0,(double)L/samples[II].amp_fac);
//     fprintf( fout_conf,"%d %g\n", 0,(double)L/samples[II].amp_fac);
//     fprintf( fout_conf,"%d %g\n", 0,(double)L/samples[II].amp_fac);
//
//     fprintf(fout_conf, "ITEM: ATOMS id type x y z radius Transparency\n");
//
//     for(ii=0; ii < N; ++ii){
//          fprintf( fout_conf,"%d %d ",ii, 1);
//
//         fprintf( fout_conf,"%g ", (samples[II].particles[ii].x )/samples[II].amp_fac);
//         fprintf( fout_conf,"%g ", (samples[II].particles[ii].y )/samples[II].amp_fac);
//         fprintf( fout_conf,"%g ", (samples[II].particles[ii].z )/samples[II].amp_fac);
//
//         fprintf( fout_conf,"%g %g\n", 0.5*samples[II].particles[ii]._2r, 0.2);
//     }

//     for(ii=0; ii < N; ++ii){
//          fprintf( fout_conf,"%d %d ",ii, 1);
//
//         fprintf( fout_conf,"%.16lf ", (samples[II].particles[ii].x + (double)L*samples[II].particles[ii].x_pbc)/samples[II].amp_fac);
//         fprintf( fout_conf,"%.16lf ", (samples[II].particles[ii].y + (double)L*samples[II].particles[ii].y_pbc)/samples[II].amp_fac);
//         fprintf( fout_conf,"%.16lf ", (samples[II].particles[ii].z + (double)L*samples[II].particles[ii].z_pbc)/samples[II].amp_fac);
//
//         fprintf( fout_conf,"%.16lf %g\n", 0.5*samples[II].particles[ii]._2r, 0.2);
//     }

//     fclose(fout_conf);

    FILE *fout_conf2;
    sprintf(file_name,"./output/output_%.2d/conf/c_%.3d_%.3d.dat", input_index, tt,label);
//     sprintf(file_name,"./output/conf/c_%.3d_%.3d.dat", tt,label);
    fout_conf2 = fopen(file_name, "w");

    for(ii=0; ii < N; ++ii)
        fprintf( fout_conf2,"%.16lf %.16lf %.16lf %.16lf %d %d %d\n", samples[II].particles[ii].x, samples[II].particles[ii].y, samples[II].particles[ii].z, samples[II].particles[ii]._2r, samples[II].particles[ii].x_pbc, samples[II].particles[ii].y_pbc, samples[II].particles[ii].z_pbc);


    fclose(fout_conf2);

    FILE *fout_conf3;
    sprintf(file_name,"./output/output_%.2d/conf/par_%.3d_%.3d.dat", input_index, tt,label);
//     sprintf(file_name,"./output/conf/par_%.3d_%.3d.dat", tt,label);
    fout_conf3 = fopen(file_name, "w");

    fprintf( fout_conf3,"%.16lf\n%g\n%g\n%.16lf", samples[II].amp_fac, samples[II].A, samples[II].B, 1./samples[II].p_measured);

    fclose(fout_conf3);


    return;
}
