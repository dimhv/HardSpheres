void neighbor_cell_list(int label_bac, int index, int list[]){

    int bac, sx, sy;

    switch (label_bac%9){

        case 0:
            list[0] = index;
            list[1] = 0;
            list[2] = 0;
            list[3] = 0;

            if( (int)(label_bac/9 )==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }
            break;

        case 1:
            bac = index - 1; //left
            sx = sy = 0;
            if(index%L_mesh==0){
                bac += L_mesh;
                sx = L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 2:
            bac = index + 1; //right
            sx = sy = 0;
            if( (index+1)%L_mesh==0 ){
                bac -= L_mesh;
                sx = -L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 3:
            bac = index - L_mesh; //bottom
            sx = sy = 0;
            if(index%L2_mesh<L_mesh){
                bac += L2_mesh;
                sy = L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 4:
            bac = index + L_mesh; //top
            sx = sy = 0;
            if( index%L2_mesh >= L_mesh*(L_mesh-1)){
                bac -= L2_mesh;
                sy = -L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 5:
            bac = index + L_mesh - 1; //top/left
            sx = sy = 0;

            if(index%L2_mesh>=L_mesh*(L_mesh-1)){
                bac -= L2_mesh;
                sy = -L;
            }

            if(index%L_mesh==0){
                bac += L_mesh;
                sx = L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }


            break;

        case 6:
            bac = index + L_mesh + 1; //top/right
            sx = sy = 0;

            if(index%L2_mesh>=L_mesh*(L_mesh-1)){
                bac -= L2_mesh;
                sy = -L;
            }

            if( (index+1)%L_mesh==0 ){
                bac -= L_mesh;
                sx = -L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 7:
            bac = index - L_mesh - 1; //bottom/left
            sx = sy = 0;

            if(index%L2_mesh<L_mesh){
                bac += L2_mesh;
                sy = L;
            }

            if(index%L_mesh==0){
                bac += L_mesh;
                sx = L;
            }


            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 8:
            bac = index - L_mesh + 1; //bottom/rigth
            sx = sy = 0;


            if(index%L2_mesh<L_mesh){
                bac += L2_mesh;
                sy = L;
            }

            if( (index+1)%L_mesh==0 ){
                bac -= L_mesh;
                sx = -L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;
    }

    return;
}

int step_check(
int II,
Sample* samples,
double x_bac,
double y_bac,
double z_bac,
int index_part,
int index_cell,
int sx,
int sy,
int sz,
double amp_fac_bac
){

    int ii, JJ;

    double delta_x, delta_y, delta_z, d, D;


    for(ii=0; ii < samples[II].cell_ind[index_cell]; ++ii){

        JJ = samples[II].cell[index_cell][ii];

        delta_x = x_bac + sx - samples[II].particles[JJ].x;
        delta_y = y_bac + sy - samples[II].particles[JJ].y;
        delta_z = z_bac + sz - samples[II].particles[JJ].z;

        d = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;
        d = sqrt(d);

        D = 0.5*amp_fac_bac*samples[II].amp_fac*(samples[II].particles[index_part]._2r + samples[II].particles[JJ]._2r);

        if( d < D && JJ!=index_part )
            return 0;
    }

    return 1;
}

int legal_configuration_stat(int II, Sample* samples, double amp_fac_bac){

    int ii, jj, kk, check;
    int list[4];

    double x_bac, y_bac, z_bac;

    for(ii=0; ii < N; ++ii){


        x_bac = samples[II].particles[ii].x;
        y_bac = samples[II].particles[ii].y;
        z_bac = samples[II].particles[ii].z;

        jj = (int)( x_bac/mesh_fac ) + L_mesh*(int)( y_bac/mesh_fac ) + L2_mesh*(int)( z_bac/mesh_fac );


        for(kk=0; kk<27; ++kk){

            neighbor_cell_list(kk, jj, list);
            check = step_check(II, samples,x_bac,y_bac,z_bac,ii,list[0],list[1],list[2], list[3], amp_fac_bac);

            if(check==0)
                return 0;
        }
    }

    return 1;
}



int legal_configuration_dyn(int thread_id, Sample* samples, Seed* seeds){

    int ii, jj, kk, nn, accept_num, ii_index_cell, check, aa, bb, jj2;
    int list[4];

    double x_bac, y_bac,z_bac, ra_bac, rb_bac;

    int II = samples[thread_id].index_list;
    int x_pbc_bac, y_pbc_bac, z_pbc_bac;

    accept_num = 0;

    //     double theta, phi;

    for(nn=0; nn < N; ++nn){

        x_pbc_bac = y_pbc_bac = z_pbc_bac = 0;

        ii = (int)( xor64(&seeds[thread_id])*N );

        x_bac = samples[II].particles[ii].x + samples[II].A*( 2.*xor64(&seeds[thread_id]) - 1. );
        y_bac = samples[II].particles[ii].y + samples[II].A*( 2.*xor64(&seeds[thread_id]) - 1. );
        z_bac = samples[II].particles[ii].z + samples[II].A*( 2.*xor64(&seeds[thread_id]) - 1. );


        while(x_bac < 0.){
            x_bac += L;
            x_pbc_bac -= 1;
        }

        while(x_bac >= L){
            x_bac -= L;
            x_pbc_bac += 1;
        }


        while(y_bac < 0.){
            y_bac += L;
            y_pbc_bac -= 1;
        }


        while(y_bac >= L){
            y_bac -= L;
            y_pbc_bac += 1;
        }


        while(z_bac < 0.){
            z_bac += L;
            z_pbc_bac -= 1;
        }


        while(z_bac >= L){
            z_bac -= L;
            z_pbc_bac += 1;
        }


        jj = (int)( x_bac/mesh_fac ) + L_mesh*(int)( y_bac/mesh_fac ) + L2_mesh*(int)( z_bac/mesh_fac );


        for(kk=0; kk<27; ++kk){

            neighbor_cell_list(kk, jj, list);
            check = step_check(II, samples,x_bac,y_bac,z_bac,ii,list[0],list[1],list[2], list[3], 1);

            if(check==0)
                break;
        }

        if(check==1){

            ii_index_cell = (int)(samples[II].particles[ii].x/mesh_fac) + L_mesh*(int)( samples[II].particles[ii].y/mesh_fac ) + L2_mesh*(int)( samples[II].particles[ii].z/mesh_fac );


            if( ii_index_cell != jj ){

                for(kk=0; kk < samples[II].cell_ind[ii_index_cell]; ++kk)
                    if( samples[II].cell[ii_index_cell][kk] == ii )
                        break;


                samples[II].cell_ind[ii_index_cell] -= 1;
                samples[II].cell[ii_index_cell][kk] = samples[II].cell[ii_index_cell][ samples[II].cell_ind[ii_index_cell] ];

                samples[II].cell[jj][ samples[II].cell_ind[jj] ] = ii;
                samples[II].cell_ind[jj] += 1;

            }

            samples[II].particles[ii].x  = x_bac;
            samples[II].particles[ii].y  = y_bac;
            samples[II].particles[ii].z  = z_bac;

             samples[II].particles[ii].x_pbc += x_pbc_bac;
             samples[II].particles[ii].y_pbc += y_pbc_bac;
             samples[II].particles[ii].z_pbc += z_pbc_bac;

            accept_num += 1;
        }

    }

    return accept_num;
}


void g_r_dist_count(int II, Sample* samples){

    int ii, jj, ig;
    double delta_x, delta_y, delta_z, d, D;

    double delg = (double)L/(2.*nhis);


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
                bins[ig] += 2.;
        }

    }

    return;
}

double g_r_dist_return(int II, Sample* samples, double norm_fact){

    int ii;
    double d, nid, vb;
    double delg = (double)L/(2.*nhis);
    double rho = (double)N/( samples[II].V );

    double bac_bins;

    int list[4];


    for(ii=0; ii < nhis; ++ii){

        d = delg * (ii + 0.5);
        vb = ( (ii + 1) * (ii + 1) * (ii + 1) - ii * ii * ii ) *delg * delg * delg;
        nid =  (4./3.)*M_PI*vb*rho;

        bins[ii] /= ((double)N * nid * norm_fact);

        if(d > 1){
            bac_bins = 1. + 0.5*(4./3.)*M_PI*rho*bins[ii];
            return bac_bins;
        }
    }

    return 0.0;
}

void MC(int thread_id, Sample* samples, Seed* seeds){

    int ii;
    double accept_ratio, accept_ratio_Vol, delta_V, delta_H, bac_double;

    accept_ratio_Vol = 0.;

    int II = samples[thread_id].index_list;


    for(ii=0; ii < t_check; ii++){

        accept_ratio = (double)legal_configuration_dyn(thread_id, samples, seeds);

        accept_ratio /= (double)N;
        samples[II].A += (accept_ratio - 0.3)*samples[II].A;

        if(samples[II].A > 0.5*L)
            samples[II].A = 0.5*L;

        samples[II].accept_mov = accept_ratio;

        delta_V = samples[II].B *( 2.*xor64(&seeds[thread_id]) - 1. );
//         delta_H = samples[thread_id].P*delta_V - ( (double)N/beta )*log( (samples[II].V+delta_V)/samples[II].V );
        delta_H = samples[thread_id].P*samples[II].V*( exp(delta_V) - 1. ) - ( ((double)N+1.)/beta )*delta_V;


        if( delta_H < 0. || xor64(&seeds[thread_id]) < exp(-beta*delta_H ) ){

//             bac_double = ( samples[II].V/(samples[II].V+delta_V) );
            bac_double = 1./exp( delta_V );

            bac_double = pow( bac_double, 1./3. ); //this instruction can not be commented

            if ( legal_configuration_stat(II, samples, bac_double) ){

//                 samples[II].V += delta_V;
                samples[II].V *= exp(delta_V);
                samples[II].amp_fac *= bac_double;

                accept_ratio_Vol += 1;
            }
        }

        samples[thread_id].P += samples[thread_id].rate;
    }

    accept_ratio_Vol /= (double)t_check;
    samples[II].B += (accept_ratio_Vol - acept_Vbac)*samples[II].B;
    samples[II].accept_V = accept_ratio_Vol;



    return;
}
