void MC_PT(Sample* samples, int time_bac, Seed* seeds){

    int ii, II, II0;
    double delta;

    if(time_bac%2==0)
        for(ii=1; ii < num_threads; ii+=2){

            II = samples[ii].index_list;
            II0 = samples[ii-1].index_list;

            delta = (samples[ii].P - samples[ii-1].P)*(samples[II].V - samples[II0].V);
            if( delta > 0. || xor64(&seeds[num_threads]) < exp(delta) ){

                samples[ii].index_list = II0;
                samples[ii-1].index_list = II;

                samples[ii-1].accept_num_PT += 1.;
            }
        }

    else
        for(ii=2; ii < num_threads; ii+=2){

            II = samples[ii].index_list;
            II0 = samples[ii-1].index_list;

            delta = (samples[ii].P - samples[ii-1].P)*(samples[II].V - samples[II0].V);
            if( delta > 0. || xor64(&seeds[num_threads]) < exp(delta) ){

                samples[ii].index_list = II0;
                samples[ii-1].index_list = II;

                samples[ii-1].accept_num_PT += 1.;
            }
        }

    return;
}
