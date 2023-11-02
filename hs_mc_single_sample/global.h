//######################################################################################
//paramters//###########################################################################

int N; //input
double amp_fac0;
int L_mesh, L2_mesh, cell_num;

double P0; //input

int t_check; //input
double P_delta;

double p_max;
double pc_fr_max;
double B_accept;

//######################################################################################
//Structure definition//###############################################################

typedef struct {
    double x0;
    double y0;
    double z0;


    double x;
    double y;
    double z;

    double _2r;

} Particle;

typedef struct {
    Particle* particles;

    int* cell_ind;
    int** cell;

    double P;
    double P_prime;
    double rate;


    double accept_num_PT;
    double accept_num_PT_bac;

    double accept_mov;
    double accept_V;

    double A;
    double B;

    double V;

    double V_hs;

    int index_list;

    double amp_fac;

    double p_measured;
    double V_bac;

} Sample;

typedef struct {
    unsigned int ux, uy;
} Seed;
//######################################################################################
