//######################################################################################
//paramters//###########################################################################

int N; //input
double amp_fac0;
int L_mesh, L2_mesh, cell_num;

int num_threads; //input
int num_cores; //input

int t_check; //input
int tobs; //input

double acept_Vbac;
double acept_Xbac;
//######################################################################################
//Structure definition//###############################################################
typedef struct {
    double volume_plot;
    double msd_plot;
    double com_plot;
    double p_plot;
    double pf_plot;

} Time;

typedef struct {
    double x0;
    double y0;
    double z0;

    double x;
    double y;
    double z;

    int x_pbc;
    int y_pbc;
    int z_pbc;

    int x0_pbc;
    int y0_pbc;
    int z0_pbc;


    double _2r;

} Particle;

typedef struct {
    Particle* particles;
    Time* time;

    int* cell_ind;
    int** cell;

    double P;
    double P_prime;

    double accept_num_PT;
    double accept_num_PT_bac;
    double accept_av_PT;

    double accept_mov;
    double accept_V;

    double A;
    double B;
    double V;

    double V_hs;

    int index_list;

    double amp_fac;

    double x_com;
    double y_com;
    double z_com;

    double x0_com;
    double y0_com;
    double z0_com;

    double COM0;
    double COM;

    double p_measured;
    int* bins;

} Sample;

typedef struct {
    unsigned int ux, uy;
} Seed;
//######################################################################################
