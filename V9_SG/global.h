//######################################################################################
//paramters//###########################################################################

int L;
int L2;
int N;

int tobs;
int t_check;
int t_rlx;
int num_threads;
int num_cores;

double beta0;
double betaf;

int fac;
int label;

//######################################################################################
//Structure definition//###############################################################

typedef struct {
    double q_plot;

} Time;

typedef struct {
    int x0;
    int x;

    int* J;

} Particle;


typedef struct {
    Particle* particles;
    Time* time;

    double beta;
    double beta_prime;

    int e;

    int index_list;
    double accept_num_PT;
    double accept_num_PT_bac;

} Sample;

typedef struct {
    unsigned ux,uy;
} Seed;

//######################################################################################
