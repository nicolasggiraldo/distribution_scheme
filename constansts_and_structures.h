////////////////////////////////////////////////////
// HEADER WITH ALL DATA SRUCTURES FOR THE PROGRAM //
////////////////////////////////////////////////////



/////////////////////////////
// PREPROCESSOR DIRECTIVES //
/////////////////////////////

/* Index preprocessor for the C-Order */
#define INDEX(i,j,k) (1L*k) + GV.NGRID * ( (1L*j) + GV.NGRID * (1L*i) )
#define X      0
#define Y      1
#define Z      2



////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES FOR INTERPOLATION IN DAUBECHIES-TYPES SCHEMES //
///////////////////////////////////////////////////////////////////
gsl_interp_accel *acc = NULL;
gsl_spline *spline = NULL;
int len_array_D20 = 0;
double *x_D20 = NULL;
double *y_D20 = NULL;
FILE *fin_D20 = NULL;



///////////////////////////////
// STRUCTURES OF THE PROGRAM //
///////////////////////////////

/* Global variales */
struct globalVariables{
  int      NGRID;          // Number of cell in each axis.
  long int NGRID3;         // Total number of cells (NGRID3 = NGRID^3)
  int      GADGET_VERSION; // GADGET version of the snapshot
  double   L;              // Lenght of the simulation in Mpc
  long int NP_TOT;         // Total number of particles in the simulation
  double   TOTAL_MASS;     // Total mass of all particles in the simulation
  double   RHO_MEAN;       // Mean density of ALL the simulation
  double   VOL_CELL;       // Volume of each cell
  double   H;              // Size of the cell
  char     *FILE_NAME;     // Path of the GADGET binary
  char     *SCHEME;        // Scheme used for grid assignation

  /* COSMOLOGICAL PARAMETERS OF THE SIMULATION */
  double OMEGA_M0;         //Omega matter at present time
  double OMEGA_L0;         //Omega Lambda at present time
  double ZRS;              //Redshift of the simulation
  double HUBBLEPARAM;      //Hubble parameter of the simulation
  
}GV;


/* Structure particle */
struct particle{
  double pos[3]; // Particles positions 
  double mass;   // Particle mass
};
struct particle *part=NULL;


/* Structure gadget_head */
struct gadget_head{
  unsigned int Npart[6];
  double       mass[6];
  double       time;
  double       redshift;
  int          flag_sfr;
  int          flag_feedback;
  int          npartTotal[6];
  int          flag_cooling;
  int          num_files;
  double       BoxSize;
  double       Omega0;
  double       OmegaLambda;
  double       HubbleParam; 
  char         fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8]; // Fills to 256 Bytes
}Header;


/* Struct cell  */
struct Cell{
  double   denCon;  // Density contrast of the cell
  double   rho;     // Density in the cell
}*cells=NULL;
