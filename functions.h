/*
 * Function:  read_parameters
 * --------------------
 * Reads the parameter file in which are the main parameters 
 * necessary to run the code.
 *
 * The information loaded are: 
 * NGRID:          Number of elements in one axis of the grid.
 * FILE_NAME:      File name path of the GADGET binary file.
 * GADGET_VERSION: Version of the GADGET snapshot, could be 1 or 2.
 * SCHEME:         Scheme which is going to be used to assign 
 *                 particles to the grid, could be Nearest Grid
 *                 point (NGP), Cloud In Cell (CIC), Triangular
 *                 Shaped Cloud (TSC), Daubechies D20 scaling
 *                 function (D20) schemes.
 * 
 * The parameter file is read with the help of the library libconfig, 
 * for more information of the library use go to the manual:
 * http://www.hyperrealm.com/libconfig/libconfig_manual.html
 *
 *
 *  param_file_name: String with the name of the parameter file.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the parameter file.
 *           -2 --> There is an error whith the settings of the 
 *                  parameter file.
 */
int read_parameters(char param_file_name[]){

  config_t cfg;            /* Returns all parameters in this structure */
  const char *str1, *str2; /* Going to be used to read strings variables */
  

  /*Initialization */
  config_init(&cfg);
  
  /* Read the file. If there is an error, report it and exit. */
  if(!config_read_file(&cfg, param_file_name)){
    printf("%s:%d - %s\n",
	   config_error_file(&cfg),
	   config_error_line(&cfg),
	   config_error_text(&cfg));
    config_destroy(&cfg);
    // Value -1 means there is an error loading the param file
    return -1;
  }
  
  /* Get the NGRID value in one axis of the cell. */
  if( config_lookup_int(&cfg, "NGRID", &(GV.NGRID) ) ){
    // Checking if NGRID value is valid, that is, NGRID > 0
    if(GV.NGRID > 0){
      printf("Number of Grids in each axis: %d\n", GV.NGRID);
    }
    else{
      printf("Invalid 'NGRID' setting in configuration file.\n");
      return -2;
    }
  }
  else{
    printf("No 'NGRID' setting in configuration file.\n");
    return -2;
  }
  
  // Calculating GV.NGRID3 from GV.NGRID
  GV.NGRID3 = GV.NGRID * GV.NGRID * GV.NGRID;
  
  
  /* Get the configuration file name. */
  if(config_lookup_string( &cfg, "FILE_NAME", &str1)){
    GV.FILE_NAME = strdup(str1);
    printf("Reading from File: %s\n", GV.FILE_NAME);
  }
  else{
    printf("No 'FILE_NAME' setting in configuration file.\n");
    return -2;
  }
  
  
  /* Get the configuration file name. */
  if(config_lookup_int(&cfg, "GADGET_VERSION", &(GV.GADGET_VERSION) )){
    // Checking if GADGET version is valid
    if(GV.GADGET_VERSION == 1 || GV.GADGET_VERSION == 2){
      printf("GADGET version of the snapshot: %d\n", GV.GADGET_VERSION);
    }
    else{
      printf("Invalid 'GADGET_VERSION' setting in configuration file.\n");
      return -2;
    }
  }
  else{
    printf("\nNo 'GADGET_VERSION' setting in configuration file.");
    return -2;
  }
  
  
  /* Get the configuration file name. */
  if(config_lookup_string(&cfg, "SCHEME", &str2)){
    char ngp[]    = "NGP";
    char cic[]    = "CIC";
    char tsc[]    = "TSC";
    char daub20[] = "D20";
    if(strcmp(str2,ngp) == 0 || strcmp(str2,cic) == 0 ||
       strcmp(str2,tsc) == 0 || strcmp(str2,daub20) == 0){
      GV.SCHEME = strdup(str2);
    printf("Scheme used: %s\n", GV.SCHEME);
    }
    else{
      printf("Invalid 'SCHEME' setting in configuration file.\n");
      return -2;
    }
  }
  else{
    printf("No 'SCHEME' setting in configuration file.\n");
    return -2;
  }
  
  
  config_destroy(&cfg);
  
  
  return 0;
}



/*
 * Function:  readGADGET1BinaryFile
 * --------------------
 * Reads a GADGET binary file and store the information in the data structure
 * variable *part* it also returns the total number of particles.
 *
 *  There are no arguments in the routiene.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the snapshot.
 *           -2 --> Structure particle could not be allocated.
 *           -3 --> Can not read properly ids
 */
int readGADGET1BinaryFile(){
  FILE *fdata = NULL;
  int i, j;
  int N_min, N_max, dummy;
  //float Maux;
  float faux[3];
  //unsigned int uintaux;
  long int N_tot = 0L;
  size_t err;
  double totalMass;

  printf("\n-----------------------------------------------\n");
  printf("Reading snapshot: %s\n", GV.FILE_NAME);
  printf("GADGET version: %d\n",  GV.GADGET_VERSION);
  
  fdata = fopen(GV.FILE_NAME,"r");
  if(fdata == NULL){
    printf("File %s cannot be open\n", GV.FILE_NAME);
    return -1;
  }

  err = fread(&dummy,  sizeof(dummy),              1, fdata);
  err = fread(&Header, sizeof(struct gadget_head), 1, fdata);
  err = fread(&dummy,  sizeof(dummy),              1, fdata);

  printf("\n-----------------------------------------------\n");
  printf("Reading snapshot:\n");
  printf("----------------------------------------\n");
  for(i = 0; i<6; i++){
    N_tot += Header.npartTotal[i];
    printf("Type %d has Npart=%12d NpartTotal=%12d with mass %16.8lf\n", i,
	   Header.Npart[i], Header.npartTotal[i], Header.mass[i]);
  }//for i

  printf("----------------------------------------\n");
  printf("There is a total %ld particles in the snapshot\n", N_tot);
  printf("----------------------------------------\n");
  printf(" * Frame's Time... %16.8f\n", Header.time);
  printf(" * Redshift...     %16.8f\n", Header.redshift);
  printf(" * Flagsfr...      %16d\n",   Header.flag_sfr);
  printf(" * Flagfed...      %16d\n",   Header.flag_feedback);
  printf(" * Flagcool...     %16d\n",   Header.flag_cooling);
  printf(" * numfiles...     %16d\n",   Header.num_files);
  printf(" * Boxsize...      %16.8f\n", Header.BoxSize);
  printf(" * Omega0...       %16.8f\n", Header.Omega0);
  printf(" * OmageLa...      %16.8f\n", Header.OmegaLambda);
  printf(" * Hubbleparam...  %16.8f\n", Header.HubbleParam);
  printf("-----------------------------------------------\n\n");

  // Memory allocation for part variable
  part = (struct particle *) calloc((size_t) N_tot, sizeof(struct particle));

  if(part == NULL){
    printf("Structure particles could not be allocated\n");
    return -2;
  }//if

  /*********************/
  /* Getting positions */
  /*********************/
  err = fread(&dummy, sizeof(dummy), 1, fdata);
  for(i=0; i<N_tot; i++){
    err = fread(&faux[0], sizeof(float), 3, fdata);
    part[i].pos[X] = faux[0];
    part[i].pos[Y] = faux[1];
    part[i].pos[Z] = faux[2];
  }//for i
  
  err = fread(&dummy, sizeof(dummy), 1, fdata);
  if(dummy != (3*N_tot*sizeof(float))){
    printf(" Can not read properly positions %d %lu\n",dummy,3*N_tot*sizeof(float));
    return -3;
  }//if

  // Closing file
  fclose(fdata);

  
  /**********************/
  /* Getting velocities */
  /**********************/
  // VELOCITIES ARE NOT NECESSARY FOR OUR WORK
  //err = fread(&dummy,sizeof(dummy),1,fdata);
  //for(i=0; i<N_tot; i++){
  //err = fread(&faux[0],sizeof(float),3,fdata);
  // SKIPPING!
  //}//for i
  
  //err = fread(&dummy, sizeof(dummy), 1, fdata);
  //if(dummy != (3*N_tot*sizeof(float))){
  //printf(" Can not read properly velocities %d %lu\n",
  //dummy, 3*N_tot*sizeof(float));
  //return -3;
  //}//if

  
  /****************/
  /* Getting ID's */
  /****************/
  // AS IDENTIFICATION METHOD,
  // WE ARE GOING TO USE THE INDEX OF THE PARTICLE ARRAY
  //err = fread(&dummy, sizeof(dummy), 1, fdata);
  //for(i=0; i<N_tot; i++){
  //err = fread(&uintaux, sizeof(unsigned int), 1, fdata);
  // SKIPPING;
  //}//for i

  //err = fread(&dummy, sizeof(dummy), 1, fdata);
  //if(dummy != (N_tot*sizeof(unsigned int))){
  //printf(" Can not read properly ids %d %lu\n",
  //dummy, N_tot*sizeof(unsigned int));
  //return -3;
  //}//if

  
  /******************/
  /* Getting masses */
  /******************/
  //err = fread(&dummy, sizeof(dummy),1,fdata);
  //N_min = N_max=0;
  //for(j=0; j<=5; j++){
  //N_max=N_max+Header.npartTotal[j];
  //if((Header.mass[j]==0) && (Header.npartTotal[j]!=0)){
  //for(i=N_min;i<N_max;i++){
  //err = fread(&Maux,sizeof(float),1,fdata);
  //part[i].mass = Maux;
  //}//for i
  //}//if
  //if((Header.mass[j]!=0) && (Header.npartTotal[j]!=0)){
  //for(i=N_min;i<N_max;i++){
  //part[i].mass = Header.mass[j];
  //}//for i
  //}//if
  //N_min=N_max;
  //}//for j

  //err = fread(&dummy,sizeof(dummy),1,fdata);

  //////////////////////////////////////////
  N_min = N_max=0;
  for(j=0; j<=5; j++){
    N_max=N_max+Header.npartTotal[j];
    
    if((Header.mass[j]!=0) && (Header.npartTotal[j]!=0)){
      for(i=N_min;i<N_max;i++){
	part[i].mass = Header.mass[j];
      }//for i
    }//if
    
    N_min=N_max;

  }//for j
  
  // Total Mass calculation
  totalMass = 0.0;
  for(j=0; j<=5; j++){
    totalMass += Header.Npart[j] * Header.mass[j];
  }

  // Closing file
  //fclose(fdata);


  /***************************************************************/
  /* Storing parameters of the simulation in the gloval structure */
  /***************************************************************/

  GV.L = Header.BoxSize;                 // Lenght of the simulation in Mpc
  GV.NP_TOT = N_tot;                     /* Total number of particles in 
					    the simulation */
  GV.TOTAL_MASS = totalMass;             /* Total mass of all particles in 
					    the simulation */
  GV.RHO_MEAN = totalMass / pow(GV.L,3); // Mean density of ALL the simulation
  GV.H = GV.L / (1.0*GV.NGRID);          // Size of the cell
  GV.VOL_CELL = GV.H * GV.H * GV.H;      // Volume of each cell

  /* Cosmological parameters  */
  GV.OMEGA_M0 = Header.Omega0;           //Omega matter at present time
  GV.OMEGA_L0 = Header.OmegaLambda;      //Omega Lambda at present time
  GV.ZRS = Header.redshift;              //Redshift of the simulation
  GV.HUBBLEPARAM = Header.HubbleParam;   //Hubble parameter of the simulation


  // To avoiding warnings
  if(err){};
  
  return 0;
}



/*
 * Function:  readGADGET2BinaryFile
 * --------------------
 * Reads a GADGET2 binary file and store the information in the data structure
 * variable *part* it also returns the total number of particles.
 *
 *  There are no arguments in the routiene.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the snapshot.
 *           -2 --> Structure particle could not be allocated.
 *           -3 --> Can not read properly ids
 */
long int readGADGET2BinaryFile(){
  FILE *fdata = NULL;
  int i, j;
  int N_min, N_max, dummy;
  float Maux, faux[3];
  unsigned int uintaux;
  long int N_tot = 0L;
  char label[4];
  size_t err;
  double totalMass;

  printf("\n-----------------------------------------------\n");
  printf("Reading snapshot: %s\n", GV.FILE_NAME);
  printf("GADGET version: %d\n",  GV.GADGET_VERSION);
  
  fdata = fopen(GV.FILE_NAME,"r");
  if(fdata == NULL){
    printf("File %s cannot be open\n", GV.FILE_NAME);
    return -1;
  }


  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&label,sizeof(char),4,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  //printf(" ** %s\n",label);

  err = fread(&dummy, sizeof(dummy), 1, fdata);
  err = fread(&Header, sizeof(struct gadget_head), 1, fdata);
  err = fread(&dummy, sizeof(dummy), 1, fdata);

  printf("\n-----------------------------------------------\n");
  printf("Reading snapshot:\n");
  printf("----------------------------------------\n");
  for(i = 0; i<6; i++){
    N_tot += Header.npartTotal[i];
    printf("Type %d has Npart=%12d NpartTotal=%12d with mass %16.8lf\n", i,
	   Header.Npart[i], Header.npartTotal[i], Header.mass[i]);
  }//for i

  printf("----------------------------------------\n");
  printf("There is a total %ld particles in the snapshot\n", N_tot);
  printf("----------------------------------------\n");
  printf(" * Frame's Time... %16.8f\n", Header.time);
  printf(" * Redshift...     %16.8f\n", Header.redshift);
  printf(" * Flagsfr...      %16d\n",     Header.flag_sfr);
  printf(" * Flagfed...      %16d\n",     Header.flag_feedback);
  printf(" * Flagcool...     %16d\n",     Header.flag_cooling);
  printf(" * numfiles...     %16d\n",     Header.num_files);
  printf(" * Boxsize...      %16.8f\n", Header.BoxSize);
  printf(" * Omega0...       %16.8f\n", Header.Omega0);
  printf(" * OmageLa...      %16.8f\n", Header.OmegaLambda);
  printf(" * Hubbleparam...  %16.8f\n", Header.HubbleParam);
  printf("-----------------------------------------------\n\n");

  // Memory allocation for part variable
  part = (struct particle *) calloc((size_t) N_tot, sizeof(struct particle));

  if(part == NULL){
    printf("Structure particles could not be allocated\n");
    return -2;
  }//if

  /*********************/
  /* Getting positions */
  /*********************/
  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&label,sizeof(char),4,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  //printf(" ** %s\n",label);


  err = fread(&dummy, sizeof(dummy), 1, fdata);
  for(i=0; i<N_tot; i++){
    err = fread(&faux[0], sizeof(float), 3, fdata);
    part[i].pos[X] = faux[0];
    part[i].pos[Y] = faux[1];
    part[i].pos[Z] = faux[2];
  }//for i
  
  err = fread(&dummy, sizeof(dummy), 1, fdata);
  if(dummy != (3*N_tot*sizeof(float))){
    printf(" Can not read properly ids %d %lu\n",
	   dummy, 3*N_tot*sizeof(float));
    return -3;
  }//if

  /**********************/
  /* Getting velocities */
  /**********************/
  // VELOCITIES ARE NOT NECESSARY FOR OUR WORK
  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&label,sizeof(char),4,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  //printf(" ** %s\n",label);

  err = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0; i<N_tot; i++){
    err = fread(&faux[0],sizeof(float),3,fdata);
    // SKIPPING
  }//for i
  
  err = fread(&dummy, sizeof(dummy), 1, fdata);
  if(dummy != (3*N_tot*sizeof(float))){
    printf(" Can not read properly ids %d %lu\n",
	   dummy, 3*N_tot*sizeof(float));
    return -3;
  }//if

  /****************/
  /* Getting ID's */
  /****************/
  // AS IDENTIFICATION METHOD,
  // WE ARE GOING TO USE THE INDEX OF THE PARTICLE ARRAY
  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&label,sizeof(char),4,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  //printf(" ** %s\n",label);

  err = fread(&dummy, sizeof(dummy), 1, fdata);
  for(i=0; i<N_tot; i++){
    err = fread(&uintaux, sizeof(unsigned int), 1, fdata);
    // SKIPPING
  }//for i

  err = fread(&dummy, sizeof(dummy), 1, fdata);
  if(dummy != (N_tot*sizeof(unsigned int))){
    printf(" Can not read properly ids %d %lu\n",
	   dummy, N_tot*sizeof(unsigned int));
    return -3;
  }//if

  /*************************************/
  /* Getting individual masses masses  */
  /*************************************/

  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&label,sizeof(char),4,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  err = fread(&dummy,sizeof(dummy),1,fdata);
  //printf(" ** %s\n",label);

  err = fread(&dummy, sizeof(dummy),1,fdata);
  N_min = N_max=0;
  for(j=0; j<=5; j++){
    N_max=N_max+Header.npartTotal[j];
    
    if((Header.mass[j]==0) && (Header.npartTotal[j]!=0)){
      printf("  Reading individual masses for particles of type %d\n",j);
      for(i=N_min;i<N_max;i++){
	err = fread(&Maux,sizeof(float),1,fdata);
	part[i].mass = Maux;
      }//for i
    }//if
    
    if((Header.mass[j]!=0) && (Header.npartTotal[j]!=0)){
      for(i=N_min;i<N_max;i++){
	part[i].mass = Header.mass[j];
      }//for i
    }//if
    N_min=N_max;
  }//for j

  // Total Mass calculation
  totalMass = 0.0;
  for(j=0; j<=5; j++){
    totalMass += Header.Npart[j] * Header.mass[j];
  }
  
  err = fread(&dummy,sizeof(dummy),1,fdata);

  // Closing file
  fclose(fdata);

  
  /***************************************************************/
  /* Storing parameters of the simulation in the gloval structure */
  /***************************************************************/
  
  
  GV.L = Header.BoxSize;                 // Lenght of the simulation in Mpc
  GV.NP_TOT = N_tot;                     /* Total number of particles in 
					    the simulation */
  GV.TOTAL_MASS = totalMass;             /* Total mass of all particles in 
					    the simulation */
  GV.RHO_MEAN = totalMass / pow(GV.L,3); // Mean density of ALL the simulation
  GV.H = GV.L / (1.0*GV.NGRID);          // Size of the cell
  GV.VOL_CELL = GV.H * GV.H * GV.H;      // Volume of each cell
  
  /* Cosmological parameters  */
  GV.OMEGA_M0 = Header.Omega0;           //Omega matter at present time
  GV.OMEGA_L0 = Header.OmegaLambda;      //Omega Lambda at present time
  GV.ZRS = Header.redshift;              //Redshift of the simulation
  GV.HUBBLEPARAM = Header.HubbleParam;   //Hubble parameter of the simulation


  // To avoiding warnings
  if(err){};
  
  return 0;
}




/*
 * Function:  compare  
 * --------------------
 * Comparation function for the use of the qsort routine.
 * This function gives the necessary rules for the sorting.
 * 
 *  eleme1 and elem2: Elements of the struct particle type 
 *                    to sort according to the z-axis value.
 *   
 *  returns: 1,  if the z value of the first particle is 
 *               greater than the second particle.
 *           -1, if the z value of the second particle is 
 *               greater than the first particle.
 *           0,  if both particles have the same z value.                                                
 */
int compare(struct particle *elem1, struct particle *elem2){
  
  // Comparition for the qsort routine
  if( elem1->pos[Z] < elem2->pos[Z] ){
    return -1;
  }else if( elem1->pos[Z] > elem2->pos[Z] ){
    return 1;
  }else{
    return 0;
  }
  
}

/* This type definition is necessary 
   for the use of the qsort routine */
typedef int (*compfn)(const void*, const void*);




/*
 * Function:  W_NGP
 * --------------------
 * Computes the window function for the Nearest GRID POINT (NGP) scheme. 
 *
 *  x, y, z: position in x, y and z.
 *
 *  H: separation of grids ( H = L / N_grid ).
 *
 *  returns: value of the window function according to the distribution 
 *           scheme.
 */
double W_NGP(double x, double y, double z, double H){
  //double Wx, Wy, Wz;

  /* Nearest Grid Point */
  // One dimensional window function in the X-axis
  /*
    if( fabs(x) < H*0.5 ){
    Wx = 1.0;
    }else if( fabs(x) == H*0.5  ){
    Wx = 0.5;
    }else{
    Wx = 0.0;
    }
  */
  
  // One dimensional window function in the Y-axis
  /*
    if( fabs(y) < H*0.5 ){
    Wy = 1.0;
    }else if( fabs(y) == H*0.5  ){
    Wy = 0.5;
    }else{
    Wy = 0.0;
    }
  */
  
  // One dimensional window function in the Z-axis
  /*
    if( fabs(z) < H*0.5 ){
    Wz = 1.0;
    }else if( fabs(z) == H*0.5  ){
    Wz = 0.5;
    }else{
    Wz = 0.0;
    }
  */

  
  //return Wx * Wy * Wz;
  /* As we use a regular cubic grid, the 
     three dimensional window function is 
     given as the multiplication of three 
     one dimensional window functions.
     WindowFunction = W(x) * W(y) * W(z).
  */
  return 1.0;
}



/*
 * Function:  W_CIC
 * --------------------
 * Computes the window function for the Cloud In Cell (CIC) scheme. 
 *
 *  x, y, z: position in x, y and z.
 *
 *  H: separation of grids ( H = L / N_grid ).
 *
 *  returns: value of the window function according to the distribution 
 *           scheme.
 */
double W_CIC(double x, double y, double z, double H){
  double Wx, Wy, Wz;

  /* Cloud In Cell */
  // One dimensional window function in the X-axis
  if( fabs(x) < H ){
    Wx = 1.0 - fabs(x)/H;
  }else{
    Wx = 0.0;
  }
  
  // One dimensional window function in the Y-axis
  if( fabs(y) < H ){
    Wy = 1.0 - fabs(y)/H;
  }else{
    Wy = 0.0;
  }
  
  // One dimensional window function in the Z-axis
  if( fabs(z) < H ){
    Wz = 1.0 - fabs(z)/H;
  }else{
    Wz = 0.0;
  }

  
  return Wx * Wy * Wz; /* As we use a regular cubic grid, the 
			  three dimensional window function is 
			  given as the multiplication of three 
			  one dimensional window functions.
			  WindowFunction = W(x) * W(y) * W(z).
		       */
}



/*
 * Function:  W_TSC
 * --------------------
 * Computes the window function for the Triangular Shaped Cloud (TSC) 
 * scheme. 
 *
 *  x, y, z: position in x, y and z.
 *
 *  H: separation of grids ( H = L / N_grid ).
 *
 *  returns: value of the window function according to the distribution 
 *           scheme.
 */
double W_TSC(double x, double y, double z, double H){
  double Wx, Wy, Wz;

  /* Triangular Shaped Cloud */
  // One dimensional window function in the X-axis
  if( fabs(x) <= H*0.5 ){
    Wx = 0.75 - ( (x*x)/(H*H) );
  }else if( H*0.5 <= fabs(x) && fabs(x) <= 1.5*H  ){
    Wx = 0.5*(1.5 - fabs(x)/H)*(1.5 - fabs(x)/H);
  }else{
    Wx = 0.0;
  }
  
  // One dimensional window function in the Y-axis
  if( fabs(y) <= H*0.5 ){
    Wy = 0.75 - ( (y*y)/(H*H) );
  }else if( H*0.5 <= fabs(y) && fabs(y) <= 1.5*H  ){
    Wy = 0.5*(1.5 - fabs(y)/H)*(1.5 - fabs(y)/H);
  }else{
    Wy = 0.0;
  }
  
  // One dimensional window function in the Z-axis
  if( fabs(z) <= H*0.5 ){
    Wz = 0.75 - ( (z*z)/(H*H) );
  }else if( H*0.5 <= fabs(z) && fabs(z) <= 1.5*H  ){
    Wz = 0.5*(1.5 - fabs(z)/H)*(1.5 - fabs(z)/H);
  }else{
    Wz = 0.0;
  }
  

  return Wx * Wy * Wz; /* As we use the regular cubic grid, the 
			  three dimensional window function is 
			  given as the multiplication of three 
			  one dimensional window functions.
			  WindowFunction = W(x) * W(y) * W(z).
		       */
}



/*
 * Function:  W_D20
 * --------------------
 * Computes the window function for the Daubechies D20 scaling 
 * function scheme. 
 *
 *  x, y, z: position in x, y and z.
 *
 *  H: separation of grids ( H = L / N_grid ).
 *
 *  returns: value of the window function according to the distribution 
 *           scheme.
 */
double W_D20(double x, double y, double z, double H){
  double Wx, Wy, Wz;

  // One dimensional window function in the X-axis
  Wx = gsl_spline_eval(spline, fabs(x)/GV.H, acc);

  // One dimensional window function in the Y-axis
  Wy = gsl_spline_eval(spline, fabs(y)/GV.H, acc);

  // One dimensional window function in the Z-axis
  Wz = gsl_spline_eval(spline, fabs(z)/GV.H, acc);

  return Wx * Wy * Wz; /* As we use the regular cubic grid, the 
			  three dimensional window function is 
			  given as the multiplication of three 
			  one dimensional window functions.
			  WindowFunction = W(x) * W(y) * W(z).
		       */
}



/*
 * Function:  mod
 * --------------------
 * Calculate the modulo operation for two numbers a and b (a%b) 
 * includding negative numbers.
 *
 *  a: Numerator of the division.
 *  b: Denominator of the division.
 *
 *  returns: The modulo a%b includding the option for negative numbers.
 */
int mod(int a, int b){
  int mod = a%b;
  while(mod<0){
    mod += b;
  }
  return mod;
}
