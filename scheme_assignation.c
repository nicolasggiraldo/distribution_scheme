#include          <stdio.h>
#include         <stdlib.h>
#include           <math.h>
#include      <libconfig.h>
#include         <string.h>
#include <gsl/gsl_spline.h>

#include "constansts_and_structures.h"
#include                 "functions.h"

int main(int argc, char *argv[]){
  
  int i, j, k; // Counter for X, Y, Z axis
  int ii, jj, kk; // Aditional counters for the scheme assignation calculation
  long int id; // Counter for the identificator of the particles
  long int index_cell, indexaux; /* counters for the cells */
  int index_start=0, index_end=0; /* Index to start and finish for the scheme 
				     assignation according to the scheme used */
  double (*W)(double,double,double,double) = NULL; /* Memory addres to 
						      the window function */
  double xp, yp, zp; // Position of the particles in X, Y, Z axis
  double xc, yc, zc; // Position of the cells in X, Y, Z axis
  // Buffer with the output path of the cell
  char buffer[200];
  // Outfile in which save the cell
  FILE *outfile = NULL;
  // Variable in which save the total mass
  double totalMass;
  
  
  //////////////////////////////
  //* READING PARAMETER FILE *//
  //////////////////////////////
  if(argc<2){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("%s: You must specify the name of the parameter file\n",argv[0]);
    printf("For example: %s pathOfFile/parameterFile.txt\n",argv[0]);
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  // Reading parameter file and verifying there is no error.
  switch( read_parameters(argv[1]) ){
  case -1 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad path to the parameter file.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  case -2 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad settings in the parameter file.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  //////////////////////////////////
  //* READING GADGET BINARY FILE *//
  //////////////////////////////////
  // According to the GADGET version, snapshot has different (subtle) format.
  switch( GV.GADGET_VERSION ){
    
  case 1 :
    
    switch( readGADGET1BinaryFile() ){
    case -1 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Problems loading the snapshot.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    case -2 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Structure particle could not be allocated.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    case -3 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Can not read properly ids.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
    
    break;
    
  case 2 :
    
    switch( readGADGET2BinaryFile() ){
    case -1 :
      printf("\n***********************************");
      printf("***********************************\n");	
      printf("Error: Problems loading the snapshot.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    case -2 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Structure particle could not be allocated.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    case -3 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Can not read properly ids.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
    
    break;
    
  }
  
  
  
  ///////////////////////////////////////////
  //* DEFINING THE WINDOW FUNCTION TO USE *//
  ///////////////////////////////////////////
  if(      strcmp(GV.SCHEME, "NGP") == 0 ){
    // NGP
    W = W_NGP;
    index_start = 0;
    index_end   = 0;
  }
  else if( strcmp(GV.SCHEME, "CIC") == 0 ){
    // CIC
    W = W_CIC;
    index_start = -1;
    index_end   =  1;
  }
  else if( strcmp(GV.SCHEME, "TSC") == 0 ){
    // TSC
    W = W_TSC;
    index_start = -2;
    index_end   =  2;
  }
  else if( strcmp(GV.SCHEME, "D20") == 0 ){
    int err;
    // D20
    
    //len_array_D20 = 192;
    //x_D20 = (double *) calloc( len_array_D20, sizeof(double));
    //y_D20 = (double *) calloc( len_array_D20, sizeof(double));
    //fin_D20 = fopen("./D20.bin", "rb");

    //fread(x_D20, sizeof(double), len_array_D20, fin_D20);
    //fread(y_D20, sizeof(double), len_array_D20, fin_D20);
    //fclose(fin_D20);

    len_array_D20 = 610;
    x_D20 = (double *) calloc( len_array_D20, sizeof(double));
    y_D20 = (double *) calloc( len_array_D20, sizeof(double));

    fin_D20 = fopen("./D20.txt", "r");

    for(i=0; i<len_array_D20; i++)
      err=fscanf(fin_D20, "%lf %lf", &x_D20[i], &y_D20[i]);
    
    fclose(fin_D20);

    // GSL interpolation allocation
    acc    = gsl_interp_accel_alloc(); // accelerator
    spline = gsl_spline_alloc(gsl_interp_cspline, len_array_D20); // spline

    // GSL init
    gsl_spline_init(spline, x_D20, y_D20, len_array_D20);
    
    W = W_D20;
    index_start = 0;
    index_end   = 10;
    //index_end   = 8;
    
    if(err){}
  }



  //////////////////////////////
  //* CELL STRUCT ALLOCATION *//
  //////////////////////////////
  /* Array of structure Cell, size NGRID^3 */
  cells = (struct Cell *) calloc( GV.NGRID3, sizeof( struct Cell) );

  // Setting values to zero at the beggining
  /*
    for(id=0L; id < GV.NGRID3; id++){
    cells[id].Np_cell = 0L;
    cells[id].denCon  = 0.0;
    cells[id].rho     = 0.0;
    }
  */
  
  
  
  //////////////////////////////
  //* FROM PARTICLES TO GRID *//
  //////////////////////////////
  printf("Locating particles into the grid\n");
  printf("-----------------------------------------------\n");

  /* Locating cells */
  printf("I am in the %3.0lf%% of the particles", 0.0 );
  fflush(stdout);
  for(id=0L; id < GV.NP_TOT; id++){
    
    if( (id%1000000) == 0 ){
      //printf("\rI am in the particle %ld of %ld", long_id, GV.NpTot);
      printf("\rI am in the %3.0lf%% of the particles", (1.0*id)/GV.NP_TOT * 100.0 );
      fflush(stdout);
    }
    
    /* Coordinates of the particle  */
    xp = part[id].pos[X];
    yp = part[id].pos[Y];
    zp = part[id].pos[Z];

    // index in the x, y and z-axis
    i = (int) floor( (xp / GV.L) * GV.NGRID );
    j = (int) floor( (yp / GV.L) * GV.NGRID );
    k = (int) floor( (zp / GV.L) * GV.NGRID );
    
    // index in C-order
    index_cell = INDEX(i,j,k);

    /* Calculating scheme assignation */
    for(     ii=index_start; ii<=index_end; ii++){
      for(   jj=index_start; jj<=index_end; jj++){
	for( kk=index_start; kk<=index_end; kk++){
	  indexaux = INDEX(mod(i+ii, GV.NGRID),
			   mod(j+jj, GV.NGRID),
			   mod(k+kk, GV.NGRID));
	  xc = GV.H * (0.5 + i+ii);
	  yc = GV.H * (0.5 + j+jj);
	  zc = GV.H * (0.5 + k+kk);
	  cells[indexaux].rho += part[id].mass * W(xc-xp, yc-yp, zc-zp, GV.H);
	}
      }
    }
  }
  printf("\n");
  
  
  
  //////////////////////////////////////////////////////
  //* TERMINATING THE CALCULATION AND SAVING IN FILE *//
  //////////////////////////////////////////////////////
  if( strcmp(GV.SCHEME, "NGP") == 0 ){
    // NGP
    sprintf(buffer, "%s_NGRID_%d_NGP.dens", GV.FILE_NAME, GV.NGRID);
  }
  else if( strcmp(GV.SCHEME, "CIC") == 0 ){
    // CIC
    sprintf (buffer, "%s_NGRID_%d_CIC.dens", GV.FILE_NAME, GV.NGRID);
  }
  else if( strcmp(GV.SCHEME, "TSC") == 0 ){
    // TSC
    sprintf (buffer, "%s_NGRID_%d_TSC.dens", GV.FILE_NAME, GV.NGRID);
  }
  else if( strcmp(GV.SCHEME, "D20") == 0 ){
    // D20
    sprintf (buffer, "%s_NGRID_%d_D20.dens", GV.FILE_NAME, GV.NGRID);
  }
  
  // Opening file for output of the cell
  outfile = fopen(buffer, "wb");

  printf("Saving data in %s\n", buffer);
  printf("-----------------------------------------------\n");
  
  // Setting total mass to zero
  totalMass = 0.0;
  
  /* Saving cosmological parameters of the simulation */
  fwrite(    &GV.OMEGA_M0, sizeof(double), 1, outfile);
  fwrite(    &GV.OMEGA_L0, sizeof(double), 1, outfile);
  fwrite(         &GV.ZRS, sizeof(double), 1, outfile);
  fwrite( &GV.HUBBLEPARAM, sizeof(double), 1, outfile);

  /* Saving simulation parameters */
  fwrite(          &GV.NGRID,       sizeof(int), 1, outfile);
  fwrite( &GV.GADGET_VERSION,       sizeof(int), 1, outfile);
  fwrite(              &GV.L,    sizeof(double), 1, outfile);
  fwrite(         &GV.NP_TOT,  sizeof(long int), 1, outfile);
  fwrite(     &GV.TOTAL_MASS,    sizeof(double), 1, outfile);
  fwrite(       &GV.RHO_MEAN,    sizeof(double), 1, outfile);
  fwrite(       &GV.VOL_CELL,    sizeof(double), 1, outfile);
  fwrite(              &GV.H,    sizeof(double), 1, outfile);
  fwrite(    &(GV.SCHEME[0]),      sizeof(char), 1, outfile);
  fwrite(    &(GV.SCHEME[1]),      sizeof(char), 1, outfile);
  fwrite(    &(GV.SCHEME[2]),      sizeof(char), 1, outfile);

  /* Saving cell data */
  for(i=0; i<GV.NGRID; i++){
    for(j=0; j<GV.NGRID; j++){
      for(k=0; k<GV.NGRID; k++){

	index_cell = INDEX(i,j,k);
	
	/* positions of the cells */
	//xc = GV.H * (0.5 + i);
	//yc = GV.H * (0.5 + j);
	//zc = GV.H * (0.5 + k);
	
	/* Calculating the final density in the cell.
	   This is made as a verification of the scheme,
	   the mass must conservate with any scheme used. */
	totalMass += cells[index_cell].rho;

	/* Calculating the density in the cell */
	cells[index_cell].rho /= GV.VOL_CELL;
	
	/* Calculating the final density contrast in the cell */
	cells[index_cell].denCon = (cells[index_cell].rho / GV.RHO_MEAN) - 1.0;
		
	// Writing density contrast
	fwrite(&cells[index_cell].denCon, sizeof(double), 1, outfile);
      }
    }
  }
  fclose(outfile);
  
  
  /* Validation of the scheme */
  if( strcmp(GV.SCHEME, "NGP") == 0 ){
    // NGP
    printf("Mass NGP        = %lf\n", totalMass);
  }
  else if( strcmp(GV.SCHEME, "CIC") == 0 ){
    // CIC
    printf("Mass CIC        = %lf\n", totalMass);
  }
  else if( strcmp(GV.SCHEME, "TSC") == 0 ){
    // TSC
    printf("Mass TSC        = %lf\n", totalMass);
  }
  else if( strcmp(GV.SCHEME, "D20") == 0 ){
    // D20
    printf("Mass D20        = %lf\n", totalMass);
    
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    free(x_D20);
    free(y_D20);
  }
  
  printf("Mass Simulation = %lf\n", GV.TOTAL_MASS);
  printf("%% Difference    = %.4e%%\n",
	 (100.0 * (totalMass-GV.TOTAL_MASS)) / GV.TOTAL_MASS);

  
  /* Freeing up memory allocation */  
  free(part);
  free(cells);
  
  return 0;
}
