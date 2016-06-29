/*
 * Filename   : cg.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : mié 29 jun 2016 17:32:19 CEST
 *
 * Author     : jatorre@fisfun.uned.es
 *
 * Purpose    : Obtain mesoscopic variables
 *              from microscopic details
 *
 */

#include "cg.h"

int main (int argc, char *argv[]) {

  PrintInitInfo();

  PrintComputingOptions();

  PrintMsg("INIT");

  char * filestr;
  if (argc == 2)
  {
    filestr = argv[1];
    PrintMsg("argv[1] provided.");
    printf("Using %s as input file list...\n", filestr);
  }
  else 
  {
    filestr = "sim";
    PrintMsg("argv[1] not provided. Generating input data and storing file list in 'sim'");
    PrintMsg("Creating snapshots from the simulation files...");
    PrepareInputFiles();
  }
  
  char str[100]; 
  memset(str,'\0',sizeof(str));
  char basename[7];

  // Read input list
  char Snapshot[NSteps][7];
  ReadInputFiles(filestr, Snapshot);
 
  // INIT OF BLOCK. Create output files
  struct OutputFiles oFile;
  
  // Create output directory if it does not exist
  struct stat status;
  if (stat("./output", &status) != 0) 
    mkdir("./output",0755);
  
  // Microscopic files
  //   strcpy (str, "./output/");
  //   strcat (str, filestr);
  //   strcat (str, ".MicrozForce.dat");
  //   oFile.MicrozForce = fopen(str, "w");
  //   
  //   strcpy (str, "./output/");
  //   strcat (str, filestr);
  //   strcat (str, ".MicroEnergy.dat");
  //   oFile.MicroEnergy = fopen(str, "w");
  //   
  //   strcpy (str, "./output/");
  //   strcat (str, filestr);
  //   strcat (str, ".MicroKinetic.dat");
  //   oFile.MicroKinetic = fopen(str, "w");
  //   
  //   strcpy (str, "./output/");
  //   strcat (str, filestr);
  //   strcat (str, ".MicroVmod.dat");
  //   oFile.MicroVmod = fopen(str, "w");

  // Mesoscopic files
  #if __COMPUTE_DENSITY__
    sprintf(str, "./output/%s.MesoDensity_0.dat", filestr);
    oFile.MesoDensity_0 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoDensity_1.dat", filestr);
    oFile.MesoDensity_1 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoDensity_2.dat", filestr);
    oFile.MesoDensity_2 = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_matrix * BinaryMesoDensity_0 = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoDensity_1 = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoDensity_2 = gsl_matrix_calloc(NSteps,NNodes);
    #endif

  #endif

  #if __COMPUTE_FORCE__
    sprintf(str, "./output/%s.MesoForce_x.dat", filestr);
    oFile.MesoxForce = fopen(str, "w");
    sprintf(str, "./output/%s.MesoForce_y.dat", filestr);
    oFile.MesoyForce = fopen(str, "w");
    sprintf(str, "./output/%s.MesoForce_z.dat", filestr);
    oFile.MesozForce = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_matrix * BinaryMesoForce_x = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoForce_y = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoForce_z = gsl_matrix_calloc(NSteps,NNodes);
    #endif

  #endif

  #if __COMPUTE_ENERGY__
    sprintf(str, "./output/%s.MesoEnergy.dat", filestr);
    oFile.MesoEnergy = fopen(str, "w");
    sprintf(str, "./output/%s.MesoKinetic.dat", filestr);
    oFile.MesoKinetic = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_matrix * BinaryMesoEnergy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoKinetic = gsl_matrix_calloc(NSteps,NNodes);
    #endif

  #endif

  #if __COMPUTE_TEMPERATURE__
    sprintf(str, "./output/%s.MesoTemp.dat", filestr);
    oFile.MesoTemp = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_matrix * BinaryMesoTemp = gsl_matrix_calloc(NSteps,NNodes);
    #endif

  #endif 

  #if __COMPUTE_STRESS__
    // Kinetic stress tensor
    sprintf(str, "./output/%s.MesoSigma1_xx.dat", filestr);
    oFile.MesoSigma1_00 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma1_xy.dat", filestr);
    oFile.MesoSigma1_01 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma1_xz.dat", filestr);
    oFile.MesoSigma1_02 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma1_yx.dat", filestr);
    oFile.MesoSigma1_10 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma1_yy.dat", filestr);
    oFile.MesoSigma1_11 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma1_yz.dat", filestr);
    oFile.MesoSigma1_12 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma1_zx.dat", filestr);
    oFile.MesoSigma1_20 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma1_zy.dat", filestr);
    oFile.MesoSigma1_21 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma1_zz.dat", filestr);
    oFile.MesoSigma1_22 = fopen(str, "w");
  
    // Virial stress tensor
    sprintf(str, "./output/%s.MesoSigma2_xx.dat", filestr);
    oFile.MesoSigma2_00 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma2_xy.dat", filestr);
    oFile.MesoSigma2_01 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma2_xz.dat", filestr);
    oFile.MesoSigma2_02 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma2_yx.dat", filestr);
    oFile.MesoSigma2_10 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma2_yy.dat", filestr);
    oFile.MesoSigma2_11 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma2_yz.dat", filestr);
    oFile.MesoSigma2_12 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma2_zx.dat", filestr);
    oFile.MesoSigma2_20 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma2_zy.dat", filestr);
    oFile.MesoSigma2_21 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma2_zz.dat", filestr);
    oFile.MesoSigma2_22 = fopen(str, "w");

    // Total stress tensor
    sprintf(str, "./output/%s.MesoSigma_xx.dat", filestr);
    oFile.MesoSigma_00 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma_xy.dat", filestr);
    oFile.MesoSigma_01 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma_xz.dat", filestr);
    oFile.MesoSigma_02 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma_yx.dat", filestr);
    oFile.MesoSigma_10 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma_yy.dat", filestr);
    oFile.MesoSigma_11 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma_yz.dat", filestr);
    oFile.MesoSigma_12 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma_zx.dat", filestr);
    oFile.MesoSigma_20 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma_zy.dat", filestr);
    oFile.MesoSigma_21 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoSigma_zz.dat", filestr);
    oFile.MesoSigma_22 = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_matrix * BinaryMesoSigma_xx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma_xy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma_xz  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma_yx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma_yy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma_yz  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma_zx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma_zy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma_zz  = gsl_matrix_calloc(NSteps,NNodes);
      
      gsl_matrix * BinaryMesoSigma1_xx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma1_xy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma1_xz  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma1_yx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma1_yy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma1_yz  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma1_zx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma1_zy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma1_zz  = gsl_matrix_calloc(NSteps,NNodes);
      
      gsl_matrix * BinaryMesoSigma2_xx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma2_xy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma2_xz  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma2_yx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma2_yy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma2_yz  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma2_zx  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma2_zy  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoSigma2_zz  = gsl_matrix_calloc(NSteps,NNodes);
    #endif

  #endif

  #if __COMPUTE_MOMENTUM__
    sprintf(str, "./output/%s.MesoMomentum_x.dat", filestr);
    oFile.MesoMomentum_x = fopen(str, "w");
    sprintf(str, "./output/%s.MesoMomentum_y.dat", filestr);
    oFile.MesoMomentum_y = fopen(str, "w");
    sprintf(str, "./output/%s.MesoMomentum_z.dat", filestr);
    oFile.MesoMomentum_z = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_matrix * BinaryMesoMomentum_x  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoMomentum_y  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoMomentum_z  = gsl_matrix_calloc(NSteps,NNodes);
    #endif

  #endif

  #if __COMPUTE_VELOCITY__
    sprintf(str, "./output/%s.MesoVelocity_x.dat", filestr);
    oFile.MesoVelocity_x = fopen(str, "w");
    sprintf(str, "./output/%s.MesoVelocity_y.dat", filestr);
    oFile.MesoVelocity_y = fopen(str, "w");
    sprintf(str, "./output/%s.MesoVelocity_z.dat", filestr);
    oFile.MesoVelocity_z = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_matrix * BinaryMesoVelocity_x  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoVelocity_y  = gsl_matrix_calloc(NSteps,NNodes);
      gsl_matrix * BinaryMesoVelocity_z  = gsl_matrix_calloc(NSteps,NNodes);
    #endif

  #endif

  #if __COMPUTE_INTERNAL_ENERGY__
    sprintf(str, "./output/%s.MesoInternalEnergy.dat", filestr);
    oFile.MesoInternalEnergy = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_matrix * BinaryMesoInternalEnergy = gsl_matrix_calloc(NSteps,NNodes);
    #endif

  #endif
    
  #if __COMPUTE_MACRO_ENERGY__
    sprintf(str, "./output/%s.MacroEnergyUpperWall.dat", filestr);
    oFile.MacroEnergyUpperWall = fopen(str, "w");
    sprintf(str, "./output/%s.MacroEnergyLowerWall.dat", filestr);
    oFile.MacroEnergyLowerWall = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_vector * BinaryMacroEnergyUpperWall = gsl_vector_calloc(NSteps);
      gsl_vector * BinaryMacroEnergyLowerWall = gsl_vector_calloc(NSteps);
    #endif

  #endif
  
  #if __COMPUTE_MACRO_MOMENTUM__
    sprintf(str, "./output/%s.MacroMomentumUpperWall.dat", filestr);
    oFile.MacroMomentumUpperWall = fopen(str, "w");
    sprintf(str, "./output/%s.MacroMomentumLowerWall.dat", filestr);
    oFile.MacroMomentumLowerWall = fopen(str, "w");

    #if __BINARY_OUTPUT__
      sprintf(str, "./output/%s.MacroMomentumUpperWallx.dat", filestr);
      oFile.MacroMomentumUpperWallx = fopen(str, "w");
      sprintf(str, "./output/%s.MacroMomentumUpperWally.dat", filestr);
      oFile.MacroMomentumUpperWally = fopen(str, "w");
      sprintf(str, "./output/%s.MacroMomentumUpperWallz.dat", filestr);
      oFile.MacroMomentumUpperWallz = fopen(str, "w");
      sprintf(str, "./output/%s.MacroMomentumLowerWallx.dat", filestr);
      oFile.MacroMomentumLowerWallx = fopen(str, "w");
      sprintf(str, "./output/%s.MacroMomentumLowerWally.dat", filestr);
      oFile.MacroMomentumLowerWally = fopen(str, "w");
      sprintf(str, "./output/%s.MacroMomentumLowerWallz.dat", filestr);
      oFile.MacroMomentumLowerWallz = fopen(str, "w");
      gsl_vector * BinaryMacroMomentumUpperWallx = gsl_vector_calloc(NSteps);
      gsl_vector * BinaryMacroMomentumUpperWally = gsl_vector_calloc(NSteps);
      gsl_vector * BinaryMacroMomentumUpperWallz = gsl_vector_calloc(NSteps);
      gsl_vector * BinaryMacroMomentumLowerWallx = gsl_vector_calloc(NSteps);
      gsl_vector * BinaryMacroMomentumLowerWally = gsl_vector_calloc(NSteps);
      gsl_vector * BinaryMacroMomentumLowerWallz = gsl_vector_calloc(NSteps);
    #endif

  #endif

  #if __COMPUTE_CENTER_OF_MASS__
    sprintf(str, "./output/%s.CenterOfMassUpperWall.dat", filestr);
    oFile.CenterOfMassUpperWall = fopen(str, "w");
    sprintf(str, "./output/%s.CenterOfMassLowerWall.dat", filestr);
    oFile.CenterOfMassLowerWall = fopen(str, "w");

    #if __BINARY_OUTPUT__
      gsl_vector * BinaryCenterOfMassUpperWall = gsl_vector_calloc(NSteps);
      gsl_vector * BinaryCenterOfMassLowerWall = gsl_vector_calloc(NSteps);
    #endif

  #endif
  
  // END OF BLOCK. All output files created

  // INIT OF BLOCK. Computing vectors and matrices that
  // are constant throughout the program

  PrintMsg("Generating node positions...");
  gsl_vector * z = gsl_vector_calloc(NNodes);
  Compute_Node_Positions(z);
  sprintf(str, "./output/%s.MesoNodes.dat", filestr);
  SaveVectorWithoutIndex(z, str);

  PrintMsg("Obtaining neighboring matrix...");
  gsl_matrix * Neighbors = gsl_matrix_calloc (Mx*My*Mz,27);
  Compute_NeighborMatrix(Neighbors);

  // END OF BLOCK. All constant quantities created

  // BEGIN OF BLOCK. Definition of needed vectors, matrices, and so on

  // Positions, velocities and momentum
  gsl_matrix * BaseData   = gsl_matrix_calloc (NParticles,8);

  gsl_matrix * Positions  = gsl_matrix_calloc (NParticles,4);
  gsl_matrix * Velocities = gsl_matrix_calloc (NParticles,3);
  gsl_matrix * Momentum   = gsl_matrix_calloc (NParticles,3);

  FILE *BaseFile;

  // Linked list
  gsl_vector * List      = gsl_vector_calloc (NParticles);
  gsl_vector * ListHead  = gsl_vector_calloc (Mx*My*Mz);

  // Microscopic variables
  gsl_matrix * Force   = gsl_matrix_calloc (NParticles,3);
  gsl_vector * Energy  = gsl_vector_calloc (NParticles);
  gsl_vector * Kinetic = gsl_vector_calloc (NParticles);

  // Mesoscopic variables
  gsl_matrix * MesoDensity      = gsl_matrix_calloc (NNodes,2);
  gsl_matrix * MesoForce        = gsl_matrix_calloc (NNodes,3);
  gsl_vector * TotalMesoDensity = gsl_vector_calloc (NNodes);
  gsl_vector * MesoEnergy       = gsl_vector_calloc (NNodes);
  gsl_vector * MesoKinetic      = gsl_vector_calloc (NNodes);
  gsl_vector * MesoTemp         = gsl_vector_calloc (NNodes);
  // Sigma matrices will be stored in the following order 
  // 00, 01, 02, 10, 11, 12, 20, 21, 22
  gsl_matrix * MesoSigma1   = gsl_matrix_calloc (NNodes,9);
  gsl_matrix * MesoSigma2   = gsl_matrix_calloc (NNodes,9);
  gsl_matrix * MesoSigma    = gsl_matrix_calloc (NNodes,9);
  
  gsl_matrix * MesoMomentum = gsl_matrix_calloc (NNodes,3);
  gsl_matrix * MesoVelocity = gsl_matrix_calloc (NNodes,3);
  
  gsl_vector * MesoInternalEnergy = gsl_vector_calloc (NNodes);

  // END OF BLOCK

  // START COMPUTATION

  for (int Step=0;Step<NSteps;Step++)
  {
    strcpy(basename,Snapshot[Step]);

    // BaseData is a matrix that stores: 
    // ID TYPE x y z vx vy vz
    // The ID of a particle corresponds to the row
    PrintMsg("Reading microscopic data");
    // REDUNDANT? gsl_matrix_set_zero(BaseData);
    sprintf(str, "./data/%s", basename);
    printf("\tInput file: %s\n", str);
    BaseFile = fopen(str, "r");

      gsl_matrix_fscanf(BaseFile, BaseData);

      gsl_vector_view Position_0      = gsl_matrix_column(BaseData,1);
      gsl_matrix_set_col(Positions,0,&Position_0.vector);
      gsl_vector_view Position_1      = gsl_matrix_column(BaseData,2);
      gsl_matrix_set_col(Positions,1,&Position_1.vector);
      gsl_vector_view Position_2      = gsl_matrix_column(BaseData,3);
      gsl_matrix_set_col(Positions,2,&Position_2.vector);
      gsl_vector_view Position_3      = gsl_matrix_column(BaseData,4);
      gsl_matrix_set_col(Positions,3,&Position_3.vector);
      gsl_vector_view Velocity_0      = gsl_matrix_column(BaseData,5);
      gsl_matrix_set_col(Velocities,0,&Velocity_0.vector);
      gsl_vector_view Velocity_1      = gsl_matrix_column(BaseData,6);
      gsl_matrix_set_col(Velocities,1,&Velocity_1.vector);
      gsl_vector_view Velocity_2      = gsl_matrix_column(BaseData,7);
      gsl_matrix_set_col(Velocities,2,&Velocity_2.vector);

    fclose(BaseFile);

    // There are some positions coordinates  that are outside the box lammps can
    // deal  with this issue without  major problems.  Here,  we use  PBC to put
    // all the atom inside the box
    PrintMsg("Fixing PBC in the positions file...");
    FixPBC(Positions);

    #if __COMPUTE_MOMENTUM__
      Compute_Momentum(Positions,Velocities,Momentum);
    #endif
    
    PrintMsg("Obtaining linked list...");
    Compute_Linked_List(Positions, List, ListHead);

    PrintMsg("Computing forces in the fluid (type 2 particles) due to the wall (type 1 particles)");
    Compute_Forces(Positions, Velocities, Neighbors, ListHead, List, 2, 1, Force, Energy, Kinetic);
    
    //  Checkpoint: Print the force exerted on type2 particles
    //              and the energy of all the particles
    // 
    //     gsl_vector_view  zPart = gsl_matrix_column(Positions,3);
    //     gsl_vector_view FzPart = gsl_matrix_column(Force,2);
    // 
    //     PrintInfo(Step, &zPart.vector,  oFile.MicrozForce);
    //     PrintInfo(Step, &FzPart.vector, oFile.MicrozForce);
    //     PrintInfo(Step, &zPart.vector,  oFile.MicroEnergy);
    //     PrintInfo(Step, Energy,         oFile.MicroEnergy);
    //     PrintInfo(Step, &zPart.vector,  oFile.MicroKinetic);
    //     PrintInfo(Step, Kinetic,        oFile.MicroKinetic);
    //     
    //     PrintMsg("Computing the module of the velocity as a estimator for the temperature...");
    //     gsl_vector * Vmod = gsl_vector_calloc(NParticles);
    //     Compute_Velocity_Module(Velocities, Vmod);
    //     PrintInfo(Step, &zPart.vector, oFile.MicroVmod);
    //     PrintInfo(Step, Vmod,          oFile.MicroVmod);
    //     gsl_vector_free(Vmod);

    // MESOSCOPIC INFORMATION
    
    #pragma omp parallel sections num_threads(NPROC)
    {
      #if __COMPUTE_DENSITY__
      #pragma omp section
      {
        PrintMsg("Obtaining node densities...");
        // Perform the calculation
        Compute_Meso_Density(Positions, z, MesoDensity);

        // Define a few views to print vectors
        gsl_vector_view MesoDensity_1 = gsl_matrix_column(MesoDensity,0);
        gsl_vector_view MesoDensity_2 = gsl_matrix_column(MesoDensity,1);

        // Sum the two types of densities
        gsl_vector_memcpy(TotalMesoDensity, &MesoDensity_1.vector);
        gsl_vector_add   (TotalMesoDensity, &MesoDensity_2.vector);

        // Print result
        #if __BINARY_OUTPUT__
          gsl_vector_view MesoDensity0Row = gsl_matrix_row(BinaryMesoDensity_0,Step);
          gsl_vector_view MesoDensity1Row = gsl_matrix_row(BinaryMesoDensity_1,Step);
          gsl_vector_view MesoDensity2Row = gsl_matrix_row(BinaryMesoDensity_2,Step);
          
          gsl_vector_memcpy(&MesoDensity0Row.vector,     TotalMesoDensity);
          gsl_vector_memcpy(&MesoDensity1Row.vector,&MesoDensity_1.vector);
          gsl_vector_memcpy(&MesoDensity2Row.vector,&MesoDensity_2.vector);
        #else
          PrintInfo(Step, &MesoDensity_1.vector, oFile.MesoDensity_1);
          PrintInfo(Step, &MesoDensity_2.vector, oFile.MesoDensity_2);
          PrintInfo(Step,      TotalMesoDensity, oFile.MesoDensity_0);
        #endif
      }
      #endif
      #if __COMPUTE_FORCE__
      #pragma omp section
      {
        PrintMsg("Obtaining node forces...");
        // Perform the calculation
        Compute_Meso_Force(Positions, Force, z, MesoForce);

        // Define a few views to print vectors
        gsl_vector_view MesoForce_x = gsl_matrix_column(MesoForce,0);
        gsl_vector_view MesoForce_y = gsl_matrix_column(MesoForce,1);
        gsl_vector_view MesoForce_z = gsl_matrix_column(MesoForce,2);

        // Print result
        #if __BINARY_OUTPUT__
          gsl_vector_view MesoForcexRow = gsl_matrix_row(BinaryMesoForce_x,Step);
          gsl_vector_view MesoForceyRow = gsl_matrix_row(BinaryMesoForce_y,Step);
          gsl_vector_view MesoForcezRow = gsl_matrix_row(BinaryMesoForce_z,Step);
          
          gsl_vector_memcpy(&MesoForcexRow.vector,&MesoForce_x.vector);
          gsl_vector_memcpy(&MesoForceyRow.vector,&MesoForce_y.vector);
          gsl_vector_memcpy(&MesoForcezRow.vector,&MesoForce_z.vector);
        #else
          PrintInfo(Step, &MesoForce_x.vector, oFile.MesoxForce);
          PrintInfo(Step, &MesoForce_y.vector, oFile.MesoyForce);
          PrintInfo(Step, &MesoForce_z.vector, oFile.MesozForce);
        #endif
      }
      #endif
      #if __COMPUTE_ENERGY__
      #pragma omp section
      {
        PrintMsg("Obtaining node energies...");
        #if __BINARY_OUTPUT__
          gsl_vector_view MesoEnergyRow = gsl_matrix_row(BinaryMesoEnergy,Step);
          Compute_Meso_Profile(Positions, Energy, z, &MesoEnergyRow.vector, 2);
        #else
          Compute_Meso_Profile(Positions, Energy, z, MesoEnergy, 2);
          PrintInfo(Step, MesoEnergy, oFile.MesoEnergy);
        #endif
      }
      #endif
      #if __COMPUTE_MOMENTUM__
      #pragma omp section
      {
        PrintMsg("Obtaining node momentum...");
        // Microscopic views
        gsl_vector_view Momentum_x      = gsl_matrix_column(Momentum,0);
        gsl_vector_view Momentum_y      = gsl_matrix_column(Momentum,1);
        gsl_vector_view Momentum_z      = gsl_matrix_column(Momentum,2);
        
        #if __BINARY_OUTPUT__
          gsl_vector_view MesoMomentumxRow = gsl_matrix_row(BinaryMesoMomentum_x,Step);
          gsl_vector_view MesoMomentumyRow = gsl_matrix_row(BinaryMesoMomentum_y,Step);
          gsl_vector_view MesoMomentumzRow = gsl_matrix_row(BinaryMesoMomentum_z,Step);

          Compute_Meso_Profile(Positions, &Momentum_x.vector, z, &MesoMomentumxRow.vector, 2);
          Compute_Meso_Profile(Positions, &Momentum_y.vector, z, &MesoMomentumyRow.vector, 2);
          Compute_Meso_Profile(Positions, &Momentum_z.vector, z, &MesoMomentumzRow.vector, 2);
        #else
          // Mesoscopic views
          gsl_vector_view MesoMomentum_x  = gsl_matrix_column(MesoMomentum,0);
          gsl_vector_view MesoMomentum_y  = gsl_matrix_column(MesoMomentum,1);
          gsl_vector_view MesoMomentum_z  = gsl_matrix_column(MesoMomentum,2);

          // Perform the computation
          Compute_Meso_Profile(Positions, &Momentum_x.vector, z, &MesoMomentum_x.vector, 2);
          Compute_Meso_Profile(Positions, &Momentum_y.vector, z, &MesoMomentum_y.vector, 2);
          Compute_Meso_Profile(Positions, &Momentum_z.vector, z, &MesoMomentum_z.vector, 2);

          // Print info
          PrintInfo(Step, &MesoMomentum_x.vector, oFile.MesoMomentum_x);
          PrintInfo(Step, &MesoMomentum_y.vector, oFile.MesoMomentum_y);
          PrintInfo(Step, &MesoMomentum_z.vector, oFile.MesoMomentum_z);
        #endif
      }
      #endif
      #if __COMPUTE_VELOCITY__
      #pragma omp section
      {
        PrintMsg("Obtaining node velocity...");
        // Mesoscopic views
        gsl_vector_view MesoDensity_2  = gsl_matrix_column( MesoDensity,1);
       
        // Perform the computation
        Compute_Meso_Velocity(MesoMomentum,&MesoDensity_2.vector,MesoVelocity);
        
        // Mesoscopic views
        gsl_vector_view MesoVelocity_x = gsl_matrix_column(MesoVelocity,0);
        gsl_vector_view MesoVelocity_y = gsl_matrix_column(MesoVelocity,1);
        gsl_vector_view MesoVelocity_z = gsl_matrix_column(MesoVelocity,2);

        #if __BINARY_OUTPUT__
          gsl_vector_view MesoVelocityxRow = gsl_matrix_row(BinaryMesoVelocity_x,Step);
          gsl_vector_view MesoVelocityyRow = gsl_matrix_row(BinaryMesoVelocity_y,Step);
          gsl_vector_view MesoVelocityzRow = gsl_matrix_row(BinaryMesoVelocity_z,Step);

          gsl_vector_memcpy(&MesoVelocityxRow.vector,&MesoVelocity_x.vector);
          gsl_vector_memcpy(&MesoVelocityyRow.vector,&MesoVelocity_y.vector);
          gsl_vector_memcpy(&MesoVelocityzRow.vector,&MesoVelocity_z.vector);
        #else
          // Print info
          PrintInfo(Step, &MesoVelocity_x.vector, oFile.MesoVelocity_x);
          PrintInfo(Step, &MesoVelocity_y.vector, oFile.MesoVelocity_y);
          PrintInfo(Step, &MesoVelocity_z.vector, oFile.MesoVelocity_z);
        #endif
      }
      #endif
      #if __COMPUTE_ENERGY__
      #pragma omp section
      {
        PrintMsg("Obtaining node kinetic energies...");
        #if __BINARY_OUTPUT__
          gsl_vector_view MesoKineticRow = gsl_matrix_row(BinaryMesoKinetic,Step);
          Compute_Meso_Profile(Positions, Kinetic, z, &MesoKineticRow.vector, 2);
        #else
          Compute_Meso_Profile(Positions, Kinetic, z, MesoKinetic, 2);
          PrintInfo(Step, MesoKinetic, oFile.MesoKinetic);
        #endif
      }
      #endif
      #if __COMPUTE_STRESS__
      #pragma omp section
      {
        PrintMsg("Obtaining node kinetic stress tensors...");
        // Perform the computation
        Compute_Meso_Sigma1(Positions, Velocities, MesoSigma1);

        // MesoSigma = MesoSigma1 + MesoSigma2
        // This is the first part
        gsl_matrix_memcpy(MesoSigma,MesoSigma1);

        // Vector views
        gsl_vector_view  MesoSigma1_00 = gsl_matrix_column(MesoSigma1,0);
        gsl_vector_view  MesoSigma1_01 = gsl_matrix_column(MesoSigma1,1);
        gsl_vector_view  MesoSigma1_02 = gsl_matrix_column(MesoSigma1,2);
        gsl_vector_view  MesoSigma1_10 = gsl_matrix_column(MesoSigma1,3);
        gsl_vector_view  MesoSigma1_11 = gsl_matrix_column(MesoSigma1,4);
        gsl_vector_view  MesoSigma1_12 = gsl_matrix_column(MesoSigma1,5);
        gsl_vector_view  MesoSigma1_20 = gsl_matrix_column(MesoSigma1,6);
        gsl_vector_view  MesoSigma1_21 = gsl_matrix_column(MesoSigma1,7);
        gsl_vector_view  MesoSigma1_22 = gsl_matrix_column(MesoSigma1,8);
     
        #if __BINARY_OUTPUT__
          gsl_vector_view MesoSigma1_00Row = gsl_matrix_row(BinaryMesoSigma1_xx,Step);
          gsl_vector_view MesoSigma1_01Row = gsl_matrix_row(BinaryMesoSigma1_xy,Step);
          gsl_vector_view MesoSigma1_02Row = gsl_matrix_row(BinaryMesoSigma1_xz,Step);
          gsl_vector_view MesoSigma1_10Row = gsl_matrix_row(BinaryMesoSigma1_yx,Step);
          gsl_vector_view MesoSigma1_11Row = gsl_matrix_row(BinaryMesoSigma1_yy,Step);
          gsl_vector_view MesoSigma1_12Row = gsl_matrix_row(BinaryMesoSigma1_yz,Step);
          gsl_vector_view MesoSigma1_20Row = gsl_matrix_row(BinaryMesoSigma1_zx,Step);
          gsl_vector_view MesoSigma1_21Row = gsl_matrix_row(BinaryMesoSigma1_zy,Step);
          gsl_vector_view MesoSigma1_22Row = gsl_matrix_row(BinaryMesoSigma1_zz,Step);

          gsl_vector_memcpy(&MesoSigma1_00Row.vector,&MesoSigma1_00.vector);
          gsl_vector_memcpy(&MesoSigma1_01Row.vector,&MesoSigma1_01.vector);
          gsl_vector_memcpy(&MesoSigma1_02Row.vector,&MesoSigma1_02.vector);
          gsl_vector_memcpy(&MesoSigma1_10Row.vector,&MesoSigma1_10.vector);
          gsl_vector_memcpy(&MesoSigma1_11Row.vector,&MesoSigma1_11.vector);
          gsl_vector_memcpy(&MesoSigma1_12Row.vector,&MesoSigma1_12.vector);
          gsl_vector_memcpy(&MesoSigma1_20Row.vector,&MesoSigma1_20.vector);
          gsl_vector_memcpy(&MesoSigma1_21Row.vector,&MesoSigma1_21.vector);
          gsl_vector_memcpy(&MesoSigma1_22Row.vector,&MesoSigma1_22.vector);
        #else
          // Print info
          PrintInfo(Step, &MesoSigma1_00.vector, oFile.MesoSigma1_00);
          PrintInfo(Step, &MesoSigma1_01.vector, oFile.MesoSigma1_01);
          PrintInfo(Step, &MesoSigma1_02.vector, oFile.MesoSigma1_02);
          PrintInfo(Step, &MesoSigma1_10.vector, oFile.MesoSigma1_10);
          PrintInfo(Step, &MesoSigma1_11.vector, oFile.MesoSigma1_11);
          PrintInfo(Step, &MesoSigma1_12.vector, oFile.MesoSigma1_12);
          PrintInfo(Step, &MesoSigma1_20.vector, oFile.MesoSigma1_20);
          PrintInfo(Step, &MesoSigma1_21.vector, oFile.MesoSigma1_21);
          PrintInfo(Step, &MesoSigma1_22.vector, oFile.MesoSigma1_22);
        #endif
      }
      #endif
    }

    #if __COMPUTE_STRESS__
      PrintMsg("Obtaining node virial stress tensor...");
      // Perform the computation
      Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, MesoSigma2, z);
      // MesoSigma = MesoSigma1 + MesoSigma2
      // This is the second computation
      gsl_matrix_add (MesoSigma, MesoSigma2);

      // Vector views
      gsl_vector_view  MesoSigma2_00 = gsl_matrix_column(MesoSigma2,0);
      gsl_vector_view  MesoSigma2_01 = gsl_matrix_column(MesoSigma2,1);
      gsl_vector_view  MesoSigma2_02 = gsl_matrix_column(MesoSigma2,2);
      gsl_vector_view  MesoSigma2_10 = gsl_matrix_column(MesoSigma2,3);
      gsl_vector_view  MesoSigma2_11 = gsl_matrix_column(MesoSigma2,4);
      gsl_vector_view  MesoSigma2_12 = gsl_matrix_column(MesoSigma2,5);
      gsl_vector_view  MesoSigma2_20 = gsl_matrix_column(MesoSigma2,6);
      gsl_vector_view  MesoSigma2_21 = gsl_matrix_column(MesoSigma2,7);
      gsl_vector_view  MesoSigma2_22 = gsl_matrix_column(MesoSigma2,8);
      
      #if __BINARY_OUTPUT__
        gsl_vector_view MesoSigma2_00Row = gsl_matrix_row(BinaryMesoSigma2_xx,Step);
        gsl_vector_view MesoSigma2_01Row = gsl_matrix_row(BinaryMesoSigma2_xy,Step);
        gsl_vector_view MesoSigma2_02Row = gsl_matrix_row(BinaryMesoSigma2_xz,Step);
        gsl_vector_view MesoSigma2_10Row = gsl_matrix_row(BinaryMesoSigma2_yx,Step);
        gsl_vector_view MesoSigma2_11Row = gsl_matrix_row(BinaryMesoSigma2_yy,Step);
        gsl_vector_view MesoSigma2_12Row = gsl_matrix_row(BinaryMesoSigma2_yz,Step);
        gsl_vector_view MesoSigma2_20Row = gsl_matrix_row(BinaryMesoSigma2_zx,Step);
        gsl_vector_view MesoSigma2_21Row = gsl_matrix_row(BinaryMesoSigma2_zy,Step);
        gsl_vector_view MesoSigma2_22Row = gsl_matrix_row(BinaryMesoSigma2_zz,Step);

        gsl_vector_memcpy(&MesoSigma2_00Row.vector,&MesoSigma2_00.vector);
        gsl_vector_memcpy(&MesoSigma2_01Row.vector,&MesoSigma2_01.vector);
        gsl_vector_memcpy(&MesoSigma2_02Row.vector,&MesoSigma2_02.vector);
        gsl_vector_memcpy(&MesoSigma2_10Row.vector,&MesoSigma2_10.vector);
        gsl_vector_memcpy(&MesoSigma2_11Row.vector,&MesoSigma2_11.vector);
        gsl_vector_memcpy(&MesoSigma2_12Row.vector,&MesoSigma2_12.vector);
        gsl_vector_memcpy(&MesoSigma2_20Row.vector,&MesoSigma2_20.vector);
        gsl_vector_memcpy(&MesoSigma2_21Row.vector,&MesoSigma2_21.vector);
        gsl_vector_memcpy(&MesoSigma2_22Row.vector,&MesoSigma2_22.vector);
      #else
        // Print info
        PrintInfo(Step, &MesoSigma2_00.vector, oFile.MesoSigma2_00);
        PrintInfo(Step, &MesoSigma2_01.vector, oFile.MesoSigma2_01);
        PrintInfo(Step, &MesoSigma2_02.vector, oFile.MesoSigma2_02);
        PrintInfo(Step, &MesoSigma2_10.vector, oFile.MesoSigma2_10);
        PrintInfo(Step, &MesoSigma2_11.vector, oFile.MesoSigma2_11);
        PrintInfo(Step, &MesoSigma2_12.vector, oFile.MesoSigma2_12);
        PrintInfo(Step, &MesoSigma2_20.vector, oFile.MesoSigma2_20);
        PrintInfo(Step, &MesoSigma2_21.vector, oFile.MesoSigma2_21);
        PrintInfo(Step, &MesoSigma2_22.vector, oFile.MesoSigma2_22);
      #endif

      PrintMsg("Saving stress tensors...");
      // More views
      gsl_vector_view  MesoSigma_00 = gsl_matrix_column(MesoSigma,0);
      gsl_vector_view  MesoSigma_01 = gsl_matrix_column(MesoSigma,1);
      gsl_vector_view  MesoSigma_02 = gsl_matrix_column(MesoSigma,2);
      gsl_vector_view  MesoSigma_10 = gsl_matrix_column(MesoSigma,3);
      gsl_vector_view  MesoSigma_11 = gsl_matrix_column(MesoSigma,4);
      gsl_vector_view  MesoSigma_12 = gsl_matrix_column(MesoSigma,5);
      gsl_vector_view  MesoSigma_20 = gsl_matrix_column(MesoSigma,6);
      gsl_vector_view  MesoSigma_21 = gsl_matrix_column(MesoSigma,7);
      gsl_vector_view  MesoSigma_22 = gsl_matrix_column(MesoSigma,8);
     
      #if __BINARY_OUTPUT__
        gsl_vector_view MesoSigma_00Row = gsl_matrix_row(BinaryMesoSigma_xx,Step);
        gsl_vector_view MesoSigma_01Row = gsl_matrix_row(BinaryMesoSigma_xy,Step);
        gsl_vector_view MesoSigma_02Row = gsl_matrix_row(BinaryMesoSigma_xz,Step);
        gsl_vector_view MesoSigma_10Row = gsl_matrix_row(BinaryMesoSigma_yx,Step);
        gsl_vector_view MesoSigma_11Row = gsl_matrix_row(BinaryMesoSigma_yy,Step);
        gsl_vector_view MesoSigma_12Row = gsl_matrix_row(BinaryMesoSigma_yz,Step);
        gsl_vector_view MesoSigma_20Row = gsl_matrix_row(BinaryMesoSigma_zx,Step);
        gsl_vector_view MesoSigma_21Row = gsl_matrix_row(BinaryMesoSigma_zy,Step);
        gsl_vector_view MesoSigma_22Row = gsl_matrix_row(BinaryMesoSigma_zz,Step);

        gsl_vector_memcpy(&MesoSigma_00Row.vector,&MesoSigma_00.vector);
        gsl_vector_memcpy(&MesoSigma_01Row.vector,&MesoSigma_01.vector);
        gsl_vector_memcpy(&MesoSigma_02Row.vector,&MesoSigma_02.vector);
        gsl_vector_memcpy(&MesoSigma_10Row.vector,&MesoSigma_10.vector);
        gsl_vector_memcpy(&MesoSigma_11Row.vector,&MesoSigma_11.vector);
        gsl_vector_memcpy(&MesoSigma_12Row.vector,&MesoSigma_12.vector);
        gsl_vector_memcpy(&MesoSigma_20Row.vector,&MesoSigma_20.vector);
        gsl_vector_memcpy(&MesoSigma_21Row.vector,&MesoSigma_21.vector);
        gsl_vector_memcpy(&MesoSigma_22Row.vector,&MesoSigma_22.vector);
      #else
        // Print info
        PrintInfo(Step, &MesoSigma_00.vector, oFile.MesoSigma_00);
        PrintInfo(Step, &MesoSigma_01.vector, oFile.MesoSigma_01);
        PrintInfo(Step, &MesoSigma_02.vector, oFile.MesoSigma_02);
        PrintInfo(Step, &MesoSigma_10.vector, oFile.MesoSigma_10);
        PrintInfo(Step, &MesoSigma_11.vector, oFile.MesoSigma_11);
        PrintInfo(Step, &MesoSigma_12.vector, oFile.MesoSigma_12);
        PrintInfo(Step, &MesoSigma_20.vector, oFile.MesoSigma_20);
        PrintInfo(Step, &MesoSigma_21.vector, oFile.MesoSigma_21);
        PrintInfo(Step, &MesoSigma_22.vector, oFile.MesoSigma_22);
      #endif
    #endif

    #if __COMPUTE_TEMPERATURE__
      PrintMsg("Obtaining node temperature...");
      gsl_vector_view MesoDensity_2 = gsl_matrix_column(MesoDensity,1);
      #if __BINARY_OUTPUT__
        gsl_vector_view MesoTempRow = gsl_matrix_row(BinaryMesoTemp,Step);
        Compute_Meso_Temp(MesoKinetic, &MesoDensity_2.vector, &MesoTempRow.vector);
      #else
        Compute_Meso_Temp(MesoKinetic, &MesoDensity_2.vector, MesoTemp);
        PrintInfo(Step, MesoTemp, oFile.MesoTemp);
      #endif
    #endif
        
    #if __COMPUTE_INTERNAL_ENERGY__
      PrintMsg("Obtaining node internal energies...");
      #if __BINARY_OUTPUT__
        gsl_vector_view MesoInternalEnergyRow = gsl_matrix_row(BinaryMesoInternalEnergy,Step);
        Compute_InternalEnergy(MesoEnergy, MesoMomentum, &MesoDensity_2.vector, &MesoInternalEnergyRow.vector);
      #else
        Compute_InternalEnergy(MesoEnergy, MesoMomentum, &MesoDensity_2.vector, MesoInternalEnergy);
        PrintInfo(Step, MesoInternalEnergy, oFile.MesoInternalEnergy);
      #endif
    #endif
    
    // MACROSCOPIC INFORMATION

    #if __COMPUTE_MACRO_ENERGY__
      double MacroEnergy;
      
      PrintMsg("Computing the energy of upper wall");
      MacroEnergy = Compute_Macro(Energy, Positions, 1, "top");
      #if __BINARY_OUTPUT__
        gsl_vector_set(BinaryMacroEnergyUpperWall,Step,MacroEnergy);
      #else
        PrintScalarWithIndex(Step, MacroEnergy, oFile.MacroEnergyUpperWall); 
      #endif

      PrintMsg("Computing the energy of lower wall");
      MacroEnergy = Compute_Macro(Energy, Positions, 1, "bottom");
      #if __BINARY_OUTPUT__
        gsl_vector_set(BinaryMacroEnergyLowerWall,Step,MacroEnergy);
      #else
        PrintScalarWithIndex(Step, MacroEnergy, oFile.MacroEnergyLowerWall); 
      #endif
    #endif

    #if __COMPUTE_MACRO_MOMENTUM__
      double TemporalMomentum;

      #if __BINARY_OUTPUT == false
        gsl_vector * MacroMomentum = gsl_vector_calloc(3);
      #endif

      gsl_vector_view Momentum_x = gsl_matrix_column(Momentum,0);
      gsl_vector_view Momentum_y = gsl_matrix_column(Momentum,1);
      gsl_vector_view Momentum_z = gsl_matrix_column(Momentum,2);
      
      PrintMsg("Computing the momentum of upper wall");
      
      TemporalMomentum = Compute_Macro(&Momentum_x.vector, Positions, 1, "top"); 
      #if __BINARY_OUTPUT__
        gsl_vector_set(BinaryMacroMomentumUpperWallx,Step,TemporalMomentum);
      #else
        gsl_vector_set(MacroMomentum,0,TemporalMomentum);
      #endif

      TemporalMomentum = Compute_Macro(&Momentum_y.vector, Positions, 1, "top"); 
      #if __BINARY_OUTPUT__
        gsl_vector_set(BinaryMacroMomentumUpperWally,Step,TemporalMomentum);
      #else
        gsl_vector_set(MacroMomentum,1,TemporalMomentum);
      #endif

      TemporalMomentum = Compute_Macro(&Momentum_z.vector, Positions, 1, "top"); 
      #if __BINARY_OUTPUT__
        gsl_vector_set(BinaryMacroMomentumUpperWallz,Step,TemporalMomentum);
      #else
        gsl_vector_set(MacroMomentum,2,TemporalMomentum);
      #endif
    
      #if __BINARY_OUTPUT == false
        PrintInfo(Step, MacroMomentum, oFile.MacroMomentumUpperWall);
      #endif
    
      PrintMsg("Computing the momentum of lower wall");

      TemporalMomentum = Compute_Macro(&Momentum_x.vector, Positions, 1, "bottom"); 
      #if __BINARY_OUTPUT__
        gsl_vector_set(BinaryMacroMomentumLowerWallx,Step,TemporalMomentum);
      #else
        gsl_vector_set(MacroMomentum,0,TemporalMomentum);
      #endif
      
      TemporalMomentum = Compute_Macro(&Momentum_y.vector, Positions, 1, "bottom"); 
      #if __BINARY_OUTPUT__
        gsl_vector_set(BinaryMacroMomentumLowerWally,Step,TemporalMomentum);
      #else
        gsl_vector_set(MacroMomentum,1,TemporalMomentum);
      #endif
      
      TemporalMomentum = Compute_Macro(&Momentum_z.vector, Positions, 1, "bottom"); 
      #if __BINARY_OUTPUT__
        gsl_vector_set(BinaryMacroMomentumLowerWallz,Step,TemporalMomentum);
      #else
        gsl_vector_set(MacroMomentum,2,TemporalMomentum);
      #endif

      #if __BINARY_OUTPUT == false
        PrintInfo(Step, MacroMomentum, oFile.MacroMomentumLowerWall);
      #endif
    
      #if __BINARY_OUTPUT == false
        gsl_vector_free(MacroMomentum);
      #endif

    #endif

    #if __COMPUTE_CENTER_OF_MASS__
      gsl_vector * CenterOfMass = gsl_vector_calloc(4);
      
      PrintMsg("Computing Center of Mass upper wall");
      Compute_CenterOfMass(Positions, 1, "top", CenterOfMass);
      PrintInfo(Step, CenterOfMass, oFile.CenterOfMassUpperWall);
      
      PrintMsg("Computing Center of Mass lower wall");
      Compute_CenterOfMass(Positions, 1, "bottom", CenterOfMass);
      PrintInfo(Step, CenterOfMass, oFile.CenterOfMassLowerWall);
      
      gsl_vector_free(CenterOfMass);
    #endif
  }

  #if __BINARY_OUTPUT__
    #if __COMPUTE_DENSITY__
      gsl_matrix_fwrite(oFile.MesoDensity_0,BinaryMesoDensity_0);
      gsl_matrix_fwrite(oFile.MesoDensity_1,BinaryMesoDensity_1);
      gsl_matrix_fwrite(oFile.MesoDensity_2,BinaryMesoDensity_2);
    #endif

    #if __COMPUTE_FORCE__
      gsl_matrix_fwrite(oFile.MesoxForce,BinaryMesoForce_x);
      gsl_matrix_fwrite(oFile.MesoyForce,BinaryMesoForce_y);
      gsl_matrix_fwrite(oFile.MesozForce,BinaryMesoForce_z);
    #endif

    #if __COMPUTE_ENERGY__
      gsl_matrix_fwrite(oFile.MesoEnergy,  BinaryMesoEnergy);
      gsl_matrix_fwrite(oFile.MesoKinetic,BinaryMesoKinetic);
    #endif

    #if __COMPUTE_TEMPERATURE__
      gsl_matrix_fwrite(oFile.MesoTemp,  BinaryMesoTemp);
    #endif

    #if __COMPUTE_STRESS__
      gsl_matrix_fwrite(oFile.MesoSigma1_00,BinaryMesoSigma1_xx);
      gsl_matrix_fwrite(oFile.MesoSigma1_01,BinaryMesoSigma1_xy);
      gsl_matrix_fwrite(oFile.MesoSigma1_02,BinaryMesoSigma1_xz);
      gsl_matrix_fwrite(oFile.MesoSigma1_10,BinaryMesoSigma1_yx);
      gsl_matrix_fwrite(oFile.MesoSigma1_11,BinaryMesoSigma1_yy);
      gsl_matrix_fwrite(oFile.MesoSigma1_12,BinaryMesoSigma1_yz);
      gsl_matrix_fwrite(oFile.MesoSigma1_20,BinaryMesoSigma1_zx);
      gsl_matrix_fwrite(oFile.MesoSigma1_21,BinaryMesoSigma1_zy);
      gsl_matrix_fwrite(oFile.MesoSigma1_22,BinaryMesoSigma1_zz);
      
      gsl_matrix_fwrite(oFile.MesoSigma2_00,BinaryMesoSigma2_xx);
      gsl_matrix_fwrite(oFile.MesoSigma2_01,BinaryMesoSigma2_xy);
      gsl_matrix_fwrite(oFile.MesoSigma2_02,BinaryMesoSigma2_xz);
      gsl_matrix_fwrite(oFile.MesoSigma2_10,BinaryMesoSigma2_yx);
      gsl_matrix_fwrite(oFile.MesoSigma2_11,BinaryMesoSigma2_yy);
      gsl_matrix_fwrite(oFile.MesoSigma2_12,BinaryMesoSigma2_yz);
      gsl_matrix_fwrite(oFile.MesoSigma2_20,BinaryMesoSigma2_zx);
      gsl_matrix_fwrite(oFile.MesoSigma2_21,BinaryMesoSigma2_zy);
      gsl_matrix_fwrite(oFile.MesoSigma2_22,BinaryMesoSigma2_zz);
      
      gsl_matrix_fwrite(oFile.MesoSigma_00,BinaryMesoSigma_xx);
      gsl_matrix_fwrite(oFile.MesoSigma_01,BinaryMesoSigma_xy);
      gsl_matrix_fwrite(oFile.MesoSigma_02,BinaryMesoSigma_xz);
      gsl_matrix_fwrite(oFile.MesoSigma_10,BinaryMesoSigma_yx);
      gsl_matrix_fwrite(oFile.MesoSigma_11,BinaryMesoSigma_yy);
      gsl_matrix_fwrite(oFile.MesoSigma_12,BinaryMesoSigma_yz);
      gsl_matrix_fwrite(oFile.MesoSigma_20,BinaryMesoSigma_zx);
      gsl_matrix_fwrite(oFile.MesoSigma_21,BinaryMesoSigma_zy);
      gsl_matrix_fwrite(oFile.MesoSigma_22,BinaryMesoSigma_zz);
    #endif

    #if __COMPUTE_MOMENTUM__
      gsl_matrix_fwrite(oFile.MesoMomentum_x,BinaryMesoMomentum_x);
      gsl_matrix_fwrite(oFile.MesoMomentum_y,BinaryMesoMomentum_y);
      gsl_matrix_fwrite(oFile.MesoMomentum_z,BinaryMesoMomentum_z);
    #endif

    #if __COMPUTE_VELOCITY__
      gsl_matrix_fwrite(oFile.MesoVelocity_x,BinaryMesoVelocity_x);
      gsl_matrix_fwrite(oFile.MesoVelocity_y,BinaryMesoVelocity_y);
      gsl_matrix_fwrite(oFile.MesoVelocity_z,BinaryMesoVelocity_z);
    #endif

    #if __COMPUTE_INTERNAL_ENERGY__
      gsl_matrix_fwrite(oFile.MesoInternalEnergy,BinaryMesoInternalEnergy);
    #endif

    #if __COMPUTE_MACRO_ENERGY__
      gsl_vector_fwrite(oFile.MacroEnergyUpperWall,BinaryMacroEnergyUpperWall);
      gsl_vector_fwrite(oFile.MacroEnergyLowerWall,BinaryMacroEnergyLowerWall);
    #endif

    #if __COMPUTE_MACRO_MOMENTUM__
      gsl_vector_fwrite(oFile.MacroMomentumUpperWallx,BinaryMacroMomentumUpperWallx);
      gsl_vector_fwrite(oFile.MacroMomentumUpperWally,BinaryMacroMomentumUpperWally);
      gsl_vector_fwrite(oFile.MacroMomentumUpperWallz,BinaryMacroMomentumUpperWallz);
      gsl_vector_fwrite(oFile.MacroMomentumLowerWallx,BinaryMacroMomentumLowerWallx);
      gsl_vector_fwrite(oFile.MacroMomentumLowerWally,BinaryMacroMomentumLowerWally);
      gsl_vector_fwrite(oFile.MacroMomentumLowerWallz,BinaryMacroMomentumLowerWallz);
    #endif

    #if __COMPUTE_CENTER_OF_MASS__
      gsl_vector_fwrite(oFile.CenterOfMassUpperWall,BinaryCenterOfMassUpperWall);
      gsl_vector_fwrite(oFile.CenterOfMassLowerWall,BinaryCenterOfMassLowerWall);
    #endif

  #endif
  // Close micro files
  // fclose(oFile.MicrozForce);
  // fclose(oFile.MicroEnergy);
  // fclose(oFile.MicroKinetic);
  // fclose(oFile.MicroVmod);

  // Close meso files
  #if __COMPUTE_DENSITY__
  fclose(oFile.MesoDensity_0);
  fclose(oFile.MesoDensity_1);
  fclose(oFile.MesoDensity_2);
  #endif

  #if __COMPUTE_FORCE__
  fclose(oFile.MesoxForce);
  fclose(oFile.MesoyForce);
  fclose(oFile.MesozForce);
  #endif

  #if __COMPUTE_ENERGY__
  fclose(oFile.MesoEnergy);
  fclose(oFile.MesoKinetic);
  #endif

  #if __COMPUTE_TEMPERATURE__
  fclose(oFile.MesoTemp);
  #endif

  #if __COMPUTE_STRESS__
  fclose(oFile.MesoSigma1_00);
  fclose(oFile.MesoSigma1_01);
  fclose(oFile.MesoSigma1_02);
  fclose(oFile.MesoSigma1_10);
  fclose(oFile.MesoSigma1_11);
  fclose(oFile.MesoSigma1_12);
  fclose(oFile.MesoSigma1_20);
  fclose(oFile.MesoSigma1_21);
  fclose(oFile.MesoSigma1_22);
  
  fclose(oFile.MesoSigma2_00);
  fclose(oFile.MesoSigma2_01);
  fclose(oFile.MesoSigma2_02);
  fclose(oFile.MesoSigma2_10);
  fclose(oFile.MesoSigma2_11);
  fclose(oFile.MesoSigma2_12);
  fclose(oFile.MesoSigma2_20);
  fclose(oFile.MesoSigma2_21);
  fclose(oFile.MesoSigma2_22);

  fclose(oFile.MesoSigma_00);
  fclose(oFile.MesoSigma_01);
  fclose(oFile.MesoSigma_02);
  fclose(oFile.MesoSigma_10);
  fclose(oFile.MesoSigma_11);
  fclose(oFile.MesoSigma_12);
  fclose(oFile.MesoSigma_20);
  fclose(oFile.MesoSigma_21);
  fclose(oFile.MesoSigma_22);
  #endif

  #if __COMPUTE_MOMENTUM__
  fclose(oFile.MesoMomentum_x);
  fclose(oFile.MesoMomentum_y);
  fclose(oFile.MesoMomentum_z);
  #endif

  #if __COMPUTE_VELOCITY__
  fclose(oFile.MesoVelocity_x);
  fclose(oFile.MesoVelocity_y);
  fclose(oFile.MesoVelocity_z);
  #endif

  #if __COMPUTE_INTERNAL_ENERGY__
  fclose(oFile.MesoInternalEnergy);
  #endif

  #if __COMPUTE_MACRO_ENERGY__
  fclose(oFile.MacroEnergyUpperWall);
  fclose(oFile.MacroEnergyLowerWall);
  #endif
  
  #if __COMPUTE_MACRO_MOMENTUM__
  fclose(oFile.MacroMomentumUpperWall);
  fclose(oFile.MacroMomentumLowerWall);
  #endif
  
  #if __COMPUTE_CENTER_OF_MASS__
  fclose(oFile.CenterOfMassUpperWall);
  fclose(oFile.CenterOfMassLowerWall);
  #endif

  // SECOND COMPUTATION. OBTAIN MEAN VALUES

  PrintMsg("Computing mean values...");
    
  #pragma omp parallel sections num_threads(NPROC)
  {
    #pragma omp section
    {
      gsl_vector * MesoAverage = gsl_vector_calloc(NNodes);
    
      #if __COMPUTE_DENSITY__
      Compute_Mean_Values(filestr, ".MesoDensity_0.dat",        MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoDensity_0.avg.dat", z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoDensity_1.dat",        MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoDensity_1.avg.dat", z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoDensity_2.dat",        MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoDensity_2.avg.dat", z, MesoAverage);
      #endif

      #if __COMPUTE_FORCE__
      Compute_Mean_Values(filestr, ".MesoForce_x.dat",           MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoForce_x.avg.dat",    z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoForce_y.dat",           MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoForce_y.avg.dat",    z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoForce_z.dat",           MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoForce_z.avg.dat",    z, MesoAverage);
      #endif

      #if __COMPUTE_ENERGY__
      Compute_Mean_Values(filestr, ".MesoKinetic.dat",          MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoKinetic.avg.dat",   z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoEnergy.dat",           MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoEnergy.avg.dat",    z, MesoAverage);
      #endif

      #if __COMPUTE_TEMPERATURE__
      Compute_Mean_Values(filestr, ".MesoTemp.dat",             MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoTemp.avg.dat",      z, MesoAverage);
      #endif
      
      gsl_vector_free(MesoAverage);
    } 
    #if __COMPUTE_STRESS__
    #pragma omp section
    {
      gsl_vector * MesoAverage = gsl_vector_calloc(NNodes);

      // Kinetic stress tensor
      Compute_Mean_Values(filestr, ".MesoSigma1_xx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_xx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma1_xy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_xy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma1_xz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_xz.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma1_yx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_yx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma1_yy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_yy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma1_yz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_yz.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma1_zx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_zx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma1_zy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_zy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma1_zz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma1_zz.avg.dat",  z, MesoAverage);
      
      gsl_vector_free(MesoAverage);
    } 
    #pragma omp section
    {
      gsl_vector * MesoAverage = gsl_vector_calloc(NNodes);
  
      // Virial stress tensor
      Compute_Mean_Values(filestr, ".MesoSigma2_xx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_xx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma2_xy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_xy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma2_xz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_xz.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma2_yx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_yx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma2_yy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_yy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma2_yz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_yz.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma2_zx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_zx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma2_zy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_zy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma2_zz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma2_zz.avg.dat",  z, MesoAverage);
      
      gsl_vector_free(MesoAverage);
    } 
    #pragma omp section
    {
      gsl_vector * MesoAverage = gsl_vector_calloc(NNodes);
  
      // Total stress tensor
      Compute_Mean_Values(filestr, ".MesoSigma_xx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_xx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma_xy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_xy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma_xz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_xz.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma_yx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_yx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma_yy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_yy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma_yz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_yz.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma_zx.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_zx.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma_zy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_zy.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoSigma_zz.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoSigma_zz.avg.dat",  z, MesoAverage);
      
      gsl_vector_free(MesoAverage);
    } 
    #endif
    #pragma omp section
    {
      gsl_vector * MesoAverage = gsl_vector_calloc(NNodes);

      #if __COMPUTE_MOMENTUM__
      Compute_Mean_Values(filestr, ".MesoMomentum_x.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoMomentum_x.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoMomentum_y.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoMomentum_y.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoMomentum_z.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoMomentum_z.avg.dat",  z, MesoAverage);
      #endif

      #if __COMPUTE_VELOCITY__
      Compute_Mean_Values(filestr, ".MesoVelocity_x.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoVelocity_x.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoVelocity_y.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoVelocity_y.avg.dat",  z, MesoAverage);
      Compute_Mean_Values(filestr, ".MesoVelocity_z.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoVelocity_z.avg.dat",  z, MesoAverage);
      #endif

      #if __COMPUTE_INTERNAL_ENERGY__
      Compute_Mean_Values(filestr, ".MesoInternalEnergy.dat",         MesoAverage);
      SaveVectorWithIndex(filestr, ".MesoInternalEnergy.avg.dat",  z, MesoAverage);
      #endif

      gsl_vector_free(MesoAverage);
    }
  }
  // END OF BLOCK. COMPUTATION DONE

  // BEGIN OF BLOCK. FREE MEM
  
  // Free micro vectors and matrices
  gsl_matrix_free(Force);
  gsl_vector_free(Energy);
  gsl_vector_free(Kinetic);
  
  // gsl_vector_free(zPart);
  gsl_vector_free(z);

  gsl_vector_free(List);
  gsl_vector_free(ListHead);
  gsl_matrix_free(Neighbors);
  
  gsl_matrix_free(Positions);
  gsl_matrix_free(Velocities);
  gsl_matrix_free(Momentum);
  gsl_matrix_free(BaseData);
   
  // Free meso vectors and matrices
  gsl_vector_free(TotalMesoDensity);
  gsl_matrix_free(MesoDensity);
  gsl_matrix_free(MesoForce);
  gsl_vector_free(MesoEnergy);
  gsl_vector_free(MesoKinetic);
  gsl_vector_free(MesoTemp);
  gsl_matrix_free(MesoSigma1);
  gsl_matrix_free(MesoSigma2);
  gsl_matrix_free(MesoSigma);
  
  gsl_matrix_free(MesoMomentum);
  gsl_matrix_free(MesoVelocity);
  gsl_vector_free(MesoInternalEnergy);

  #if __BINARY_OUTPUT__
      gsl_matrix_free(BinaryMesoDensity_0); 
      gsl_matrix_free(BinaryMesoDensity_1); 
      gsl_matrix_free(BinaryMesoDensity_2); 
      gsl_matrix_free(BinaryMesoForce_x); 
      gsl_matrix_free(BinaryMesoForce_y); 
      gsl_matrix_free(BinaryMesoForce_z); 
      gsl_matrix_free(BinaryMesoEnergy); 
      gsl_matrix_free(BinaryMesoKinetic); 
      gsl_matrix_free(BinaryMesoTemp); 
      gsl_matrix_free(BinaryMesoSigma_xx); 
      gsl_matrix_free(BinaryMesoSigma_xy); 
      gsl_matrix_free(BinaryMesoSigma_xz); 
      gsl_matrix_free(BinaryMesoSigma_yx); 
      gsl_matrix_free(BinaryMesoSigma_yy); 
      gsl_matrix_free(BinaryMesoSigma_yz); 
      gsl_matrix_free(BinaryMesoSigma_zx); 
      gsl_matrix_free(BinaryMesoSigma_zy); 
      gsl_matrix_free(BinaryMesoSigma_zz); 
      gsl_matrix_free(BinaryMesoSigma1_xx); 
      gsl_matrix_free(BinaryMesoSigma1_xy); 
      gsl_matrix_free(BinaryMesoSigma1_xz); 
      gsl_matrix_free(BinaryMesoSigma1_yx); 
      gsl_matrix_free(BinaryMesoSigma1_yy); 
      gsl_matrix_free(BinaryMesoSigma1_yz); 
      gsl_matrix_free(BinaryMesoSigma1_zx); 
      gsl_matrix_free(BinaryMesoSigma1_zy); 
      gsl_matrix_free(BinaryMesoSigma1_zz); 
      gsl_matrix_free(BinaryMesoSigma2_xx); 
      gsl_matrix_free(BinaryMesoSigma2_xy); 
      gsl_matrix_free(BinaryMesoSigma2_xz); 
      gsl_matrix_free(BinaryMesoSigma2_yx); 
      gsl_matrix_free(BinaryMesoSigma2_yy); 
      gsl_matrix_free(BinaryMesoSigma2_yz); 
      gsl_matrix_free(BinaryMesoSigma2_zx); 
      gsl_matrix_free(BinaryMesoSigma2_zy); 
      gsl_matrix_free(BinaryMesoSigma2_zz); 
      gsl_matrix_free(BinaryMesoMomentum_x);  
      gsl_matrix_free(BinaryMesoMomentum_y); 
      gsl_matrix_free(BinaryMesoMomentum_z); 
      gsl_matrix_free(BinaryMesoVelocity_x); 
      gsl_matrix_free(BinaryMesoVelocity_y); 
      gsl_matrix_free(BinaryMesoVelocity_z); 
      gsl_matrix_free(BinaryMesoInternalEnergy); 
      gsl_vector_free(BinaryMacroEnergyUpperWall); 
      gsl_vector_free(BinaryMacroEnergyLowerWall); 
      gsl_vector_free(BinaryMacroMomentumUpperWallx); 
      gsl_vector_free(BinaryMacroMomentumUpperWally); 
      gsl_vector_free(BinaryMacroMomentumUpperWallz); 
      gsl_vector_free(BinaryMacroMomentumLowerWallx); 
      gsl_vector_free(BinaryMacroMomentumLowerWally); 
      gsl_vector_free(BinaryMacroMomentumLowerWallz); 
      gsl_vector_free(BinaryCenterOfMassUpperWall); 
      gsl_vector_free(BinaryCenterOfMassLowerWall); 
  #endif

  // END OF BLOCK. MEM FREE

  // LoadInfo();
  
  PrintMsg("EOF. Have a nice day.");
  return 0;
}
