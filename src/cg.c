/*
 * Filename   : cg.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : mié 15 jun 2016 17:26:57 CEST
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
  // system("if [ ! -d output ]; then mkdir output; fi");
  struct stat status;
  if (stat("./output", &status) != 0) // && S_ISDIR(status.st_mode)))
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
  #endif

  #if __COMPUTE_FORCE__
    sprintf(str, "./output/%s.MesoForce_x.dat", filestr);
    oFile.MesoxForce = fopen(str, "w");
    sprintf(str, "./output/%s.MesoForce_y.dat", filestr);
    oFile.MesoyForce = fopen(str, "w");
    sprintf(str, "./output/%s.MesoForce_z.dat", filestr);
    oFile.MesozForce = fopen(str, "w");
  #endif

  #if __COMPUTE_ENERGY__
    sprintf(str, "./output/%s.MesoEnergy.dat", filestr);
    oFile.MesoEnergy = fopen(str, "w");
    sprintf(str, "./output/%s.MesoKinetic.dat", filestr);
    oFile.MesoKinetic = fopen(str, "w");
  #endif

  #if __COMPUTE_TEMPERATURE__
    sprintf(str, "./output/%s.MesoTemp.dat", filestr);
    oFile.MesoTemp = fopen(str, "w");
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
  #endif

  #if __COMPUTE_Q__
    // First part
    sprintf(str, "./output/%s.MesoQ1_x.dat", filestr);
    oFile.MesoQ1_0 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoQ1_y.dat", filestr);
    oFile.MesoQ1_1 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoQ1_z.dat", filestr);
    oFile.MesoQ1_2 = fopen(str, "w");
    // Second part
    sprintf(str, "./output/%s.MesoQ2_x.dat", filestr);
    oFile.MesoQ2_0 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoQ2_y.dat", filestr);
    oFile.MesoQ2_1 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoQ2_z.dat", filestr);
    oFile.MesoQ2_2 = fopen(str, "w");
    // Total heat flux fluid-fluid
    sprintf(str, "./output/%s.MesoQ_x.dat", filestr);
    oFile.MesoQ_0 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoQ_y.dat", filestr);
    oFile.MesoQ_1 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoQ_z.dat", filestr);
    oFile.MesoQ_2 = fopen(str, "w");
  #endif

  #if __COMPUTE_PI__
    sprintf(str, "./output/%s.MesoPi_x.dat", filestr);
    oFile.MesoPi_0 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoPi_y.dat", filestr);
    oFile.MesoPi_1 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoPi_z.dat", filestr);
    oFile.MesoPi_2 = fopen(str, "w");
  #endif
  
  #if __COMPUTE_MOMENTUM__
    sprintf(str, "./output/%s.MesoMomentum_x.dat", filestr);
    oFile.MesoMomentum_0 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoMomentum_y.dat", filestr);
    oFile.MesoMomentum_1 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoMomentum_z.dat", filestr);
    oFile.MesoMomentum_2 = fopen(str, "w");
  #endif

  #if __COMPUTE_VELOCITY__
    sprintf(str, "./output/%s.MesoVelocity_x.dat", filestr);
    oFile.MesoVelocity_0 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoVelocity_y.dat", filestr);
    oFile.MesoVelocity_1 = fopen(str, "w");
    sprintf(str, "./output/%s.MesoVelocity_z.dat", filestr);
    oFile.MesoVelocity_2 = fopen(str, "w");
  #endif

  #if __COMPUTE_INTERNAL_ENERGY__
    sprintf(str, "./output/%s.MesoInternalEnergy.dat", filestr);
    oFile.MesoInternalEnergy = fopen(str, "w");
  #endif
    
  #if __COMPUTE_MACRO_ENERGY__
    sprintf(str, "./output/%s.MacroEnergyUpperWall.dat", filestr);
    oFile.MacroEnergyUpperWall = fopen(str, "w");
    sprintf(str, "./output/%s.MacroEnergyLowerWall.dat", filestr);
    oFile.MacroEnergyLowerWall = fopen(str, "w");
  #endif

  #if __COMPUTE_MACRO_INTERNAL_ENERGY__
    sprintf(str, "./output/%s.MacroInternalEnergyUpperWall.dat", filestr);
    oFile.MacroInternalEnergyUpperWall = fopen(str, "w");
    sprintf(str, "./output/%s.MacroInternalEnergyLowerWall.dat", filestr);
    oFile.MacroInternalEnergyLowerWall = fopen(str, "w");
  #endif

  #if __COMPUTE_MACRO_MOMENTUM__
    sprintf(str, "./output/%s.MacroMomentumUpperWall.dat", filestr);
    oFile.MacroMomentumUpperWall = fopen(str, "w");
    sprintf(str, "./output/%s.MacroMomentumLowerWall.dat", filestr);
    oFile.MacroMomentumLowerWall = fopen(str, "w");
  #endif

  #if __COMPUTE_CENTER_OF_MASS__
    sprintf(str, "./output/%s.CenterOfMassUpperWall.dat", filestr);
    oFile.CenterOfMassUpperWall = fopen(str, "w");
    sprintf(str, "./output/%s.CenterOfMassLowerWall.dat", filestr);
    oFile.CenterOfMassLowerWall = fopen(str, "w");
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
  gsl_matrix * PositionsBase  = gsl_matrix_calloc (NParticles,5);
  gsl_matrix * VelocitiesBase = gsl_matrix_calloc (NParticles,5);

  gsl_matrix * Positions  = gsl_matrix_calloc (NParticles,4);

  gsl_matrix * Velocities = gsl_matrix_calloc (NParticles,3);
  gsl_matrix * Momentum   = gsl_matrix_calloc (NParticles,3);

  FILE *PositionsFile;
  FILE *VelocitiesFile;

  // Linked list
  gsl_vector * List      = gsl_vector_calloc (NParticles);
  gsl_vector * ListHead  = gsl_vector_calloc (Mx*My*Mz);

  // Microscopic variables
  gsl_matrix * Force   = gsl_matrix_calloc (NParticles,3);
  gsl_vector * Energy  = gsl_vector_calloc (NParticles);
  gsl_vector * Kinetic = gsl_vector_calloc (NParticles);

  // Mesoscopic variables
  gsl_matrix * MesoForce     = gsl_matrix_calloc (NNodes,3);
  gsl_vector * MesoDensity_0 = gsl_vector_calloc (NNodes);
  gsl_vector * MesoDensity_1 = gsl_vector_calloc (NNodes);
  gsl_vector * MesoDensity_2 = gsl_vector_calloc (NNodes);
  gsl_vector * MesoEnergy    = gsl_vector_calloc (NNodes);
  gsl_vector * MesoKinetic   = gsl_vector_calloc (NNodes);
  gsl_vector * MesoTemp      = gsl_vector_calloc (NNodes);
  // Sigma matrices will be stored in the following order 
  // 00, 01, 02, 10, 11, 12, 20, 21, 22
  gsl_matrix * MesoSigma1   = gsl_matrix_calloc (NNodes,9);
  gsl_matrix * MesoSigma2   = gsl_matrix_calloc (NNodes,9);
  gsl_matrix * MesoSigma    = gsl_matrix_calloc (NNodes,9);
 
  gsl_matrix * MesoQ1       = gsl_matrix_calloc (NNodes,3);
  gsl_matrix * MesoQ2       = gsl_matrix_calloc (NNodes,3);
  gsl_matrix * MesoQ        = gsl_matrix_calloc (NNodes,3);
  
  gsl_matrix * MesoPi       = gsl_matrix_calloc (NNodes,3);

  gsl_matrix * MesoMomentum = gsl_matrix_calloc (NNodes,3);
  gsl_matrix * MesoVelocity = gsl_matrix_calloc (NNodes,3);
  
  gsl_vector * MesoInternalEnergy = gsl_vector_calloc (NNodes);
  
  // Macroscopic variables
  gsl_vector * MacroMomentumUpper = gsl_vector_calloc(3);
  gsl_vector * MacroMomentumLower = gsl_vector_calloc(3);

  // END OF BLOCK

  // START COMPUTATION

  for (int Step=0;Step<NSteps;Step++)
  {
    strcpy(basename,Snapshot[Step]);

    //pragma omp parallel sections num_threads(2)
    {
      //pragma omp section
      {
        // Positions is a matrix that stores: 
        // TYPE x y z 
        // The ID of a particle corresponds to the row
        PrintMsg("Reading microscopic positions");
        gsl_matrix_set_zero(PositionsBase);
        sprintf(str, "./data/positions/%s", basename);
        printf("\tInput file: %s\n", str);
        PositionsFile = fopen(str, "r");
        gsl_matrix_fscanf(PositionsFile, PositionsBase);

        gsl_vector_view Position_0      = gsl_matrix_column(PositionsBase,1);
        gsl_matrix_set_col(Positions,0,&Position_0.vector);
        gsl_vector_view Position_1      = gsl_matrix_column(PositionsBase,2);
        gsl_matrix_set_col(Positions,1,&Position_1.vector);
        gsl_vector_view Position_2      = gsl_matrix_column(PositionsBase,3);
        gsl_matrix_set_col(Positions,2,&Position_2.vector);
        gsl_vector_view Position_3      = gsl_matrix_column(PositionsBase,4);
        gsl_matrix_set_col(Positions,3,&Position_3.vector);

        fclose(PositionsFile);

        // There are some positions coordinates  that are outside the box lammps can
        // deal  with this issue without  major problems.  Here,  we use  PBC to put
        // all the atom inside the box
        PrintMsg("Fixing PBC in the positions file...");
        FixPBC(Positions);
      }
      //pragma omp section
      {
        // Velocities is a matrix that stores:
        // vx vy vz
        // The ID of a particle corresponds to the row
        PrintMsg("Reading microscopic velocities");
        gsl_matrix_set_zero(VelocitiesBase);
        sprintf(str, "./data/velocities/%s", basename);
        printf("\tInput file: %s\n", str);
        VelocitiesFile = fopen(str, "r");
        gsl_matrix_fscanf(VelocitiesFile, VelocitiesBase);

        gsl_vector_view Velocity_0      = gsl_matrix_column(VelocitiesBase,2);
        gsl_matrix_set_col(Velocities,0,&Velocity_0.vector);
        gsl_vector_view Velocity_1      = gsl_matrix_column(VelocitiesBase,3);
        gsl_matrix_set_col(Velocities,1,&Velocity_1.vector);
        gsl_vector_view Velocity_2      = gsl_matrix_column(VelocitiesBase,4);
        gsl_matrix_set_col(Velocities,2,&Velocity_2.vector);

        fclose(VelocitiesFile);

        // Compute microscopic momentum
        Compute_Momentum(Positions,Velocities,Momentum);
      }
    }

    PrintMsg("Obtaining linked list...");
    Compute_Linked_List(Positions, List, ListHead);

    // Checkpoint: Compute neighboring cells of a TestCell
    //     int TestCell = 60;
    //     printf("Compute neighbors of cell %d:\n", TestCell);
    //     gsl_vector * neighbors = gsl_vector_calloc (27);
    //     Compute_NeighborCells(TestCell, neighbors, Mx, My, Mz);
    //     for (int i=0; i<27; i++)
    //         printf("%f, ", gsl_vector_get(neighbors,i));
    //     printf("\n");
    
    // Checkpoint: Find particles that belong to TestCell
    //
    //     int j = gsl_vector_get(ListHead,TestCell);
    //     while (j >= 0)
    //     {
    //      printf("Particle %d is in cell %d\n", j, TestCell);
    //      j = gsl_vector_get(List,j);
    //     } 
    
    PrintMsg("Computing forces in the fluid (type 2 particles) due to the wall (type 1 particles)");
    Compute_Forces(Positions, Velocities, Neighbors, ListHead, List, 2, 1, Force, Energy, Kinetic);
    
    // Checkpoint: Compare velocities and momentum
    //     gsl_vector_view  gx = gsl_matrix_column(Momentum,0);
    //     gsl_vector_view  vx = gsl_matrix_column(Velocities,0);
    //     PrintInfo(Step, &gx.vector, oFile.MicroG);
    //     PrintInfo(Step, &vx.vector, oFile.MicroV);
    
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

    // Checkpoint: Draw temperature in povray
    //     PrintMsg("Drawing the temperature of the particles...");
    //     gsl_vector * vr = RescaleVector (Vmod);
    //     strcpy (str, "./povray/");
    //     strcat (str, basename);
    //     strcat (str, ".Temperature.inc");
    //     DrawTemperature (Positions,vr,str);
    //     gsl_vector_free (vr);

    // Checkpoint: Find the neighboring cells of the cell in which a TestParticle is into
    //    int TestParticle = 14412;
    //    printf("TESTING PARTICLE %d (type %d) at (%f,%f,%f)\n", TestParticle, ((int) gsl_matrix_get(Positions,TestParticle,0)), 
    //        gsl_matrix_get(Positions,TestParticle,1), gsl_matrix_get(Positions,TestParticle,2), gsl_matrix_get(Positions,TestParticle,3));
    //    int TestCell = FindParticle(Positions,TestParticle);
    //    gsl_vector * NeighboringCells = gsl_vector_calloc(27);
    //    gsl_matrix_get_row(NeighboringCells, Neighbors, TestCell);
    //    printf("Particle %d is in Cell %d\n", TestParticle, TestCell);
    //    printf("Neighboring cells of cell %d are (", TestCell);
    //    for (int i=0;i<27;i++)
    //        printf("%d, ",((int) gsl_vector_get(NeighboringCells,i)));
    //    printf(")\n");
 
    // Checkpoint: Find the Verlet list of TestParticle
    //    int *Verlet;
    //     
    //    Verlet = malloc(27 * NParticles * sizeof(int) / (Mx*My*Mz) );
    //    int NumberOfNeighbors = Compute_VerletList(Positions, TestParticle, NeighboringCells, TestCell, ListHead, List, Verlet);
    //    Verlet = realloc(Verlet, NumberOfNeighbors * sizeof(int));
    //     
    //    printf("Particle %d has %d neighbors\n", TestParticle, NumberOfNeighbors);
    //    for (int i=0;i<NumberOfNeighbors;i++)
    //        printf("%d (type %d) at (%f,%f,%f)\n", Verlet[i], ((int) gsl_matrix_get(Positions,Verlet[i],0)),
    //              gsl_matrix_get(Positions,Verlet[i],1), gsl_matrix_get(Positions,Verlet[i],2), gsl_matrix_get(Positions,Verlet[i],3));
    //    printf(")\n");

    // Checkpoint: Draw a povray script to visualize neighboring cells of particle and its Verlet list.
    //    DrawSim(Micro, TestParticle, TestCell, NeighboringCells, Verlet, NumberOfNeighbors);

    // Checkpoint: Compute Verlet list of a given particle
    //    printf("[%s]\tCompute Verlet list...\n",__TIME__);
    //    //int TestParticle = 172;
    //    int TestParticle = 533;
    //    // int TestParticle = 1068;
    //    int TestCell = FindParticle(Micro,TestParticle);
    //    gsl_vector * NeighboringCells = gsl_vector_calloc(27);
    //        
    //    gsl_matrix_get_row(NeighboringCells, Neighbors, TestCell);
    //    
    //    printf("Particle %d is in Cell %d\n", TestParticle, TestCell);
    //    printf("Neighboring cells of cell %d are (", TestCell);
    //    for (int i=0;i<27;i++)
    //        printf("%d, ",((int) gsl_vector_get(NeighboringCells,i)));
    //    printf(")\n");
    //  
    //    int *Verlet;
    //    
    //    Verlet = malloc(27 * NParticles * sizeof(int) / (Mx*My*Mz) );
    //    int NumberOfNeighbors = Compute_VerletList(Micro, TestParticle, NeighboringCells, TestCell, ListHead, List, Verlet);
    //    Verlet = realloc(Verlet, NumberOfNeighbors * sizeof(int));
    //    
    //    printf("Particle %d has %d neighbors\n", TestParticle, NumberOfNeighbors);
    //    for (int i=0;i<NumberOfNeighbors;i++)
    //        printf("%d, ", Verlet[i]);
    //    printf(")\n");
    // 
    //    DrawSim(Micro, TestParticle, TestCell, NeighboringCells, Verlet, NumberOfNeighbors);

    // MESOSCOPIC INFORMATION
    
    #pragma omp parallel sections
    {
      #if __COMPUTE_DENSITY__
      #pragma omp section
      {
        PrintMsg("Obtaining node densities...");
        Compute_Meso_Density(Positions, z, 1, MesoDensity_1);
        PrintInfo(Step, MesoDensity_1, oFile.MesoDensity_1);
        Compute_Meso_Density(Positions, z, 2, MesoDensity_2);
        PrintInfo(Step, MesoDensity_2, oFile.MesoDensity_2);
        gsl_vector_memcpy(MesoDensity_0, MesoDensity_1);
        gsl_vector_add(MesoDensity_0, MesoDensity_2);
        PrintInfo(Step, MesoDensity_0, oFile.MesoDensity_0);
      }
      #endif
      #if __COMPUTE_FORCE__
      #pragma omp section
      {
        PrintMsg("Obtaining node forces...");
        Compute_Meso_Force(Positions, Force, z, MesoForce);
        gsl_vector_view MesoxForce = gsl_matrix_column(MesoForce,0);
        PrintInfo(Step, &MesoxForce.vector, oFile.MesoxForce);
        gsl_vector_view MesoyForce = gsl_matrix_column(MesoForce,1);
        PrintInfo(Step, &MesoyForce.vector, oFile.MesoyForce);
        gsl_vector_view MesozForce = gsl_matrix_column(MesoForce,2);
        PrintInfo(Step, &MesozForce.vector, oFile.MesozForce);
      }
      #endif
      #if __COMPUTE_ENERGY__
      #pragma omp section
      {
        PrintMsg("Obtaining node energies...");
        Compute_Meso_Profile(Positions, Energy, z, MesoEnergy, 2);
        PrintInfo(Step, MesoEnergy, oFile.MesoEnergy);
      }
      #endif
      #if __COMPUTE_MOMENTUM__
      #pragma omp section
      {
        PrintMsg("Obtaining node momentum...");

        gsl_vector_view MesoMomentum_0  = gsl_matrix_column(MesoMomentum,0);
        gsl_vector_view Momentum_0      = gsl_matrix_column(Momentum,0);
        Compute_Meso_Profile(Positions, &Momentum_0.vector, z, &MesoMomentum_0.vector, 2);
        PrintInfo(Step, &MesoMomentum_0.vector, oFile.MesoMomentum_0);
        
        gsl_vector_view MesoMomentum_1  = gsl_matrix_column(MesoMomentum,1);
        gsl_vector_view Momentum_1      = gsl_matrix_column(Momentum,1);
        Compute_Meso_Profile(Positions, &Momentum_1.vector, z, &MesoMomentum_1.vector, 2);
        PrintInfo(Step, &MesoMomentum_1.vector, oFile.MesoMomentum_1);
        
        gsl_vector_view MesoMomentum_2  = gsl_matrix_column(MesoMomentum,2);
        gsl_vector_view Momentum_2      = gsl_matrix_column(Momentum,2);
        Compute_Meso_Profile(Positions, &Momentum_2.vector, z, &MesoMomentum_2.vector, 2);
        PrintInfo(Step, &MesoMomentum_2.vector, oFile.MesoMomentum_2);
      }
      #endif
      #if __COMPUTE_VELOCITY__
      #pragma omp section
      {
        PrintMsg("Obtaining node velocity...");
        Compute_Meso_Velocity(MesoMomentum,MesoDensity_2,MesoVelocity);

        gsl_vector_view MesoVelocity_0 = gsl_matrix_column(MesoVelocity,0);
        PrintInfo(Step, &MesoVelocity_0.vector, oFile.MesoVelocity_0);
        gsl_vector_view MesoVelocity_1 = gsl_matrix_column(MesoVelocity,1);
        PrintInfo(Step, &MesoVelocity_1.vector, oFile.MesoVelocity_1);
        gsl_vector_view MesoVelocity_2 = gsl_matrix_column(MesoVelocity,2);
        PrintInfo(Step, &MesoVelocity_2.vector, oFile.MesoVelocity_2);
      }
      #endif
      #if __COMPUTE_ENERGY__
      #pragma omp section
      {
        PrintMsg("Obtaining node kinetic energies...");
        Compute_Meso_Profile(Positions, Kinetic, z, MesoKinetic, 2);
        PrintInfo(Step, MesoKinetic, oFile.MesoKinetic);
      }
      #endif
      #if __COMPUTE_STRESS__
      #pragma omp section
      {
        PrintMsg("Obtaining node kinetic stress tensors...");
        Compute_Meso_Sigma1(Positions, Velocities, MesoSigma1);
        gsl_matrix_memcpy(MesoSigma,MesoSigma1);

        gsl_vector_view  MesoSigma1_00 = gsl_matrix_column(MesoSigma1,0);
        PrintInfo(Step, &MesoSigma1_00.vector, oFile.MesoSigma1_00);
        gsl_vector_view  MesoSigma1_01 = gsl_matrix_column(MesoSigma1,1);
        PrintInfo(Step, &MesoSigma1_01.vector, oFile.MesoSigma1_01);
        gsl_vector_view  MesoSigma1_02 = gsl_matrix_column(MesoSigma1,2);
        PrintInfo(Step, &MesoSigma1_02.vector, oFile.MesoSigma1_02);
        gsl_vector_view  MesoSigma1_10 = gsl_matrix_column(MesoSigma1,3);
        PrintInfo(Step, &MesoSigma1_10.vector, oFile.MesoSigma1_10);
        gsl_vector_view  MesoSigma1_11 = gsl_matrix_column(MesoSigma1,4);
        PrintInfo(Step, &MesoSigma1_11.vector, oFile.MesoSigma1_11);
        gsl_vector_view  MesoSigma1_12 = gsl_matrix_column(MesoSigma1,5);
        PrintInfo(Step, &MesoSigma1_12.vector, oFile.MesoSigma1_12);
        gsl_vector_view  MesoSigma1_20 = gsl_matrix_column(MesoSigma1,6);
        PrintInfo(Step, &MesoSigma1_20.vector, oFile.MesoSigma1_20);
        gsl_vector_view  MesoSigma1_21 = gsl_matrix_column(MesoSigma1,7);
        PrintInfo(Step, &MesoSigma1_21.vector, oFile.MesoSigma1_21);
        gsl_vector_view  MesoSigma1_22 = gsl_matrix_column(MesoSigma1,8);
        PrintInfo(Step, &MesoSigma1_22.vector, oFile.MesoSigma1_22);
      }
      #endif
    }

    #if __COMPUTE_STRESS__
      PrintMsg("Obtaining node virial stress tensor...");

      Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, MesoSigma2, z);
      gsl_matrix_add (MesoSigma, MesoSigma2);

      gsl_vector_view  MesoSigma2_00 = gsl_matrix_column(MesoSigma2,0);
      PrintInfo(Step, &MesoSigma2_00.vector, oFile.MesoSigma2_00);
      gsl_vector_view  MesoSigma2_01 = gsl_matrix_column(MesoSigma2,1);
      PrintInfo(Step, &MesoSigma2_01.vector, oFile.MesoSigma2_01);
      gsl_vector_view  MesoSigma2_02 = gsl_matrix_column(MesoSigma2,2);
      PrintInfo(Step, &MesoSigma2_02.vector, oFile.MesoSigma2_02);
      gsl_vector_view  MesoSigma2_10 = gsl_matrix_column(MesoSigma2,3);
      PrintInfo(Step, &MesoSigma2_10.vector, oFile.MesoSigma2_10);
      gsl_vector_view  MesoSigma2_11 = gsl_matrix_column(MesoSigma2,4);
      PrintInfo(Step, &MesoSigma2_11.vector, oFile.MesoSigma2_11);
      gsl_vector_view  MesoSigma2_12 = gsl_matrix_column(MesoSigma2,5);
      PrintInfo(Step, &MesoSigma2_12.vector, oFile.MesoSigma2_12);
      gsl_vector_view  MesoSigma2_20 = gsl_matrix_column(MesoSigma2,6);
      PrintInfo(Step, &MesoSigma2_20.vector, oFile.MesoSigma2_20);
      gsl_vector_view  MesoSigma2_21 = gsl_matrix_column(MesoSigma2,7);
      PrintInfo(Step, &MesoSigma2_21.vector, oFile.MesoSigma2_21);
      gsl_vector_view  MesoSigma2_22 = gsl_matrix_column(MesoSigma2,8);
      PrintInfo(Step, &MesoSigma2_22.vector, oFile.MesoSigma2_22);

      PrintMsg("Saving stress tensors...");

      gsl_vector_view  MesoSigma_00 = gsl_matrix_column(MesoSigma,0);
      PrintInfo(Step, &MesoSigma_00.vector, oFile.MesoSigma_00);
      gsl_vector_view  MesoSigma_01 = gsl_matrix_column(MesoSigma,1);
      PrintInfo(Step, &MesoSigma_01.vector, oFile.MesoSigma_01);
      gsl_vector_view  MesoSigma_02 = gsl_matrix_column(MesoSigma,2);
      PrintInfo(Step, &MesoSigma_02.vector, oFile.MesoSigma_02);
      gsl_vector_view  MesoSigma_10 = gsl_matrix_column(MesoSigma,3);
      PrintInfo(Step, &MesoSigma_10.vector, oFile.MesoSigma_10);
      gsl_vector_view  MesoSigma_11 = gsl_matrix_column(MesoSigma,4);
      PrintInfo(Step, &MesoSigma_11.vector, oFile.MesoSigma_11);
      gsl_vector_view  MesoSigma_12 = gsl_matrix_column(MesoSigma,5);
      PrintInfo(Step, &MesoSigma_12.vector, oFile.MesoSigma_12);
      gsl_vector_view  MesoSigma_20 = gsl_matrix_column(MesoSigma,6);
      PrintInfo(Step, &MesoSigma_20.vector, oFile.MesoSigma_20);
      gsl_vector_view  MesoSigma_21 = gsl_matrix_column(MesoSigma,7);
      PrintInfo(Step, &MesoSigma_21.vector, oFile.MesoSigma_21);
      gsl_vector_view  MesoSigma_22 = gsl_matrix_column(MesoSigma,8);
      PrintInfo(Step, &MesoSigma_22.vector, oFile.MesoSigma_22);
    #endif

    #if __COMPUTE_Q__
      PrintMsg("Obtaining node heat flux fluid-fluid (Q1)...");
      
      Compute_Meso_Q1(Positions, Velocities, Energy, MesoQ1);
      gsl_matrix_memcpy(MesoQ,MesoQ1);
      
      gsl_vector_view  MesoQ1_0 = gsl_matrix_column(MesoQ1,0);
      PrintInfo(Step, &MesoQ1_0.vector, oFile.MesoQ1_0);
      gsl_vector_view  MesoQ1_1 = gsl_matrix_column(MesoQ1,1);
      PrintInfo(Step, &MesoQ1_1.vector, oFile.MesoQ1_1);
      gsl_vector_view  MesoQ1_2 = gsl_matrix_column(MesoQ1,2);
      PrintInfo(Step, &MesoQ1_2.vector, oFile.MesoQ1_2);
      
      PrintMsg("Obtaining node heat flux fluid-fluid (Q2)...");

      Compute_Meso_Q2(Positions, Velocities, Neighbors, ListHead, List, MesoQ2, z);
      gsl_matrix_add (MesoQ, MesoQ2);

      gsl_vector_view  MesoQ2_0 = gsl_matrix_column(MesoQ2,0);
      PrintInfo(Step, &MesoQ2_0.vector, oFile.MesoQ2_0);
      gsl_vector_view  MesoQ2_1 = gsl_matrix_column(MesoQ2,1);
      PrintInfo(Step, &MesoQ2_1.vector, oFile.MesoQ2_1);
      gsl_vector_view  MesoQ2_2 = gsl_matrix_column(MesoQ2,2);
      PrintInfo(Step, &MesoQ2_2.vector, oFile.MesoQ2_2);
      
      PrintMsg("Saving heat flux fluid-fluid...");

      gsl_vector_view  MesoQ_0 = gsl_matrix_column(MesoQ,0);
      PrintInfo(Step, &MesoQ_0.vector, oFile.MesoQ_0);
      gsl_vector_view  MesoQ_1 = gsl_matrix_column(MesoQ,1);
      PrintInfo(Step, &MesoQ_1.vector, oFile.MesoQ_1);
      gsl_vector_view  MesoQ_2 = gsl_matrix_column(MesoQ,2);
      PrintInfo(Step, &MesoQ_2.vector, oFile.MesoQ_2);
    #endif

    #if __COMPUTE_PI__
      PrintMsg("Obtaining and saving node heat flux solid-fluid (Pi)...");
      
      Compute_Meso_Pi(Positions, Velocities, Neighbors, ListHead, List, MesoPi, z);
      
      gsl_vector_view  MesoPi_0 = gsl_matrix_column(MesoPi,0);
      PrintInfo(Step, &MesoPi_0.vector, oFile.MesoPi_0);
      gsl_vector_view  MesoPi_1 = gsl_matrix_column(MesoPi,1);
      PrintInfo(Step, &MesoPi_1.vector, oFile.MesoPi_1);
      gsl_vector_view  MesoPi_2 = gsl_matrix_column(MesoPi,2);
      PrintInfo(Step, &MesoPi_2.vector, oFile.MesoPi_2);
    #endif

    #if __COMPUTE_TEMPERATURE__
      PrintMsg("Obtaining node temperature...");
      Compute_Meso_Temp(MesoKinetic, MesoDensity_2, MesoTemp);
      PrintInfo(Step, MesoTemp, oFile.MesoTemp);
    #endif
        
    #if __COMPUTE_INTERNAL_ENERGY__
      PrintMsg("Obtaining node internal energies...");
      Compute_InternalEnergy(MesoEnergy, MesoMomentum, MesoDensity_2, MesoInternalEnergy);
      PrintInfo(Step, MesoInternalEnergy, oFile.MesoInternalEnergy);
    #endif
    
    // MACROSCOPIC INFORMATION

    #if __COMPUTE_MACRO_ENERGY__
      double MacroEnergyUpper;
      double MacroEnergyLower;

      PrintMsg("Computing the energy of upper wall");
      MacroEnergyUpper = Compute_Macro(Energy, Positions, 1, "top");
      PrintScalarWithIndex(Step, MacroEnergyUpper, oFile.MacroEnergyUpperWall); 

      PrintMsg("Computing the energy of lower wall");
      MacroEnergyLower = Compute_Macro(Energy, Positions, 1, "bottom");
      PrintScalarWithIndex(Step, MacroEnergyLower, oFile.MacroEnergyLowerWall); 
    #endif

    #if __COMPUTE_MACRO_MOMENTUM__
      double TemporalMomentum;

      gsl_vector_view Momentum_0 = gsl_matrix_column(Momentum,0);
      gsl_vector_view Momentum_1 = gsl_matrix_column(Momentum,1);
      gsl_vector_view Momentum_2 = gsl_matrix_column(Momentum,2);
      
      PrintMsg("Computing the momentum of upper wall");
      
      TemporalMomentum = Compute_Macro(&Momentum_0.vector, Positions, 1, "top"); 
      gsl_vector_set(MacroMomentumUpper,0,TemporalMomentum);
      TemporalMomentum = Compute_Macro(&Momentum_1.vector, Positions, 1, "top"); 
      gsl_vector_set(MacroMomentumUpper,1,TemporalMomentum);
      TemporalMomentum = Compute_Macro(&Momentum_2.vector, Positions, 1, "top"); 
      gsl_vector_set(MacroMomentumUpper,2,TemporalMomentum);
      
      PrintInfo(Step, MacroMomentumUpper, oFile.MacroMomentumUpperWall);
    
      PrintMsg("Computing the momentum of lower wall");

      TemporalMomentum = Compute_Macro(&Momentum_0.vector, Positions, 1, "bottom"); 
      gsl_vector_set(MacroMomentumLower,0,TemporalMomentum);
      TemporalMomentum = Compute_Macro(&Momentum_1.vector, Positions, 1, "bottom"); 
      gsl_vector_set(MacroMomentumLower,1,TemporalMomentum);
      TemporalMomentum = Compute_Macro(&Momentum_2.vector, Positions, 1, "bottom"); 
      gsl_vector_set(MacroMomentumLower,2,TemporalMomentum);

      PrintInfo(Step, MacroMomentumLower, oFile.MacroMomentumLowerWall);
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


    #if __COMPUTE_MACRO_INTERNAL_ENERGY__
      double MacroInternalEnergy;
      double TotalMass;

      PrintMsg("Computing Total Mass upper wall");
      TotalMass = Compute_TotalMass(Positions, 1, "top");

      PrintMsg("Computing the internal energy of upper wall");
      MacroInternalEnergy = Compute_MacroInternalEnergy(MacroEnergyUpper, MacroMomentumUpper, TotalMass);
      PrintScalarWithIndex(Step, MacroInternalEnergy, oFile.MacroInternalEnergyUpperWall); 

      PrintMsg("Computing Total Mass lower wall");
      TotalMass = Compute_TotalMass(Positions, 1, "bottom");
      
      PrintMsg("Computing the internal energy of lower wall");
      MacroInternalEnergy = Compute_MacroInternalEnergy(MacroEnergyLower, MacroMomentumLower, TotalMass);
      PrintScalarWithIndex(Step, MacroInternalEnergy, oFile.MacroInternalEnergyLowerWall); 
    #endif
  }
 
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

  #if __COMPUTE_Q_
  fclose(oFile.MesoQ1_0);
  fclose(oFile.MesoQ1_1);
  fclose(oFile.MesoQ1_2);
  
  fclose(oFile.MesoQ2_0);
  fclose(oFile.MesoQ2_1);
  fclose(oFile.MesoQ2_2);
  
  fclose(oFile.MesoQ_0);
  fclose(oFile.MesoQ_1);
  fclose(oFile.MesoQ_2);
  #endif
  
  #if __COMPUTE_PI_
  fclose(oFile.MesoPi_0);
  fclose(oFile.MesoPi_1);
  fclose(oFile.MesoPi_2);
  #endif

  #if __COMPUTE_MOMENTUM__
  fclose(oFile.MesoMomentum_0);
  fclose(oFile.MesoMomentum_1);
  fclose(oFile.MesoMomentum_2);
  #endif

  #if __COMPUTE_VELOCITY__
  fclose(oFile.MesoVelocity_0);
  fclose(oFile.MesoVelocity_1);
  fclose(oFile.MesoVelocity_2);
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
  
  #if __COMPUTE_MACRO_INTERNAL_ENERGY__
  fclose(oFile.MacroInternalEnergyUpperWall);
  fclose(oFile.MacroInternalEnergyLowerWall);
  #endif

  // SECOND COMPUTATION. OBTAIN MEAN VALUES

  PrintMsg("Computing mean values...");
    
  #pragma omp parallel sections
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
      SaveVectorWithIndex(filestr, ".MesoSigma1_z.avg.dat",  z, MesoAverage);
      
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
  gsl_matrix_free(PositionsBase);
  gsl_matrix_free(VelocitiesBase);
  gsl_matrix_free(Momentum);
   
  // Free meso vectors and matrices
  gsl_vector_free(MesoDensity_0);
  gsl_vector_free(MesoDensity_1);
  gsl_vector_free(MesoDensity_2);
  gsl_matrix_free(MesoForce);
  gsl_vector_free(MesoEnergy);
  gsl_vector_free(MesoKinetic);
  gsl_vector_free(MesoTemp);
  gsl_matrix_free(MesoSigma1);
  gsl_matrix_free(MesoSigma2);
  gsl_matrix_free(MesoSigma);
  gsl_matrix_free(MesoQ1);
  gsl_matrix_free(MesoQ2);
  gsl_matrix_free(MesoQ);
  gsl_matrix_free(MesoPi);
  gsl_matrix_free(MesoMomentum);
  gsl_matrix_free(MesoVelocity);
  gsl_vector_free(MesoInternalEnergy);
  
  // Free macro vectors
  gsl_vector_free(MacroMomentumUpper);
  gsl_vector_free(MacroMomentumLower);
  // END OF BLOCK. MEM FREE
  
  PrintMsg("EOF. Have a nice day.");
  return 0;
  }
