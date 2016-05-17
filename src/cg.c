/*
 * Filename   : cg.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : mar 17 may 2016 19:38:43 CEST
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
  system("if [ ! -d output ]; then mkdir output; fi");
  
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
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoDensity_0.dat");
  oFile.MesoDensity_0 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoDensity_1.dat");
  oFile.MesoDensity_1 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoDensity_2.dat");
  oFile.MesoDensity_2 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoxForce.dat");
  oFile.MesoxForce = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoyForce.dat");
  oFile.MesoyForce = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesozForce.dat");
  oFile.MesozForce = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoEnergy.dat");
  oFile.MesoEnergy = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoKinetic.dat");
  oFile.MesoKinetic = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoTemp.dat");
  oFile.MesoTemp = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_00.dat");
  oFile.MesoSigma1_00 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_01.dat");
  oFile.MesoSigma1_01 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_02.dat");
  oFile.MesoSigma1_02 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_10.dat");
  oFile.MesoSigma1_10 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_11.dat");
  oFile.MesoSigma1_11 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_12.dat");
  oFile.MesoSigma1_12 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_20.dat");
  oFile.MesoSigma1_20 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_21.dat");
  oFile.MesoSigma1_21 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_22.dat");
  oFile.MesoSigma1_22 = fopen(str, "w");
  
  // Virial stress tensor
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_00.dat");
  oFile.MesoSigma2_00 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_01.dat");
  oFile.MesoSigma2_01 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_02.dat");
  oFile.MesoSigma2_02 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_10.dat");
  oFile.MesoSigma2_10 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_11.dat");
  oFile.MesoSigma2_11 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_12.dat");
  oFile.MesoSigma2_12 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_20.dat");
  oFile.MesoSigma2_20 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_21.dat");
  oFile.MesoSigma2_21 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_22.dat");
  oFile.MesoSigma2_22 = fopen(str, "w");

  // Total stress tensor
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_00.dat");
  oFile.MesoSigma_00 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_01.dat");
  oFile.MesoSigma_01 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_02.dat");
  oFile.MesoSigma_02 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_10.dat");
  oFile.MesoSigma_10 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_11.dat");
  oFile.MesoSigma_11 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_12.dat");
  oFile.MesoSigma_12 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_20.dat");
  oFile.MesoSigma_20 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_21.dat");
  oFile.MesoSigma_21 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma_22.dat");
  oFile.MesoSigma_22 = fopen(str, "w");

  // END OF BLOCK. All output files created

  // INIT OF BLOCK. Computing vectors and matrices that
  // are constant throughout the program

  PrintMsg("Generating node positions...");
  gsl_vector * z = gsl_vector_calloc(NNodes);
  Compute_Node_Positions(z);
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoNodes.dat");
  SaveVectorWithoutIndex(z, str);

  PrintMsg("Obtaining neighboring matrix...");
  gsl_matrix * Neighbors = gsl_matrix_calloc (Mx*My*Mz,27);
  Compute_NeighborMatrix(Neighbors);

  // END OF BLOCK. All constant quantities created

  // BEGIN OF BLOCK. Definition of needed vectors, matrices, and so on

  // Positions and velocities
  gsl_matrix * Positions  = gsl_matrix_calloc (NParticles,4);
  gsl_matrix * Velocities = gsl_matrix_calloc (NParticles,3);

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

  // END OF BLOCK

  // START COMPUTATION

  for (int Step=0;Step<NSteps;Step++)
  {
    strcpy(basename,Snapshot[Step]);

    #pragma omp parallel sections num_threads(2)
    {
      #pragma omp section
      {
        // Positions is a matrix that stores: 
        // TYPE x y z 
        // The ID of a particle corresponds to the row
        PrintMsg("Reading microscopic positions");
        gsl_matrix_set_zero(Positions);
        strcpy (str, "./data/positions/");
        strcat (str, basename);
        strcat (str, ".pos");
        printf("\tInput file: %s\n", str);
        PositionsFile = fopen(str, "r");
          gsl_matrix_fscanf(PositionsFile, Positions);
        fclose(PositionsFile);

        // There are some positions coordinates  that are outside the box lammps can
        // deal  with this issue without  major problems.  Here,  we use  PBC to put
        // all the atom inside the box
        PrintMsg("Fixing PBC in the positions file...");
        FixPBC(Positions);
      }
      #pragma omp section
      {
        // Velocities is a matrix that stores:
        // vx vy vz
        // The ID of a particle corresponds to the row
        PrintMsg("Reading microscopic velocities");
        gsl_matrix_set_zero(Velocities);
        strcpy (str, "./data/velocities/");
        strcat (str, basename);
        strcat (str, ".vel");
        printf("\tInput file: %s\n", str);
        VelocitiesFile = fopen(str, "r");
          gsl_matrix_fscanf(VelocitiesFile, Velocities);
        fclose(VelocitiesFile);
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
    
    #pragma omp parallel sections num_threads(6)
    {
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
      #pragma omp section
      {
        PrintMsg("Obtaining node energies...");
        Compute_Meso_Energy(Positions, Energy, z, MesoEnergy);
        PrintInfo(Step, MesoEnergy, oFile.MesoEnergy);
      }
      #pragma omp section
      {
        PrintMsg("Obtaining node kinetic energies...");
        Compute_Meso_Energy(Positions, Kinetic, z, MesoKinetic);
        PrintInfo(Step, MesoKinetic, oFile.MesoKinetic);
      }
      #pragma omp section
      {
        PrintMsg("Obtaining node temperature...");
        Compute_Meso_Temp(MesoKinetic, MesoDensity_0, MesoTemp);
        PrintInfo(Step, MesoTemp, oFile.MesoTemp);
      }
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
    }

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

  }
 
  // Close micro files
  // fclose(oFile.MicrozForce);
  // fclose(oFile.MicroEnergy);
  // fclose(oFile.MicroKinetic);
  // fclose(oFile.MicroVmod);

  // Close meso files
  fclose(oFile.MesoDensity_0);
  fclose(oFile.MesoDensity_1);
  fclose(oFile.MesoDensity_2);
  fclose(oFile.MesoxForce);
  fclose(oFile.MesoyForce);
  fclose(oFile.MesozForce);
  fclose(oFile.MesoEnergy);
  fclose(oFile.MesoKinetic);
  
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

  // SECOND COMPUTATION. OBTAIN MEAN VALUES

  PrintMsg("Computing mean values...");

  gsl_vector * MesoAverage = gsl_vector_calloc(NNodes);

  Compute_Mean_Values(filestr, ".MesoDensity_0.dat",        MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoDensity_0.avg.dat", z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoDensity_1.dat",        MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoDensity_1.avg.dat", z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoDensity_2.dat",        MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoDensity_2.avg.dat", z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoxForce.dat",           MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoxForce.avg.dat",    z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoyForce.dat",           MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoyForce.avg.dat",    z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesozForce.dat",           MesoAverage);
  SaveVectorWithIndex(filestr, ".MesozForce.avg.dat",    z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoKinetic.dat",          MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoKinetic.avg.dat",   z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoEnergy.dat",           MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoEnergy.avg.dat",    z, MesoAverage);

  // Stress tensor averages

  // Kinetic stress tensor
  Compute_Mean_Values(filestr, ".MesoSigma1_00.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_00.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma1_01.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_01.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma1_02.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_02.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma1_10.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_10.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma1_11.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_11.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma1_12.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_12.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma1_20.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_20.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma1_21.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_21.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma1_22.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma1_22.avg.dat",  z, MesoAverage);
  
  // Virial stress tensor
  Compute_Mean_Values(filestr, ".MesoSigma2_00.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_00.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma2_01.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_01.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma2_02.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_02.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma2_10.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_10.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma2_11.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_11.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma2_12.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_12.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma2_20.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_20.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma2_21.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_21.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma2_22.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma2_22.avg.dat",  z, MesoAverage);
  
  // Total stress tensor
  Compute_Mean_Values(filestr, ".MesoSigma_00.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_00.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma_01.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_01.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma_02.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_02.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma_10.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_10.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma_11.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_11.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma_12.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_12.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma_20.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_20.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma_21.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_21.avg.dat",  z, MesoAverage);
  Compute_Mean_Values(filestr, ".MesoSigma_22.dat",         MesoAverage);
  SaveVectorWithIndex(filestr, ".MesoSigma_22.avg.dat",  z, MesoAverage);
  
  gsl_vector_free(MesoAverage);

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

  // END OF BLOCK. MEM FREE
  
  PrintMsg("EOF. Have a nice day.");
  return 0;
}
