/*
 * Filename   : cg.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : mar 03 may 2016 18:36:06 CEST
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

  if (argc != 2)
  {
    PrintMsg("ERROR");
    PrintMsg("ERROR: CG needs an input basename. Exiting now...");
    PrintMsg("ERROR");
    return 1;
  }

  PrintMsg("INIT");
  
  // Read the list of all snapshots and create the output files

  // Define variables
  char str[100]; 
  memset(str,'\0',sizeof(str));
  char * filestr = argv[1];
  char basename[6];

  // Read input list
  char Snapshot[NSteps][6];
  ReadInputFiles(filestr, Snapshot);
 
  // INIT OF BLOCK. Create output files
  struct OutputFiles oFile;

  // Microscopic files
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MicrozForce.dat");
  oFile.MicrozForce = fopen(str, "w");
  
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MicroEnergy.dat");
  oFile.MicroEnergy = fopen(str, "w");
  
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MicroKinetic.dat");
  oFile.MicroKinetic = fopen(str, "w");
  
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MicroVmod.dat");
  oFile.MicroVmod = fopen(str, "w");

  // Mesoscopic files
  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoDensity.dat");
  oFile.MesoDensity = fopen(str, "w");

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
  strcat (str, ".MesoSigma1_11.dat");
  oFile.MesoSigma1_11 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_22.dat");
  oFile.MesoSigma1_22 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_01.dat");
  oFile.MesoSigma1_01 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma1_12.dat");
  oFile.MesoSigma1_12 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_00.dat");
  oFile.MesoSigma2_00 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_11.dat");
  oFile.MesoSigma2_11 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_22.dat");
  oFile.MesoSigma2_22 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_01.dat");
  oFile.MesoSigma2_01 = fopen(str, "w");

  strcpy (str, "./output/");
  strcat (str, filestr);
  strcat (str, ".MesoSigma2_12.dat");
  oFile.MesoSigma2_12 = fopen(str, "w");

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
  // gsl_vector * zPart   = gsl_vector_calloc (NParticles);
  // gsl_vector * FzPart  = gsl_vector_calloc (NParticles);

  // Mesoscopic variables
  gsl_matrix * MesoForce   = gsl_matrix_calloc (NNodes,3);
  gsl_vector * MesoDensity = gsl_vector_calloc (NNodes);
  gsl_vector * MesoEnergy  = gsl_vector_calloc (NNodes);
  gsl_vector * MesoKinetic = gsl_vector_calloc (NNodes);
  gsl_vector * MesoTemp    = gsl_vector_calloc (NNodes);
  gsl_vector * MesoSigma1  = gsl_vector_calloc (NNodes);
  gsl_vector * MesoSigma2  = gsl_vector_calloc (NNodes);

  // END OF BLOCK

  // START COMPUTATION

  for (int Step=0;Step<NSteps;Step++)
  {
    strcpy(basename,Snapshot[Step]);
    
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

    PrintMsg("Obtaining linked list...");
    Compute_Linked_List(Positions, List, ListHead);
    
    // Checkpoint: Compute neighbors of a TestCell
    //     int TestCell = 60;
    //     printf("Compute neighbors of cell %d:\n", TestCell);
    //     gsl_vector * neighbors = gsl_vector_calloc (27);
    //     Compute_NeighborCells(TestCell, neighbors, Mx, My, Mz);
    //     for (int i=0; i<27; i++)
    //         printf("%f, ", gsl_vector_get(neighbors,i));
    //     printf("\n");
    
    // Checkpoint: Find particles that belong to TestCell
    //
    //     // jatorre@12apr16
    //     // In FORTRAN the loop is done while j != 0. Here, there exists
    //     // a particle labelled with a zero index, so we need to perform
    //     // the loop including that particle
    //     
    //     int j = gsl_vector_get(ListHead,TestCell);
    //     while (j >= 0)
    //     {
    //      printf("Particle %d is in cell %d\n", j, TestCell);
    //      j = gsl_vector_get(List,j);
    //     } 
    //
    
    // clock_t t1, t2;

    // t1 = clock();
    PrintMsg("Computing forces in the fluid (type 2 particles) due to the wall (type 1 particles)");
    printf("\tNote that the energy of all the particles is also computed, regardless of the type of particle\n");
    
    Compute_Forces(Positions, Velocities, Neighbors, ListHead, List, 2, 1, Force, Energy, Kinetic);
    // t2 = clock();
    // printf("\tTime elapsed computing forces and energies: %ld ms\n", timediff(t1,t2));
    
    //  Checkpoint: Print the force exerted on type2 particles
    //              and the energy of all the particles
    // TODO: See if vector_view has a better performance than matrix_get_col
    // gsl_matrix_get_col( zPart, Positions, 3);

    gsl_vector_view FzPart = gsl_matrix_column(Force,2);
    PrintInfo(Step, &FzPart.vector,  oFile.MicrozForce);
    PrintInfo(Step, Energy,  oFile.MicroEnergy);
    PrintInfo(Step, Kinetic, oFile.MicroKinetic);
    
    PrintMsg("Computing the module of the velocity as a estimator for the temperature...");
    gsl_vector * Vmod = Compute_Velocity_Module(Velocities);
    PrintInfo(Step, Vmod, oFile.MicroVmod);
    gsl_vector_free (Vmod);

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

    PrintMsg("Obtaining node densities...");
    Compute_Meso_Density(Positions,z,MesoDensity);
    PrintInfo(Step, MesoDensity, oFile.MesoDensity);

    PrintMsg("Obtaining node forces...");
    Compute_Meso_Force(Positions, Force, z, MesoForce);
    gsl_vector_view MesoxForce = gsl_matrix_column(MesoForce,0);
    PrintInfo(Step, &MesoxForce.vector, oFile.MesoxForce);
    gsl_vector_view MesoyForce = gsl_matrix_column(MesoForce,1);
    PrintInfo(Step, &MesoyForce.vector, oFile.MesoyForce);
    gsl_vector_view MesozForce = gsl_matrix_column(MesoForce,2);
    PrintInfo(Step, &MesozForce.vector, oFile.MesozForce);

    PrintMsg("Obtaining node energies...");
    Compute_Meso_Energy(Positions, Energy, z, MesoEnergy);
    PrintInfo(Step, MesoEnergy, oFile.MesoEnergy);

    PrintMsg("Obtaining node kinetic energies...");
    Compute_Meso_Energy(Positions, Kinetic, z, MesoKinetic);
    PrintInfo(Step, MesoKinetic, oFile.MesoKinetic);

    PrintMsg("Obtaining node temperature...");
    Compute_Meso_Temp(MesoKinetic, MesoDensity, MesoTemp);
    PrintInfo(Step, MesoTemp, oFile.MesoTemp);

    // Remember that velocities are stored as following:
    // VX VY VZ
    // So
    // K_\mu^{xz} \propto v^x v^z = 0 2
   
    PrintMsg("Obtaining node kinetic stress tensors...");

    Compute_Meso_Sigma1(Positions, Velocities, 0, 0, MesoSigma1);
    PrintInfo(Step, MesoSigma1, oFile.MesoSigma1_00);
    Compute_Meso_Sigma1(Positions, Velocities, 1, 1, MesoSigma1);
    PrintInfo(Step, MesoSigma1, oFile.MesoSigma1_11);
    Compute_Meso_Sigma1(Positions, Velocities, 2, 2, MesoSigma1);
    PrintInfo(Step, MesoSigma1, oFile.MesoSigma1_22);
    Compute_Meso_Sigma1(Positions, Velocities, 0, 1, MesoSigma1);
    PrintInfo(Step, MesoSigma1, oFile.MesoSigma1_01);
    Compute_Meso_Sigma1(Positions, Velocities, 1, 2, MesoSigma1);
    PrintInfo(Step, MesoSigma1, oFile.MesoSigma1_12);

    PrintMsg("Obtaining node virial stress tensor...");

    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 0, 0, MesoSigma2, z);
    PrintInfo(Step, MesoSigma2, oFile.MesoSigma2_00);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 1, 1, MesoSigma2, z);
    PrintInfo(Step, MesoSigma2, oFile.MesoSigma2_11);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 2, 2, MesoSigma2, z);
    PrintInfo(Step, MesoSigma2, oFile.MesoSigma2_22);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 0, 1, MesoSigma2, z);
    PrintInfo(Step, MesoSigma2, oFile.MesoSigma2_01);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 1, 2, MesoSigma2, z);
    PrintInfo(Step, MesoSigma2, oFile.MesoSigma2_12);
  
  }
 
  // END OF BLOCK. COMPUTATION DONE

  // BEGIN OF BLOCK. FREE MEM
  
  // Close micro files
  fclose(oFile.MicrozForce);
  fclose(oFile.MicroEnergy);
  fclose(oFile.MicroKinetic);
  fclose(oFile.MicroVmod);

  // Close meso files
  fclose(oFile.MesoDensity);
  fclose(oFile.MesoxForce);
  fclose(oFile.MesoyForce);
  fclose(oFile.MesozForce);
  fclose(oFile.MesoEnergy);
  fclose(oFile.MesoKinetic);
 
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
  gsl_vector_free(MesoDensity);
  gsl_matrix_free(MesoForce);
  gsl_vector_free(MesoEnergy);
  gsl_vector_free(MesoKinetic);
  gsl_vector_free(MesoTemp);
  gsl_vector_free(MesoSigma1);
  gsl_vector_free(MesoSigma2);

  // END OF BLOCK. MEM FREE
  
  PrintMsg("EOF. Have a nice day.");
  return 0;
}
