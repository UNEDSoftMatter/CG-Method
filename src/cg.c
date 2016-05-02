/*
 * Filename   : cg.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : lun 02 may 2016 22:52:36 CEST
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

  // Use argv[1] as input basename to write outputs
  char str[100]; 
  memset(str,'\0',sizeof(str));

  PrintMsg("INIT");

  char * filestr = argv[1];
  char basename[6];
  char Snapshot[NSteps][6];
  ReadInputFiles(filestr, Snapshot);

  for (int Step=0;Step<NSteps;Step++)
  {
    strcpy(basename,Snapshot[Step]);
    
    // Positions is a matrix that stores: 
    // TYPE x y z 
    // The ID of a particle corresponds to the row
    PrintMsg("Reading microscopic positions");
    strcpy (str, "./data/positions/");
    strcat (str, basename);
    strcat (str, ".pos");
    printf("\tInput file: %s\n", str);
    gsl_matrix * Positions = gsl_matrix_calloc (NParticles,4);
    FILE *iFile;
    iFile = fopen(str, "r");
      gsl_matrix_fscanf(iFile, Positions);
    fclose(iFile);

    // Velocities is a matrix that stores:
    // vx vy vz
    // The ID of a particle corresponds to the row
    PrintMsg("Reading microscopic velocities");
    strcpy (str, "./data/velocities/");
    strcat (str, basename);
    strcat (str, ".vel");
    printf("\tInput file: %s\n", str);
    gsl_matrix * Velocities = gsl_matrix_calloc (NParticles,3);
    FILE *iFile2;
    iFile2 = fopen(str, "r");
      gsl_matrix_fscanf(iFile2, Velocities);
    fclose(iFile2);

    PrintMsg("Obtaining linked list...");

    // As C arrays goes from 0 to N-1,  we specify the end of a linked list with
    // -1
    gsl_vector * List     = gsl_vector_calloc (NParticles);
    gsl_vector * ListHead = gsl_vector_calloc (Mx*My*Mz);
    
    gsl_vector_add_constant(List,-1.0);
    gsl_vector_add_constant(ListHead,-1.0);

    printf("\tComponents of List:     \t%zu\n", List->size);
    printf("\tComponents of ListHead: \t%zu\n", ListHead->size);
    
    Compute_Linked_List(Positions, List, ListHead);
    
    PrintMsg("Obtaining neighboring matrix...");
    gsl_matrix * Neighbors = gsl_matrix_calloc (Mx*My*Mz,27);
    Compute_NeighborMatrix(Neighbors);

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
    
    clock_t t1, t2;

    t1 = clock();
    PrintMsg("Computing forces in the fluid (type 2 particles) due to the wall (type 1 particles)");
    printf("\tNote that the energy of all the particles is also computed, regardless of the type of particle\n");
    gsl_matrix * Force   = gsl_matrix_calloc (NParticles,3);
    gsl_vector * Energy  = gsl_vector_calloc (NParticles);
    gsl_vector * Kinetic = gsl_vector_calloc (NParticles);
    Compute_Forces(Positions, Velocities, Neighbors, ListHead, List, 2, 1, Force, Energy, Kinetic);
    t2 = clock();
    printf("\tTime elapsed computing forces and energies: %ld ms\n", timediff(t1,t2));
    
    //  Checkpoint: Print the force exerted on type2 particles
    //              and the energy of all the particles
    gsl_vector * zPart  = gsl_vector_calloc(NParticles);
    gsl_vector * FzPart = gsl_vector_calloc(NParticles);
    // TODO: See if vector_view has a better performance than matrix_get_col
    gsl_matrix_get_col( zPart, Positions, 3);
    gsl_matrix_get_col(FzPart, Force, 2);
    
    // Print microscopic information
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MicrozForce.dat");
    SaveVectorWithIndex(zPart, FzPart, NParticles, str);
    
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MicroEnergy.dat");
    SaveVectorWithIndex(zPart, Energy, NParticles, str);
    
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MicroKinetic.dat");
    SaveVectorWithIndex(zPart, Kinetic, NParticles, str);
    
    PrintMsg("Computing the module of the velocity as a estimator for the temperature...");
    gsl_vector * Vmod = Compute_Velocity_Module(Velocities);

    // Print more microscopic information
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MicroVmodule.dat");
    SaveVectorWithIndex(zPart, Vmod, NParticles, str); 
    
    PrintMsg("Drawing the temperature of the particles...");
    gsl_vector * vr = RescaleVector (Vmod);
    strcpy (str, "./povray/");
    strcat (str, basename);
    strcat (str, ".Temperature.inc");
    DrawTemperature (Positions,vr,str);

    gsl_vector_free (Vmod);
    gsl_vector_free (vr);
    gsl_vector_free (FzPart);
    gsl_vector_free (zPart);

    // Checkpoint: Find the neighboring cells of the cell in which a TestParticle is into
    //
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

    PrintMsg("Generating node positions...");
    gsl_vector * z     = gsl_vector_calloc(NNodes);
    Compute_Node_Positions(z);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoNodes.dat");
    SaveVectorWithoutIndex(z, str);

    PrintMsg("Obtaining node densities...");
    gsl_vector * MesoDensity = gsl_vector_calloc (NNodes);
    Compute_Meso_Density(Positions,z,MesoDensity);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoDensity.dat");
    SaveVectorWithIndex(z, MesoDensity, NNodes, str);
    
    PrintMsg("Obtaining node forces...");
    gsl_matrix * MesoForce = gsl_matrix_calloc (NNodes,3);
    Compute_Meso_Force(Positions, Force, z, MesoForce);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoForce.dat");
    SaveMatrixWithIndex(z, MesoForce, str);
     
    gsl_matrix_free(MesoForce);
    gsl_matrix_free(Force);
    
    PrintMsg("Obtaining node energies...");
    gsl_vector * MesoEnergy = gsl_vector_calloc (NNodes);
    Compute_Meso_Energy(Positions, Energy, z, MesoEnergy);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoEnergy.dat");
    SaveVectorWithIndex(z, MesoEnergy, NNodes, str);
    
    PrintMsg("Obtaining node kinetic energies...");
    gsl_vector * MesoKinetic   = gsl_vector_calloc (NNodes);
    Compute_Meso_Energy(Positions, Kinetic, z, MesoKinetic);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoKinetic.dat");
    SaveVectorWithIndex(z, MesoKinetic, NNodes, str);
    
    PrintMsg("Obtaining node temperature...");
    gsl_vector * MesoTemp      = gsl_vector_calloc (NNodes);
    Compute_Meso_Temp(MesoKinetic, MesoDensity, MesoTemp);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoTemp.dat");
    SaveVectorWithIndex(z, MesoTemp, NNodes, str);
    
    gsl_vector_free(MesoTemp);
   
    gsl_vector_free(MesoEnergy);
    gsl_vector_free(Energy);
    
    gsl_vector_free(MesoKinetic);
    gsl_vector_free(Kinetic);

    // Remember that velocities are stored as following:
    // VX VY VZ
    // So
    // K_\mu^{xz} \propto v^x v^z = 0 2
   
    PrintMsg("Obtaining node kinetic stress tensors...");
    
    gsl_vector * MesoSigma1_xz = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma1(Positions, Velocities, 0, 2, MesoSigma1_xz);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoKxz.dat");
    SaveVectorWithIndex(z, MesoSigma1_xz, NNodes, str); 
    gsl_vector_free(MesoSigma1_xz);
    
    gsl_vector * MesoSigma1_xx = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma1(Positions, Velocities, 0, 0, MesoSigma1_xx);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoKxx.dat");
    SaveVectorWithIndex(z, MesoSigma1_xx, NNodes, str); 
    gsl_vector_free(MesoSigma1_xx);
    
    gsl_vector * MesoSigma1_yy = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma1(Positions, Velocities, 1, 1, MesoSigma1_yy);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoKyy.dat");
    SaveVectorWithIndex(z, MesoSigma1_yy, NNodes, str); 
    gsl_vector_free(MesoSigma1_yy);
    
    gsl_vector * MesoSigma1_zz = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma1(Positions, Velocities, 2, 2, MesoSigma1_zz);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoKzz.dat");
    SaveVectorWithIndex(z, MesoSigma1_zz, NNodes, str); 
    gsl_vector_free(MesoSigma1_zz);
    
    gsl_vector * MesoSigma1_xy = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma1(Positions, Velocities, 0, 1, MesoSigma1_xy);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoKxy.dat");
    SaveVectorWithIndex(z, MesoSigma1_xy, NNodes, str); 
    gsl_vector_free(MesoSigma1_xy);

    PrintMsg("Obtaining node virial stress tensor...");
    
    gsl_vector * MesoSigma2_xx = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 0, 0, MesoSigma2_xx, z);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoPxx.dat");
    SaveVectorWithIndex(z, MesoSigma2_xx, NNodes, str);
    gsl_vector_free(MesoSigma2_xx);
    
    gsl_vector * MesoSigma2_yy = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 1, 1, MesoSigma2_yy, z);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoPyy.dat");
    SaveVectorWithIndex(z, MesoSigma2_yy, NNodes, str);
    gsl_vector_free(MesoSigma2_yy);
    
    gsl_vector * MesoSigma2_zz = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 2, 2, MesoSigma2_zz, z);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoPzz.dat");
    SaveVectorWithIndex(z, MesoSigma2_zz, NNodes, str);
    gsl_vector_free(MesoSigma2_zz);
    
    gsl_vector * MesoSigma2_xy = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 0, 1, MesoSigma2_xy, z);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoPxy.dat");
    SaveVectorWithIndex(z, MesoSigma2_xy, NNodes, str);
    gsl_vector_free(MesoSigma2_xy);
    
    gsl_vector * MesoSigma2_xz = gsl_vector_calloc (NNodes);
    Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, 0, 2, MesoSigma2_xz, z);
    strcpy (str, "./output/");
    strcat (str, basename);
    strcat (str, ".MesoPxz.dat");
    SaveVectorWithIndex(z, MesoSigma2_xz, NNodes, str);
    gsl_vector_free(MesoSigma2_xz);

    gsl_vector_free(z);
    gsl_matrix_free(Positions);
    gsl_matrix_free(Velocities);
    gsl_vector_free(MesoDensity);

    gsl_vector_free(List);
    gsl_vector_free(ListHead);
    gsl_matrix_free(Neighbors);
  
  }

  PrintMsg("EOF. Have a nice day.");
  return 0;
}
