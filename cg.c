/*
 * Filename   : cg.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : mié 27 abr 2016 13:50:20 CEST
 *
 * Author     : jatorre@fisfun.uned.es
 *
 * Purpose    : Obtain mesoscopic variables
 *              from microscopic details
 *
 */
#include "cg.h"
int main (int argc, char *argv[]) {
  
    if (argc != 2)
    {
      PrintMsg("ERROR");
      PrintMsg("ERROR: CG needs an input basename. Exiting now...");
      PrintMsg("ERROR");
      return 1;
    }

    char * basename = argv[1];
    char * str = strcpy(str, basename);
   
    PrintInitInfo();

    PrintMsg("INIT");

    PrintMsg("Reading microscopic positions");
    str = strcpy (str, "./data/positions/");
    str = strcat (str, basename);
    str = strcat (str, ".pos");
    printf("\tInput file: %s\n", str);

    gsl_matrix * Positions = gsl_matrix_calloc (NParticles,4);
    FILE *iFile;
    iFile = fopen(str, "r");
    gsl_matrix_fscanf(iFile, Positions);
    fclose(iFile);

    PrintMsg("Reading microscopic velocities");
    str = strcpy (str, "./data/velocities/");
    str = strcat (str, basename);
    str = strcat (str, ".vel");
    printf("\tInput file: %s\n", str);
    gsl_matrix * Velocities = gsl_matrix_calloc (NParticles,3);
    FILE *iFile2;
    iFile2 = fopen(str, "r");
    gsl_matrix_fscanf(iFile2, Velocities);
    fclose(iFile2);

    PrintMsg("Obtaining linked list...");

    // jatorre@12apr16
    // As C arrays begin at 0, we specify the end of a linked list with -1
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
  //  long elapsed;

    t1 = clock();
    PrintMsg("Computing forces due to the wall and energy of each particle...");
    gsl_matrix * Force   = gsl_matrix_calloc (NParticles,3);
    gsl_vector * Energy  = gsl_vector_calloc (NParticles);
    gsl_vector * Kinetic = gsl_vector_calloc (NParticles);
    Compute_Forces(Positions, Velocities, Neighbors, ListHead, List, 1, 2, Force, Energy, Kinetic);
    t2 = clock();
//    elapsed = timediff(t1, t2);
    printf("Time elapsed computing forces and energies: %ld ms\n", timediff(t1,t2));
    
    gsl_vector_free(List);
    gsl_vector_free(ListHead);
    gsl_matrix_free(Neighbors);

    //  Checkpoint: Print the force exerted on type1 particles
    //              and the energy of all the particles
    gsl_vector * zPart  = gsl_vector_calloc(NParticles);
    gsl_vector * FzPart = gsl_vector_calloc(NParticles);
    gsl_matrix_get_col( zPart, Positions, 3);
    gsl_matrix_get_col(FzPart, Force, 2);
    
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MicrozForce.dat");
    SaveVectorWithIndex(zPart, FzPart, NParticles, str);
    
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MicroEnergy.dat");
    SaveVectorWithIndex(zPart, Energy, NParticles, str);
    
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MicroKinetic.dat");
    SaveVectorWithIndex(zPart, Kinetic, NParticles, str);
    
    PrintMsg("Computing the module of the velocity as a estimator for the temperature...");
    gsl_vector * Vmod = Compute_Velocity_Module(Velocities);
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MicroVmodule.dat");
    SaveVectorWithIndex(zPart, Vmod, NParticles, str); 
    
    PrintMsg("Drawing the temperature of the particles...");
    gsl_vector * vr = RescaleVector (Vmod);
    str = strcpy (str, "./povray/");
    str = strcat (str, basename);
    str = strcat (str, ".Temperature.inc");
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

    PrintMsg("Generating node positions...");
    gsl_vector * z     = gsl_vector_calloc(NNodes);
    Compute_Node_Positions(z);
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MesoNodes.dat");
    SaveVectorWithoutIndex(z, str);

    PrintMsg("Obtaining node densities...");
    gsl_vector * MesoDensity = gsl_vector_calloc (NNodes);
    Compute_Meso_Density(Positions,z,MesoDensity);
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MesoDensity.dat");
    SaveVectorWithIndex(z, MesoDensity, NNodes, str);
    
    gsl_vector_free(MesoDensity);

    PrintMsg("Obtaining node forces...");
    gsl_matrix * MesoForce = gsl_matrix_calloc (NNodes,3);
    Compute_Meso_Force(Positions, Force, z, MesoForce);
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MesoForce.dat");
    SaveMatrixWithIndex(z, MesoForce, str);
     
    gsl_matrix_free(MesoForce);
    gsl_matrix_free(Force);
    
    PrintMsg("Obtaining node energies...");
    gsl_vector * MesoEnergy = gsl_vector_calloc (NNodes);
    Compute_Meso_Energy(Positions, Energy, z, MesoEnergy);
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MesoEnergy.dat");
    SaveVectorWithIndex(z, MesoEnergy, NNodes, str);
    
    PrintMsg("Obtaining node kinetic energies...");
    gsl_vector * MesoKinetic   = gsl_vector_calloc (NNodes);
    Compute_Meso_Energy(Positions, Kinetic, z, MesoKinetic);
    str = strcpy (str, "./output/");
    str = strcat (str, basename);
    str = strcat (str, ".MesoKinetic.dat");
    SaveVectorWithIndex(z, MesoKinetic, NNodes, str);
    
    gsl_vector_free(MesoEnergy);
    gsl_vector_free(Energy);
    
    gsl_vector_free(MesoKinetic);
    gsl_vector_free(Kinetic);

    // PrintMsg("Obtaining node kinetic stress tensor...");
    // gsl_matrix * MesoSigma1 = gsl_matrix_calloc (NNodes,3);
    // Compute_Meso_Sigma1(Positions, Velocities, MesoSigma1);
    // SaveMatrixWithIndex(z, MesoSigma1, "MesoKineticStress.dat");
    // 
    // PrintMsg("Obtaining node virial stress tensor...");
    // gsl_matrix * MesoSigma2 = gsl_matrix_calloc (NNodes,3);
    // Compute_Meso_Sigma2(Positions, Neighbors, ListHead, List, MesoSigma2);
    // SaveMatrixWithIndex(z, MesoSigma2, "MesoVirialStress.dat");
 
    gsl_vector_free(z);
    gsl_matrix_free(Positions);
    gsl_matrix_free(Velocities);
    
    // gsl_matrix_free(MesoSigma1);
    // gsl_matrix_free(MesoSigma2);

    PrintMsg("EOF. Have a nice day.");
    return 0;
}
