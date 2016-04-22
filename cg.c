/*
 * Filename   : cg.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : vie 22 abr 2016 13:09:43 CEST
 *
 * Author     : jatorre@fisfun.uned.es
 *
 * Purpose    : Obtain mesoscopic variables
 *              from microscopic details
 *
 */
#include "cg.h"
int main (void) {

    printf("This program computes mesoscopic variables from microscopic configurations\n");
    PrintMsg("INIT");

    PrintMsg("Reading microscopic positions");
    printf("\tInput file: %s\n", iFilePosStr);

    gsl_matrix * Positions = gsl_matrix_calloc (NParticles,4);
    FILE *iFile;
    iFile = fopen(iFilePosStr, "r");
    gsl_matrix_fscanf(iFile, Positions);
    fclose(iFile);

    PrintMsg("Reading microscopic velocities");
    printf("\tInput file: %s\n", iFileVelStr);
    gsl_matrix * Velocities = gsl_matrix_calloc (NParticles,3);
    FILE *iFile2;
    iFile2 = fopen(iFileVelStr, "r");
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
    long elapsed;

    t1 = clock();
    PrintMsg("Computing forces...");
    gsl_matrix * Forces = gsl_matrix_calloc (NParticles,3);
    Compute_Forces(Positions, Neighbors, ListHead, List, 1, 2, Forces);
    t2 = clock();
    
    elapsed = timediff(t1, t2);
    printf("Time elapsed computing forces: %ld ms\n", elapsed);

    //  Checkpoint: Print the force exerted on type1 particles 
    gsl_vector * zPart  = gsl_vector_calloc(NParticles);
    gsl_vector * FzPart = gsl_vector_calloc(NParticles);
    gsl_matrix_get_col( zPart, Positions, 3);
    gsl_matrix_get_col(FzPart, Forces, 2);
    SaveVectorWithIndex(zPart, FzPart, "MicrozForce.dat");
    gsl_vector_free (zPart);
    gsl_vector_free (FzPart);

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

    // Checkpoint
    for (int i=0;i<NNodes;i++)
      printf("z(%d) = %f\n", i, gsl_vector_get(z,i));

    PrintMsg("Obtaining node densities...");
    gsl_vector * MesoDensity = gsl_vector_calloc (NNodes);
    Compute_Meso_Density(Positions,z,MesoDensity);
    SaveVectorWithIndex(z, MesoDensity, "MesoDensity.dat");

    PrintMsg("Obtaining node forces...");
    gsl_matrix * MesoForce = gsl_matrix_calloc (NNodes,3);
    Compute_Meso_Force(Positions, Forces, z, MesoForce);
    SaveMatrixWithIndex(z, MesoForce, "MesoForce.dat");
    
    PrintMsg("Obtaining node stress tensor...");
    gsl_matrix * MesoSigma1 = gsl_matrix_calloc (NNodes,3);
    Compute_Meso_Sigma1(Positions, Velocities, MesoSigma1);
    SaveMatrixWithIndex(z, MesoSigma1, "MesoSigma1.dat");
 
    gsl_vector_free(List);
    gsl_vector_free(ListHead);
    gsl_vector_free(MesoDensity);
    gsl_vector_free(z);
    gsl_matrix_free(Positions);
    // gsl_matrix_free(Velocities);
    gsl_matrix_free(Neighbors);
    gsl_matrix_free(Forces);
    gsl_matrix_free(MesoForce);

    PrintMsg("EOF. Have a nice day.");
    return 0;
}
