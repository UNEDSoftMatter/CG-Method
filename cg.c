/*
 * Filename   : cg.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : lun 18 abr 2016 21:10:27 CEST
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

    printf("[%s]\tReading microscopic positions from file %s...\n", __TIME__, iFilePosStr);
    gsl_matrix * Positions = gsl_matrix_calloc (NParticles,4);
    FILE *iFile;
    iFile = fopen(iFilePosStr, "r");
    gsl_matrix_fscanf(iFile, Positions);
    fclose(iFile);

//     printf("[%s]\tReading microscopic velocities from file %s...\n", __TIME__, iFileVelStr);
//     gsl_matrix * Velocities = gsl_matrix_calloc (NParticles,3);
//     FILE *iFile2;
//     iFile2 = fopen(iFileVelStr, "r");
//     gsl_matrix_fscanf(iFile2, Velocities);
//     fclose(iFile2);

    printf("[%s]\tObtaining linked list...\n", __TIME__);

    // jatorre@12apr16
    // As C arrays begin at 0, we specify the end of a linked list with -1
    gsl_vector * List     = gsl_vector_calloc (NParticles);
    gsl_vector * ListHead = gsl_vector_calloc (Mx*My*Mz);
    
    gsl_vector_add_constant(List,-1.0);
    gsl_vector_add_constant(ListHead,-1.0);

    printf("Components of List:     \t%zu\n", List->size);
    printf("Components of ListHead: \t%zu\n", ListHead->size);
    
    Compute_Linked_List(Positions, List, ListHead);
    
    printf("[%s]\tObtaining Neighboring Matrix...\n", __TIME__);
    gsl_matrix * Neighbors = gsl_matrix_calloc (Mx*My*Mz,27);
    Compute_NeighborMatrix(Neighbors);

    // Checkpoint
    //     printf("Compute neighbors of cell 60:\n");
    //     gsl_vector * neighbors = gsl_vector_calloc (27);
    //     gsl_matrix_get_row (neighbors, Neighbors, 60);
    //     for (int i=0; i<27; i++)
    //         printf("%f, ", gsl_vector_get(neighbors,i));
    //     printf("\n");
    //     gsl_vector_free(neighbors);

    //     // Checkpoint
    //     for (int i=0;i<NParticles;i++)
    //       printf("List(%d) = %f\n", i, gsl_vector_get(List,i));
    // 
    //     // Checkpoint
    //     printf("Compute neighbors of cell %d:\n", TestCell);
    //     int TestCell = 60;
    //     gsl_vector * neighbors = gsl_vector_calloc (27);
    //     Compute_NeighborCells(TestCell, neighbors, Mx, My, Mz);
    //     for (int i=0; i<27; i++)
    //         printf("%f, ", gsl_vector_get(neighbors,i));
    //     printf("\n");
    //     // jatorre@12apr16
    //     // In FORTRAN the loop is done while j != 0. Here, there exists
    //     // a particle labelled with a zero index, so we need to perform
    //     // the loop including that particle
    //     //
    //     int j = gsl_vector_get(ListHead,TestCell);
    //     while (j >= 0)
    //     {
    //      printf("Particle %d is in cell %d\n", j, TestCell);
    //      j = gsl_vector_get(List,j);
    //     } 
    //

    printf("[%s]\tComputing forces...\n", __TIME__);
    gsl_matrix * Forces = gsl_matrix_calloc (NParticles,3);
    Compute_Forces(Positions, Neighbors, ListHead, List, 1, 2, Forces);

    for (int i=0;i<NParticles;i++)
    {
      if (gsl_matrix_get(Positions,i,0) == 1.0) 
      {
        double distance1 = gsl_matrix_get(Positions,i,3)-2.0;
        double distance2 = 18.0-gsl_matrix_get(Positions,i,3);
        double distance = distance1;
        (distance > distance2) ? distance = distance2 : distance;
        printf("The wall exerts on particle %d (at %f to the nearest wall)) the following force: (%f,%f,%f)\n", i, distance,
        gsl_matrix_get(Forces,i,0), gsl_matrix_get(Forces,i,1), gsl_matrix_get(Forces,i,2));
      }
    }
//     printf("[%s]\tCompute Verlet list...\n",__TIME__);
//     //int TestParticle = 172;
//     int TestParticle = 533;
//     // int TestParticle = 1068;
//     int TestCell = FindParticle(Micro,TestParticle);
//     gsl_vector * NeighboringCells = gsl_vector_calloc(27);
//         
//     gsl_matrix_get_row(NeighboringCells, Neighbors, TestCell);
//     
//     printf("Particle %d is in Cell %d\n", TestParticle, TestCell);
//     printf("Neighboring cells of cell %d are (", TestCell);
//     for (int i=0;i<27;i++)
//         printf("%d, ",((int) gsl_vector_get(NeighboringCells,i)));
//     printf(")\n");
//  
//     int *Verlet;
//     
//     Verlet = malloc(27 * NParticles * sizeof(int) / (Mx*My*Mz) );
//     int NumberOfNeighbors = Compute_VerletList(Micro, TestParticle, NeighboringCells, TestCell, ListHead, List, Verlet);
//     Verlet = realloc(Verlet, NumberOfNeighbors * sizeof(int));
//     
//     printf("Particle %d has %d neighbors\n", TestParticle, NumberOfNeighbors);
//     for (int i=0;i<NumberOfNeighbors;i++)
//         printf("%d, ", Verlet[i]);
//     printf(")\n");
// 
//     DrawSim(Micro, TestParticle, TestCell, NeighboringCells, Verlet, NumberOfNeighbors);

//     printf("[%s]\tGenerating node positions...\n",__TIME__);
//     gsl_vector * z     = gsl_vector_calloc(NNodes);
//     Compute_Node_Positions(z);

    // Checkpoint
    //   for (int i=0;i<NNodes;i++)
    //     printf("z(%d) = %f\n", i, gsl_vector_get(z,i));

//      printf("[%s]\tObtaining node densities...\n",__TIME__);
//      gsl_vector * n = gsl_vector_calloc (NNodes);
//     Compute_Meso_Density(Micro,z,n);
 
    // Checkpoint
    //   for (int i=0;i<NNodes;i++)
    //      printf("n(%d) = %f, at z = %f\n", i, gsl_vector_get(n,i), gsl_vector_get(z,i));

    gsl_vector_free(List);
    gsl_vector_free(ListHead);
//    gsl_vector_free(n);
    gsl_matrix_free(Positions);
//     gsl_matrix_free(Velocities);
    gsl_matrix_free(Neighbors);
    gsl_matrix_free(Forces);

    PrintMsg("EOF");
    return 0;
}
