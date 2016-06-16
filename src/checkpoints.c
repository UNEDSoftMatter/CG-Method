/*
 * Filename   : checkpoints.c
 *
 * Created    : 16.06.2016
 *
 * Modified   : jue 16 jun 2016 16:47:51 CEST
 *
 * Author     : jatorre
 *
 * Purpose    :  
 *
 */
#include "cg.h"

void Check_Neighboring_Cells(int iCell)
{
  printf("Computing neighbors of cell %d:\n", iCell);
  gsl_vector * neighbors = gsl_vector_calloc (27);
  Compute_NeighborCells(iCell, neighbors, Mx, My, Mz);
  for (int i=0; i<27; i++)
    printf("%f, ", gsl_vector_get(neighbors,i));
  printf("\n");
}
    
// Checkpoint: Find particles that belong to TestCell
void Check_Particles_in_Cell(int iCell, gsl_vector * ListHead, gsl_vector * List)
{
  int j = gsl_vector_get(ListHead,iCell);
  while (j >= 0)
  {
    printf("Particle %d is in cell %d\n", j, iCell);
    j = gsl_vector_get(List,j);
  } 
}
    
// Checkpoint: Draw temperature in povray
void Check_Draw_Temperature(gsl_vector * Vmod, gsl_matrix * Positions, char * basename)
{
  PrintMsg("Drawing the temperature of the particles...");
  
  char str[100]; 
  memset(str,'\0',sizeof(str));
  gsl_vector * vr = RescaleVector (Vmod);
    sprintf(str, "./povray/%s.Temperature.inc", basename);
    DrawTemperature (Positions,vr,str);
  gsl_vector_free (vr);
}

// Checkpoint: Find the Verlet list of TestParticle
void Check_Neighbors_of_a_Particle(gsl_matrix * Positions, int iParticle, gsl_matrix * NeighboringCells, 
                                   int iCell, gsl_vector * ListHead, gsl_vector * List)
{
  int * Verlet = malloc(27 * NParticles * sizeof(int) / (Mx*My*Mz) );
  int iCell (Positions, iParticle);
  int NumberOfNeighbors = Compute_VerletList(Positions, iParticle, NeighboringCells, iCell, ListHead, List, Verlet);
    
  printf("Particle %d has %d neighbors\n", TestParticle, NumberOfNeighbors);
  for (int i=0;i<NumberOfNeighbors;i++)
    printf("%d (type %d) at (%f,%f,%f)\n", Verlet[i], ((int) gsl_matrix_get(Positions,Verlet[i],0)), 
         gsl_matrix_get(Positions,Verlet[i],1), gsl_matrix_get(Positions,Verlet[i],2), gsl_matrix_get(Positions,Verlet[i],3));
  printf(")\n");

  free(Verlet);
}
