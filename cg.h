#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#define NParticles 23406
#define NNodes        64
#define Lx            40
#define Ly            40
#define Lz            20
#define Rcut           2.5
#define Mx     (int) (Lx/Rcut)
#define My     (int) (Ly/Rcut)
#define Mz     (int) (Lz/Rcut)

#define iFilePosStr  "Positions.sort.3col.dat"
#define iFileVelStr "Velocities.sort.3col.dat"

void Compute_Meso_Density(gsl_matrix * Micro, gsl_vector * z, gsl_vector * n);
void Compute_Node_Positions(gsl_vector * z);
void Compute_Linked_List(gsl_matrix * Micro, gsl_vector * List, gsl_vector * ListHead);
void Compute_NeighborCells(int cell, gsl_vector * neighbors);
void Compute_NeighborMatrix(gsl_matrix * Neighbors);
int FindParticle(gsl_matrix * Micro, int TestParticle);
int Compute_VerletList(gsl_matrix * Micro, int TestParticle, gsl_vector * NeighboringCells,
                       int TestCell, gsl_vector * LinkedHead, gsl_vector * LinkedList, int * Verlet);

void DrawSim(gsl_matrix * Micro, int TestParticle, int TestCell, gsl_vector * NeighboringCells, int * Verlet, int NumberOfNeighbors);
