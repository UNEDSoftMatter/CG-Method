#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
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

#define e1             5.2895
#define e2             1.0
#define e12            2.2998
#define s1             1.1205 
#define s2             1.0
#define s12            1.0602

#define m1             1.34
#define m2             1.0

#define iFilePosStr  "data/Positions.sort.3col.dat"
#define iFileVelStr "data/Velocities.sort.3col.dat"

void Compute_Meso_Density(gsl_matrix * Micro, gsl_vector * z, gsl_vector * n);
void Compute_Node_Positions(gsl_vector * z);
void Compute_Linked_List(gsl_matrix * Micro, gsl_vector * List, gsl_vector * ListHead);
void Compute_NeighborCells(int cell, gsl_vector * neighbors);
void Compute_NeighborMatrix(gsl_matrix * Neighbors);
int FindParticle(gsl_matrix * Micro, int TestParticle);
int Compute_VerletList(gsl_matrix * Micro, int TestParticle, gsl_vector * NeighboringCells, int TestCell, gsl_vector * LinkedHead, gsl_vector * LinkedList, int * Verlet);

void DrawSim(gsl_matrix * Micro, int TestParticle, int TestCell, gsl_vector * NeighboringCells, int * Verlet, int NumberOfNeighbors);
void Compute_Forces(gsl_matrix * Positions, gsl_matrix * Neighbors, gsl_vector * ListHead, 
                    gsl_vector * List, int type1, int type2, gsl_matrix * Forces);
double GetLJsigma(int type1, int type2);
double GetLJepsilon(int type1, int type2);
void PrintMsg(char *msg);
