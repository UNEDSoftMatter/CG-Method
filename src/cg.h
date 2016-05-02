#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* #############################################################################
#  Parameters of the simulation 
############################################################################# */

#define NParticles 12908
#define NSteps        10
#define NNodes        32
#define Lx            16.9154
#define Ly            16.9154
#define Lz           153.776 
#define Rcut           2.5
#define Mx     (int) (Lx/Rcut)
#define My     (int) (Ly/Rcut)
#define Mz     (int) (Lz/Rcut)
#define RealLz       153.776 

#define e1             1000.0
#define e2             1.0
#define e12            1.0
#define s1             1.0 
#define s2             1.0
#define s12            1.0
#define ecut1          ((double)  4.0*e1*(pow( s1/Rcut,12)-pow( s1/Rcut,6)))
#define ecut2          ((double)  4.0*e2*(pow( s2/Rcut,12)-pow( s2/Rcut,6)))
#define ecut12         ((double) 4.0*e12*(pow(s12/Rcut,12)-pow(s12/Rcut,6)))

#define m1             1.0
#define m2             1.0

/* #############################################################################
#  Microscopic functions 
############################################################################# */

//  Compute the  force that  particle j  (type2) exerts  on particle  i (type1).
// This function only computes the force  if particle j is of type2 and particle
// i is of type1.  It will compute  the force between all the particles if type1
// = type2 = 0 The energy will be always computed between all the particles.
//
// Note that,  provided a VerletList, the best performance is obtained computing
// the interaction between i and VerletList[j].

double Compute_Force_ij (gsl_matrix * Positions, int i, int j, int type1, 
                         int type2, double * fij);

void Compute_Forces (gsl_matrix * Positions, gsl_matrix * Velocities, 
                     gsl_matrix * Neighbors, gsl_vector * ListHead, 
                     gsl_vector * List, int type1, int type2, 
                     gsl_matrix * Forces, gsl_vector * Energy,
                     gsl_vector * Kinetic);

double * GetLJParams (double type1, double type2);

double GetLJsigma (int type1, int type2);

double GetLJepsilon (int type1, int type2);

double KineticEnergy (gsl_vector * i, int type);

gsl_vector * Compute_Velocity_Module (gsl_matrix * Velocities);

/* #############################################################################
#  Auxiliary functions that appear in aux.c 
############################################################################# */

double       MaxVector     (gsl_vector * v); // Obtain the maximum of a vector
double       MinVector     (gsl_vector * v); // Obtain the minimum of a vector
gsl_vector * RescaleVector (gsl_vector * v); // Rescale a vector between 0 and 1


/* #############################################################################
#  Draw functions for povray
############################################################################# */

//  Draw the  position of  each particle  in Micro  and color  from blue  to red
// according to its Velocity

void DrawTemperature (gsl_matrix * Micro, gsl_vector * Velocity, char * str); 

// Draw all the particles in Micro,  focus on TestParticle and draw its Verlet's
// list.  Also draw the neighboring cells

void DrawSim (gsl_matrix * Micro, int TestParticle, int TestCell, 
              gsl_vector * NeighboringCells, int * Verlet, 
              int NumberOfNeighbors);

/* #############################################################################
#  Microscopic functions that are used to build a Verlet list (in verlet.c) 
############################################################################# */

// Obtain a linked list (count the neighboring particles of all the particles )

void Compute_Linked_List (gsl_matrix * Micro, gsl_vector * List, 
                          gsl_vector * ListHead);

// Identify the neighbors of a given cell

void Compute_NeighborCells (int cell, gsl_vector * neighbors);

// Group all the neighboring cells into a matrix

void Compute_NeighborMatrix (gsl_matrix * Neighbors);

// Find the cell in which a particle is located. The function returns the index
// of the cell

int FindParticle (gsl_matrix * Micro, int TestParticle);

// Compute  the list of particles  that are neighbors of  a given particle.  The
// function gives the list of neighboring  particles as an array (Verlet) and it
// returns the total number of neighbors found.

int Compute_VerletList (gsl_matrix * Micro, int TestParticle, 
                        gsl_vector * NeighboringCells, int TestCell, 
                        gsl_vector * LinkedHead, gsl_vector * LinkedList, 
                        int * Verlet);

/* #############################################################################
#  IO functions that appears in io.c 
############################################################################# */

// Store into "File" the matrix Matrix, with the first column given by vector z

void SaveMatrixWithIndex (gsl_vector * z, gsl_matrix * Matrix, char * File);

// Store into "File" the vector Vector, with the first column given by vector z

void SaveVectorWithIndex (gsl_vector * z, gsl_vector * Vector, int size, 
                          char * File);

// Store into "File" the vector z

void SaveVectorWithoutIndex (gsl_vector * z, char * File);

// Print an info message into stdout

void PrintMsg (char *msg);

// Compute time differences

long timediff (clock_t t1, clock_t t2);

// Show initial info about the program

void PrintInitInfo(void);

// Read a txt file and import as an array

void ReadInputFiles(char * iFileStr, char (*iFiles)[6]);

/* #############################################################################
#  Mesoscopic functions in functions.c 
############################################################################# */

void Compute_Node_Positions (gsl_vector * z);

void Compute_Meso_Energy (gsl_matrix * Micro, gsl_vector * MicroEnergy, 
                          gsl_vector * z, gsl_vector * MesoEnergy);

void Compute_Meso_Density (gsl_matrix * Positions, gsl_vector * z, 
                           gsl_vector * MesoDensity);

void Compute_Meso_Force (gsl_matrix * Positions, gsl_matrix * Forces, 
                         gsl_vector * n, gsl_matrix * MesoForce);

void Compute_Meso_Temp(gsl_vector * MesoKinetic, gsl_vector * MesoDensity, 
                       gsl_vector * MesoTemp);

void Compute_Meso_Sigma1 (gsl_matrix * Positions, gsl_matrix * Velocities,
                          int idx1, int idx2, gsl_vector * MesoSigma1);

void Compute_Meso_Sigma2 (gsl_matrix * Positions, gsl_matrix * Neighbors, 
                          gsl_vector * ListHead,  gsl_vector * List, 
                          int idx1, int idx2, gsl_vector * MesoSigma2, 
                          gsl_vector * z);
