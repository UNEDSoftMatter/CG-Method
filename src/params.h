/* #############################################################################
#  Parameters of the simulation 
############################################################################# */

// Relative path (from ~/CG-Method) where positions and velocites can be found
#define PositionsFileStr  "../simulation/production.positions"
#define VelocitiesFileStr "../simulation/production.velocities"

// Total number of particles (type1+type2) in the simulation
#define NParticles 11025

// Number of snapshots in dump files
#define NSteps      5

// Number of nodes you want to create
#define NNodes        100

// Size of the simulation box
#define Lx            17.3162
#define Ly            17.3162
#define Lz            34.6325

// Interaction parameters of LJ
#define m1             5.0
#define m2             1.0
#define Rcut           2.5
#define e1             50.0
#define e2             1.0
#define e12            1.0
#define s1             0.73 
#define s2             1.0
#define s12            0.84

// Number of cells in which the simulation box is splitted 
// (Do not change)
#define Mx     (int) (Lx/Rcut)
#define My     (int) (Ly/Rcut)
#define Mz     (int) (Lz/Rcut)

// Cut-energy in the LJ potential
// (used if shift yes is specified in lammps)
// (Do not change)
#define ecut1          ((double)  4.0*e1*(pow( s1/Rcut,12)-pow( s1/Rcut,6)))
#define ecut2          ((double)  4.0*e2*(pow( s2/Rcut,12)-pow( s2/Rcut,6)))
#define ecut12         ((double) 4.0*e12*(pow(s12/Rcut,12)-pow(s12/Rcut,6)))

// Computations to be done
#define true  1
#define false 0

#define __COMPUTE_DENSITY__                 true
#define __COMPUTE_FORCE__                   false
#define __COMPUTE_ENERGY__                  true
#define __COMPUTE_TEMPERATURE__             false
#define __COMPUTE_STRESS__                  false
#define __COMPUTE_MOMENTUM__                true
#define __COMPUTE_VELOCITY__                false
#define __COMPUTE_INTERNAL_ENERGY__         true
#define __COMPUTE_MACRO_ENERGY__            true
#define __COMPUTE_MACRO_MOMENTUM__          true
#define __COMPUTE_CENTER_OF_MASS__          true
#define __COMPUTE_MACRO_INTERNAL_ENERGY__   true
#define __COMPUTE_Q__                       true
#define __COMPUTE_PI__                      true

