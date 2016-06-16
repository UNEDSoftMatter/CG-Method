/* #############################################################################
#  Parameters of the simulation 
############################################################################# */

// Relative path (from ~/CG-Method) where positions and velocites can be found
#define FileStr  "./output.data"

// Total number of particles (type1+type2) in the simulation
#define NParticles 18086

// Number of snapshots in dump files
#define NSteps      1000

// Number of nodes you want to create
#define NNodes        64

// Size of the simulation box
#define Lx            40.0  
#define Ly            40.0
#define Lz            16.0

// Interaction parameters of LJ
#define m1             1.0
#define m2             1.34
#define Rcut           2.5
#define e1             1.0
#define e2             5.2895
#define e12            2.2998
#define s1             1.0 
#define s2             1.1205
#define s12            1.0602

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

#define __COMPUTE_DENSITY__         true
#define __COMPUTE_FORCE__           true
#define __COMPUTE_ENERGY__          true
#define __COMPUTE_TEMPERATURE__     true
#define __COMPUTE_STRESS__          true
#define __COMPUTE_MOMENTUM__        true
#define __COMPUTE_VELOCITY__        true
#define __COMPUTE_INTERNAL_ENERGY__ true
#define __COMPUTE_MACRO_ENERGY__    true
#define __COMPUTE_MACRO_MOMENTUM__  true
#define __COMPUTE_CENTER_OF_MASS__  true

