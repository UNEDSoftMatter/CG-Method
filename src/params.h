/* #############################################################################
#  Parameters of the simulation 
############################################################################# */

// Relative path (from ~/CG-Method) where positions and velocites can be found
#define PositionsFileStr  "../output_simulations/production.positions.53500"
#define VelocitiesFileStr "../output_simulations/production.velocities.53500"

// Total number of particles (type1+type2) in the simulation
#define NParticles 11475

// Number of snapshots in dump files
#define NSteps 200

// Number of nodes you want to create
#define NNodes        64

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

