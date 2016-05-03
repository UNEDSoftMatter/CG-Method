/* #############################################################################
#  Parameters of the simulation 
############################################################################# */

#define NParticles 12908
#define NSteps         5
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


