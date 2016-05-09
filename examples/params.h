/* #############################################################################
#  Parameters of the simulation 
############################################################################# */

#define NParticles 18086
#define NSteps       501
#define NNodes        64
#define Lx            40.0  
#define Ly            40.0
#define Lz            16.0001 
#define Rcut           2.5
#define Mx     (int) (Lx/Rcut)
#define My     (int) (Ly/Rcut)
#define Mz     (int) (Lz/Rcut)
#define RealLz        16.0001

#define e1             1.0
#define e2             5.2895
#define e12            2.2998
#define s1             1.0 
#define s2             1.1205
#define s12            1.0602
#define ecut1          ((double)  4.0*e1*(pow( s1/Rcut,12)-pow( s1/Rcut,6)))
#define ecut2          ((double)  4.0*e2*(pow( s2/Rcut,12)-pow( s2/Rcut,6)))
#define ecut12         ((double) 4.0*e12*(pow(s12/Rcut,12)-pow(s12/Rcut,6)))

#define m1             1.0
#define m2             1.34


