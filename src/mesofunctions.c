/*
 * Filename   : mesofunctions.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : jue 19 may 2016 17:24:48 CEST
 *
 * Author     : jatorre
 *
 * Purpose    : Mesoscopic functions for cg.c 
 *
 */
#include "cg.h"

void Compute_Node_Positions(gsl_vector * z)
{
  // Valid only for a regular lattice
  // TODO: Consider irregular lattices
  for (int mu=0;mu<NNodes;mu++)
    gsl_vector_set(z,mu,(double) (mu+1)*Lz/NNodes);
}

void Compute_Meso_Density(gsl_matrix * Micro, gsl_vector * z, int type, 
                          gsl_vector * n)
{

  // RESET vector
  gsl_vector_set_zero(n);

  // Valid for PBC,  this function obtains the  density of a slab of volume Lx *
  //  Ly *  dz The  slab is  build as  a  finite  element  based  on  a Delaunay
  // tessellation
  double dv = ((float) Lx * Ly * Lz) / NNodes;
  double dz = ((float) Lz) / NNodes;
  double zi;
  int muRight, muLeft;

  for (int i=0;i<NParticles;i++)
  {
    zi = gsl_matrix_get(Micro,i,3);
    if (((int) gsl_matrix_get(Micro,i,0) == 0)||((int) gsl_matrix_get(Micro,i,0) == type)) 
    {
      muRight = (int) floor(zi*NNodes/Lz);        
      muLeft  = muRight-1;
      if (muLeft < 0) 
      {
        n->data[muRight*n->stride] += zi/dz;
        n->data[ NNodes*n->stride] += (gsl_vector_get(z,muRight) - zi)/dz;
      } 
      else if (muRight == NNodes)
      {
        n->data[ NNodes*n->stride] += 1.0/dz;
      }
      else 
      {
        n->data[muRight*n->stride] += (zi -  gsl_vector_get(z,muLeft))/dz;
        n->data[ muLeft*n->stride] += (gsl_vector_get(z,muRight) - zi)/dz;
      }
    }
  }
  gsl_vector_scale(n,1.0/dv);
}

void Compute_Meso_Force(gsl_matrix * Positions, gsl_matrix * Forces, 
                        gsl_vector * z, gsl_matrix * MesoForce)
{
  double zi, fx, fy, fz;
  int muRight, muLeft;
  double dv = ((float) Lx * Ly * Lz) / NNodes;
  double dz = ((float) Lz) / NNodes;

  for (int i=0;i<NParticles;i++)
  {
    zi      = gsl_matrix_get(Positions,i,3);
    fx      = gsl_matrix_get(Forces,i,0);
    fy      = gsl_matrix_get(Forces,i,1);
    fz      = gsl_matrix_get(Forces,i,2);
    
    muRight = (int) floor(zi*NNodes/Lz);        
    muLeft  = muRight-1;
    
    if (muLeft < 0) 
    {
      MesoForce->data[muRight*MesoForce->tda+0] += fx * zi/dz;
      MesoForce->data[ muLeft*MesoForce->tda+0] += fx * (gsl_vector_get(z,muRight) - zi)/dz;
      MesoForce->data[muRight*MesoForce->tda+1] += fy * zi/dz;
      MesoForce->data[ muLeft*MesoForce->tda+1] += fy * (gsl_vector_get(z,muRight) - zi)/dz;
      MesoForce->data[muRight*MesoForce->tda+2] += fz * zi/dz;
      MesoForce->data[ muLeft*MesoForce->tda+2] += fz * (gsl_vector_get(z,muRight) - zi)/dz;
    } 
    else if (muRight == NNodes)
    {
      MesoForce->data[ NNodes*MesoForce->tda+0] += fx/dz;
      MesoForce->data[ NNodes*MesoForce->tda+1] += fy/dz;
      MesoForce->data[ NNodes*MesoForce->tda+2] += fz/dz;
    }
    else 
    {
      MesoForce->data[muRight*MesoForce->tda+0] += fx * (zi -  gsl_vector_get(z,muLeft))/dz;
      MesoForce->data[ muLeft*MesoForce->tda+0] += fx * (gsl_vector_get(z,muRight) - zi)/dz;
      MesoForce->data[muRight*MesoForce->tda+1] += fy * (zi -  gsl_vector_get(z,muLeft))/dz;
      MesoForce->data[ muLeft*MesoForce->tda+1] += fy * (gsl_vector_get(z,muRight) - zi)/dz;
      MesoForce->data[muRight*MesoForce->tda+2] += fz * (zi -  gsl_vector_get(z,muLeft))/dz;
      MesoForce->data[ muLeft*MesoForce->tda+2] += fz * (gsl_vector_get(z,muRight) - zi)/dz;
    }
  }
  gsl_matrix_scale(MesoForce,1.0/dv);
}

void Compute_Meso_Sigma1 (gsl_matrix * Positions, gsl_matrix * Velocities, 
                          gsl_matrix * MesoSigma1)
{
  int mu = 0;
  double mass = 0.0;
  double dv = ((float) Lx * Ly * Lz) / NNodes;
  
  gsl_matrix_set_zero(MesoSigma1);
  
  // Loop over all i-particles
  for (int i=0;i<NParticles;i++)
  {
    // Consider only type2 (fluid) particles 
    if ((int) gsl_matrix_get(Positions,i,0) == 2)
    {
      // Obtain the bin to where the i-particle belongs to
      mu = floor(gsl_matrix_get(Positions,i,3)*NNodes/Lz) - 1;
      // PBC: If mu == -1, the particle belongs to the upper bin
      ( mu == -1 ) ? mu = NNodes-1 : mu ;

      // mass = ( gsl_matrix_get(Positions,i,0) == 1 ? m1 : m2 );
      // If type of atom == 1 then it is a wall particle, so it does not contribute to the
      // stress tensor
      mass = ( (int) gsl_matrix_get(Positions,i,0) == 1 ? 0.0 : m2 );

      double * sigma1 = malloc(9*sizeof(double));
      sigma1[0] = gsl_matrix_get(Velocities,i,0) * gsl_matrix_get(Velocities,i,0);
      sigma1[1] = gsl_matrix_get(Velocities,i,0) * gsl_matrix_get(Velocities,i,1);
      sigma1[2] = gsl_matrix_get(Velocities,i,0) * gsl_matrix_get(Velocities,i,2);
      sigma1[3] = gsl_matrix_get(Velocities,i,1) * gsl_matrix_get(Velocities,i,0);
      sigma1[4] = gsl_matrix_get(Velocities,i,1) * gsl_matrix_get(Velocities,i,1);
      sigma1[5] = gsl_matrix_get(Velocities,i,1) * gsl_matrix_get(Velocities,i,2);
      sigma1[6] = gsl_matrix_get(Velocities,i,2) * gsl_matrix_get(Velocities,i,0);
      sigma1[7] = gsl_matrix_get(Velocities,i,2) * gsl_matrix_get(Velocities,i,1);
      sigma1[8] = gsl_matrix_get(Velocities,i,2) * gsl_matrix_get(Velocities,i,2);
     
      for (int j=0;j<9;j++)
        MesoSigma1->data[mu*MesoSigma1->tda+j] += mass * sigma1[j];

      free(sigma1);
    }
  }
  gsl_matrix_scale(MesoSigma1,1.0/dv);
}

void Compute_Meso_Sigma2 (gsl_matrix * Positions, gsl_matrix * Neighbors, gsl_vector * ListHead,
                          gsl_vector * List, gsl_matrix * MesoSigma2, gsl_vector * z)
{

  gsl_matrix_set_zero(MesoSigma2);

  double dv = ((float) Lx * Ly * Lz) / NNodes;
  double dz = ((float) Lz) / NNodes;
  
  int omp_get_max_threads();
  int chunks = NParticles / omp_get_max_threads();
  
  #pragma omp parallel
  {
    #pragma omp for schedule (dynamic,chunks) 
    // Forall i particles
    for (int i=0;i<NParticles;i++)
    {
      // Only for fluid (type 2) particle
      if ((int) gsl_matrix_get(Positions,i,0) == 2)
      {
        // Find the bin mu to which the particle i belongs
        int mu = floor(gsl_matrix_get(Positions,i,3)*NNodes/Lz) - 1;
        // NEVER APPLIED (bc there is no type2 particles in bin NNodes)
        // Checkpoint
        if (mu == -1) 
          printf("ERROR! Fluid particle %d in bin %d!\n", i, mu);
        // ( mu == -1 ) ? mu = NNodes-1 : mu ;

        // Find the cell to which the particle i belongs and all its neighboring cells
        int iCell = FindParticle(Positions,i);
        gsl_vector_view NeighboringCells = gsl_matrix_row(Neighbors, iCell);
        
        // Find the neighbors of particle i
        int * Verlet = malloc(27 * NParticles * sizeof(int) / (Mx*My*Mz) );
        int NNeighbors = Compute_VerletList(Positions, i, &NeighboringCells.vector, iCell, ListHead, List, Verlet);
        
        // Forall Verlet[j] neighboring particles
        for (int j=0;j<NNeighbors;j++)
        {
          // Only for fluid (type 2) particles
          if ((int) gsl_matrix_get(Positions,Verlet[j],0) == 2)
          {
            // Find the bin nu to which the particle Verlet[j] belongs
            int nu = floor(gsl_matrix_get(Positions,Verlet[j],3)*NNodes/Lz) - 1;
            // NEVER APPLIED (bc there is no type2 particles in bin NNodes)
            // Checkpoint
            if (nu == -1) 
              printf("ERROR! Fluid particle %d in bin %d!\n", Verlet[j], nu);
            // ( nu == -1 ) ? nu = NNodes-1 : nu ;

            // Compute only the force between  particles of type 2 and particle of
            // type 2 (fluid-fluid interaction)
            // This  portion is  redundant,  the program does  not enter  into the
            // loop for particles different from type 2
            double * fij = malloc(3*sizeof(double));
            
            double eij = Compute_Force_ij (Positions, i, Verlet[j], 2, 2, fij);
    
            double * rij = malloc(3*sizeof(double));

            rij[0]  = gsl_matrix_get(Positions,i,1) - gsl_matrix_get(Positions,Verlet[j],1);
            rij[0] -= Lx*round(rij[0]/Lx);
            rij[1]  = gsl_matrix_get(Positions,i,2) - gsl_matrix_get(Positions,Verlet[j],2);
            rij[1] -= Ly*round(rij[1]/Ly);
            rij[2]  = gsl_matrix_get(Positions,i,3) - gsl_matrix_get(Positions,Verlet[j],3);
            rij[2] -= Lz*round(rij[2]/Lz);

            double * sigma2 = malloc(9*sizeof(double));

            sigma2[0] = rij[0]*fij[0];
            sigma2[1] = rij[0]*fij[1];
            sigma2[2] = rij[0]*fij[2];
            sigma2[3] = rij[1]*fij[0];
            sigma2[4] = rij[1]*fij[1];
            sigma2[5] = rij[1]*fij[2];
            sigma2[6] = rij[2]*fij[0];
            sigma2[7] = rij[2]*fij[1];
            sigma2[8] = rij[2]*fij[2];

            if (mu == nu)
            {
              for (int k=0;k<9;k++)
                MesoSigma2->data[mu*MesoSigma2->tda+k] += sigma2[k];
            }
            else if (mu > nu)
            {
              // {
              // NEVER APPLIED (bc there are no interacting type2 particles between the walls)
              if (mu - nu > NNodes/2)
                printf("ERROR! Boundary conditions ij applied! between %d and %d!\n", i, Verlet[j]);

              //   // z is out of range for mu = NNodes-1
              //   // zmu = (gsl_vector_get(z,mu+1)-gsl_matrix_get(Positions,i,3))/zij;
              //   int bin = ( mu+1 > NNodes - 1 ? (int) Lz : mu+1 );
              //   zmu = (gsl_vector_get(z,bin)-gsl_matrix_get(Positions,i,3))/zij;
              //   MesoSigma2->data[mu*MesoSigma2->stride] += val*zmu;
              //   // out of range issue
              //   // for (int sigma=mu+1;sigma<=NNodes-1;sigma++)
              //   for (int sigma=bin;sigma<=NNodes-1;sigma++)
              //   {
              //     zsigma = (NNodes-1 == bin ? 0.0 : dz/zij);
              //     MesoSigma2->data[sigma*MesoSigma2->stride] += val*zsigma;
              //   }
              //   for (int sigma=0;sigma<=nu-1;sigma++)
              //   {
              //     zsigma = (nu-1 == 0 ? 0.0 : dz/zij);
              //     MesoSigma2->data[sigma*MesoSigma2->stride] += val*zsigma;
              //   }
              //   znu = (gsl_matrix_get(Positions,Verlet[j],3)-gsl_vector_get(z,nu))/zij;
              //   MesoSigma2->data[nu*MesoSigma2->stride] += val*znu;
              // }
              // else 
              // {
              
              // We should consider mu == NNodes
              //  zmu = (gsl_matrix_get(Positions,i,3)-gsl_vector_get(z,mu))/zij;
              double zmu = (mu == NNodes ? gsl_matrix_get(Positions,i,3)/rij[2] : (gsl_matrix_get(Positions,i,3)-gsl_vector_get(z,mu))/rij[2]);
              for (int k=0;k<9;k++)
                MesoSigma2->data[mu*MesoSigma2->tda+k] += sigma2[k]*zmu;
              for (int sigma=mu-1;sigma>nu;sigma--)
              {
                double zsigma = (nu == mu-1 ? 0.0 : dz/rij[2]);
                for (int k=0;k<9;k++)
                  MesoSigma2->data[mu*MesoSigma2->tda+k] += sigma2[k]*zsigma;
              }
              double znu = (gsl_vector_get(z,nu+1)-gsl_matrix_get(Positions,Verlet[j],3))/rij[2];
              for (int k=0;k<9;k++)
                MesoSigma2->data[mu*MesoSigma2->tda+k] += sigma2[k]*znu;

              // }
            } 
            else 
            {
              // NEVER APPLIED (bc there are no interacting type2 particles between the walls)
              if (nu-mu > NNodes/2)
              // {
                printf("ERROR! Boundary conditions ji applied! between %d and %d!\n", i, Verlet[j]);
              //   // z is out of range for nu = NNodes-1
              //   // znu = (gsl_vector_get(z,nu+1)-gsl_matrix_get(Positions,Verlet[j],3))/rij[2];
              //   int bin = ( nu+1 > NNodes - 1 ? (int) Lz : nu+1 );
              //   znu = (gsl_vector_get(z,bin)-gsl_matrix_get(Positions,Verlet[j],3))/rij[2];
              //   MesoSigma2->data[nu*MesoSigma2->stride] += val*znu;
              //   //for (int sigma=nu+1;sigma<=NNodes-1;sigma++)
              //   for (int sigma=bin;sigma<=NNodes-1;sigma++)
              //   {
              //     zsigma = (NNodes-1 == bin ? 0.0 : dz/rij[2]);
              //     MesoSigma2->data[sigma*MesoSigma2->stride] += val*zsigma;
              //   }
              //   for (int sigma=0;sigma<=mu-1;sigma++)
              //   {
              //     zsigma = (mu-1 == 0 ? 0.0 : dz/rij[2]);
              //     MesoSigma2->data[sigma*MesoSigma2->stride] += val*zsigma;
              //   }
              //   zmu = (gsl_matrix_get(Positions,i,3)-gsl_vector_get(z,mu))/rij[2];
              //   MesoSigma2->data[mu*MesoSigma2->stride] += val*zmu;
              // }
              // else
              // {

              // Note that nu > mu implies that rij[2] < 0

              // We should consider nu == NNodes
              // znu = (gsl_matrix_get(Positions,Verlet[j],3)-gsl_vector_get(z,nu))/rij[2];
              double znu = (nu == NNodes ? gsl_matrix_get(Positions,Verlet[j],3)/fabs(rij[2]) : (gsl_matrix_get(Positions,Verlet[j],3)-gsl_vector_get(z,nu))/fabs(rij[2]));
              for (int k=0;k<9;k++)
                MesoSigma2->data[mu*MesoSigma2->tda+k] += sigma2[k]*znu;
              for (int sigma=nu-1;sigma>mu;sigma--)
              {
                double zsigma = (mu == nu-1 ? 0.0 : dz/fabs(rij[2]));
                for (int k=0;k<9;k++)
                  MesoSigma2->data[mu*MesoSigma2->tda+k] += sigma2[k]*zsigma;
              }
              double zmu = (gsl_vector_get(z,mu+1)-gsl_matrix_get(Positions,i,3))/fabs(rij[2]);
              for (int k=0;k<9;k++)
                MesoSigma2->data[mu*MesoSigma2->tda+k] += sigma2[k]*zmu;
              // }
            }
            free(fij);
            free(rij);
            free(sigma2);
          }
        }
        free(Verlet);
      }
    }
  }
  gsl_matrix_scale(MesoSigma2,0.5/dv);
}
          
void Compute_Meso_Energy(gsl_matrix * Micro, gsl_vector * MicroEnergy, gsl_vector * z, gsl_vector * MesoEnergy)
{
  double dv = ((float) Lx * Ly * Lz) / NNodes;
  double dz = ((float) Lz) / NNodes;
  double zi, ei;
  int muRight, muLeft;

  // RESET vector
  gsl_vector_set_zero(MesoEnergy);

  for (int i=0;i<NParticles;i++)
  {
    zi      = gsl_matrix_get(Micro,i,3);
    ei      = gsl_vector_get(MicroEnergy,i);
    muRight = (int) floor(zi*NNodes/Lz);        
    muLeft  = muRight-1;
    if (muLeft < 0) 
    {
      MesoEnergy->data[muRight*MesoEnergy->stride] += ei * zi/dz;
      MesoEnergy->data[ NNodes*MesoEnergy->stride] += ei * (gsl_vector_get(z,muRight) - zi)/dz;
    } 
    else if (muRight == NNodes)
    {
      MesoEnergy->data[ NNodes*MesoEnergy->stride] += ei * 1.0/dz;
    }
    else 
    {
      MesoEnergy->data[muRight*MesoEnergy->stride] += ei * (zi -  gsl_vector_get(z,muLeft))/dz;
      MesoEnergy->data[ muLeft*MesoEnergy->stride] += ei * (gsl_vector_get(z,muRight) - zi)/dz;
    }
  }
  gsl_vector_scale(MesoEnergy,1.0/dv);
}

void Compute_Meso_Temp(gsl_vector * MesoKinetic, gsl_vector * MesoDensity, gsl_vector * MesoTemp)
{
  double val;

  for (int mu=0;mu<NNodes;mu++)
  {
    if (gsl_vector_get(MesoDensity,mu) != 0)
    {
      val = gsl_vector_get(MesoKinetic,mu) / gsl_vector_get(MesoDensity,mu);
    }
    else
    {
      val = 0.0;
    }
    gsl_vector_set(MesoTemp,mu,val);
  }
  gsl_vector_scale(MesoTemp,2.0/3.0);
}
  
void Compute_Mean_Values(char * basename, char * filename, gsl_vector * MeanValues)
{
  gsl_matrix * InputMatrix = gsl_matrix_calloc (NSteps,NNodes+1);     
  
  FILE *iFile;
  char str[100]; 
  strcpy (str, "./output/");
  strcat (str, basename);
  strcat (str, filename);
  iFile = fopen(str, "r");
    gsl_matrix_fscanf(iFile, InputMatrix);
  fclose(iFile);

  gsl_vector_set_zero(MeanValues);

  for (int i=0;i<NSteps;i++)
  {
    for (int mu=1;mu<NNodes+1;mu++)
    {
      MeanValues->data[(mu-1)*MeanValues->stride] += gsl_matrix_get(InputMatrix,i,mu); 
    }
  }

  gsl_vector_scale(MeanValues,1.0/NSteps);
}

void Compute_Meso_Velocity(gsl_matrix * MesoMomentum, gsl_vector * MesoDensity, gsl_matrix * MesoVelocity)
{
  gsl_matrix_memcpy(MesoVelocity,MesoMomentum);
  
  gsl_vector_view vx = gsl_matrix_column(MesoVelocity,0);
  gsl_vector_div(&vx.vector,MesoDensity);
  gsl_vector_view vy = gsl_matrix_column(MesoVelocity,1);
  gsl_vector_div(&vy.vector,MesoDensity);
  gsl_vector_view vz = gsl_matrix_column(MesoVelocity,2);
  gsl_vector_div(&vz.vector,MesoDensity);

  for (int mu=0;mu<NNodes;mu++)
  {
    gsl_vector_view vmu = gsl_matrix_row(MesoVelocity,mu); 
    if (fabs(gsl_vector_get(MesoDensity,mu)) < 10E-8)
      gsl_vector_set_zero(&vmu.vector);
  }
}
          
void Compute_Meso_Profile(gsl_matrix * Positions, gsl_vector * Micro, gsl_vector * z, gsl_vector * Meso)
{
  double dv = ((float) Lx * Ly * Lz) / NNodes;
  double dz = ((float) Lz) / NNodes;
  double zi, ei;
  int muRight, muLeft;

  // RESET vector
  gsl_vector_set_zero(Meso);

  for (int i=0;i<NParticles;i++)
  {
    zi      = gsl_matrix_get(Positions,i,3);
    ei      = gsl_vector_get(Micro,i);
    muRight = (int) floor(zi*NNodes/Lz);        
    muLeft  = muRight-1;
    if (muLeft < 0) 
    {
      Meso->data[muRight*Meso->stride] += ei * zi/dz;
      Meso->data[ NNodes*Meso->stride] += ei * (gsl_vector_get(z,muRight) - zi)/dz;
    } 
    else if (muRight == NNodes)
    {
      Meso->data[ NNodes*Meso->stride] += ei * 1.0/dz;
    }
    else 
    {
      Meso->data[muRight*Meso->stride] += ei * (zi -  gsl_vector_get(z,muLeft))/dz;
      Meso->data[ muLeft*Meso->stride] += ei * (gsl_vector_get(z,muRight) - zi)/dz;
    }
  }
  gsl_vector_scale(Meso,1.0/dv);
}
        
void Compute_InternalEnergy(gsl_vector * MesoEnergy, gsl_matrix * MesoMomentum, 
                            gsl_vector * MesoDensity, gsl_vector * InternalEnergy)
{
  gsl_vector * GMod2 = gsl_vector_calloc(NNodes);

  for (int mu=0;mu<NNodes;mu++)
  {
    double gmu2 = pow(gsl_matrix_get(MesoMomentum,mu,0),2) + pow(gsl_matrix_get(MesoMomentum,mu,1),2) + pow(gsl_matrix_get(MesoMomentum,mu,2),2);
    gsl_vector_set(GMod2,mu,gmu2);
  }

  gsl_vector_memcpy(InternalEnergy,GMod2);
  gsl_vector_scale(InternalEnergy,-1.0);
  gsl_vector_div(InternalEnergy,MesoDensity);
  
  for (int mu=0;mu<NNodes;mu++)
  {
    if (fabs(gsl_vector_get(MesoDensity,mu)) < 10E-8)
      gsl_vector_set(InternalEnergy,mu,0.0);
  }
  
  gsl_vector_add(InternalEnergy,MesoEnergy);

  gsl_vector_free(GMod2);
}

