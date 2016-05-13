/*
 * Filename   : mesofunctions.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : vie 13 may 2016 14:28:01 CEST
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
    if (gsl_matrix_get(Micro,i,0) == 2) 
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

void Compute_Meso_Sigma1 (gsl_matrix * Positions, gsl_matrix * Velocities, int idx1, int idx2, 
                          gsl_vector * MesoSigma1)
{
  int mu = 0;
  double mass = 0.0;
  double dv = ((float) Lx * Ly * Lz) / NNodes;
  
  gsl_vector_set_zero(MesoSigma1);
  
  for (int i=0;i<NParticles;i++)
  {
    if (gsl_matrix_get(Positions,i,0) == 2)
    {
      mu = floor(gsl_matrix_get(Positions,i,3)*NNodes/Lz) - 1;
      ( mu == -1 ) ? mu = NNodes-1 : mu ;
      // mass = ( gsl_matrix_get(Positions,i,0) == 1 ? m1 : m2 );
      // If type of atom == 1 then it is a wall particle, so it does not contribute to the
      // stress tensor
      mass = ( gsl_matrix_get(Positions,i,0) == 1 ? 0.0 : m2 );

      MesoSigma1->data[mu*MesoSigma1->stride] += mass * gsl_matrix_get(Velocities,i,idx1) * gsl_matrix_get(Velocities,i,idx2);
    }
  }
  gsl_vector_scale(MesoSigma1,1.0/dv);
}

void Compute_Meso_Sigma2 (gsl_matrix * Positions, gsl_matrix * Neighbors, gsl_vector * ListHead,
                          gsl_vector * List, int idx1, int idx2, gsl_vector * MesoSigma2, gsl_vector * z)
{

  gsl_vector_set_zero(MesoSigma2);

  double dv = ((float) Lx * Ly * Lz) / NNodes;
  double dz = ((float) Lz) / NNodes;
  int mu,nu;
  double eij, val, zij, zmu, zsigma, znu;
  double * fij = malloc(3*sizeof(double));


  // Forall i particles
  for (int i=0;i<NParticles;i++)
  {
    // Only for fluid (type 2) particle
    if (gsl_matrix_get(Positions,i,0) == 2)
    {
      // Find the bin mu to which the particle i belongs
      mu = floor(gsl_matrix_get(Positions,i,3)*NNodes/Lz) - 1;
      // NEVER APPLIED (bc there is no type2 particles in bin NNodes)
      // Checkpoint
      if (mu == -1) 
        printf("ERROR! Fluid particle %d in bin %d!\n", i, mu);
      // ( mu == -1 ) ? mu = NNodes-1 : mu ;

      // Find the cell to which the particle i belongs and all its neighboring cells
      int iCell = FindParticle(Positions,i);
      // Code optimization: changed get_row to vector_view
      // gsl_vector * NeighboringCells = gsl_vector_calloc(27);
      // gsl_matrix_get_row(NeighboringCells, Neighbors, iCell);
      gsl_vector_view NeighboringCells = gsl_matrix_row(Neighbors, iCell);
      
      // Find the neighbors of particle i
      int * Verlet = malloc(27 * NParticles * sizeof(int) / (Mx*My*Mz) );
      // Code optimization: changed get_row to vector_view
      // int NNeighbors = Compute_VerletList(Positions, i, NeighboringCells, iCell, ListHead, List, Verlet);
      int NNeighbors = Compute_VerletList(Positions, i, &NeighboringCells.vector, iCell, ListHead, List, Verlet);
      // Not needed
      // Verlet = realloc(Verlet, NNeighbors * sizeof(int));
      
      // Forall Verlet[j] neighboring particles
      for (int j=0;j<NNeighbors;j++)
      {
        // Only for fluid (type 2) particles
        if (gsl_matrix_get(Positions,Verlet[j],0) == 2)
        {
          // Find the bin nu to which the particle Verlet[j] belongs
          nu = floor(gsl_matrix_get(Positions,Verlet[j],3)*NNodes/Lz) - 1;
          // NEVER APPLIED (bc there is no type2 particles in bin NNodes)
          // Checkpoint
          if (nu == -1) 
            printf("ERROR! Fluid particle %d in bin %d!\n", Verlet[j], nu);
          // ( nu == -1 ) ? nu = NNodes-1 : nu ;

          // Compute only the force between  particles of type 2 and particle of
          // type 2 (fluid-fluid interaction)
          // This  portion is  redundant,  the program does  not enter  into the
          // loop for particles different from type 2
          eij = Compute_Force_ij (Positions, i, Verlet[j], 2, 2, fij);
   
          double distance  = gsl_matrix_get(Positions,i,idx1+1) - gsl_matrix_get(Positions,Verlet[j],idx1+1);

          switch (idx1)
          {
            case 0:
              distance -= Lx*round(distance/Lx);
              break;
            case 1:
              distance -= Ly*round(distance/Ly);
              break;
            case 2:
              distance -= Lz*round(distance/Lz);
              break;
          }
          
          val = distance*fij[idx2];

          zij  = gsl_matrix_get(Positions,i,3) - gsl_matrix_get(Positions,Verlet[j],3);
          zij -= Lz*round(zij/Lz);
          
          if (mu == nu)
          {
            MesoSigma2->data[mu*MesoSigma2->stride] += val;
          }
          else if (mu > nu)
          {
            // NEVER APPLIED (bc there are no interacting type2 particles between the walls)
            if (mu - nu > NNodes/2)
            // {
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
            zmu = (mu == NNodes ? gsl_matrix_get(Positions,i,3)/zij : (gsl_matrix_get(Positions,i,3)-gsl_vector_get(z,mu))/zij);
            MesoSigma2->data[mu*MesoSigma2->stride] += val*zmu;
            for (int sigma=mu-1;sigma>nu;sigma--)
            {
              zsigma = (nu == mu-1 ? 0.0 : dz/zij);
              MesoSigma2->data[sigma*MesoSigma2->stride] += val*zsigma;
            }
            znu = (gsl_vector_get(z,nu+1)-gsl_matrix_get(Positions,Verlet[j],3))/zij;
            MesoSigma2->data[nu*MesoSigma2->stride] += val*znu;

            // }
          } 
          else 
          {
            // NEVER APPLIED (bc there are no interacting type2 particles between the walls)
            if (nu-mu > NNodes/2)
            // {
              printf("ERROR! Boundary conditions ji applied! between %d and %d!\n", i, Verlet[j]);
            //   // z is out of range for nu = NNodes-1
            //   // znu = (gsl_vector_get(z,nu+1)-gsl_matrix_get(Positions,Verlet[j],3))/zij;
            //   int bin = ( nu+1 > NNodes - 1 ? (int) Lz : nu+1 );
            //   znu = (gsl_vector_get(z,bin)-gsl_matrix_get(Positions,Verlet[j],3))/zij;
            //   MesoSigma2->data[nu*MesoSigma2->stride] += val*znu;
            //   //for (int sigma=nu+1;sigma<=NNodes-1;sigma++)
            //   for (int sigma=bin;sigma<=NNodes-1;sigma++)
            //   {
            //     zsigma = (NNodes-1 == bin ? 0.0 : dz/zij);
            //     MesoSigma2->data[sigma*MesoSigma2->stride] += val*zsigma;
            //   }
            //   for (int sigma=0;sigma<=mu-1;sigma++)
            //   {
            //     zsigma = (mu-1 == 0 ? 0.0 : dz/zij);
            //     MesoSigma2->data[sigma*MesoSigma2->stride] += val*zsigma;
            //   }
            //   zmu = (gsl_matrix_get(Positions,i,3)-gsl_vector_get(z,mu))/zij;
            //   MesoSigma2->data[mu*MesoSigma2->stride] += val*zmu;
            // }
            // else
            // {

            // Note that nu > mu implies that zij < 0

            // We should consider nu == NNodes
            // znu = (gsl_matrix_get(Positions,Verlet[j],3)-gsl_vector_get(z,nu))/zij;
            znu = (nu == NNodes ? gsl_matrix_get(Positions,Verlet[j],3)/fabs(zij) : (gsl_matrix_get(Positions,Verlet[j],3)-gsl_vector_get(z,nu))/fabs(zij));
            MesoSigma2->data[nu*MesoSigma2->stride] += val*znu;
            for (int sigma=nu-1;sigma>mu;sigma--)
            {
              zsigma = (mu == nu-1 ? 0.0 : dz/fabs(zij));
              MesoSigma2->data[sigma*MesoSigma2->stride] += val*zsigma;
            }
            zmu = (gsl_vector_get(z,mu+1)-gsl_matrix_get(Positions,i,3))/fabs(zij);
            MesoSigma2->data[mu*MesoSigma2->stride] += val*zmu;
            // }
          }
        }
      }
      free(Verlet);
    }
  }
  free(fij);
  gsl_vector_scale(MesoSigma2,0.5/dv);
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
  
void Compute_Mean_Values(FILE * MatrixFile, gsl_vector * MeanValues)
{
  gsl_matrix * InputMatrix = gsl_matrix_calloc (NSteps,NNodes+1);     
  
  gsl_matrix_fscanf(MatrixFile, InputMatrix);

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
