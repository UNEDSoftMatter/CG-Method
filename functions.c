/*
 * Filename   : functions.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : mar 26 abr 2016 12:09:42 CEST
 *
 * Author     : jatorre
 *
 * Purpose    : Functions for cg.c 
 *
 */
#include "cg.h"

void Compute_Node_Positions(gsl_vector * z)
{
    for (int mu=0;mu<NNodes;mu++)
        gsl_vector_set(z,mu,(double) (mu+1)*RealLz/NNodes);
}

void Compute_Meso_Density(gsl_matrix * Micro, gsl_vector * z, gsl_vector * n)
{
    double dv = ((float) Lx * Ly * RealLz) / NNodes;
    double dz = ((float) RealLz) / NNodes;
    double zi;
    int muRight, muLeft;

    for (int i=0;i<NParticles;i++)
    {
        zi      = gsl_matrix_get(Micro,i,3);
        muRight = (int) floor(zi*NNodes/RealLz);        
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
    gsl_vector_scale(n,1.0/dv);
}

void Compute_Forces(gsl_matrix * Positions, gsl_matrix * Velocities, gsl_matrix * Neighbors, gsl_vector * ListHead, 
                    gsl_vector * List, int type1, int type2, gsl_matrix * Forces, gsl_vector * Energy)
{

  gsl_matrix_scale(Forces,0.0);
  gsl_vector_scale(Energy,0.0);

  // Begin of parallel region
  
  int omp_get_thread_num();
  int omp_get_max_threads();
  int chunks = NParticles / omp_get_max_threads();

  #pragma omp parallel
  {
    #pragma omp for schedule (dynamic,chunks)
    for (int i=0;i<NParticles;i++)
    {
      gsl_vector * vi = gsl_vector_calloc(3);
      gsl_matrix_get_row(vi, Velocities,i);
      double * fij = malloc(3*sizeof(double));
      double   eij = KineticEnergy(vi, (int) gsl_matrix_get(Positions,i,0));
      
      gsl_vector * NeighboringCells = gsl_vector_calloc(27);
      int * Verlet = malloc(27 * NParticles * sizeof(int) / (Mx*My*Mz));
      int iCell    = FindParticle(Positions,i);
      gsl_matrix_get_row(NeighboringCells, Neighbors, iCell);
           
      int NNeighbors = Compute_VerletList(Positions, i, NeighboringCells, iCell, ListHead, List, Verlet);
      gsl_vector_free(NeighboringCells);

      Verlet = realloc(Verlet, NNeighbors * sizeof(int));
 
      for (int j=0;j<NNeighbors;j++)
      {
        if (i != Verlet[j])
        {
          eij = Compute_Force_ij(Positions, i, Verlet[j], fij);
          Forces->data[i*Forces->tda + 0] += fij[0];
          Forces->data[i*Forces->tda + 1] += fij[1];
          Forces->data[i*Forces->tda + 2] += fij[2];
          Energy->data[i*Energy->stride]  += eij;
        }
      }
      Verlet = realloc(Verlet, 27 * NParticles * sizeof(int) / (Mx * My * Mz));
      gsl_vector_free(vi);
      free(fij);
    }
  }

  // End of parallel region

  gsl_vector_scale(Energy,0.5);
  
}

double *GetLJParams(double type1, double type2)
{
  double *lj = malloc(3*sizeof(double));

  if (type1 == type2)
  {
    if (type1 == 1.0) 
    {
      lj[0] = e1;
      lj[1] = s1;
      lj[2] = ecut1;
    }
    else
    {
      lj[0] = e2;
      lj[1] = s2;
      lj[2] = ecut2;
    }
  } 
  else
  {
      lj[0] = e12;
      lj[1] = s12;
      lj[2] = ecut12;
  }
  return lj;
}

double GetLJepsilon(int type1, int type2)
{
  double epsilon;
  if (type1 == type2)
  {
    if (type1 == 1.0) 
    {
      epsilon = e1;
    }
    else
    {
      epsilon = e2;
    }
  } 
  else
  {
    epsilon = e12;
  }
  return epsilon;
}

double GetLJsigma(int type1, int type2)
{
  double sigma;
  if (type1 == type2)
  {
    if (type1 == 1.0) 
    {
      sigma   = s1;
    }
    else
    {
      sigma = s2;
    }
  } 
  else
  {
    sigma = s12;
  }
  return sigma;
}

void Compute_Meso_Force(gsl_matrix * Positions, gsl_matrix * Forces, gsl_vector * z, gsl_matrix * MesoForce)
{
  double zi, fx, fy, fz;
  int muRight, muLeft;
  double dv = ((float) Lx * Ly * RealLz) / NNodes;
  double dz = ((float) RealLz) / NNodes;

  for (int i=0;i<NParticles;i++)
  {
      zi      = gsl_matrix_get(Positions,i,3);
      fx      = gsl_matrix_get(Forces,i,0);
      fy      = gsl_matrix_get(Forces,i,1);
      fz      = gsl_matrix_get(Forces,i,2);
      
      muRight = (int) floor(zi*NNodes/RealLz);        
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

void Compute_Meso_Sigma1 (gsl_matrix * Positions, gsl_matrix * Velocities, gsl_matrix * MesoSigma1)
{
  int mu = 0;
  double mass = 0.0;
  gsl_matrix_scale(MesoSigma1, 0.0);
  double dv = ((float) Lx * Ly * RealLz) / NNodes;
  for (int i=0;i<NParticles;i++)
  {
    mu = floor(gsl_matrix_get(Positions,i,3)*NNodes/RealLz) - 1;
    ( mu == -1 ) ? mu = NNodes-1 : mu ;
    mass = ( gsl_matrix_get(Positions,i,0) == 1 ? m1 : m2 );

    MesoSigma1->data[mu*MesoSigma1->tda+0] += mass * pow(gsl_matrix_get(Velocities,i,0),2);
    MesoSigma1->data[mu*MesoSigma1->tda+1] += mass * pow(gsl_matrix_get(Velocities,i,1),2);
    MesoSigma1->data[mu*MesoSigma1->tda+2] += mass * pow(gsl_matrix_get(Velocities,i,2),2);
  
  }
  gsl_matrix_scale(MesoSigma1,1.0/dv);

}

void Compute_Meso_Sigma2 (gsl_matrix * Positions, gsl_matrix * Neighbors, gsl_vector * ListHead, gsl_vector * List, gsl_matrix * MesoSigma2)
{

  double dv = ((float) Lx * Ly * RealLz) / NNodes;

  gsl_matrix_scale(MesoSigma2,0.5/dv);
}
          
double Compute_Force_ij (gsl_matrix * Positions, int i, int j, double * fij)
{
   double r2i, r6i, ff;
   double eij = 0.0;

   double deltax  = Positions->data[i*Positions->tda + 1] - Positions->data[j*Positions->tda + 1];
          deltax -= Lx*round(deltax/Lx);
   double deltay  = Positions->data[i*Positions->tda + 2] - Positions->data[j*Positions->tda + 2];
          deltay -= Ly*round(deltax/Ly);
   double deltaz  = Positions->data[i*Positions->tda + 3] - Positions->data[j*Positions->tda + 3];
          deltaz -= RealLz*round(deltaz/RealLz);

   double r2      = deltax*deltax + deltay*deltay + deltaz*deltaz;

   double * lj = malloc(3*sizeof(float));
   lj = GetLJParams(gsl_matrix_get(Positions,i,0), gsl_matrix_get(Positions,j,0));
   double epsilon = lj[0];
   double sigma   = lj[1];
   double ecut    = lj[2];
   free (lj);

   if (r2 <= Rcut*Rcut)
   {
     r2i      = sigma*sigma/r2;
     r6i      = pow(r2i,3);
     ff       = 48.0*epsilon*r2i*r6i*(r6i-0.5);
     if (( ((int) gsl_matrix_get(Positions,i,0)) == 1)&&(((int) gsl_matrix_get(Positions,j,0) == 2)))
     {
       fij[0] = ff*deltax;  
       fij[1] = ff*deltay;  
       fij[2] = ff*deltaz;  
     }
     eij = 4.0*epsilon*r6i*(r6i-1.0)-ecut;
   }

   return eij;
}

double KineticEnergy (gsl_vector * i, int type)
{
  double K = pow(gsl_vector_get(i,0),2) + pow(gsl_vector_get(i,1),2) + pow(gsl_vector_get(i,2),2);
  
  if (type == 1)
  {
    K *= m1;
  }
  else
  {
    K *= m2;
  }
  return K;
   
}

void Compute_Meso_Energy(gsl_matrix * Micro, gsl_vector * MicroEnergy, gsl_vector * z, gsl_vector * MesoEnergy)
{
    double dv = ((float) Lx * Ly * RealLz) / NNodes;
    double dz = ((float) RealLz) / NNodes;
    double zi, ei;
    int muRight, muLeft;

    for (int i=0;i<NParticles;i++)
    {
        zi      = gsl_matrix_get(Micro,i,3);
        ei      = gsl_vector_get(MicroEnergy,i);
        muRight = (int) floor(zi*NNodes/RealLz);        
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

gsl_vector * Compute_Velocity_Module (gsl_matrix * Velocities)
{
  gsl_vector * vel = gsl_vector_calloc (NParticles);

  double vx, vy, vz;
  for (int i=0;i<NParticles;i++)
  {
    vx = gsl_matrix_get(Velocities,i,0);
    vy = gsl_matrix_get(Velocities,i,1);
    vz = gsl_matrix_get(Velocities,i,2);

    gsl_vector_set(vel,i,sqrt(vx*vx+vy*vy+vz*vz));
  }

  return vel;
}
