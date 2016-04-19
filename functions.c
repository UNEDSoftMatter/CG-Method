/*
 * Filename   : functions.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : mar 19 abr 2016 17:22:43 CEST
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
    
void Compute_Linked_List(gsl_matrix * Micro, gsl_vector * List, gsl_vector * ListHead)
{
    unsigned int    iCell;
    double xi, yi, zi;

    // PARTICLES GO FROM 0 TO N-1
    for (int i=0;i<NParticles;i++)
    {
        xi = gsl_matrix_get(Micro,i,1);
        yi = gsl_matrix_get(Micro,i,2);
        zi = gsl_matrix_get(Micro,i,3);
      
        // CELLS GO FROM 0 TO Mx*My*Mz-1
        // LIST IS FILLED FROM 0 TO N-1
        //
        // This is an example for boxes between 1 to MxMyMz and List from 1 to N
        // iCell = 1 + floor(xi*Mx/Lx) + floor(yi*My/Ly)*Mx + floor(zi*Mz/Lz)*Mx*My; 
        // printf("Particle %d  at (%f, %f, %f) assigned to iCell %d\n", i+1, xi, yi, zi, iCell);
        // List->data[(i+1)*List->stride] = ListHead->data[iCell*ListHead->stride];
        // ListHead->data[iCell*ListHead->stride] = (i+1);

        iCell = floor(xi*Mx/Lx) + floor(yi*My/Ly)*Mx + floor(zi*Mz/Lz)*Mx*My; 
        // printf("Particle %d  at (%f, %f, %f) assigned to iCell %d\n", i, xi, yi, zi, iCell);
        List->data[i*List->stride] = ListHead->data[iCell*ListHead->stride];
        ListHead->data[iCell*ListHead->stride] = i;
    }
}

void Compute_NeighborCells(int cell, gsl_vector * neighbors)
{

    // This is an example in a 4x4x4 box
    // int i = (cell%16)%4;
    // int j = floor((cell%16)/4);
    // int k = floor(cell/16);
    // 
    // int xip =   1;
    // int xim =  -1;
    // int yip =   4;
    // int yim =  -4;
    // int zip =  16;
    // int zim = -16;
    // 
    // if (i == 0) xim =   3;
    // if (i == 3) xip =  -3;
    // if (j == 0) yim =  12;
    // if (j == 3) yip = -12;
    // if (k == 0) zim =  48;
    // if (k == 3) zip = -48;
    
    int i = (cell%(My*Mz))%Mx;
    int j = floor((cell%(Mx*My))/My);
    int k = floor(cell/(Mx*My));
    
    int xip =   1;
    int xim =  -1;
    int yip =   Mx;
    int yim =  -Mx;
    int zip =  Mx*My;
    int zim = -Mx*My;

    if (i == 0)    xim =        (Mx-1);
    if (i == Mx-1) xip =       -(Mx-1);
    if (j == 0)    yim =     Mx*(My-1);
    if (j == My-1) yip =    -Mx*(My-1);
    if (k == 0)    zim =  Mx*My*(Mz-1);
    if (k == Mz-1) zip = -Mx*My*(Mz-1);

    //  Checkpoint. Print index of iCell
    //  printf("i = %d, j = %d, k = %d. CELL = %d\n", i, j, k, cell);

    gsl_vector_set(neighbors, 0,((float) cell));
    gsl_vector_set(neighbors, 1,((float) cell + xip));
    gsl_vector_set(neighbors, 2,((float) cell + xim));

    gsl_vector_set(neighbors, 3,((float) cell       + yip));
    gsl_vector_set(neighbors, 4,((float) cell + xip + yip));
    gsl_vector_set(neighbors, 5,((float) cell + xim + yip));

    gsl_vector_set(neighbors, 6,((float) cell       + yim));
    gsl_vector_set(neighbors, 7,((float) cell + xip + yim));
    gsl_vector_set(neighbors, 8,((float) cell + xim + yim));
    
    gsl_vector_set(neighbors, 9,((float) cell       + zip));
    gsl_vector_set(neighbors,10,((float) cell + xip + zip));
    gsl_vector_set(neighbors,11,((float) cell + xim + zip));

    gsl_vector_set(neighbors,12,((float) cell       + yip + zip));
    gsl_vector_set(neighbors,13,((float) cell + xip + yip + zip));
    gsl_vector_set(neighbors,14,((float) cell + xim + yip + zip));

    gsl_vector_set(neighbors,15,((float) cell       + yim + zip));
    gsl_vector_set(neighbors,16,((float) cell + xip + yim + zip));
    gsl_vector_set(neighbors,17,((float) cell + xim + yim + zip));
    
    gsl_vector_set(neighbors,18,((float) cell       + zim));
    gsl_vector_set(neighbors,19,((float) cell + xip + zim));
    gsl_vector_set(neighbors,20,((float) cell + xim + zim));

    gsl_vector_set(neighbors,21,((float) cell       + yip + zim));
    gsl_vector_set(neighbors,22,((float) cell + xip + yip + zim));
    gsl_vector_set(neighbors,23,((float) cell + xim + yip + zim));

    gsl_vector_set(neighbors,24,((float) cell       + yim + zim));
    gsl_vector_set(neighbors,25,((float) cell + xip + yim + zim));
    gsl_vector_set(neighbors,26,((float) cell + xim + yim + zim));
}

void Compute_NeighborMatrix(gsl_matrix * Neighbors)
{
    gsl_vector * neighborVector = gsl_vector_calloc (27);
    for (int cell=0;cell<Mx*My*Mz;cell++)
    {
        Compute_NeighborCells(cell, neighborVector);
        gsl_matrix_set_row(Neighbors, cell, neighborVector);
    }
    gsl_vector_free(neighborVector);
}
    
int FindParticle(gsl_matrix * Micro,int TestParticle)
{
    int Cell;
    double xi = gsl_matrix_get(Micro,TestParticle,1);
    double yi = gsl_matrix_get(Micro,TestParticle,2);
    double zi = gsl_matrix_get(Micro,TestParticle,3);
    
    Cell = floor(xi*Mx/Lx) + floor(yi*My/Ly)*Mx + floor(zi*Mz/Lz)*Mx*My; 
    
    return Cell;
}

int Compute_VerletList(gsl_matrix * Micro, int TestParticle, gsl_vector * NeighboringCells, 
                       int TestCell, gsl_vector * ListHead, gsl_vector * LinkedList, int * Verlet)
{
    int NumberOfNeighbors = 0;
    int    icell;
    int    j;
    double x0,y0,z0,x1,y1,z1,r;
    double dx, dy, dz;

    x0 = gsl_matrix_get(Micro,TestParticle,1);
    y0 = gsl_matrix_get(Micro,TestParticle,2);
    z0 = gsl_matrix_get(Micro,TestParticle,3);

    for (int i=0;i<27;i++)
    {
        icell = gsl_vector_get(NeighboringCells,i);
        j     = gsl_vector_get(ListHead,icell);
        while (j >= 0)
        {
            x1 = gsl_matrix_get(Micro,j,1);
            dx = x1 - x0;
            (dx > Lx/2) ? dx-= Lx : dx;

            y1 = gsl_matrix_get(Micro,j,2);
            dy = y1 - y0;
            (dy > Ly/2) ? dy-= Ly : dy;

            z1 = gsl_matrix_get(Micro,j,3);
            dz = z1 - z0;
            (dz > Lz/2) ? dz-= Lz : dz;

            r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r <= Rcut)
            { 
                Verlet[NumberOfNeighbors] = j;
                NumberOfNeighbors++;
            }
            j = gsl_vector_get(LinkedList,j);
            // printf("Testing cell %d, particle %d\n", icell, j);
        }
    }

    return NumberOfNeighbors;
}
    
void Compute_Forces(gsl_matrix * Positions, gsl_matrix * Neighbors, gsl_vector * ListHead, 
                    gsl_vector * List, int type1, int type2, gsl_matrix * Forces)
{
  double epsilon = GetLJepsilon(type1, type2);
  double sigma   = GetLJsigma(type1, type2);

  // printf("epsilon = %f, sigma = %f\n", epsilon, sigma);

  gsl_matrix_scale(Forces,0.0);

  #pragma omp parallel
  {
    #pragma omp for
    for (int i=0;i<NParticles;i++)
    {
      // TEST to check the performance of a parallel thread
      // nanosleep((const struct timespec[]){{0, 200000L}}, NULL);

      double xi, yi, zi, xj, yj, zj, dx, dy, dz, r2, r2i, r6i, ff;
      gsl_vector * NeighboringCells = gsl_vector_calloc(27);
      int * Verlet = malloc(27 * NParticles * sizeof(int) / (Mx*My*Mz));

      // printf("[%s]\tComputing Verlet list for particle %d...\n",__TIME__,i);
      int iCell = FindParticle(Positions,i);
      gsl_matrix_get_row(NeighboringCells, Neighbors, iCell);
           
      xi = gsl_matrix_get(Positions,i,1);    
      yi = gsl_matrix_get(Positions,i,2);    
      zi = gsl_matrix_get(Positions,i,3);    
 
      // Checkpoint: Print particles inside a cell and neighboring cells
      //   printf("Particle %d (type %d) at (%f, %f, %f) is in Cell %d\n", i, ((int) gsl_matrix_get(Positions,i,0)), xi, yi, zi, iCell);
      //   printf("Neighboring cells of cell %d are (", iCell);
      //   for (int nu=0;nu<27;nu++)
      //       printf("%d, ",((int) gsl_vector_get(NeighboringCells,nu)));
      //   printf(")\n");
      
      Verlet = realloc(Verlet, 27 * NParticles * sizeof(int) / (Mx * My * Mz));
      int NumberOfNeighbors = Compute_VerletList(Positions, i, NeighboringCells, iCell, ListHead, List, Verlet);
      Verlet = realloc(Verlet, NumberOfNeighbors * sizeof(int));
 
      // Checkpoint: Print Verlet list of a particle
      //   printf("Particle %d has %d neighbors\n", i, NumberOfNeighbors);
      //   for (int i=0;i<NumberOfNeighbors;i++)
      //     printf("%d, ", Verlet[i]);
      //   printf(")\n");

      for (int j=0;j<NumberOfNeighbors;j++)
      {
         if (( ((int) gsl_matrix_get(Positions,i,0)) == type1)&&(((int) gsl_matrix_get(Positions,Verlet[j],0) == type2)))
         {
           xj = gsl_matrix_get(Positions,Verlet[j],1);    
           yj = gsl_matrix_get(Positions,Verlet[j],2);    
           zj = gsl_matrix_get(Positions,Verlet[j],3);    

           dx = (xi - xj);
           dx -= Lx*round(dx/Lx);
         
           dy = (yi - yj);
           dy -= Ly*round(dy/Ly);
         
           dz = (zi - zj);
           dz -= Lz*round(dz/Lz);
 
           r2 = dx*dx + dy*dy + dz*dz;

           if (r2 <= Rcut*Rcut)
           {
             r2i = sigma*sigma/r2;
             r6i = pow(r2i,3);
             ff  = 48.0*epsilon*r2i*r6i*(r6i-0.5);
             
             // printf("Force calculation for particles %d and %d at distance %f\n", i, j, sqrt(r2)); 
             // printf("Force calculation gives r2i %f, r6i %f, ff %f\n", r2i, r6i, ff);
           
             Forces->data[i*Forces->tda + 0] += ff*dx;
             Forces->data[i*Forces->tda + 1] += ff*dy;
             Forces->data[i*Forces->tda + 2] += ff*dz;
         
             //Forces->data[Verlet[j]*Forces->tda + 0] -= ff*dx;
             //Forces->data[Verlet[j]*Forces->tda + 1] -= ff*dy;
             //Forces->data[Verlet[j]*Forces->tda + 2] -= ff*dz;

             //  printf("Force calculation gives fx %f, fy %f, fz %f\n", ff*dx, ff*dy, ff*dz);
             //  printf("Force calculation on %d gives fx %f, fy %f, fz %f\n", Verlet[j], 
             //      gsl_matrix_get(Forces,Verlet[j],0), gsl_matrix_get(Forces,Verlet[j],1), gsl_matrix_get(Forces,Verlet[j],2));
             //  printf("Force calculation on %d gives fx %f, fy %f, fz %f\n", i, 
             //      gsl_matrix_get(Forces,i,0), gsl_matrix_get(Forces,i,1), gsl_matrix_get(Forces,i,2));
           }
         }
      }
    }
  }
}

double GetLJepsilon(int type1, int type2)
{
  double epsilon;
  if (type1 == type2)
  {
    if (type1 == 1) 
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
    if (type1 == 1) 
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
  double dv = ((float) Lx * Ly * RealLz) / NNodes ;
  double dz = ((float) RealLz) / NNodes;
  double zi, fx, fy, fz;
  int muRight, muLeft;

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
