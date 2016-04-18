/*
 * Filename   : functions.c
 *
 * Created    : 07.04.2016
 *
 * Modified   : vie 15 abr 2016 18:23:00 CEST
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
        gsl_vector_set(z,mu,(double) (mu+1)*Lz/NNodes);
}

void Compute_Meso_Density(gsl_matrix * Micro, gsl_vector * z, gsl_vector * n)
{
    double dv = Lx * Ly * Lz / NNodes;
    double zi;
    int muRight, muLeft;

    for (int i=0;i<NParticles;i++)
    {
        zi      = gsl_matrix_get(Micro,i,2);
        muRight = (int) floor(zi*NNodes/Lz);        
        muLeft  = muRight-1;
        if (muLeft < 0) 
        {
            n->data[muRight*n->stride] += zi/dv;
            n->data[ NNodes*n->stride] += (gsl_vector_get(z,muRight) - zi)/dv;
        
        } else {
            n->data[muRight*n->stride] += (zi -  gsl_vector_get(z,muLeft))/dv;
            n->data[ muLeft*n->stride] += (gsl_vector_get(z,muRight) - zi)/dv;
        }
    }
    gsl_vector_scale(n,dv);
}
    
void Compute_Linked_List(gsl_matrix * Micro, gsl_vector * List, gsl_vector * ListHead)
{
    unsigned int    iCell;
    double xi, yi, zi;

    // PARTICLES GO FROM 0 TO N-1
    for (int i=0;i<NParticles;i++)
    {
        xi = gsl_matrix_get(Micro,i,0);
        yi = gsl_matrix_get(Micro,i,1);
        zi = gsl_matrix_get(Micro,i,2);
      
        // CELLS GO FROM 0 TO Mx*My*Mz-1
        // LIST IS FILLED FROM 0 TO N-1
        //
        // This is an example for boxes between 1 to MxMyMz and List from 1 to N
        // iCell = 1 + floor(xi*Mx/Lx) + floor(yi*My/Ly)*Mx + floor(zi*Mz/Lz)*Mx*My; 
        // printf("Particle %d  at (%f, %f, %f) assigned to iCell %d\n", i+1, xi, yi, zi, iCell);
        // List->data[(i+1)*List->stride] = ListHead->data[iCell*ListHead->stride];
        // ListHead->data[iCell*ListHead->stride] = (i+1);

        iCell = floor(xi*Mx/Lx) + floor(yi*My/Ly)*Mx + floor(zi*Mz/Lz)*Mx*My; 
        printf("Particle %d  at (%f, %f, %f) assigned to iCell %d\n", i, xi, yi, zi, iCell);
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
    int j = floor((cell%(Mx*Mz))/My);
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
    //  printf("i = %d, j = %d, k = %d\n", i, j, k);

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
    double xi = gsl_matrix_get(Micro,TestParticle,0);
    double yi = gsl_matrix_get(Micro,TestParticle,1);
    double zi = gsl_matrix_get(Micro,TestParticle,2);
    
    Cell = floor(xi*Mx/Lx) + floor(yi*My/Ly)*Mx + floor(zi*Mz/Lz)*Mx*My; 
    
    return Cell;
}

int Compute_VerletList(gsl_matrix * Micro, int TestParticle, gsl_vector * NeighboringCells, int TestCell, gsl_vector * ListHead, gsl_vector * LinkedList, int * Verlet)
{
    int NumberOfNeighbors = 0;
    int    icell;
    int    j;
    double x0,y0,z0,x1,y1,z1,r;
    double dx, dy, dz;

    x0 = gsl_matrix_get(Micro,TestParticle,0);
    y0 = gsl_matrix_get(Micro,TestParticle,1);
    z0 = gsl_matrix_get(Micro,TestParticle,2);

    for (int i=0;i<27;i++)
    {
        icell = gsl_vector_get(NeighboringCells,i);
        j     = gsl_vector_get(ListHead,icell);
        while (j >= 0)
        {
            x1 = gsl_matrix_get(Micro,j,0);
            dx = x1 - x0;
            (dx > Lx/2) ? dx-= Lx : dx;

            y1 = gsl_matrix_get(Micro,j,1);
            dy = y1 - y0;
            (dy > Ly/2) ? dy-= Ly : dy;

            z1 = gsl_matrix_get(Micro,j,2);
            dz = z1 - z0;
            (dz > Lz/2) ? dz-= Lz : dz;

            r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r <= Rcut)
            { 
                Verlet[NumberOfNeighbors] = j;
                NumberOfNeighbors++;
            }
            j = gsl_vector_get(LinkedList,j);
        }
    }

    return NumberOfNeighbors;
}
