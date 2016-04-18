/*
 * Filename   : draw.c
 *
 * Created    : 13.04.2016
 *
 * Modified   : mié 13 abr 2016 19:49:59 CEST
 *
 * Author     : jatorre
 *
 * Purpose    : povray files 
 *
 */
#include "cg.h"

void DrawSim(gsl_matrix * Micro, int TestParticle, int TestCell, gsl_vector * NeighboringCells, int * Verlet, int NumberOfNeighbors)
{
  double xi = gsl_matrix_get(Micro,TestParticle,0);
  double yi = gsl_matrix_get(Micro,TestParticle,1);
  double zi = gsl_matrix_get(Micro,TestParticle,2);

  FILE *iFile;
  iFile = fopen("povray/objects.inc", "w");

  fprintf(iFile, "#declare Rcut = %f;\n\n", Rcut);
  
  fprintf(iFile, "#declare Verlet = \n");
  fprintf(iFile, "  sphere{<%f,%f,%f>,%f\n", xi, yi, zi, Rcut);
  fprintf(iFile, "    interior{ior 1.01}\n");
  fprintf(iFile, "    texture{\n");
  fprintf(iFile, "      pigment {color rgbf <0.9,0.9,0.9,.99>}\n");
  fprintf(iFile, "      finish { phong 1 }\n");
  fprintf(iFile, "    }\n");
  fprintf(iFile, "  }\n\n");

  fprintf(iFile, "#declare vertex = union {\n");
  fprintf(iFile, "    cylinder {<0.0,0.0,0.0>,<0.0,Rcut,0.0>,0.01}\n");
  fprintf(iFile, "    cylinder {<Rcut,0.0,0.0>,<Rcut,Rcut,0.0>,0.01}\n");
  fprintf(iFile, "    cylinder {<0.0,0.0,Rcut>,<0.0,Rcut,Rcut>,0.01}\n");
  fprintf(iFile, "    cylinder {<Rcut,0.0,Rcut>,<Rcut,Rcut,Rcut>,0.01}\n");
  fprintf(iFile, "\n");
  fprintf(iFile, "    cylinder {<0.0,0.0,0.0>,<Rcut,0.0,0.0>,0.01}\n");
  fprintf(iFile, "    cylinder {<Rcut,0.0,0.0>,<Rcut,0.0,Rcut>,0.01}\n");
  fprintf(iFile, "    cylinder {<Rcut,0.0,Rcut>,<0.0,0.0,Rcut>,0.01}\n");
  fprintf(iFile, "    cylinder {<0.0,0.0,Rcut>,<0.0,0.0,0.0>,0.01}\n");
  fprintf(iFile, "\n");
  fprintf(iFile, "    cylinder {<0.0,Rcut,0.0>,<Rcut,Rcut,0.0>,0.01}\n");
  fprintf(iFile, "    cylinder {<Rcut,Rcut,0.0>,<Rcut,Rcut,Rcut>,0.01}\n");
  fprintf(iFile, "    cylinder {<Rcut,Rcut,Rcut>,<0.0,Rcut,Rcut>,0.01}\n");
  fprintf(iFile, "    cylinder {<0.0,Rcut,Rcut>,<0.0,Rcut,0.0>,0.01}\n");
  fprintf(iFile, "}\n\n");

  fprintf(iFile, "#declare Boxes = \n");
  fprintf(iFile, "  union{\n");
              
  int    iCell = 0;
  int    CellColor = 0;
  double CellFilter = 0.9;
  int PointBreak = 0;

  for (int i=0;i<Mx;i++)
  {
      for (int j=0;j<My;j++)
      {
          for (int k=0;k<Mz;k++)
          {
              iCell = i + j*Mx + k*Mx*My;
              CellColor = 0;
              CellFilter = 0.9;
              PointBreak = 0;

              for (int l=0;l<27;l++)
              {
                if (iCell == ((int) gsl_vector_get(NeighboringCells,l)))
                {
                    PointBreak = 1;
                    CellFilter = 0.9;
                    if (iCell == TestCell)
                        CellColor = 1;
                } else 
                {
                    if (PointBreak == 0)
                        CellFilter = 0.5;
                }
              }
              fprintf(iFile,"  object{box {<0,0,0>, <%f,%f,%f>} translate <%f,%f,%f> ", Rcut, Rcut, Rcut, i*Rcut, j*Rcut, k*Rcut);
              if (CellColor == 0) 
              {
                fprintf(iFile,"texture{pigment{color YellowGreen filter %f }finish{phong .8}}}\n", CellFilter);
              } else
              {
                fprintf(iFile,"texture{pigment{color Yellow      filter %f }finish{phong .8}}}\n", CellFilter);
              }
              fprintf(iFile,"object {vertex translate <%f,%f,%f> }\n", i*Rcut, j*Rcut, k*Rcut); 
          }
      }
  }
  fprintf(iFile, "  }\n\n");

fprintf(iFile, "#declare Particles = \n");
fprintf(iFile, "union{\n");

    double ix, iy, iz;
    int ColorParticle;

    for (int iPart=0;iPart<NParticles;iPart++)
    {
        ix = gsl_matrix_get(Micro,iPart,0);
        iy = gsl_matrix_get(Micro,iPart,1);
        iz = gsl_matrix_get(Micro,iPart,2);

        ColorParticle = 0;
        for (int j=0;j<NumberOfNeighbors;j++)
            if (Verlet[j] == iPart) ColorParticle = 1;

        if (iPart == TestParticle) ColorParticle = 2;

        if (ColorParticle == 0) 
        {
            fprintf(iFile, "text { ttf \"garamond.ttf\" \"%d\" 0.05, 0 pigment { Black } ", iPart);
        } else if (ColorParticle == 1)
        {
            fprintf(iFile, "text { ttf \"garamond.ttf\" \"%d\" 0.05, 0 pigment { White } ", iPart);
        } else 
        {
            fprintf(iFile, "text { ttf \"garamond.ttf\" \"%d\" 0.05, 0 pigment { Red } ", iPart);
        }
        fprintf(iFile, "scale 0.5 rotate <0,180,0> translate <%f, %f, %f>}\n", ix, iy, iz);
    }
  fprintf(iFile, "  }\n\n");

  fclose(iFile);

}

