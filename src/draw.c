/*
 * Filename   : draw.c
 *
 * Created    : 13.04.2016
 *
 * Modified   : mar 26 abr 2016 14:29:57 CEST
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

void DrawTemperature (gsl_matrix * Micro, gsl_vector * Velocity, char * str)
{

  FILE *iFile;
  iFile = fopen(str, "w");

  fprintf(iFile, "#declare TempParticles = \n");
  fprintf(iFile, "union{\n");

    double ix, iy, iz, iv;

    for (int i=0;i<NParticles;i++)
    {
      ix = gsl_matrix_get(Micro,i,1);
      iy = gsl_matrix_get(Micro,i,2);
      iz = gsl_matrix_get(Micro,i,3);
      iv = gsl_vector_get(Velocity,i);

      fprintf(iFile, "  sphere{<%f,%f,%f>,%f texture{pigment {color rgb <%.2f,0,%.2f> filter 0.8 }finish{phong .7}}}\n", iy, iz, ix, 0.5, iv, (1.0-iv));
    }
 
  fprintf(iFile, "  }\n\n");
  
  fclose(iFile);
}

