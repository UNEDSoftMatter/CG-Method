/*
 * Filename   : io.c
 *
 * Created    : 19.04.2016
 *
 * Modified   : mar 26 abr 2016 16:16:02 CEST
 *
 * Author     : jatorre
 *
 * Purpose    : IO functions for CG-Method 
 *
 */
#include "cg.h"

void PrintMsg(char *msg)
{
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  char * timestamp = asctime(timeinfo);
  timestamp[strlen(timestamp)-1] = '\0';
  printf("[%s]\t%s\n", timestamp, msg);
}
    
void SaveMatrixWithIndex(gsl_vector * z, gsl_matrix * Matrix, char * File)
{
    FILE *iFile;
    iFile = fopen(File, "w");
    int Nrows = Matrix->size1;
    int Ncols = Matrix->size2;
    for (int i=0;i<Nrows;i++)
    {
      fprintf(iFile, "%8.6f",gsl_vector_get(z,i));
      for (int j=0;j<Ncols;j++)
        fprintf(iFile, "\t%20.14f",gsl_matrix_get(Matrix,i,j));
      fprintf(iFile, "\n");
    }
    fclose(iFile);
}

void SaveVectorWithIndex(gsl_vector * z1, gsl_vector * z2, int Nrows, char * File)
{
  FILE *iFile;
  iFile = fopen(File, "w");
  for (int i=0;i<Nrows;i++)
  {
    fprintf(iFile, "%8.6f\t %20.15f\n",gsl_vector_get(z1,i), gsl_vector_get(z2,i));
  }
  fclose(iFile);
}

void SaveVectorWithoutIndex(gsl_vector * z, char * File)
{
    FILE *iFile;
    iFile = fopen(File, "w");
    int Nrows = z->size;
    for (int i=0;i<Nrows;i++)
      fprintf(iFile, "%d\t %8.6f\n",i, gsl_vector_get(z,i));

    fclose(iFile);
}

long timediff(clock_t t1, clock_t t2)
{
  long elapsed;
  elapsed = ((double)t2 -t1) / CLOCKS_PER_SEC * 1000;
  return elapsed;
}
    
void PrintInitInfo(void) 
{
  printf("##############################################################################\n");
  printf("#                                                                            #\n");
  printf("# CG-METHOD                                                                  #\n");
  printf("#                                                                            #\n");
  printf("# This program computes mesoscopic variables from microscopic configurations #\n");
  printf("#                                                                            #\n");
  printf("##############################################################################\n");
  printf("#                                                                            #\n");
  printf("# The positions of the particles should be in data/$1.pos in matrixform      #\n");
  printf("#                                                                            #\n");
  printf("# TYPE  x   y   z                                                            #\n");
  printf("#                                                                            #\n");
  printf("# The velocities of the particles should be in data/$1.vel in matrixform     #\n");
  printf("#                                                                            #\n");
  printf("# vx    vy  vz                                                               #\n");
  printf("#                                                                            #\n");
  printf("# (See preprocess.sh to know how to format a lammps dump file)               #\n");
  printf("#                                                                            #\n");
  printf("##############################################################################\n");
  printf("#                                                                            #\n");
  printf("# Mesoscopic profiles will be stored in ./output/                            #\n");
  printf("#                                                                            #\n");
  printf("# Microscopic visualization will be stored in ./povray/                      #\n");
  printf("#                                                                            #\n");
  printf("# Log files will be stored in ./log/                                         #\n");
  printf("#                                                                            #\n");
  printf("##############################################################################\n\n");
}