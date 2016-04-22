/*
 * Filename   : io.c
 *
 * Created    : 19.04.2016
 *
 * Modified   : vie 22 abr 2016 12:10:45 CEST
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
        fprintf(iFile, "\t%8.6f",gsl_matrix_get(Matrix,i,j));
      fprintf(iFile, "\n");
    }
}

void SaveVectorWithIndex(gsl_vector * z, gsl_vector * Vector, char * File)
{
    FILE *iFile;
    iFile = fopen(File, "w");
    int Nrows = Vector->size;
    for (int i=0;i<Nrows;i++)
      fprintf(iFile, "%8.6f\t %8.6f\n",gsl_vector_get(z,i), gsl_vector_get(Vector,i));

    fclose(iFile);
}

long timediff(clock_t t1, clock_t t2)
{
  long elapsed;
  elapsed = ((double)t2 -t1) / CLOCKS_PER_SEC * 1000;
  return elapsed;
}
