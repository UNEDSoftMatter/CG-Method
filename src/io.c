/*
 * Filename   : io.c
 *
 * Created    : 19.04.2016
 *
 * Modified   : jue 16 jun 2016 15:29:16 CEST
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

void SaveVectorWithIndex(char * basename, char * filename, gsl_vector * z1, gsl_vector * z2)
{
  
  char str[100]; 
  strcpy (str, "./output/");
  strcat (str, basename);
  strcat (str, filename);
  
  int NRows = z1->size;

  FILE * iFile = fopen(str, "w");
  for (int i=0;i<NRows;i++)
  {
    fprintf(iFile, "%8.6e\t %8.6e\n",gsl_vector_get(z1,i), gsl_vector_get(z2,i));
  }
  fclose(iFile);
}

void SaveVectorWithoutIndex(gsl_vector * z, char * File)
{
    FILE *iFile;
    iFile = fopen(File, "w");
    int NRows = z->size;
    for (int i=0;i<NRows;i++)
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
  printf("# Lastest stable version can be found on                                     #\n");
  printf("#    https://github.com/UNEDSoftMatter/CG-Method.git                         #\n");
  printf("#                                                                            #\n");
  printf("##############################################################################\n");
  printf("#                                                                            #\n");
  printf("# If called without arguments, CG-Method will create snapshots from lammps   #\n");
  printf("# output files (see README.md to find information about the correct format)  #\n");
  printf("#                                                                            #\n");
  printf("# If called with an argument, CG-Method will use the argument as an input    #\n");
  printf("# file which contains the list of snapshots                                  #\n");
  printf("#                                                                            #\n");
  printf("# In both cases, the snapshots are (or will be) located in                   #\n");
  printf("# data/positions/x????? with the format 'ID TYPE x y z' for the positions,   #\n");
  printf("# and in data/velocities/*.vel with the format 'ID TYPE vx vy vz' for the    #\n");
  printf("# velocities                                                                 #\n");
  printf("#                                                                            #\n");
  printf("# (See README.md to know how to format lammps dump files)                    #\n");
  printf("#                                                                            #\n");
  printf("##############################################################################\n");
  printf("#                                                                            #\n");
  printf("# Mesoscopic profiles will be stored in ./output/                            #\n");
  printf("#                                                                            #\n");
  printf("# Microscopic visualization (if any) will be stored in ./povray/             #\n");
  printf("#                                                                            #\n");
  printf("# Log files (if any) will be stored in ./log/                                #\n");
  printf("#                                                                            #\n");
  printf("##############################################################################\n\n");
}

void ReadInputFiles(char * iFileStr, char iFiles[][7])
{
    FILE *iFile;
    iFile = fopen(iFileStr, "r");
    if (!iFile)
    {
        PrintMsg("Error reading input file. Exiting now...");
        printf("\tThe input file was: %s\n", iFileStr);
    }

    int line=NSteps;
    int i=0;

    while (line--)
    {
        fscanf(iFile,"%s",iFiles[i]);
        i++;
    }

    fclose(iFile);
}

void PrintInfo(int Step, gsl_vector * vector, FILE* fileptr)
{
  fprintf(fileptr, "%10d", Step);

  for (int i=0;i<vector->size;i++)
    fprintf(fileptr, "\t%8.6e", gsl_vector_get(vector,i));

  fprintf(fileptr,"\n");
}

void PrepareInputFiles(void)
{
  struct stat status;
 
  if (stat("./data", &status) != 0) 
    mkdir("./data",0755);
  
  if (stat("./data/positions",  &status) == 0)
    system("rm -r ./data/positions");
  mkdir("./data/positions",0755);
  
  if (stat("./data/velocities", &status) == 0)
    system("rm -r ./data/velocities");
  mkdir("./data/velocities",0755);

  Split_File("data/positions",  PositionsFileStr);
  Split_File("data/velocities", VelocitiesFileStr);

  // Create snapshot list
  PrintMsg("Creating snapshot list in file 'sim'...");
  system("for i in $(ls data/positions/); do echo $i; done > sim");
}

void Split_File(char *directory, char *iFile)
{
  // (via) https://www.codingunit.com/c-tutorial-splitting-a-text-file-into-multiple-files
  FILE *ptr_iFile;
  FILE *ptr_oFile;
  char line[128];
  char oFileName[100];
  int  filecounter=0; 
  int  linecounter=0;
  int  SizeOfChunk = NParticles+9;

  ptr_iFile = fopen(iFile, "r");
  sprintf(oFileName, "%s/x%05d", directory, filecounter);
  ptr_oFile = fopen(oFileName, "w");
  while (fgets(line, sizeof(line), ptr_iFile) != NULL)
  {

    if (linecounter == SizeOfChunk)
    {
      fclose(ptr_oFile);
      linecounter = 0;
      filecounter++;
      sprintf(oFileName, "%s/x%05d", directory, filecounter);
      ptr_oFile = fopen(oFileName, "w");
    }
    
    if (linecounter > 8 )
      fprintf(ptr_oFile,"%s", line);

    linecounter++;
  }
  fclose(ptr_iFile);
}

void PrintComputingOptions(void)
{
  PrintMsg("Computations to be done:");
  
  printf("\tMesoscopic densities:\t\t\t\t");
  #if __COMPUTE_DENSITY__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMesoscopic forces exerted by walls:\t\t");
  #if __COMPUTE_FORCE__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMesoscopic kinetic and total energies:\t\t");
  #if __COMPUTE_FORCE__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMesoscopic temperature:\t\t\t\t");
  #if __COMPUTE_TEMPERATURE__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMesoscopic kinetic and virial stress tensor:\t");
  #if __COMPUTE_STRESS__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMesoscopic momentum:\t\t\t\t");
  #if __COMPUTE_MOMENTUM__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMesoscopic velocities:\t\t\t\t");
  #if __COMPUTE_VELOCITY__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMesoscopic internal energy:\t\t\t");
  #if __COMPUTE_INTERNAL_ENERGY__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMacroscopic energy:\t\t\t\t");
  #if __COMPUTE_MACRO_ENERGY__
    printf("true\n");
  #else
    printf("false\n");
  #endif
  
  printf("\tMacroscopic momentum:\t\t\t\t");
  #if __COMPUTE_MACRO_MOMENTUM__
    printf("true\n");
  #else
    printf("false\n");
  #endif

  printf("\tCenter of mass:\t\t\t\t");
  #if __COMPUTE_CENTER_OF_MASS__
    printf("true\n");
  #else
    printf("false\n");
  #endif

}

void PrintScalarWithIndex(int Step, double Value, FILE*fileptr)
{
    fprintf(fileptr, "%10d\t%8.10e\n", Step, Value);
}
