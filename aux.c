/*
 * Filename   : aux.c
 *
 * Created    : 25.04.2016
 *
 * Modified   : lun 25 abr 2016 12:39:05 CEST
 *
 * Author     : jatorre
 *
 * Purpose    : Auxiliar functions 
 *
 */

#include "cg.h"

double MaxVector (gsl_vector * v)
{
  double max = 0.0;

  for (int i=0;i<v->size;i++)
    (gsl_vector_get(v,i) > max ) ? max = gsl_vector_get(v,i) : max ;

  return max;
}

double MinVector (gsl_vector * v)
{
  double min = 0.0;

  for (int i=0;i<v->size;i++)
    (gsl_vector_get(v,i) < min ) ? min = gsl_vector_get(v,i) : min ;

  return min;
}

gsl_vector * RescaleVector (gsl_vector * v)
{
  gsl_vector * vrescaled = gsl_vector_calloc (v->size);

  double max = MaxVector(v);
  double min = MinVector(v);
  double vi;

  for (int i=0;i<v->size;i++)
  {
    vi = (gsl_vector_get(v,i)-min)/(max-min);
    gsl_vector_set(vrescaled,i,vi);
  }

  return vrescaled; 

}
