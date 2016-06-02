/*
 * Filename   : macrofunctions.c
 *
 * Created    : 01.05.2016
 *
 * Modified   : jue 1 may 2016 18:47:47 CEST
 *
 * Author     : DiegoDZ
 *
 * Purpose    : Macroscopic functions for cg.c 
 *
 */
#include "cg.h"



void Compute_MacroEnergy(gsl_vector * MicroEnergy, gsl_matrix * Positions, int type, double MacroEnergy)
{
    MacroEnergy = 0.0;
    for(int i=0;i<NParticles;i++)
    {
        if ((int)gsl_matrix_get (Positions, i, 0) == type)
            MacroEnergy += gsl_vector_get(MicroEnergy, i);
    }
}


void Compute_MacroMomentum(gsl_vector * Mmod, gsl_matrix * Positions, int type, double MacroMomentum)
{
    MacroMomentum = 0.0;
    for(int i=0;i<NParticles;i++)
    {
        if ((int)gsl_matrix_get (Positions, i, 0) == type)
            MacroMomentum += gsl_vector_get(Mmod, i);
    }
}


void Compute_CenterOfMass(gsl_matrix * Positions, int type, double * CenterOfMass)
{

}
