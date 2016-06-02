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

// ¿Cómo sé la masa que tienen los átomos? ¿Y TotalMass?
// Estoy operando con CenterOfMass como si fuera un vector, aunque es un array ¿?¿?¿?
void Compute_CenterOfMass(gsl_matrix * Positions, int type, gsl_vector * Position, double * CenterOfMass)
{
    gsl_vector_zero(CenterOfMass);
    gsl_vector_zero(Position);
    for(int i=0;i<NParticles;i++)
    {
        if ((int)gsl_matrix_get (Positions, i, 0) == type)
            gsl_matrix_get_row(gsl_vector * Position, gls_matrix * Positions, i)
            gsl_vector_scale(gsl_vector * Position , mass)
    
    gsl_vector_add(gsl_vector * CenterOfMass, gsl_vector * Position)
    gsl_vector_scale(gsl_vector * CenterOfMass, 1 / TotalMass)
    }

}
