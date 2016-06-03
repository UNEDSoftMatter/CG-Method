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



double Compute_Macro(gsl_vector * Micro, gsl_matrix * Positions, int type, char *str)
{
    double Macro = 0.0;
    for(int i=0;i<NParticles;i++)
    {
        if ((int)gsl_matrix_get (Positions, i, 0) == type)
        {
            if (strcmp(str,"top") == 0)
            {
                if (gsl_matrix_get(Positions,i,3) >= Lz/2.0)
                {
                    //printf("Particle %d at %f >= %f\n", i, gsl_matrix_get(Positions,i,3),Lz/2.0);
                    Macro += gsl_vector_get(Micro, i);
                }
            }
            else if (strcmp(str,"bottom") == 0)
            {
                if (gsl_matrix_get(Positions,i,3) < Lz/2.0)
                {
                   // printf("Particle %d at %f < %f\n", i, gsl_matrix_get(Positions,i,3),Lz/2.0);
                    Macro += gsl_vector_get(Micro, i);
                }
            }
            else
            {
                    Macro += gsl_vector_get(Micro, i);
            }
    
        }
    }
    return Macro;
}


void Compute_CenterOfMass(gsl_matrix * Positions, int type, char *str, gsl_vector * CenterOfMass)
{
    gsl_vector_set_zero(CenterOfMass);
    gsl_vector * Position = gsl_vector_calloc(4);
    double TotalMass = 0.0;
    for(int i=0;i<NParticles;i++)
    {
        if ((int)gsl_matrix_get (Positions, i, 0) == type)
        {
            if (strcmp(str,"top") == 0)
            {
                if (gsl_matrix_get(Positions,i,3) >= Lz/2.0)
                {
                    gsl_matrix_get_row(Position, Positions, i);
                    if ((int)gsl_matrix_get (Positions, i, 0) == 1)
                    {
                        gsl_vector_scale(Position , m1);
                        TotalMass += m1;
                    }
                    else 
                    {
                        gsl_vector_scale(Position , m2);
                        TotalMass += m2;
                    }

                    gsl_vector_add(CenterOfMass, Position);
                }
            }
            else if (strcmp(str, "bottom") == 0)
            {
                if (gsl_matrix_get(Positions,i,3) < Lz/2.0)
                {
                    gsl_matrix_get_row(Position, Positions, i);
                    if ((int)gsl_matrix_get (Positions, i, 0) == 1)
                    {
                        gsl_vector_scale(Position , m1);
                        TotalMass += m1;
                    }
                    else 
                    {
                        gsl_vector_scale(Position , m2);
                        TotalMass += m2;
                    }
                    gsl_vector_add(CenterOfMass, Position);
                }
            }
            else
            {
                gsl_matrix_get_row(Position, Positions, i);
                if ((int)gsl_matrix_get (Positions, i, 0) == 1)
                {
                    gsl_vector_scale(Position , m1);
                    TotalMass += m1;
                }
                else 
                {
                    gsl_vector_scale(Position , m2);
                    TotalMass += m2;
                }
                gsl_vector_add(CenterOfMass, Position);
            }
        }
    }
    gsl_vector_free(Position);
    gsl_vector_scale(CenterOfMass, 1.0 / TotalMass);
}
