#include "../globals.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_vector.h>

// Units: kpc, Msun, km/s
// Dims:  2D/3D
// Ref:   Hu et al. (2014)

#ifndef TWO_DIM
    int NNpart = 1; //64;
#else
    int NNpart = 16;
#endif
const float Unit_Mass_in_g = 1.989e33;
const float Unit_Velocity_in_cms = 1.0e5;

const float  Supernova_Energy = 1.0e51 / (Unit_Mass_in_g * Unit_Velocity_in_cms * Unit_Velocity_in_cms); // Total injection energy 
// const float  Supernova_Energy = 6.78e53 / (Unit_Mass_in_g * Unit_Velocity_in_cms * Unit_Velocity_in_cms); // Total injection energy 

void setup_Sedov_Blast()
{
    Problem.Boxsize[0] = 0.2;//3;
    Problem.Boxsize[1] = 0.2;//3;
    Problem.Boxsize[2] = 0.2;//3;

    sprintf ( Problem.Name, "IC_SedovBlast_000" );

    printf("UnitMass_in_g:        %.3e\n", Unit_Mass_in_g);
    printf("UnitVelocity_in_cm/s: %.3e\n", Unit_Velocity_in_cms);
    printf("Supernova_Energy:     %.3e\n", Supernova_Energy);
    printf("NNpart:               %d\n", NNpart);

    const double rho = 1.24E7;

    Problem.Rho_Max = rho;

    Density_Func_Ptr = &Sedov_Blast_Density;
    U_Func_Ptr = &Sedov_Blast_U;
    PostProcessing_Func_Ptr = &Sedov_Blast_PostProcessing;
}

float Sedov_Blast_Density ( const int ipart , const double bias )
{
    return 1.24E7;
}

// This function calculates the distance from each particle to the center of the Box. It returns the distance to the NNpart times furthest neighbour of zero.

float Sedov_Blast_abs ()
{
    gsl_vector *abs_of_zero = gsl_vector_alloc ( Param.Npart );

    for ( int i = 0; i < Param.Npart; i++ ) {

        double x = P[i].Pos[0] - 0.5 * Problem.Boxsize[0];
        double y = P[i].Pos[1] - 0.5 * Problem.Boxsize[1];
#ifndef TWO_DIM
        double z = P[i].Pos[2] - 0.5 * Problem.Boxsize[2];
#else
        double z = 0.0;
#endif

        double r = sqrt (x*x + y*y + z*z);
        gsl_vector_set ( abs_of_zero, i, r );
    }

    gsl_sort_vector ( abs_of_zero );

    double dist = gsl_vector_get ( abs_of_zero, NNpart - 1 );

    gsl_vector_free ( abs_of_zero );

    return dist;

}

#ifdef KINETIC_SEDOV
float Sedov_Blast_kinetic ()
{

    int NNpart = 32;
    Radius =
        double partpos[3][32];

    float maxDistance = Sedov_Blast_abs ();

    if ( Radius <= maxDistance ) {

        for ( int i =; i < 3; i++ ) {

            for ( int j = 0; j < NNpart, j++ ) {

                double pos = P[i].Pos[j];
                partpos[i][j] = pos;

            }

        }

    }

}
#endif // KINETIC_SEDOV


float Sedov_Blast_U ( const int ipart )
{
    float AmbientTemperature = 10.0;
    float HydrogenMass_in_g  = 1.67e-24;
    float Boltzmann_in_cgs   = 1.38e-16;
    float gamma              = 5.0/3.0;
    float u_in_cms2          = (Boltzmann_in_cgs * AmbientTemperature) / (HydrogenMass_in_g * (gamma - 1.0));
    return u_in_cms2 / Unit_Velocity_in_cms / Unit_Velocity_in_cms;
}


//! @todo improvement: use gsl_sort_vector_index
void Sedov_Blast_PostProcessing ()
{
    /*const double u_sn = 4.18971E5;
    int sn_count = pow ( Param.Npart / 3200., 3.0 );
    sn_count = min ( sn_count, 1 );
    //! @todo assign u_sn to sn_count innermost particles*/

    float maxDistance = Sedov_Blast_abs ();

    int count_inj = 0;

    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {
        const double x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
        const double y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
#ifndef TWO_DIM
        const double z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;
        double Radius = sqrt ( x * x + y * y + z * z );
#else
        double Radius = sqrt ( x * x + y * y );
#endif
        if ( Radius <=  maxDistance ) {
            // SphP[ipart].U = Supernova_Energy / (double) NNpart / Problem.Mpart;
            count_inj++;
        }
    }

    printf(" inject energy to %d particles... ", count_inj);
    
    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {
        const double x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
        const double y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
#ifndef TWO_DIM
        const double z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;
        double Radius = sqrt ( x * x + y * y + z * z );
#else
        double Radius = sqrt ( x * x + y * y );
#endif
        if ( Radius <=  maxDistance ) {
            SphP[ipart].U = Supernova_Energy / (float) count_inj / Problem.Mpart;
            // printf("U %.3e, m %.3e\n", Supernova_Energy / (float) count_inj / Problem.Mpart, Problem.Mpart);
        }
    }

}
