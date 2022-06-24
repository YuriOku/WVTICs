#include "../globals.h"

/* Keplerian disk IC by Hopkins (2015) 
 * This problem requires an external gravitational potential of
 * phi = -(r^2 + e^2)^(-1/2),
 * where e = 0.25 is a softening length.
 */

void setup_Keplerian_Ring ()
{
    Problem.Boxsize[0] = 4.0;
    Problem.Boxsize[1] = 4.0;
    Problem.Boxsize[2] = 0.4; /* two-dimentional problem */

    sprintf ( Problem.Name, "IC_KeplerianRing" );

    const double rho = 0.01 + ( 1 / 0.5 ) * ( 1 / 0.5 );

    Problem.Rho_Max = rho;

    Density_Func_Ptr = &Keplerian_Ring_Density;
    U_Func_Ptr = &Keplerian_Ring_U;
    Velocity_Func_Ptr = &Keplerian_Ring_Velocity;
}

/* At first we set up a constant density in the Box */
float Keplerian_Ring_Density ( const int ipart , const double bias )
{

    const double x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    const double y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    const double Radius = sqrt ( x * x + y * y );

    if ( Radius < 0.5 ) {
        return 0.01 + pow(Radius / 0.5, 3);
    } else if ( Radius >= 0.5 && Radius <= 2.0 ) {
        return 0.01 + 1;
    } else {
        return 0.01 + 1 / pow ( ( 1 + ( Radius - 2 ) ) / 0.1, 3 );
    }

}


/* We set up the internal energy via the pressure profile and ideal equation of state */
float Keplerian_Ring_U ( const int ipart )
{
    const double gamma = 5.0 / 3.0;
    const double Pressure = 1e-6;

    return Pressure / SphP[ipart].Rho / ( gamma - 1 );

}

void Keplerian_Ring_Velocity ( const int ipart, double out[3] )
{
    const double x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    const double y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    const double Radius = sqrt ( x * x + y * y );
    const double epsilon = 0.25;
    double Vc = Radius * pow(Radius * Radius + epsilon * epsilon, -3.0/4.0);

    out[0] = -y / Radius * Vc;
    out[1] =  x / Radius * Vc;
    out[2] = 0.0;

}
