#include "../globals.h"

void setup_Blob() {
  Problem.Boxsize[0] = 40.0;
  Problem.Boxsize[1] = 10.0;
  Problem.Boxsize[2] = 10.0;

  sprintf(Problem.Name, "IC_Blob");

  const double rho = 10;

  Problem.Rho_Max = rho;

  Density_Func_Ptr = &Blob_Density;
  U_Func_Ptr = &Blob_U;
  Velocity_Func_Ptr = &Blob_Velocity;
}

/* At first we set up a constant density in the Box */
float Blob_Density(const int ipart, const double bias) {
  // cloud is centred at (5, 5, 5) with radius 1.
  double const x = P[ipart].Pos[0] - 5.0;
  double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
  double const z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;

  double Radius = sqrt(x * x + y * y + z * z);

  if (Radius < 1.0) {
    return 10;
  } else {
    return 1;
  }
}

void Blob_Velocity(const int ipart, double out[3]) {
  // cloud is centred at (5, 5, 5) with radius 1.
  double const x = P[ipart].Pos[0] - 5.0;
  double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
  double const z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;

  double Radius = sqrt(x * x + y * y + z * z);

  if (Radius < 1.0) {
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
  } else {
    out[0] = 27.0; // corresponds to Mach number 2.7
    out[1] = 0.0;
    out[2] = 0.0;
  }
}

/* We set up the internal energy via the pressure profile and ideal equation of
 * state */

float Blob_U(const int ipart) {
  // cloud is centred at (5, 5, 5) with radius 1.
  double const x = P[ipart].Pos[0] - 5.0;
  double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
  double const z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;

  double Radius = sqrt(x * x + y * y + z * z);

  // pressure equilibrium
  if (Radius < 1.0) {
    return 9;
  } else {
    return 90;
  }
}
