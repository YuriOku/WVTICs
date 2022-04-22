#include "../globals.h"

// This file sets up the IC for Kelvin Helmholtz intstability test as defined by
// McNally et al (2012).

void setup_Kelvin_Helmholtz_Instability() {
  Problem.Boxsize[0] = 1.;
  Problem.Boxsize[1] = 1.;
  Problem.Boxsize[2] = 0.25;

  sprintf(Problem.Name, "IC_KelvinHelmholtz_000");

  Problem.Rho_Max = 2.0;

  Density_Func_Ptr = &Kelvin_Helmholtz_Instability_Density;
  U_Func_Ptr = &Kelvin_Helmholtz_Instability_U;
  Velocity_Func_Ptr = &Kelvin_Helmholtz_Instability_Velocity;
}

float Kelvin_Helmholtz_Instability_Density(const int ipart, const double bias) {
  double rho1 = 1.0;
  double rho2 = 2.0;
  double rhom = 0.5 * (rho1 - rho2);
  double Lsmooth = 0.025;
  double y = P[ipart].Pos[1] / Problem.Boxsize[1];
  if (y < 0.25) {
    return rho1 - rhom * exp((y - 0.25) / Lsmooth);
    // return rho1;
  } else if (y < 0.5) {
    return rho2 + rhom * exp((-y + 0.25) / Lsmooth);
    // return rho2;
  } else if (y < 0.75) {
    return rho2 + rhom * exp((y - 0.75) / Lsmooth);
    // return rho2;
  } else {
    return rho1 - rhom * exp((-y + 0.75) / Lsmooth);
    // return rho1;
  }
}

float Kelvin_Helmholtz_Instability_U(const int ipart) {
  double pressure = 2.5;
  double rho = Kelvin_Helmholtz_Instability_Density(ipart, 0.0);
  double gamma = 5.0 / 3;
  return pressure / rho * (gamma - 1);
}

void Kelvin_Helmholtz_Instability_Velocity(const int ipart, double out[3]) {
  double v1 = 0.5;
  double v2 = -0.5;
  double vm = 0.5 * (v1 - v2);
  double Lsmooth = 0.025;
  double y = P[ipart].Pos[1] / Problem.Boxsize[1];
  if (y < 0.25) {
    out[0] = v1 - vm * exp((y - 0.25) / Lsmooth);
  } else if (y < 0.5) {
    out[0] = v2 + vm * exp((-y + 0.25) / Lsmooth);
  } else if (y < 0.75) {
    out[0] = v2 + vm * exp((y - 0.75) / Lsmooth);
  } else {
    out[0] = v1 - vm * exp((-y + 0.75) / Lsmooth);
  }

  const double deltaVy = 0.01;
  double x = P[ipart].Pos[0] / Problem.Boxsize[0];
  out[1] = deltaVy * sin(4.0 * M_PI * x);

  out[2] = 0.0;
}
