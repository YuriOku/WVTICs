#include "../globals.h"

/* PARAMETERS */
const double BoxSize = 64;  // size of the simulation box in unit length
const double zinit = 99;          // redshift of the initial condition
const double zcoll = 1;           // redshift of collapse

/* UNITS */
const double UnitLength_in_cm = 3.085678e24;  // Mpc
const double UnitMass_in_g = 1.989e43;        // 1e10 Msun
const double UnitVelocity_in_cm_per_s = 1e5;  // km/s

/* CONSTANTS */
const double GRAVITY = 6.6738e-8;     // gravitational constant in cgs
const double HUBBLE = 3.2407789e-18;  // 100 km/s/Mpc in 1/s

const double UnitDensity_in_cgs =
    UnitMass_in_g / (UnitLength_in_cm * UnitLength_in_cm * UnitLength_in_cm);
const double UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;

// critical density in unit density (independent of hubble parameter)
const double rhocrit =
    3 * HUBBLE * HUBBLE / 8 / pi / GRAVITY / UnitDensity_in_cgs;

void setup_Zeldovich_Pancake() {
  Problem.Boxsize[0] = BoxSize;
  Problem.Boxsize[1] = BoxSize;
  Problem.Boxsize[2] = BoxSize;

  sprintf(Problem.Name, "IC_Zeldovich_Pancake");

  Problem.Rho_Max = rhocrit * 1.1;

  Density_Func_Ptr = &Zeldovich_Pancake_Density;
  U_Func_Ptr = &Zeldovich_Pancake_U;
  Velocity_Func_Ptr = &Zeldovich_Pancake_Velocity;
  PostProcessing_Func_Ptr = &Zeldovich_Pancake_PostProcessing;
}

float Zeldovich_Pancake_Density(const int ipart, const double bias) {
  /* At first we set up a constant density in the Box */
  return rhocrit;
}

void Zeldovich_Pancake_Velocity(const int ipart, double out[3]) {
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;
}

float Zeldovich_Pancake_U(const int ipart) {
  /* internal energy will be set by InitGasTemp */
  return 0.0;
}

void Zeldovich_Pancake_PostProcessing() {
  /* Fourier mode of boxsize */
  double k = 2 * pi / BoxSize;

  /* collapse at x = Lbox/2.
   * note that this set up is independent of the Hubble parameter due to h in length units. */

  /* compute velocity before displacement */
  for (int i = 0; i < Param.Npart; i++) {
    P[i].Vel[0] =
        100.0 * (1 + zcoll) / sqrt(1 + zinit) * sin(k * P[i].Pos[0]) / k;

    /* convert to gadget's velocity */
    double ainit = 1 / (1 + zinit);
    P[i].Vel[0] /= sqrt(ainit);
  }

  /* displace particle using Zeldovich approximation */
  for (int i = 0; i < Param.Npart; i++) {
    P[i].Pos[0] += (1 + zcoll) / (1 + zinit) * sin(k * P[i].Pos[0]) / k;
  }
}
