double sumdr, dsumdr;
double sumdr_prev = 0;
int iter = 0;
do {
  sumdr = 0.0;
#pragma omp parallel for shared(hsml, P) \
    schedule(dynamic, nPart / Omp.NThreads / 256) reduction(+:sumdr)
  for (int ipart = 0; ipart < nPart; ipart++) {
    int ngblist[NGBMAX] = {0};
    int ngbcnt = Find_ngb(ipart, SphP[ipart].Hsml, ngblist);

    double dr[3] = {0.0};
    double sumwk;
    for (int i = 0; i < ngbcnt; i++) {  // neighbour loop

      int jpart = ngblist[i];

      if (ipart == jpart) {
        continue;
      }

      double d[3];
      double r2 = 0.0;

      for (int p = 0; p < 3; ++p) {
        d[p] = P[ipart].Pos[p] - P[jpart].Pos[p];

        if (Problem.Periodic[p]) {
          while (d[p] > boxhalf[p]) {  // find closest image
            d[p] -= boxsize[p];
          }
          while (d[p] < -boxhalf[p]) {
            d[p] += boxsize[p];
          }
        }

        r2 += d[p] * d[p];
      }

      Assert(r2 > 0,
             "Found two particles %d & %d at the same location. "
             "Consider increasing the space between your density field"
             " and the box boundaries.",
             ipart, jpart);

      float h = SphP[ipart].Hsml;

      if (r2 > p2(h)) {
        continue;
      }

      float r = sqrt(r2);
      // * p3(h) is a for legacy reasons - at some point retune the code to work
      // without it norm_hsml also plays a minor role with that
#ifdef TWO_DIM
      double kernel_fac = p2(h);
#else
      double kernel_fac = p3(h);
#endif
      float wk = sph_kernel(r, h) * kernel_fac;

      dr[0] += d[0] * wk;
      dr[1] += d[1] * wk;
#ifndef TWO_DIM
      dr[2] += d[2] * wk;
#endif
      sumwk += wk;
    }
    if (fabs(sumwk - 1.0) > 0.1) {
      P[ipart].Pos[0] += 0.1 * dr[0];
      P[ipart].Pos[1] += 0.1 * dr[1];
#ifndef TWO_DIM
      P[ipart].Pos[2] += 0.1 * dr[2];
#endif

            while ( P[ipart].Pos[0] < 0 ) { // keep it in the box
                P[ipart].Pos[0] += boxsize[0];
            }

            while ( P[ipart].Pos[0] > boxsize[0] ) {
                P[ipart].Pos[0] -= boxsize[0];
            }

            while ( P[ipart].Pos[1] < 0 ) {
                P[ipart].Pos[1] += boxsize[1];
            }

            while ( P[ipart].Pos[1] > boxsize[1] ) {
                P[ipart].Pos[1] -= boxsize[1];
            }

            while ( P[ipart].Pos[2] < 0 ) {
                P[ipart].Pos[2] += boxsize[2];
            }

            while ( P[ipart].Pos[2] > boxsize[2] ) {
                P[ipart].Pos[2] -= boxsize[2];
            }
    }
    sumdr += sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
  }
  iter++;
  printf("iter %d sumdr = %g sumdr_prev = %g\n", iter, sumdr, sumdr_prev);
  dsumdr = fabs(sumdr - sumdr_prev) / sumdr;
  sumdr_prev = sumdr;

  Find_sph_quantities();
} while (dsumdr > 1.0e-2);