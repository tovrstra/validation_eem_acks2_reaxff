

#include <cmath>


double taper(double distance, double rcut) {
  double x = distance/rcut;
  double factor = x*x;
  factor *= factor;
  double result = 1.0;
  result -= 35.0*factor;
  factor *= x;
  result += 84.0*factor;
  factor *= x;
  result -= 70.0*factor;
  factor *= x;
  result += 20.0*factor;
  return result;
}

void _set_physics_eem(double* A, double* gammai2, double* atpos, double rcut,
    int natom, int nstride) {
  for (int iatom0=0; iatom0 < natom; iatom0++) {
    for (int iatom1=0; iatom0 < iatom1; iatom1++) {
      double dx = atpos[3*iatom0] - atpos[3*iatom1];
      double dy = atpos[3*iatom0+1] - atpos[3*iatom1+1];
      double dz = atpos[3*iatom0+2] - atpos[3*iatom1+2];
      double distance = sqrt(dx*dx + dy*dy + dz*dz);
      double gammaimix = gammai2[iatom0]*gammai2[iatom1];
      double coulomb = pow(distance*distance*distance + gammaimix*gammaimix*gammaimix, 1.0/3.0);
      coulomb *= taper(distance, rcut);
      A[iatom0*nstride + iatom1] += coulomb;
      A[iatom1*nstride + iatom0] += coulomb;
    }
  }
}
