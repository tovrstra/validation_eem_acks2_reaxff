

cdef extern from "low.h":
    double taper(double distance, double rcut);

    void _set_physics_eem(double* A, double* gammai2, double* atpos, double rcut,
        int natom, int nstride);
