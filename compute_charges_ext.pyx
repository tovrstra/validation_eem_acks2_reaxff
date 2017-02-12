"""Extensions used by compute_charges."""


import numpy as np
cimport numpy as np
np.import_array()

cimport low

__all__ = ['_set_physics_eem']


def _set_physics_eem(np.ndarray[double, ndim=2] A not None, np.ndarray[double, ndim=1] B not None, atsymbols,
                     np.ndarray[double, ndim=2] atpositions not None, np.ndarray cellvecs,
                     np.ndarray recivecs, np.ndarray repeats, parameters):
    # Some checks
    assert A.flags['C_CONTIGUOUS']
    assert B.flags['C_CONTIGUOUS']
    assert atpositions.flags['C_CONTIGUOUS']
    if cellvecs is not None:
        assert cellvecs.flags['C_CONTIGUOUS']
        assert recivecs.flags['C_CONTIGUOUS']
        assert repeats.flags['C_CONTIGUOUS']

    # Reset matrices
    natom = atpositions.shape[0]
    A[:natom, :natom] = 0.0
    B[:natom] = 0.0

    # Prepare arrays with atomic parameters
    cdef np.ndarray[double, ndim=1] gammai2 = np.zeros(natom)
    for iatom, symbol in enumerate(atsymbols):
        A[iatom, iatom] = 2*parameters['etas'][symbol]
        B[iatom] = -parameters['chis'][symbol]
        gammai2[iatom] = parameters['gammas'][symbol]**(-0.5)

    # Call C++ function
    if cellvecs is None:
        low._set_physics_eem(&A[0, 0], &gammai2[0], &atpositions[0, 0], parameters['rcut'], natom, A.shape[0])
    else:
        raise NotImplementedError
