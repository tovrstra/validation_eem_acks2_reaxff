#!/usr/bin/env python3
"""Compute EEM or ACKS2 charges as in ReaxFF. Library functions.

Unlike ReaxFF, this script internally works in atomic units, because that is sooo much
easier. All parameters and geometries are converted to atomic units when they are read.
Results are converted from atomic units when they are printed on screen.
"""

import shlex

import numpy as np

try:
    from compute_charges_ext import _set_physics_eem
except ImportError:
    _set_physics_eem = None
_set_physics_acks2 = None
#_set_physics_eem = None


__all__ = ['angstrom', 'electronvolt', 'kcalmol', 'load_structure', 'load_structure_xyz',
           'load_constraints', 'EEMModel', 'ACKS2Model']


# Unit conversion in this script works as follows: internally atomic units are used. Input
# is converted into atomic units as early as possible, i.e. when reading from files.
# Output is converted from atomic units as late as possible, i.e. when printing on screen.

# To convert to atomic units, multiply by the unit name:
#     d = 5*angstrom   # Assign a distance of 5 angstrom (in atomic units) to d.
# To convert from atomic units, divide by the unit name
#     print(e/kcalmol) # print the energy e in k cal mol^-1.

# Unit conversion constants taken from NIST, which every program would ideally use.
# The calorie is just defined with a four-digit precision relative to the Joule.
angstrom = 1.0/0.52917721067
electronvolt = 1.6021766208e-19/4.359744650e-18
kcalmol = (4184/4.359744650e-18/6.022140857e23)

# Overwrite energy units to be exactly compatible with ReaxFF.
# If you are not comfortable with atomic units, uncomment the following line.
# angstrom = 1.0  # It should give exactly the same result.
kcalmol = 1.0/(angstrom*332.0638)
electronvolt = 1.0/(angstrom*14.40)


def load_structure(fn_struct):
    """Load molecular structure from the file fn_struct."""
    if fn_struct.lower().endswith('.xyz'):
        return load_structure_xyz(fn_struct)
    else:
        raise NotImplementedError


def load_structure_xyz(fn_struct):
    """Load molecular structure from the file fn_struct, in XYZ format."""
    with open(fn_struct) as f:
        natom = int(next(f))
        title = next(f)
        cellvecs = None
        for field in shlex.split(title):
            if field.count('=') == 1:
                key, value = field.split('=')
                if key == 'Lattice':
                    cellvecs = np.array([float(word) for word in value.split()]).reshape(3, 3)*angstrom
        atsymbols = []
        atpositions = np.zeros((natom, 3), dtype=float)
        for iatom in range(natom):
            words = next(f).split()
            atsymbols.append(words[0].lower())
            atpositions[iatom] = [
                float(words[1])*angstrom,
                float(words[2])*angstrom,
                float(words[3])*angstrom,
            ]
    return atsymbols, atpositions, cellvecs


def load_constraints(fn_constraints, qtot, natom):
    """Load constraints from file, if any."""
    result = [(qtot, np.arange(natom))]
    if fn_constraints is None:
        return result
    with open(fn_constraints) as f:
        for line in f:
            line = line[:line.find('#')].strip()
            if len(line) == 0:
                continue
            words = line.split()
            if len(words) == 1:
                raise IOError('Each line in the constraints file must contain at least two "words".')
            charge = float(words[0])
            indexes = np.array([int(word)-1 for word in words[1:]])
            assert (indexes >= 0).all()
            assert (indexes < natom).all()
            result.append((charge, indexes))
    return result


class EEMModel(object):
    """Compute EEM charges as in ReaxFF."""
    _set_physics_fast = _set_physics_eem

    def __init__(self, rcut=None):
        """Initialize the EEMModel object."""
        self.rcut = rcut
        self.gammas = {}
        self.chis = {}
        self.etas = {}
        self.rmin = 0.5*angstrom

    def load_parameters(self, ffield, verbose=False):
        """Load relevant parameters from a ReaxFF parameter file."""
        with open(ffield) as f:
            lines = f.readlines()
            # Parse global parameters
            self.title = lines[0].strip()
            npar_general = int(lines[1].split()[0])
            self.general_pars = []
            for line in lines[2:2+npar_general]:
                part0, sep, part1 = line.partition('!')
                self.general_pars.append((float(part0), part1))
            npar_element = int(lines[2+npar_general].split()[0])
            self.atom_pars = {}
            for ielement in range(npar_element):
                row = []
                words = lines[6+npar_general+4*ielement].split()
                symbol = words[0].lower()
                row.extend([float(word) for word in words[1:]])
                words = lines[7+npar_general+4*ielement].split()
                row.extend([float(word) for word in words])
                words = lines[8+npar_general+4*ielement].split()
                row.extend([float(word) for word in words])
                words = lines[9+npar_general+4*ielement].split()
                row.extend([float(word) for word in words])
                self.atom_pars[symbol] = np.array(row)
            self._extract_parameters(verbose)
            self._check_parameters(sorted(self.atom_pars.keys()), verbose=True)


    def _extract_parameters(self, verbose=False):
        """Extract parameters from a list of lines, loaded from a ReaxFF parameter file."""
        self.rcut = self.general_pars[12][0]*angstrom
        assert self.rcut > 0
        if verbose:
            print('Coulomb cutoff [Å]: {:10.5f}'.format(self.rcut/angstrom))
        for symbol, values in self.atom_pars.items():
            self.gammas[symbol] = values[5]*(1/angstrom)
            self.chis[symbol] = values[13]*electronvolt
            self.etas[symbol] = values[14]*electronvolt

    def _check_parameters(self, symbols, verbose=False):
        result = True
        for symbol in symbols:
            if self.etas[symbol] < 0:
                result = False
                if verbose:
                    print('    eta[{}] = {:.4f} Å^-1 < 0'.format(
                        symbol, self.etas[symbol]*angstrom))
            if self.gammas[symbol] < 0:
                result = False
                if verbose:
                    print('    gamma[{}] = {:.4f} eV < 0'.format(
                        symbol, self.gammas[symbol]/electronvolt))
            if 2*self.etas[symbol] < self.gammas[symbol]:
                result = False
                if verbose:
                    print('WARNING: polarization catastrophe safety check failed for {}'.format(symbol))
                    print('    2*eta[{symbol}]*4*pi*epsilon_0 = {:.4f} Å^-1 < gamma[{symbol}] = {:.4f} Å^-1'.format(
                        2*self.etas[symbol]*angstrom, self.gammas[symbol]*angstrom, symbol=symbol))
                    print('    eta[{symbol}] = {:.4f} eV e^-2 < gamma[{symbol}]/(8*pi*epsilon_0) = {:.4f} ev e^-2'.format(
                        self.etas[symbol]/electronvolt, self.gammas[symbol]/electronvolt/2, symbol=symbol))
        return result

    def _process_cellvecs(self, cellvecs, verbose=False):
        """Compute reciprocal cell vecs and the number neigboring images within cutoff.

        Parameters
        ----------
        cellvecs : np.ndarray, shape=(3, 3), dtype=float
            The matrix with cell vectors, each row is one vector.

        Returns
        -------
        recivecs : np.ndarray, shape=(3, 3), dtype=float
            The inverse of the matrix cellvecs, each column is one vector.
        repeats : np.ndarray, shape(3,), dtype=int
            The number of repetitions of periodic images to consider, i.e. from -repeats
            to +repeats.
        """
        if cellvecs is None:
            recivecs = None
            repeats = np.zeros(3, int)
        else:
            recivecs = np.linalg.inv(cellvecs)
            # Compute the number of images to be considered to include everything within
            # the cutoff sphere.
            spacings = (recivecs**2).sum(axis=0)**(-0.5)
            if verbose:
                print('Crystal plane spacings [Å]: {:10.5f} {:10.5f} {:10.5f}'.format(*(spacings/angstrom)))
            repeats = np.floor(self.rcut/spacings + 0.5).astype(int)
            if verbose:
                print('Supercell for electrostatics: {} {} {}'.format(*(2*repeats+1)))
        return recivecs, repeats

    @staticmethod
    def _set_constraint(A, B, index_con, target, variable_indexes):
        """Impose a (charge) constraint.

        Parameters
        ----------
        A : np.ndarray, shape=(neq, neq), dtype=float
            Matrix with linear coefficients
        B : np.ndarray, shape=(neq,), dtype=float
            Right-hand side
        index_con : int
            Row and column where constraint should be set.
        target : float
            The target value of the constraints.
        variable_indexes: np.ndarray, dtype=int
            Indexes of all variables to be constrained.
        """
        A[index_con, variable_indexes] = 1.0
        A[variable_indexes, index_con] = 1.0
        B[index_con] = target

    @staticmethod
    def _solve_constraints_least_norm(coeffs, targets, verbose=False):
        solution = np.linalg.lstsq(coeffs, targets, rcond=1e-10)[0]
        # check if the equations are solvable.
        residual = np.linalg.norm(np.dot(coeffs, solution) - targets)**2
        if verbose and residual > 1e-5:
            print('WARNING: constraints are inconsistent.')
        return solution

    def _set_physics(self, A, B, atsymbols, atpositions, cellvecs, recivecs, repeats, safe):
        if self._set_physics_fast is None:
            self._set_physics_slow(A, B, atsymbols, atpositions, cellvecs, recivecs, repeats, safe)
        else:
            self._set_physics_fast(A, B, atsymbols, atpositions, cellvecs, recivecs, repeats,
                {'rcut': self.rcut,
                 'gammas': self.gammas,
                 'chis': self.chis,
                 'etas': self.etas,
                }
            )

    def _set_physics_slow(self, A, B, atsymbols, atpositions, cellvecs, recivecs, repeats, safe):
        """Fill matrix elements in A and B due to the physics of the EEM or ACKS2 models.

        Parameters
        ----------
        A : np.ndarray, shape=(N, N), dtype=float
            The matrix with linear coefficients defining the charge equations (xmortr in
            ReaxFF).
        B : np.ndarray, shape=(N,), dtype=float
            The vector with the right-hand-sides of the equations (elcvec in ReaxFF).
        atsymbols : list of str
            List with atomic symbols.
        atpositions : np.ndarray, shape=(N, 3), dtype=float
            The atomic positions.
        cellvecs : np.ndarray, shape=(3, 3), dtype=float
            The matrix with cell vectors, each row is one vector.
        recivecs : np.ndarray, shape=(3, 3), dtype=float
            The inverse of the matrix cellvecs, each column is one vector.
        repeats : np.ndarray, shape(3,), dtype=int
            The number of repetitions of periodic images to consider, i.e. from -repeats
            to +repeats.
        safe : bool
            Set to True when the parameters are guaranteed to be safe. The positive
            definiteness of the A matrix is not checked when safe==True.
        """
        natom = len(atsymbols)
        for iatom0 in range(natom):
            self._set_physics_atom(A, B, atsymbols, iatom0, natom)
            for iatom1 in range(natom):
                # Get the (naive) minimum image convention
                delta0 = atpositions[iatom0] - atpositions[iatom1]
                # Apply minimum image convention
                if recivecs is not None:
                    delta0_frac = np.dot(recivecs.T, delta0)
                    delta0_frac -= np.round(delta0_frac)
                    delta0 = np.dot(cellvecs.T, delta0_frac)
                # Loop over all relevant neighboring cells
                for image_a in range(-repeats[0], repeats[0]+1):
                    for image_b in range(-repeats[1], repeats[1]+1):
                        for image_c in range(-repeats[2], repeats[2]+1):
                            central = image_a == 0 and image_b == 0 and image_c == 0
                            if central and iatom0 == iatom1:
                                continue
                            if central:
                                delta = delta0
                            else:
                                # Add linear combination of cell vectors
                                delta = delta0 + np.dot(cellvecs.T, [image_a, image_b, image_c])
                            distance = np.linalg.norm(delta)
                            if distance >= self.rcut:
                                continue
                            assert distance > self.rmin
                            self._set_physics_atom_pair(A, B, atsymbols, iatom0, iatom1, natom, distance)
        # Check eigenvalues of the hardness matrix
        if not safe and A is not None:
            evals = np.linalg.eigvalsh(A[:natom, :natom])
            if evals.min() <= 0:
                print('Hardness matrix has non-positive eigenvalues, which is not good!')
                print('(This is the polarization catastrophe.)')
                for e in evals:
                    print('{:10.5f}'.format(e))

    def _set_physics_atom(self, A, B, atsymbols, iatom, natom):
        """Set matrix elements related to individual atoms."""
        B[iatom] = -self.chis[atsymbols[iatom]]
        A[iatom, iatom] = 2*self.etas[atsymbols[iatom]]

    def _set_physics_atom_pair(self, A, B, atsymbols, iatom0, iatom1, natom, distance):
        """Set matrix elements related to atom pairs."""
        gamma0 = self.gammas[atsymbols[iatom0]]
        gamma1 = self.gammas[atsymbols[iatom1]]
        # In atomic units, 1/(4*pi*epsilon_0) is numerically equal
        # to one, which is nice!
        coulomb = (distance**3 + (gamma0*gamma1)**(-3.0/2.0))**(-1.0/3.0)
        coulomb *= self.taper(distance)
        A[iatom0, iatom1] += coulomb

    def compute_charges(self, atsymbols, atpositions, cellvecs, constraints,
                        reduce_constraints, verbose=False, safe=False):
        """Compute atomic charges for a single molecule or crystal.

        Parameters
        ----------
        atsymbols : list of str
            Atomic symbols, should correspond to symbols in ffield.
        atpositions : np.ndarray, dtype=float, shape=(natom, 1)
            Atomic positions in atomic units.
        cellvecs : np.ndarray, dtype=float, shape=(3, 3)
            Rows of this matrix are cell vectors. This argument may also be None.
        constraints : list of tuples
            Each element is a tuple of a charge and a list of integer indexes defining the
            group whose charge is to be constrained.
        reduce_constraints : bool
            Ignored.
        verbose : bool
            When True, intermediate results are printed.
        safe : bool
            Set to True when the parameters are guaranteed to be safe. The positive
            definiteness of the A matrix is not checked when safe==True.
        """
        natom = len(atsymbols)
        ncon = len(constraints)

        recivecs, repeats = self._process_cellvecs(cellvecs, verbose)

        # Build up the equations to be solved
        B = np.zeros(natom + ncon, float)
        A = np.zeros((natom + ncon, natom + ncon), float)
        for icon, (charge, indexes) in enumerate(constraints):
            self._set_constraint(A, B, natom + icon, charge, indexes)
        self._set_physics(A, B, atsymbols, atpositions, cellvecs, recivecs, repeats, safe)

        # Check the consistency of the constraints
        self._solve_constraints_least_norm(A[natom:natom+ncon, :natom], B[natom:natom+ncon], verbose)

        # Solve the charges
        try:
            solution = np.linalg.lstsq(A, B, rcond=1e-10)[0]
        except ValueError:
            print(A)
            raise
        charges = solution[:natom]

        # Compute the energy contributions
        terms = [
            -np.dot(B[:natom], charges),
            +0.5*np.dot(np.diag(A[:natom,:natom]), charges*charges),
            0.5*(
                np.dot(charges, np.dot(A[:natom,:natom], charges))
                - np.dot(np.diag(A[:natom,:natom]), charges*charges)
            )
        ]

        # Print out intermediate result.
        if verbose:
            # A
            A_copy = A.copy()
            A_copy[:natom, :natom] /= electronvolt
            print('A (xmortr) [mixed ReaxFF units, based on electronvolt]')
            print(A_copy)
            # B
            B_copy = B.copy()
            B_copy[:natom] /= electronvolt
            print('B (elcvec) [mixed ReaxFF units, based on electronvolt]')
            print(B_copy)
            # solution
            solution_copy = solution.copy()
            solution_copy[-1] /= electronvolt
            print('solution [mixed ReaxFF units, based on electronvolt]')
            print(solution_copy)
            # residual
            residual = np.dot(A, solution) - B
            residual[:natom] /= electronvolt
            print('residual [mixed ReaxFF units, based on electronvolt]')
            print(residual)

            # Print the energy contributions
            print('Energy 0 eneg          [k cal mol^-1]: {:10.5f}'.format(terms[0]/kcalmol))
            print('Energy 1 hard          [k cal mol^-1]: {:10.5f}'.format(terms[1]/kcalmol))
            print('Energy 2 coul          [k cal mol^-1]: {:10.5f}'.format(terms[2]/kcalmol))
            print('Following ReaxFF conventions...')
            print('Energy 0+1 Charges     [k cal mol^-1]: {:10.5f}'.format(
                (terms[0]+terms[1])/kcalmol))
            print('Energy 2 Coulomb       [k cal mol^-1]: {:10.5f}'.format(
                (terms[2])/kcalmol))

        # The total energy
        energy = sum(terms)

        return energy, charges

    def taper(self, distance):
        """Taper correction as in ReaxFF."""
        TAP7 = 20.0/self.rcut**7
        TAP6 = -70.0/self.rcut**6
        TAP5 = 84.0/self.rcut**5
        TAP4 = -35.0/self.rcut**4
        return 1.0 + TAP4*distance**4 + TAP5*distance**5 + TAP6*distance**6 + TAP7*distance**7


class ACKS2Model(EEMModel):
    """Compute ACKS2 charges as in ReaxFF."""
    _set_physics_fast = _set_physics_acks2

    def __init__(self, rcut=None):
        """Initialize the ACKS2Model object."""
        EEMModel.__init__(self, rcut)
        self.bsoft_amp = None
        self.bsoft_radii = {}

    def _extract_parameters(self, verbose):
        """Extract parameters from a list of lines, loaded from a ReaxFF parameter file."""
        EEMModel._extract_parameters(self)
        self.bsoft_amp = self.general_pars[34][0]/electronvolt
        for symbol, values in self.atom_pars.items():
            self.bsoft_radii[symbol] = values[22]*angstrom

    def _check_parameters(self, symbols, verbose):
        """Extract parameters from a list of lines, loaded from a ReaxFF parameter file."""
        result = EEMModel._check_parameters(self, symbols, verbose)
        if self.bsoft_amp <= 0:
            result = False
            if verbose:
                print('The bond softness amplitude is not positive, while it should.')
        for symbol in symbols:
            if self.bsoft_radii[symbol] < 0:
                result = False
                if verbose:
                    print('    bsoft_radius[{}] = {:.4f} Å < 0'.format(
                        symbol, self.bsoft_radii[symbol]/angstrom))
        return result

    def _set_physics_atom(self, A, B, atsymbols, iatom, natom):
        """Set matrix elements related to individual atoms."""
        EEMModel._set_physics_atom(self, A, B, atsymbols, iatom, natom)
        A[iatom, iatom+natom] = 1.0
        A[iatom+natom, iatom] = 1.0

    def _set_physics_atom_pair(self, A, B, atsymbols, iatom0, iatom1, natom, distance):
        """Set matrix elements related to atom pairs."""
        EEMModel._set_physics_atom_pair(self, A, B, atsymbols, iatom0, iatom1, natom, distance)
        bsoft_rcut = (self.bsoft_radii[atsymbols[iatom0]] + self.bsoft_radii[atsymbols[iatom1]])/2
        assert bsoft_rcut < self.rcut
        if distance < bsoft_rcut:
            x = distance/bsoft_rcut
            bsoft = self.bsoft_amp*x**3*(1-x)**6
            bsoft *= 0.5  # compensate for double counting
            A[iatom0 + natom, iatom0 + natom] -= bsoft
            A[iatom1 + natom, iatom1 + natom] -= bsoft
            A[iatom0 + natom, iatom1 + natom] += bsoft
            A[iatom1 + natom, iatom0 + natom] += bsoft

    def _reduce_constraints(self, A, B, natom, ncon):
        """Implement charge constraints in ACKS2 by zero-ing X matrix elements."""
        # A) Rewrite the first (total-charge) constraint to no longer overlap with the
        # remaining ones
        A_con = A[2*natom: 2*natom+ncon, :natom].copy()
        for icon in range(1, ncon):
            if (A_con[0] >= A_con[icon]).all():
                A_con[0] -= A_con[icon]
            else:
                raise RuntimeError('Could not reduce constraints.')
        # B) Find elements of X matrix that can remain non-zero
        mask = abs(np.dot(A_con.T, A_con)) > 1e-5
        # C) Apply mask
        A[natom:2*natom, natom:2*natom] *= mask
        for iatom in range(natom):
            A[natom + iatom, natom + iatom] = 0.0
            A[natom + iatom, natom + iatom] -= A[natom + iatom, natom:2*natom].sum()
        # D) Reduce size of equations
        if False:
            return A, B
        else:
            A[2*natom+1] = A[-1]
            A[:, 2*natom+1] = A[:, -1]
            B[2*natom+1] = B[-1]
            return A[:2*natom+2, :2*natom+2], B[:2*natom+2]

    def compute_charges(self, atsymbols, atpositions, cellvecs, constraints, reduce_constraints, verbose, safe=False):
        """Compute atomic charges for a single molecule or crystal.

        Parameters
        ----------
        atsymbols : list of str
            Atomic symbols, should correspond to symbols in ffield.
        atpositions : np.ndarray, dtype=float, shape=(natom, 1)
            Atomic positions in atomic units.
        cellvecs : np.ndarray, dtype=float, shape=(3, 3)
            Rows of this matrix are cell vectors. This argument may also be None.
        constraints : list of tuples
            Each element is a tuple of a charge and a list of integer indexes defining the
            group whose charge is to be constrained.
        reduce_constraints : bool
            Try to eliminate the charge constraints.
        verbose : bool
            When True, intermediate results are printed.
        safe : bool
            Set to True when the parameters are guaranteed to be safe. The positive
            definiteness of the A matrix is not checked when safe==True.
        """
        natom = len(atsymbols)
        ncon = len(constraints)

        recivecs, repeats = self._process_cellvecs(cellvecs)

        # Build up the equations to be solveds
        B = np.zeros(2*natom + ncon + 1, float)
        A = np.zeros((2*natom + ncon + 1, 2*natom + ncon + 1), float)
        # Constraints
        for icon, (charge, indexes) in enumerate(constraints):
            self._set_constraint(A, B, 2*natom + icon, charge, indexes)
        # Set reference charges to least-norm solution that satisfies the constraints
        B[natom:2*natom] = self._solve_constraints_least_norm(
            A[2*natom:2*natom+ncon, :natom],
            B[2*natom:2*natom+ncon], verbose)

        # Single constraint for the effective potentials
        self._set_constraint(A, B, 2*natom + ncon, 0.0, np.arange(natom, 2*natom))
        # Physics
        self._set_physics(A, B, atsymbols, atpositions, cellvecs, recivecs, repeats, safe)

        # Optionally simplify the equations
        if reduce_constraints:
            A, B = self._reduce_constraints(A, B, natom, ncon)

        # Solve
        try:
            solution = np.linalg.lstsq(A, B, rcond=1e-10)[0]
        except ValueError:
            print(A)
            raise
        charges = solution[:natom]
        potentials = solution[natom:2*natom]

        # Compute the energy contributions
        terms = [
            -np.dot(B[:natom], charges),
            0.5*np.dot(np.diag(A[:natom,:natom]), charges*charges),
            0.5*(
                np.dot(charges, np.dot(A[:natom,:natom], charges))-
                np.dot(np.diag(A[:natom,:natom]), charges*charges)
            ),
            np.dot(charges - B[natom:2*natom], potentials),
            0.5*np.dot(potentials, np.dot(A[natom:2*natom,natom:2*natom], potentials)),
        ]

        # Print out intermediate result.
        if verbose:
            # A
            A_copy = A.copy()
            A_copy[:natom, :natom] /= electronvolt
            A_copy[natom:2*natom, natom:2*natom] *= electronvolt
            print('A (xmortr) [mixed ReaxFF units, based on electronvolt]')
            print(A_copy)
            # B
            B_copy = B.copy()
            B_copy[:natom] /= electronvolt
            print('B (elcvec) [mixed ReaxFF units, based on electronvolt]')
            print(B_copy)
            # solution
            solution_copy = solution.copy()
            solution_copy[natom:2*natom] /= electronvolt
            solution_copy[-2] /= electronvolt
            print('solution [mixed ReaxFF units, based on electronvolt]')
            print(solution)
            # residual
            residual = np.dot(A, solution) - B
            residual[:natom] /= electronvolt
            print('residual [mixed ReaxFF units, based on electronvolt]')
            print(residual)

            # Print energy contributions
            print('Energy 0 eneg          [k cal mol^-1]: {:10.5f}'.format(terms[0]/kcalmol))
            print('Energy 1 hard          [k cal mol^-1]: {:10.5f}'.format(terms[1]/kcalmol))
            print('Energy 2 coul          [k cal mol^-1]: {:10.5f}'.format(terms[2]/kcalmol))
            print('Energy 3 coup          [k cal mol^-1]: {:10.5f}'.format(terms[3]/kcalmol))
            print('Energy 4 soft          [k cal mol^-1]: {:10.5f}'.format(terms[4]/kcalmol))
            # Print them like ReaxFF
            print('Following NEW ReaxFF conventions...')
            print('Energy 0+1+3+4 Charges [k cal mol^-1]: {:10.5f}'.format(
                (terms[0]+terms[1]+terms[3]+terms[4])/kcalmol))
            print('Energy 2 Coulomb       [k cal mol^-1]: {:10.5f}'.format(
                (terms[2])/kcalmol))
            print('Following OLD ReaxFF conventions...')
            print('Energy 0+1+3 Charges   [k cal mol^-1]: {:10.5f}'.format(
                (terms[0]+terms[1]+terms[3])/kcalmol))
            print('Energy 2+4 Coulomb     [k cal mol^-1]: {:10.5f}'.format(
                (terms[2]+terms[4])/kcalmol))

        # The total energy
        energy = sum(terms)

        return energy, charges


def test_eem_read_parameters():
    model = EEMModel()
    model.load_parameters('ffield_eem')
    np.testing.assert_equal(model.rcut, 10.0*angstrom)
    np.testing.assert_equal(model.gammas['c'], 1.0000/angstrom)
    np.testing.assert_equal(model.chis['c'], 4.8000*electronvolt)
    np.testing.assert_equal(model.etas['c'], 6.7000*electronvolt)

    np.testing.assert_equal(model.gammas['h'], 0.8203/angstrom)
    np.testing.assert_equal(model.chis['h'], 3.7248*electronvolt)
    np.testing.assert_equal(model.etas['h'], 9.6093*electronvolt)

    np.testing.assert_equal(model.gammas['o'], 1.0898/angstrom)
    np.testing.assert_equal(model.chis['o'], 8.5000*electronvolt)
    np.testing.assert_equal(model.etas['o'], 8.3122*electronvolt)


def test_acks2_read_parameters():
    model = ACKS2Model()
    model.load_parameters('ffield_acks2')
    np.testing.assert_equal(model.rcut, 10.0*angstrom)
    np.testing.assert_equal(model.bsoft_amp, 548.6451/electronvolt)

    np.testing.assert_equal(model.gammas['c'], 0.3500/angstrom)
    np.testing.assert_equal(model.chis['c'], 5.3422*electronvolt)
    np.testing.assert_equal(model.etas['c'], 4.5000*electronvolt)
    np.testing.assert_equal(model.bsoft_radii['c'], 3.1838*angstrom)

    np.testing.assert_equal(model.gammas['h'], 0.6683/angstrom)
    np.testing.assert_equal(model.chis['h'], 4.9673*electronvolt)
    np.testing.assert_equal(model.etas['h'], 6.2079*electronvolt)
    np.testing.assert_equal(model.bsoft_radii['h'], 3.4114*angstrom)

    np.testing.assert_equal(model.gammas['o'], 0.5500/angstrom)
    np.testing.assert_equal(model.chis['o'], 8.5000*electronvolt)
    np.testing.assert_equal(model.etas['o'], 7.9071*electronvolt)
    np.testing.assert_equal(model.bsoft_radii['o'], 5.4479*angstrom)


def test_cutoff():
    # Make sure that sufficient neigboring images are taken into account to find all
    # pairs below a given cutoff.
    class TestModel(EEMModel):
        def __init__(self):
            EEMModel.__init__(self)
            self.rcut = 5.0
            self.rmin = 0.0

        def _set_physics_atom(self, A, B, atsymbols, iatom, natom):
            pass

        def _set_physics_atom_pair(self, A, B, atsymbols, iatom0, iatom1, natom, distance):
           if distance < self.rcut:
                self.distances.append(distance)

    model = TestModel()
    for i in range(10):
        # Get a reasonable cell
        while True:
            cellvecs = np.random.normal(0, 10, (3, 3))
            if abs(np.linalg.det(cellvecs)) > 8*8*8:
                break
        # Add some atoms, well beyond wrapping convention, to make it "more difficult".
        natom = 12
        atpositions = np.dot(np.random.uniform(-1, 2, (natom, 3)), cellvecs)
        atsymbols = ['H']*natom
        # Request the settings for periodic images
        recivecs, repeats = model._process_cellvecs(cellvecs)
        # Compute distances in the regular fashion
        model.distances = []
        model._set_physics(None, None, atsymbols, atpositions, cellvecs, recivecs, repeats)
        distances1 = model.distances
        # Compute the distances again, using a larger repeat vector.
        model.distances = []
        model._set_physics(None, None, atsymbols, atpositions, cellvecs, recivecs, repeats+1)
        distances2 = model.distances
        # No new distances should be found with larger repeat vector.
        assert len(distances1) == len(distances2)


def test_reaxff_units_kcalmol():
    # Coulomb energy in k cal mol^-1 of two protons at distance of one angstrom,
    # where the number 332.0638 is taken from the ReaxFF code. This is an approximation
    # of 1/(4*pi*epsilon_0) in those units. According to the NIST database of physical
    # constants, this should be 332.0637130025968.
    # (http://physics.nist.gov/cuu/Constants/index.html)
    e_coul_reaxff = 332.0638
    # The same thing in atomic units, with conversion constants for input and output.
    e_coul_script = (1*1/(1*angstrom))/kcalmol
    # Compare
    np.testing.assert_equal(e_coul_reaxff, e_coul_script)


def test_reaxff_units_electronvolt():
    # Coulomb energy in electronvolt of two protons at distance of one angstrom,
    # where the number 14.40 is taken from the ReaxFF code. This is an approximation
    # of 1/(4*pi*epsilon_0) in those units. According to the NIST database of physical
    # constants, this should be 14.399645352261372.
    # (http://physics.nist.gov/cuu/Constants/index.html)
    e_coul_reaxff = 14.40
    # The same thing in atomic units, with conversion constants for input and output.
    e_coul_script = (1*1/(1*angstrom))/electronvolt
    # Compare
    np.testing.assert_equal(e_coul_reaxff, e_coul_script)
