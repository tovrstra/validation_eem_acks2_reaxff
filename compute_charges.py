#!/usr/bin/env python3
"""Compute EEM or ACKS2 charges as in ReaxFF.

This is a reasonably self-contained script. It only requires a working Python and Numpy
installation.

Unlike ReaxFF, this script internally works in atomic units, because that is sooo much
easier. All parameters and geometries are converted to atomic units when they are read.
"""

import argparse
import shlex

import numpy as np
np.set_printoptions(precision=4, suppress=True, linewidth=5000)


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
kcalmol = 1.0/(angstrom*332.0638)
electronvolt = 1.0/(angstrom*14.40)


def main():
    """Main program."""
    args = parse_args()
    atsymbols, atpositions, cellvecs = load_structure(args.struct)
    constraints = load_constraints(args.constrain, args.qtot, len(atsymbols))
    model = {
        'eem': EEMModel,
        'acks2': ACKS2Model,
    }[args.model]()
    model.load_parameters(args.ffield)

    energy, atcharges = model.compute_charges(atsymbols, atpositions, cellvecs,
                                              constraints, args.reduce, args.verbose)
    print('Energy [k cal mol^-1] = {:.5f}'.format(energy/kcalmol))
    print('Charges [e]:')
    for q in atcharges:
        print('{:10.5f}'.format(q))
    print('Constraints:')
    for charge, indexes in constraints:
        print('{:10.5f} {:10.5f} {}'.format(charge, atcharges[indexes].sum(), ','.join([str(i) for i in indexes])))


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser('Compute EEM or ACKS2 charges as in ReaxFF.')
    parser.add_argument(
        'model', choices=['eem', 'acks2'],
        help='The model with which to compute the charges.')
    parser.add_argument('ffield', help='A ReaxFF parameter file.')
    parser.add_argument('struct', help='An XYZ structure files.')
    parser.add_argument('-q', '--qtot', type=float, default=0.0, help='The total charge.')
    parser.add_argument(
        '-c', '--constrain',
        help='Constrain charges on groups of atoms to have a certain value. This option '
             'requires one argument: a file containing all the constraints. Non-empty '
             'lines in this file contain a float (the target charge) followed by a list '
             'of atom indexes defining the group to constrain. Atom indexes start '
             'counting from one, like in FORTRAN.')
    parser.add_argument(
        '-r', '--reduce', default=False, action='store_true',
        help='Try to eliminate the charge-constraints from the ACKS2 equations. This '
             'will only work for constraints that do not overlap.')
    parser.add_argument(
        '-v', '--verbose', default=False, action='store_true',
        help='Print out intermediate results for debugging.')
    return parser.parse_args()


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

    def __init__(self):
        """Initialize the EEMModel object."""
        self.rcut = None
        self.gammas = {}
        self.chis = {}
        self.etas = {}
        self.rmin = 0.5*angstrom

    def load_parameters(self, ffield):
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
            self._extract_parameters()

    def _extract_parameters(self):
        """Extract parameters from a list of lines, loaded from a ReaxFF parameter file."""
        self.rcut = self.general_pars[12][0]*angstrom
        assert self.rcut > 0
        print('Coulomb cutoff [Å]: {:10.5f}'.format(self.rcut/angstrom))
        for symbol, values in self.atom_pars.items():
            self.gammas[symbol] = values[5]*(1/angstrom)
            self.chis[symbol] = values[13]*electronvolt
            self.etas[symbol] = values[14]*electronvolt
            if 2*self.etas[symbol] < self.gammas[symbol]:
                print('Polarization catastrophe safety check failed for {}'.format(symbol))
                print('    2*eta[{symbol}]*4*pi*epsilon_0 = {:.4f} Å^-1 < gamma[{symbol}] = {:.4f} Å^-1'.format(
                    2*self.etas[symbol]*angstrom, self.gammas[symbol]*angstrom, symbol=symbol))
                print('    eta[{symbol}] = {:.4f} eV e^-2 < gamma[{symbol}]/(8*pi*epsilon_0) = {:.4f} ev e^-2'.format(
                    self.etas[symbol]/electronvolt, self.gammas[symbol]/electronvolt/2, symbol=symbol))

    def _process_cellvecs(self, cellvecs):
        if cellvecs is None:
            recivecs = None
            repeats = np.zeros(3, int)
        else:
            recivecs = np.linalg.inv(cellvecs)
            # Compute the number of images to be considered to include everything within
            # the cutoff sphere.
            spacings = (recivecs**2).sum(axis=0)**(-0.5)
            print('Crystal plane spacings [Å]: {:10.5f} {:10.5f} {:10.5f}'.format(*(spacings/angstrom)))
            repeats = np.floor(self.rcut/spacings + 0.5).astype(int)
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

    def _set_physics(self, A, B, atsymbols, atpositions, cellvecs, recivecs, repeats):
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
        if A is not None:
            evals = np.linalg.eigvalsh(A[:natom, :natom])
            if evals.min() <= 0:
                print('Hardness matrix has non-positive eigenvalues, which is not good!')
                print('(This is the polarization catastrophe.)')
                for e in evals:
                    print('{:10.5f}'.format(e))

    def _set_physics_atom(self, A, B, atsymbols, iatom, natom):
        B[iatom] = -self.chis[atsymbols[iatom]]
        A[iatom, iatom] = 2*self.etas[atsymbols[iatom]]

    def _set_physics_atom_pair(self, A, B, atsymbols, iatom0, iatom1, natom, distance):
        gamma0 = self.gammas[atsymbols[iatom0]]
        gamma1 = self.gammas[atsymbols[iatom1]]
        # In atomic units, 1/(4*pi*epsilon_0) is numerically equal
        # to one, which is nice!
        coulomb = (distance**3 + (gamma0*gamma1)**(-3.0/2.0))**(-1.0/3.0)
        coulomb *= self.taper(distance)
        A[iatom0, iatom1] += coulomb

    def compute_charges(self, atsymbols, atpositions, cellvecs, constraints, reduce_constraints, verbose):
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
        """
        natom = len(atsymbols)
        ncon = len(constraints)

        recivecs, repeats = self._process_cellvecs(cellvecs)

        # Build up the equations to be solveds
        B = np.zeros(natom + ncon, float)
        A = np.zeros((natom + ncon, natom + ncon), float)
        for icon, (charge, indexes) in enumerate(constraints):
            self._set_constraint(A, B, natom + icon, charge, indexes)
        self._set_physics(A, B, atsymbols, atpositions, cellvecs, recivecs, repeats)

        # Print out intermediate result.
        if verbose:
            print('A')
            print(A)
            print('B')
            print(B)

        # Solve the charges
        charges = np.linalg.solve(A, B)[:natom]

        # Compute the energy
        energy = np.dot(B[:natom], charges) + 0.5*np.dot(charges, np.dot(A[:natom,:natom], charges))

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

    def __init__(self):
        """Initialize the ACKS2Model object."""
        EEMModel.__init__(self)
        self.bsoft_amp = None
        self.bsoft_radii = {}

    def _extract_parameters(self):
        """Extract parameters from a list of lines, loaded from a ReaxFF parameter file."""
        EEMModel._extract_parameters(self)
        self.bsoft_amp = self.general_pars[34][0]/electronvolt
        if self.bsoft_amp <= 0:
            print('The bond softness amplitude is not positive, while it should.')
        for symbol, values in self.atom_pars.items():
            self.bsoft_radii[symbol] = values[22]*angstrom

    def _set_physics_atom(self, A, B, atsymbols, iatom, natom):
        EEMModel._set_physics_atom(self, A, B, atsymbols, iatom, natom)
        A[iatom, iatom+natom] = 1.0
        A[iatom+natom, iatom] = 1.0

    def _set_physics_atom_pair(self, A, B, atsymbols, iatom0, iatom1, natom, distance):
        EEMModel._set_physics_atom_pair(self, A, B, atsymbols, iatom0, iatom1, natom, distance)
        bsoft_rcut = (self.bsoft_radii[atsymbols[iatom0]] + self.bsoft_radii[atsymbols[iatom1]])/2
        assert bsoft_rcut < self.rcut
        if distance < bsoft_rcut:
            x = distance/bsoft_rcut
            bsoft = self.bsoft_amp*x**3*(1-x**6)
            A[iatom0 + natom, iatom0 + natom] -= bsoft
            A[iatom1 + natom, iatom1 + natom] -= bsoft
            A[iatom0 + natom, iatom1 + natom] += bsoft
            A[iatom1 + natom, iatom0 + natom] += bsoft

    def _reduce_constraints(self, A, B, natom, ncon):
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

    def compute_charges(self, atsymbols, atpositions, cellvecs, constraints, reduce_constraints, verbose):
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
        B[natom:2*natom] = np.linalg.lstsq(
            A[2*natom:2*natom+ncon, :natom],
            B[2*natom:2*natom+ncon], rcond=1e-10)[0]

        # Physics
        self._set_constraint(A, B, 2*natom + ncon, 0.0, np.arange(natom, 2*natom))
        self._set_physics(A, B, atsymbols, atpositions, cellvecs, recivecs, repeats)

        # Optionally simplify the equations
        if reduce_constraints:
            A, B = self._reduce_constraints(A, B, natom, ncon)

        # Print out intermediate result.
        if verbose:
            print('A')
            print(A)
            print('B')
            print(B)

        # Solve
        solution = np.linalg.solve(A, B)
        charges = solution[:natom]
        potentials = solution[natom:2*natom]

        # Compute the energy
        terms = [
            np.dot(B[:natom], charges),
            0.5*np.dot(charges, np.dot(A[:natom,:natom], charges)),
            np.dot(charges - B[natom:2*natom], potentials),
            0.5*np.dot(potentials, np.dot(A[natom:2*natom,natom:2*natom], potentials)),
        ]
        #print('energy 1 eneg [k cal mol^-1]: {:10.5f}'.format(terms[0]/kcalmol))
        #print('energy 2 hard [k cal mol^-1]: {:10.5f}'.format(terms[1]/kcalmol))
        #print('energy 3 coup [k cal mol^-1]: {:10.5f}'.format(terms[2]/kcalmol))
        #print('energy 4 soft [k cal mol^-1]: {:10.5f}'.format(terms[3]/kcalmol))
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


if __name__ == '__main__':
    main()
