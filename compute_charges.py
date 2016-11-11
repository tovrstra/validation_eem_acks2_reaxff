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


# For now, taken from NIST. Should be replaced by (less accurate) units used in ReaxFF.
# The calorie is just defined with a four-digit precision relative to the Joule.
angstrom = 1.0/0.52917721067
electronvolt = 1.6021766208e-19/4.359744650e-18
kcalmol = (4184/4.359744650e-18/6.022140857e23)


def main():
    """Main program."""
    args = parse_args()
    model = {
        'eem': EEMModel,
    }[args.model]()
    model.load_parameters(args.ffield)
    atsymbols, atpositions, cellvecs = load_structure(args.struct)

    energy, atcharges = model.compute_charges(atsymbols, atpositions, cellvecs)
    print('Energy [k cal mol^-1] = {:.5f}'.format(energy/kcalmol))
    print('Charges [e]:')
    for q in atcharges:
        print('{:10.5f}'.format(q))


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser('Compute EEM or ACKS2 charges as in ReaxFF.')
    parser.add_argument('model', choices=['eem', 'acks'],
                        help='The model with which to compute the charges.')
    parser.add_argument('ffield', help='A ReaxFF parameter file.')
    parser.add_argument('struct', help='An XYZ structure files.')
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


class EEMModel(object):
    """Compute EEM charges as in ReaxFF."""

    def __init__(self):
        """Initialize the EEMModel object."""
        self.rcut = None
        self.gammas = {}
        self.chis = {}
        self.etas = {}

    def load_parameters(self, ffield):
        """Load relevant parameters from a ReaxFF parameter file."""
        with open(ffield) as f:
            lines = f.readlines()
            self._extract_parameters(lines)

    def _extract_parameters(self, lines):
        """Extract parameters from a list of lines, loaded from a ReaxFF parameter file."""
        self.rcut = float(lines[14].split()[0])*angstrom
        print('Coulomb cutoff [Å]: {:10.5f}'.format(self.rcut/angstrom))
        nelement = int(lines[41].split()[0])
        for ielement in range(nelement):
            symbol = lines[45+ielement*4].split()[0].lower()
            self.gammas[symbol] = float(lines[45+ielement*4].split()[6])*(1/angstrom)
            self.chis[symbol] = float(lines[46+ielement*4].split()[5])*electronvolt
            self.etas[symbol] = float(lines[46+ielement*4].split()[6])*electronvolt

    def compute_charges(self, atsymbols, atpositions, cellvecs):
        """Compute atomic charges for a single molecules.

        Parameters
        ----------
        atsymbols : list of str
            Atomic symbols, should correspond to symbols in ffield.
        atpositions : np.ndarray, dtype=float, shape=(natom, 1)
            Atomic positions in atomic units.
        cellvecs : np.ndarray, dtype=float, shape=(3, 3)
            Rows of this matrix are cell vectors. This argument may also be None.
        """
        natom = len(atsymbols)
        if cellvecs is None:
            recivecs = None
            repeat_a = 0
            repeat_b = 0
            repeat_c = 0
        else:
            recivecs = np.linalg.inv(cellvecs)
            # Compute the number of images to be considered to include everything within
            # the cutoff sphere.
            spacing_a = 1.0/np.linalg.norm(recivecs[:,0])
            spacing_b = 1.0/np.linalg.norm(recivecs[:,1])
            spacing_c = 1.0/np.linalg.norm(recivecs[:,2])
            print('Crystal plane spacings [Å]: {:10.5f} {:10.5f} {:10.5f}'.format(
                spacing_a/angstrom, spacing_b/angstrom, spacing_c/angstrom))
            print(self.rcut*np.linalg.norm(recivecs[:,0]))
            repeat_a = int(np.floor(self.rcut/spacing_a + 0.5))
            repeat_b = int(np.floor(self.rcut/spacing_b + 0.5))
            repeat_c = int(np.floor(self.rcut/spacing_c + 0.5))
            print('Supercell for electrostatics: {} {} {}'.format(2*repeat_a+1, 2*repeat_b+1, 2*repeat_c+1))

        # Build up the electronegativity vector and hardness matrix, extended with
        # elements for the Lagrange multiplier, assuming the total charge is zero.
        chi_vector = np.zeros(natom + 1, float)
        eta_matrix = np.zeros((natom+1, natom+1), float)
        for iatom0 in range(natom):
            chi_vector[iatom0] = self.chis[atsymbols[iatom0]]
            eta_matrix[iatom0, iatom0] = self.etas[atsymbols[iatom0]]
            eta_matrix[natom, :natom] = 1
            eta_matrix[:natom, natom] = 1
            gamma0 = self.gammas[atsymbols[iatom0]]
            for iatom1 in range(natom):
                gamma1 = self.gammas[atsymbols[iatom1]]
                # Get the (naive) minimum image convention
                delta0 = atpositions[iatom0] - atpositions[iatom1]
                # Apply minimum image convention
                if recivecs is not None:
                    delta0_frac = np.dot(recivecs.T, delta0)
                    delta0_frac -= np.round(delta0_frac)
                    delta0 = np.dot(cellvecs.T, delta0_frac)
                # Loop over all relevant neighboring cells
                for image_a in range(-repeat_a, repeat_a+1):
                    for image_b in range(-repeat_b, repeat_b+1):
                        for image_c in range(-repeat_c, repeat_c+1):
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
                            assert distance > 1.0
                            # In atomic units, 1/(4*pi*epsilon_0) is numerically equal
                            # to one, which is nice!
                            coulomb = (distance**3 + (gamma0*gamma1)**(-3.0/2.0))**(-1.0/3.0)
                            coulomb *= self.taper(distance)
                            eta_matrix[iatom0, iatom1] += coulomb

        # Solve the charges
        charges = np.linalg.solve(eta_matrix, chi_vector)[:-1]

        # Compute the energy
        energy = np.dot(chi_vector[:-1], charges) + 0.5*np.dot(charges, np.dot(eta_matrix[:-1,:-1], charges))

        return energy, charges

    def taper(self, distance):
        """Taper correction as in ReaxFF."""
        TAP7 = 20.0/self.rcut**7
        TAP6 = -70.0/self.rcut**6
        TAP5 = 84.0/self.rcut**5
        TAP4 = -35.0/self.rcut**4
        return 1.0 + TAP4*distance**4 + TAP5*distance**5 + TAP6*distance**6 + TAP7*distance**7


if __name__ == '__main__':
    main()
