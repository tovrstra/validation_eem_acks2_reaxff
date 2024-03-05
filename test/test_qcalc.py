import numpy as np

from qcalc import EEMModel, ACKS2Model, angstrom, electronvolt, kcalmol


def test_eem_read_parameters():
    model = EEMModel()
    model.load_parameters("ffield_eem")
    np.testing.assert_equal(model.rcut, 10.0 * angstrom)
    np.testing.assert_equal(model.gammas["c"], 1.0000 / angstrom)
    np.testing.assert_equal(model.chis["c"], 4.8000 * electronvolt)
    np.testing.assert_equal(model.etas["c"], 6.7000 * electronvolt)

    np.testing.assert_equal(model.gammas["h"], 0.8203 / angstrom)
    np.testing.assert_equal(model.chis["h"], 3.7248 * electronvolt)
    np.testing.assert_equal(model.etas["h"], 9.6093 * electronvolt)

    np.testing.assert_equal(model.gammas["o"], 1.0898 / angstrom)
    np.testing.assert_equal(model.chis["o"], 8.5000 * electronvolt)
    np.testing.assert_equal(model.etas["o"], 8.3122 * electronvolt)


def test_acks2_read_parameters():
    model = ACKS2Model()
    model.load_parameters("ffield_acks2")
    np.testing.assert_equal(model.rcut, 10.0 * angstrom)
    np.testing.assert_equal(model.bsoft_amp, 548.6451 / electronvolt)

    np.testing.assert_equal(model.gammas["c"], 0.3500 / angstrom)
    np.testing.assert_equal(model.chis["c"], 5.3422 * electronvolt)
    np.testing.assert_equal(model.etas["c"], 4.5000 * electronvolt)
    np.testing.assert_equal(model.bsoft_radii["c"], 3.1838 * angstrom)

    np.testing.assert_equal(model.gammas["h"], 0.6683 / angstrom)
    np.testing.assert_equal(model.chis["h"], 4.9673 * electronvolt)
    np.testing.assert_equal(model.etas["h"], 6.2079 * electronvolt)
    np.testing.assert_equal(model.bsoft_radii["h"], 3.4114 * angstrom)

    np.testing.assert_equal(model.gammas["o"], 0.5500 / angstrom)
    np.testing.assert_equal(model.chis["o"], 8.5000 * electronvolt)
    np.testing.assert_equal(model.etas["o"], 7.9071 * electronvolt)
    np.testing.assert_equal(model.bsoft_radii["o"], 5.4479 * angstrom)


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
            if abs(np.linalg.det(cellvecs)) > 8 * 8 * 8:
                break
        # Add some atoms, well beyond wrapping convention, to make it "more difficult".
        natom = 12
        atpositions = np.dot(np.random.uniform(-1, 2, (natom, 3)), cellvecs)
        atsymbols = ["H"] * natom
        # Request the settings for periodic images
        recivecs, repeats = model._process_cellvecs(cellvecs)
        # Compute distances in the regular fashion
        model.distances = []
        model._set_physics(
            None, None, atsymbols, atpositions, cellvecs, recivecs, repeats, safe=True
        )
        distances1 = model.distances
        # Compute the distances again, using a larger repeat vector.
        model.distances = []
        model._set_physics(
            None, None, atsymbols, atpositions, cellvecs, recivecs, repeats + 1, safe=True
        )
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
    e_coul_script = (1 * 1 / (1 * angstrom)) / kcalmol
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
    e_coul_script = (1 * 1 / (1 * angstrom)) / electronvolt
    # Compare
    np.testing.assert_equal(e_coul_reaxff, e_coul_script)
