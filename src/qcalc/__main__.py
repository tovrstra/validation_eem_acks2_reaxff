"""Compute EEM or ACKS2 charges as in ReaxFF."""

import argparse
import numpy as np

from . import EEMModel, ACKS2Model, electronvolt, kcalmol, load_structure, load_constraints


def main():
    """Main program."""
    print("The number 23.02 [eV (k cal mol^-1)^-1]: {:.10f} ;)".format(electronvolt / kcalmol))
    np.set_printoptions(precision=4, suppress=True, linewidth=5000, threshold=1000000)

    args = parse_args()
    atsymbols, atpositions, cellvecs = load_structure(args.struct)
    constraints = load_constraints(args.constrain, args.qtot, len(atsymbols))
    model = {
        "eem": EEMModel,
        "acks2": ACKS2Model,
    }[args.model]()
    model.load_parameters(args.ffield, args.verbose)

    energy, atcharges = model.compute_charges(
        atsymbols, atpositions, cellvecs, constraints, args.reduce, args.verbose
    )
    print("Energy {:15s} [k cal mol^-1]: {:10.5f}".format(args.model.upper(), energy / kcalmol))
    print("Charges [e]:")
    for s, q in zip(atsymbols, atcharges):
        print("{:>3s} {:10.5f}".format(s.upper(), q))
    print("Constraints:")
    for charge, indexes in constraints:
        print(
            "{:10.5f} {:10.5f} {}".format(
                charge, atcharges[indexes].sum(), ",".join([str(i) for i in indexes])
            )
        )


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser("Compute EEM or ACKS2 charges as in ReaxFF.")
    parser.add_argument(
        "model", choices=["eem", "acks2"], help="The model with which to compute the charges."
    )
    parser.add_argument("ffield", help="A ReaxFF parameter file.")
    parser.add_argument("struct", help="An XYZ structure files.")
    parser.add_argument("-q", "--qtot", type=float, default=0.0, help="The total charge.")
    parser.add_argument(
        "-c",
        "--constrain",
        help="Constrain charges on groups of atoms to have a certain value. This option "
        "requires one argument: a file containing all the constraints. Non-empty "
        "lines in this file contain a float (the target charge) followed by a list "
        "of atom indexes defining the group to constrain. Atom indexes start "
        "counting from one, like in FORTRAN.",
    )
    parser.add_argument(
        "-r",
        "--reduce",
        default=False,
        action="store_true",
        help="Try to eliminate the charge-constraints from the ACKS2 equations. This "
        "will only work for constraints that do not overlap.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="Print out intermediate results for debugging.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
