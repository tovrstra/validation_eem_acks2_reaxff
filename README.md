# EEM and ACKS2 ReaxFF validation

## Disclaimer

This is work in progress! Do not trust the results yet. Things to be done:

- Figure out why ACKS2 results are slightly different from ReaxFF:
    - Compare `xmortr` and `elcvec` arrays.
    - Also compare EEM results
- Figure out why constraints on ACKS2 do not give the same result as disabling charge
  transfer between fragments. (This could be the right behavior.)


## Description

This is a completely self-contained reference implementation of the EEM and ACKS2 models
as they are implemented in the ReaxFF model. (In a broader context, EEM and ACKS2 may take
different forms.) Many structures are provided to generate reference data that can be used
to validate the charges obtained with full-blown ReaxFF implementations.


## Requirements

* Python >= 3.0
* Numpy >= 1.0
* Pytest (for testing only)


## Installation

### Minimal setup

```bash
git clone git@github.com:tovrstra/validation_eem_acks2_reaxff.git
cd validation_eem_acks2_reaxff
pip install .
```

### Development setup

```bash
git clone git@github.com:tovrstra/validation_eem_acks2_reaxff.git
cd validation_eem_acks2_reaxff
python -m venv venv
echo "source venv/bin/activate" > .envrc
direnv allow
pip install -e .
pip install pytest
hash -r
pytest -v
```


## How to run an example

Just type the following on the command line:

```bash
qcalc {eem,acks2} reaxff_parameter_file xyz_file
```

Energy (due to charges) and charges are just printed on screen. XYZ files for molecules
and periodic structures are included in the repository. Two dummy ReaxFF parameter files
are provided: `ffield_eem` and `ffield_acks2`.


# References

Molecular data sets taken from

* S66: http://begdb.com
* Crystals from COD: http://www.crystallography.net/cod/

The mathematical form of the EEM and ACKS2 models in ReaxFF was recently summarized in the
following paper:

* 1. Islam, M. M., Kolesov, G., Verstraelen, T., Kaxiras, E. & van Duin, A. C. T. eReaxFF:
  A Pseudoclassical Treatment of Explicit Electrons within Reactive Force Field
  Simulations. J. Chem. Theory Comput. 12, 3463–3472 (2016).
  http://dx.doi.org/10.1021/acs.jctc.6b00432

The ACKS2 model (in its more general form) is described in the following two papers:

* Verstraelen, T., Ayers, P. W., Van Speybroeck, V. & Waroquier, M. ACKS2: atom-condensed
  Kohn-Sham DFT approximated to second order. J. Chem. Phys. 138, 074108 (2013).
  http://dx.doi.org/10.1063/1.4791569

* Verstraelen, T., Vandenbrande, S. & Ayers, P. Direct computation of parameters for
  accurate polarizable force fields. J. Chem. Phys. 141, 194114 (2014).
  http://dx.doi.org/10.1063/1.4901513
