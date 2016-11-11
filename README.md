Description
-----------

This is a completely self-contained reference implementation of the EEM and ACKS2 models
as they are implemented in the ReaxFF model. Many examples are provided to generate
reference data that can be used to validate the charges obtained with full-blown ReaxFF
implementations.


Requirements
------------

* Python >= 2.7
* Numpy >= 1.0


How to run an example
---------------------

Just type the following on the command line:

    ./compute_charges.py {eem,acks2} reaxff_parameter_file {XYZ,PDB}

Energy (due to charges) and charges are just printed on screen.


References
----------

Molecular data sets taken from

* S66: http://begdb.com
* Crystals from COD

Description of EEM and ACKS2 models in ReaxFF was recently summarized in the following
paper:

http://dx.doi.org/10.1021/acs.jctc.6b00432

"eReaxFF: A Pseudoclassical Treatment of Explicit Electrons within Reactive Force Field Simulations"
Md Mahbubul Islam, Grigory Kolesov, Toon Verstraelen, Efthimios Kaxiras, and Adri C. T. van Duin
J. Chem. Theory Comput., 2016, 12 (8), pp 3463–3472
