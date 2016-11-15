Disclaimer
----------

This is work in progress! Do not trust the results yet. Things to be done:

- Figure out why constraints on ACKS2 do not give the same result as disabling charge
  transfer between fragments. (This could be the right behavior.)


Description
-----------

This is a completely self-contained reference implementation of the EEM and ACKS2 models
as they are implemented in the ReaxFF model. (In a broader context, EEM and ACKS2 may take
different forms.) Many structures are provided to generate reference data that can be used
to validate the charges obtained with full-blown ReaxFF implementations.


Requirements
------------

* Python >= 3.0
* Numpy >= 1.0


How to run an example
---------------------

Just type the following on the command line:

    ./compute_charges.py {eem,acks2} reaxff_parameter_file xyz_file

Energy (due to charges) and charges are just printed on screen. XYZ files for molecules
and periodic structures are included in the repository. Two dummy ReaxFF parameter files
are provided: ffield_eem and ffield_acks2. Unit tests can be executed as follows:

    nosetests-3.* ./compute_charges.py


References
----------

Molecular data sets taken from

* S66: http://begdb.com
* Crystals from COD: http://www.crystallography.net/cod/

The mathematical form of the EEM and ACKS2 models in ReaxFF was recently summarized in the
following paper:

http://dx.doi.org/10.1021/acs.jctc.6b00432

"eReaxFF: A Pseudoclassical Treatment of Explicit Electrons within Reactive Force Field Simulations"
Md Mahbubul Islam, Grigory Kolesov, Toon Verstraelen, Efthimios Kaxiras, and Adri C. T. van Duin
J. Chem. Theory Comput., 2016, 12 (8), pp 3463–3472
