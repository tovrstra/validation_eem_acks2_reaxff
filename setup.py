#!/usr/bin/env python3
"""Installation script."""


import numpy as np

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name='Compute charges',
    version='0.0.0',
    description='Compute EEM or ACKS2 charges as in ReaxFF.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    scripts=['compute_charges'],
    cmdclass = {
        'build_ext': build_ext,
    },
    py_modules=['compute_charges_lib'],
    ext_modules=[
        Extension("compute_charges_ext",
            sources=['compute_charges_ext.pyx', 'low.cpp'],
            depends=['low.h', 'low.pxd'],
            include_dirs=[np.get_include(), '.'],
            extra_compile_args=['-std=c++11'],
            language="c++"),
    ],
)
