.. -*- mode: rst -*-

.. image:: https://travis-ci.com/ThomasHelfer/BosonStar.svg?token=UsxKKr9jzGhcSzDspCJA&branch=master
    :target: https://travis-ci.com/ThomasHelfer/BosonStar



BosonStar - A general solver for bosonic stars 
===================================================================================

This solver gives solutions for complex scalar boson stars with a massive
potential in arbitrary dimensions and a cosmological constant as well as for a
selfinteracting proca star in four dimensions. The aim of this repo is that
code I used for some project can maybe be reused and doesn't get lost.

I used these papers to write this repo, please cite them since they did all of
the hard work.

+------------------------------------------+-------------+
| Matter type                              | arxiv id    |
+==========================================+=============+
| Scalar Complex                           | 0309131     |
+------------------------------------------+-------------+
| Proca Complex                            | 1609.01735  |
+------------------------------------------+-------------+
|                                          | 1508.05395  |
+------------------------------------------+-------------+
| complex proca star with selfinteractions | 1805.09867  |
+------------------------------------------+-------------+


Installation 
============

To install this code, run ``python setup.py install --user`` in the root directory.


Testing
============

For all modules available on this repo there I check consistency with constraint
equations automatically with every change of the code. This doesn't guarantee
that the code is bug-free, but is not completely unreasonable.

Examples
========

There are two scripts in the **example** folder, *scalar_star_solver.py* and
*complex_vector_star_selfinteracting_solver.py*. They contain examples on how to
use the code, especially some good guesses for the shooting that produce
reasonable results.
