.. -*- mode: rst -*-

.. image:: https://travis-ci.com/ThomasHelfer/BosonStar.svg?token=UsxKKr9jzGhcSzDspCJA&branch=master
    :target: https://travis-ci.com/ThomasHelfer/BosonStar



BosonStar - A general solver for bosonic stars 
===================================================================================

This solver gives solutions for complex scalar boson stars with a massive
potential in arbitrary dimensions and a cosmological constant as well as for a
selfinteracting proca star in four dimensions. The hope of this repo is that
code I used for some project can maybe be reused and doesn't get lost.

Please cite these papers when using this library for a publication

+----------------+-------------+
| Matter type    | arxiv id    |
+================+=============+
| Scalar Complex | 0309131     |
+----------------+-------------+
| Vector Complex | 1609.01735  |
+----------------+-------------+
|                | 1508.05395  |
+----------------+-------------+
|                | 1805.09867  |
+----------------+-------------+

Installation 
============

To install this code, run ``python setup.py install --user`` in the root directory.


Examples
========

There are two scripts in the **example** folder, *scalar_star_solver.py* and *complex_vector_star_selfinteracting_solver.py*.
