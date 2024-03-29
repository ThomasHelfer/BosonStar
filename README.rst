.. -*- mode: rst -*-



BosonStar - A general solver for bosonic stars 
===================================================================================

This solver gives solutions for the following Boson Stars

+------------------------------------------+-------------+----------------------------+------------------------------------------------+
| Matter type                              | dimension   | cosmological constant      | potential                                      |
+==========================================+=============+============================+================================================+
| Scalar Complex                           | arbitrary   | yes                        | 1/2 mu^2 phi^2                                 |
+------------------------------------------+-------------+----------------------------+------------------------------------------------+
| complex proca star with selfinteractions | 4           | no                         | 1/2 mu^2 A^nu A_nu + 1/4 c4 (A^nu A_nu)^2      |
+------------------------------------------+-------------+----------------------------+------------------------------------------------+

The aim of this repo is that code I used for some project can maybe be reused
and doesn't get lost and to spread the boson star fun. I used these papers to
write this repo, please cite them since they did all of the hard work.

+------------------------------------------+-------------+
| Matter type                              | arXiv id    |
+==========================================+=============+
| Scalar Complex                           | 0309131     |
+------------------------------------------+-------------+
| complex proca star with selfinteractions | 1805.09867  |
+------------------------------------------+-------------+
| complex proca star                       | 1609.01735  |
+------------------------------------------+-------------+
|                                          | 1508.05395  |
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
