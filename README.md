# DMA
Diagonal Matrix solving

Simple banded matrix class and solver for the equation A x = b where A is a banded matrix of low bandwidth, specifically 3 or 5.

The main solver is designed for the case of a dynamics problem such as balls joined by springs. In this case there is a large matrix of low bandwidth, possibly with special conditions on the ends. The matrix may change size slowly, and a large number of solves are done, meaning that preserving intermediate values is helpful for performance. An onthefly implementation is also given, which requires about 2-3 times as many operations per solve, with concomitant performance cost.

# Citations

Thomas algorithm for Tridiagonal matrices: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

L.U decomposition for a PentaDiagonal matrix: https://doi.org/10.1016/j.amc.2008.03.004
(Theorem 2.1 and Remark 2.5) used to create the equivalent algorithm