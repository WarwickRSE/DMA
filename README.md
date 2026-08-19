[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.zenodo.22011711-blue.svg)](https://doi.org/10.5281/zenodo.22011711)
# DMA
Diagonal Matrix solving

Simple banded matrix class and solver for the equation A x = b where A is a banded matrix of low bandwidth, specifically 3 or 5.

It is well know that a Gaussian Elimination approach in these cases can be simplified and sped up a lot by taking advantage of the fact that most terms are zero, and this is the basis of the extremely well known Thomas algorithm for a tri-diagonal matrix. 

It is also well-known that this approach can be duplicated for other low bandwidth cases. However, in spite of this being well-known, it proved difficult to find written down. This repo. hosts a working example, and a copy of the relevant equations fully written out. 

The main solver here is designed for the case of a dynamics problem such as balls joined by springs. In this case there is a large matrix of low bandwidth, with special conditions on the ends. The matrix may change size slowly, and a large number of solves are done, meaning that preserving intermediate values is helpful for performance. An on-the-fly implementation is also given, which requires about 2-3 times as many operations per solve, with concomitant performance cost.

# Citations

Thomas algorithm for Tridiagonal matrices: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

L.U decomposition for a PentaDiagonal matrix: https://doi.org/10.1016/j.amc.2008.03.004
(Theorem 2.1 and Remark 2.5) used to create the algorithm in this code

# How to Use this Software

Since matrix solving is deeply embedded into numerical codes, it is unlikely anybody except me would want to use this solver as written. However, if you benefit from this example when writing your own equivalent solver, be kind and drop us a citation or acknowledgement. TIA. 
