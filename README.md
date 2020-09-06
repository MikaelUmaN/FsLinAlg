# FsLinAlg
Linear algebra for FSharp.

**(Not in a useful state as of this commit)**

## Functionality
- Matrix, vectors and basic operators

## Factorizations
- PA=LU
- A=QR
- Cholesky, A=R.T R (also expressed A=LL.T)

## Solving Equations
- Backsubstitution for upper-triangular systems.
- Ax=b solved by QR with implicit Q.T * b formation.

## Eigenvalues

## Caveats
Nothing here is meant to be quick. Performance is best in dedicated libraries such as BLAS, LAPACK, CudaToolkit, ...

The plan is to add integrations to those backends and let F# be the nice, front interface.

## References
- Lloyd N. Trefethen, David Bau, III; Numerical Linear Algebra


