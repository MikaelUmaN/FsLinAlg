# FsLinAlg
Linear algebra for FSharp.

**(Not in a useful state as of this commit)**

## Functionality
- Matrix, vectors and basic operators
- PA=LU
- A=QR (Q formed incorrectly...)

## Caveats
Nothing here is meant to be quick. Performance is best in dedicated libraries such as BLAS, LAPACK, CudaToolkit, ...

The plan is to add integrations to those backends and let F# be the nice, front interface.