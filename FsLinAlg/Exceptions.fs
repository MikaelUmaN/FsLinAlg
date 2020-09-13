namespace FsLinAlg

[<AutoOpen>]
module Exceptions =
    exception LinearDependenceException of string
    let linDep = LinearDependenceException("Columns/Rows are linearly dependent")

    exception InvalidDimensionException of string
    let invDim = InvalidDimensionException("Dimensions are invalid")
    let invDimMsg msg = InvalidDimensionException(msg)

    exception PositiveDefiniteException of string
    let notPosDef = PositiveDefiniteException("Matrix is not positive definite")

    exception MaxIterationsException of string
    let maxIter k = MaxIterationsException(sprintf "Maximum number of iterations exceeded: %i" k)