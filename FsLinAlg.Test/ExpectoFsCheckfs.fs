namespace FsLinAlg.Test

open System
open Expecto
open FsCheck

open FsLinAlg

[<AutoOpen>]
module ExpectoFsCheck =

    type SquareMatrix = SquareMatrix of Matrix
    type RectangularMatrix = RectangularMatrix of Matrix

    /// A tall and thin matrix has the property that m > n.
    type TallThinMatrix = TallThinMatrix of Matrix

    let minDim = 1
    let maxDim = 10

    let minTallThinN = 1
    let maxTallThinN = 5
    let tallThinMFactor = 2
    
    let matrixDim minD maxD = gen {
        return! Gen.choose (minD, maxD)
    }
    let smallMatrixDim = gen {
        return! matrixDim minDim maxDim
    }
    let rectangularMatrixDim minD maxD = gen {
        return! matrixDim minD maxD |> Gen.two
    }
    let smallRectangularMatrixDim = gen {
        return! rectangularMatrixDim minDim maxDim
    }
    let tallThinMatrixDim minN maxN mFactor = gen {
        let! n = matrixDim minN maxN
        let minM = n + 1
        let maxM = minM * mFactor
        let! m = Gen.choose (minM, maxM)
        return (m, n)
    }
    let smallTallThinMatrixDim = gen {
        return! tallThinMatrixDim minTallThinN maxTallThinN tallThinMFactor
    }

    let nonNegInt = Arb.generate<NonNegativeInt>
    let nonNanFloat = Arb.generate<NormalFloat>

    type MatrixGens =
        static member Matrix() =
            { new Arbitrary<Matrix>() with
                override x.Generator =
                    gen {
                        let! dim = smallRectangularMatrixDim
                        let! vals = nonNanFloat |> Gen.array2DOfDim dim
                        let fvals = vals |> Array2D.map (fun v -> v.Get)
                        return Matrix(fvals)
                    } }
        static member SquareMatrix() =
            { new Arbitrary<SquareMatrix>() with
                override x.Generator =
                    gen {
                        let! dim = smallMatrixDim
                        let! vals = nonNanFloat |> Gen.array2DOfDim (dim, dim)
                        let fvals = vals |> Array2D.map (fun v -> v.Get)
                        return SquareMatrix(Matrix(fvals))
                    } }
        static member RectangularMatrix() =
            { new Arbitrary<RectangularMatrix>() with
                override x.Generator =
                    gen {
                        let! dim = smallRectangularMatrixDim
                        let! vals = nonNanFloat |> Gen.array2DOfDim dim
                        let fvals = vals |> Array2D.map (fun v -> v.Get)
                        return RectangularMatrix(Matrix(fvals))
                    } }
        static member TallThinMatrix() =
            { new Arbitrary<TallThinMatrix>() with
                override x.Generator =
                    gen {
                        let! dim = smallTallThinMatrixDim
                        let! vals = nonNanFloat |> Gen.array2DOfDim dim
                        let fvals = vals |> Array2D.map (fun v -> v.Get)
                        return TallThinMatrix(Matrix(fvals))
                    } }

    let private config = { 
        FsCheckConfig.defaultConfig with 
            maxTest = 1000
            arbitrary = [typeof<MatrixGens>] }
    
    let testProp name = testPropertyWithConfig config name