namespace FsLinAlg.Test

open System
open Expecto
open FsCheck

open MathNet.Numerics.Distributions

open FsLinAlg

[<AutoOpen>]
module ExpectoFsCheck =

    let cdfn x = Normal.CDF(0., 1., x)

    type SquareMatrix = SquareMatrix of Matrix
    type SquareMatrixSystem = SquareMatrixSystem of Matrix * Vector
    type RectangularMatrix = RectangularMatrix of Matrix

    /// A tall and thin matrix has the property that m > n.
    type TallThinMatrix = TallThinMatrix of Matrix
    type TallThinMatrixSystem = TallThinMatrixSystem of Matrix * Vector

    type SymmetricPositiveDefiniteMatrix = SymmetricPositiveDefiniteMatrix of Matrix
    type SymmetricPositiveDefiniteMatrixSystem = SymmetricPositiveDefiniteMatrixSystem of Matrix * Vector

    let minDim = 2
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

    let filterSmallFloatsExclZero f = abs f > machEps || f = 0.
    let filterSmallFloats f = abs f > machEps

    let nonNegInt = Arb.generate<NonNegativeInt>

    type NonNanNonZeroFloat = NonNanNonZeroFloat of float with
        member x.Get = match x with NonNanNonZeroFloat f -> f
        override x.ToString() = x.Get.ToString()
        static member op_Explicit(NormalFloat f) = f

    let nonNanFloat =
        Arb.generate<NormalFloat>
        |> Gen.filter (fun (NormalFloat f) -> filterSmallFloatsExclZero f)

    let nonNanNonZeroFloat = 
        Arb.generate<NormalFloat>
        |> Gen.filter (fun (NormalFloat f) -> filterSmallFloats f)
        |> Gen.map (fun (NormalFloat f) -> NonNanNonZeroFloat(f))

    type FloatGens =
        static member NonNanNonZeroFloat() =
            { new Arbitrary<NonNanNonZeroFloat>() with
                override x.Generator = nonNanNonZeroFloat
            }

    type MatrixGens =
        static member Vector() =
            { new Arbitrary<Vector>() with
                override x.Generator =
                    gen {
                        let! dim = smallMatrixDim
                        let! vals = nonNanFloat |> Gen.arrayOfLength dim
                        let fvals = vals |> Array.map (fun v -> v.Get)
                        return Vector(fvals)
                    } }
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
        static member SquareMatrixSystem() =
            { new Arbitrary<SquareMatrixSystem>() with
                override x.Generator =
                    gen {
                        let! As = Arb.generate<SquareMatrix>
                        let (SquareMatrix A) = As
                        let! vals = nonNanFloat |> Gen.arrayOfLength A.M
                        let fvals = vals |> Array.map (fun v -> v.Get)
                        return SquareMatrixSystem(A, Vector(fvals))
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
        static member TallThinMatrixSystem() =
            { new Arbitrary<TallThinMatrixSystem>() with
                override x.Generator =
                    gen {
                        let! As = Arb.generate<TallThinMatrix>
                        let (TallThinMatrix A) = As
                        let! vals = nonNanFloat |> Gen.arrayOfLength A.M
                        let fvals = vals |> Array.map (fun v -> v.Get)
                        return TallThinMatrixSystem(A, Vector(fvals))
                    } }
        static member SymmetricPositiveDefiniteMatrix() =
            { new Arbitrary<SymmetricPositiveDefiniteMatrix>() with
            override x.Generator =
                gen {
                    let! dim = smallMatrixDim
                    let! vals = 
                        nonNanFloat 
                        |> Gen.map (fun v -> cdfn v.Get) 
                        |> Gen.filter filterSmallFloatsExclZero
                        |> Gen.array2DOfDim (dim, dim)
                    
                    // Start from a random matrix.
                    let A = Matrix(vals)

                    // Make it symmetric and diagonally dominant
                    // -> Positive Symmetric Definite
                    let At = A.T
                    let AAt = 0.5*(A + At) + float(dim)*(Matrix.I dim)

                    return SymmetricPositiveDefiniteMatrix(AAt)
                } }
        static member SymmetricPositiveDefiniteMatrixSystem() =
            { new Arbitrary<SymmetricPositiveDefiniteMatrixSystem>() with
            override x.Generator =
                gen {
                    let! As = Arb.generate<SymmetricPositiveDefiniteMatrix>
                    let (SymmetricPositiveDefiniteMatrix A) = As
                    let! vals = nonNanFloat |> Gen.arrayOfLength A.M
                    let fvals = vals |> Array.map (fun v -> v.Get)

                    return SymmetricPositiveDefiniteMatrixSystem(A, Vector(fvals))
                } }

    let private config = { 
        FsCheckConfig.defaultConfig with 
            maxTest = 1000
            arbitrary = [typeof<MatrixGens>; typeof<FloatGens>]}
    
    let testProp name = testPropertyWithConfig config name