namespace FsLinAlg.Test

open System
open Expecto
open FsCheck
open FsLinAlg.DataStructure

[<AutoOpen>]
module ExpectoFsCheck =

    type SquareMatrix = SquareMatrix of Matrix
    type RectangularMatrix = RectangularMatrix of Matrix

    let minDim = 1
    let maxDim = 10
    
    let smallMatrixDim = gen {
        return! Gen.choose (minDim, maxDim)
    }
    let smallRectangularMatrixDim = gen {
        return! smallMatrixDim |> Gen.two
    }

    let nonNegInt = Arb.generate<NonNegativeInt>
    let nonNanFloat = Arb.generate<NormalFloat>

    type MatrixGens =
        static member Matrix() =
            { new Arbitrary<Matrix>() with
                override x.Generator =
                    let dim = smallRectangularMatrixDim |> Gen.sample 0 1 |> List.head
                    let vals = nonNanFloat |> Gen.array2DOfDim dim |> Gen.sample 0 1 |> List.head
                    let fvals = vals |> Array2D.map (fun v -> v.Get)
                    Matrix(fvals) |> Gen.constant }
        static member SquareMatrix() =
            { new Arbitrary<SquareMatrix>() with
                override x.Generator =
                    let dim = smallMatrixDim |> Gen.sample 0 1 |> List.head
                    let vals = nonNanFloat |> Gen.array2DOfDim (dim, dim) |> Gen.sample 0 1 |> List.head
                    let fvals = vals |> Array2D.map (fun v -> v.Get)
                    SquareMatrix(Matrix(fvals)) |> Gen.constant }
        static member RectangularMatrix() =
            { new Arbitrary<RectangularMatrix>() with
                override x.Generator =
                    let dim = smallRectangularMatrixDim |> Gen.sample 0 1 |> List.head
                    let vals = nonNanFloat |> Gen.array2DOfDim dim |> Gen.sample 0 1 |> List.head
                    let fvals = vals |> Array2D.map (fun v -> v.Get)
                    RectangularMatrix(Matrix(fvals)) |> Gen.constant }

    let private config = { 
        FsCheckConfig.defaultConfig with 
            maxTest = 1000
            arbitrary = [typeof<MatrixGens>] }
    
    let testProp name = testPropertyWithConfig config name