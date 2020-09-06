namespace FsLinAlg.Test

open Expecto
open FsCheck

open FsLinAlg

module MatrixType =

    //let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let properties =
        testList "Matrix type tests" [
            testList "Square Matrix" [
                testProp "Square" <| fun (As: SquareMatrix) ->
                    let (SquareMatrix A) = As
                    Expect.equal A.M A.N "Dimensions of a square matrix should be equal"
            ]

            testList "Positive Symmetric Definitiveness" [
                testProp "Square" <| fun (As: SymmetricPositiveDefiniteMatrix) ->
                    let (SymmetricPositiveDefiniteMatrix A) = As
                    Expect.isTrue A.IsSquare "SPD matrix should be square"

                testProp "Symmetric" <| fun (As: SymmetricPositiveDefiniteMatrix) ->
                    let (SymmetricPositiveDefiniteMatrix A) = As
                    Expect.equal A A.T "SPD matrix is not equal to its transpose"

                testProp "Energy Property" <| fun (As: SymmetricPositiveDefiniteMatrixSystem) ->
                    let (SymmetricPositiveDefiniteMatrixSystem (A, x)) = As
                    let energy = x.T*A*x |> Matrix.toScalar
                    Expect.isTrue (energy > 0.) "Energy property of positive definitiveness not fullfilled"
            ]
        ]