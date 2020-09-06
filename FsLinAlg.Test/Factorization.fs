namespace FsLinAlg.Test

open Expecto

open FsLinAlg

module Factorization =

    let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let factorizationTests =
        testList "Factorization" [
            testList "A=PLU" [
                testProp "Matrix Multiplication Equivalency" <| fun (As: SquareMatrix) ->
                    let (SquareMatrix A) = As
                    try
                        let P, L, U = A.LU
                        let ``PA*`` = L * U
                        let PA = P * A
                        Expect.equal ``PA*`` PA "PA did not equal PLU"
                    with
                    | :? LinearDependenceException -> () // Expected not to work on linearly dependent matrices.
            ]

            testList "A=QR" [
                testProp "Matrix Multiplication Equivalency" <| fun (As: TallThinMatrix) ->
                    let (TallThinMatrix A) = As
                    let Q, R = A.QR
                    let ``A*`` = Q * R
                    Expect.equal ``A*`` A "Matrices are not equal"
            ]

            testList "QRb" [
                testProp "Qtb, R equivalency with A=QR" <| fun (Abs: TallThinMatrixSystem) ->
                    let (TallThinMatrixSystem (A, b)) = Abs
                    let _, ``R*`` = A.QRb b
                    let _, R = A.QR
                    Expect.equal ``R*`` R "Matrices are not equal"
            ]

            testList "Cholesky" [
                testProp "A=R.T R multiplication equivalency" <| fun (As: SymmetricPositiveDefiniteMatrix) ->
                    let (SymmetricPositiveDefiniteMatrix A) = As
                    let R = A.Cholesky
                    let ``A*`` = R.T * R
                    Expect.equal ``A*`` A "Matrices are not equal"
            ]

            testList "Hessenberg" [
                testProp "Upper Hessenberg result" <| fun (As: SquareMatrix) ->
                    let (SquareMatrix A) = As
                    let H, _ = A.Hessenberg
                    Expect.isTrue (H.IsHessenberg()) "Matrix is not of hessenberg form"
            ]
        ]