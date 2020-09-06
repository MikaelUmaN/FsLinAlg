namespace FsLinAlg.Test

open Expecto

open FsLinAlg
open FsLinAlg.Equation

module Equation =

    let config =
        { FsCheckConfig.defaultConfig with
              maxTest = 10000 }

    [<Tests>]
    let equationTests =
        testList
            "Equation"
            [ testList "A=PLU" []

              testList
                  "Manual QRb"
                  [ test "scip.linalg.qr reference test" {

                        (*
                            From Python. There is a sign difference in R, but it is reflected in Q.

                            b = A @ x
                            Q, R = alg.qr(A, overwrite_a=False, lwork=None, mode='full', pivoting=False, check_finite=True)

                            R
                            array([[-6.40312424, -5.43484691,  1.48365074],
                                   [ 0.        ,  9.54056807,  1.09672868],
                                   [ 0.        ,  0.        , -1.35866357]])


                            Q
                            array([[-0.15617376,  0.95919021,  0.23571994],
                                   [-0.31234752,  0.1784421 , -0.93305811],
                                   [-0.93704257, -0.21934574,  0.27173271]])

                            b
                            array([16.74182727, 11.73402543, -2.71732713])
                        *)
                        let A =
                            [ [| 1.; 10.; 0.5 |] |> Vector
                              [| 2.; 3.4; 1. |] |> Vector
                              [| 6.; 3.; -2. |] |> Vector ]
                            |> Matrix.FromRowVectors
                        let x = [| -3.; 1.; 2. |] |> Vector
                        let b = (A * x).AsVector

                        let Qtb, R = A.QRb b
                        let ``x*`` = backSubstitute R Qtb
                        let ``b*`` = A * x |> Matrix.toVector

                        Expect.equal ``x*`` x "x vectors are not equal"
                        Expect.equal ``b*`` b "b vectors are not equal"
                    } ]

              testList
                  "QRb"
                  [ testProp "Rx=Qtb solves Ax=b"
                    <| fun (Abs: SquareMatrixSystem) ->
                        let (SquareMatrixSystem (A, b)) = Abs
                        let Qtb, R = A.QRb b
                        let x = backSubstitute R Qtb
                        let ``b*`` = A * x |> Matrix.toVector
                        Expect.equal ``b*`` b "b vectors are not equal" ] ]
