namespace FsLinAlg.Test

open Expecto

open FsLinAlg

module Eigenvalue =

    let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let eigenvalueTests =
        testList "Eigenvalue" [
            testList "Shifted QR" [
                test "np.linalg.eig reference test" {

                    (*
                    
                        array([6.54193618, 3.35146473, 4.5       , 4.60659909])
                    *)

                    let A =
                        [[|5.0; 0.5; 0.5; 0.5|] |> Vector
                         [|0.5; 5.0; 0.5; 0.5|] |> Vector
                         [|0.5; 0.5; 5.0; 1.0|] |> Vector
                         [|0.5; 0.5; 1.0; 4.0|] |> Vector]
                        |> Matrix.FromRowVectors

                    let D, _ = A.Hessenberg // To tridiagonal form
                    let eigs = ShiftedQR D None // To diagonal form

                    let npeigs = [6.54193618; 3.35146473; 4.5; 4.60659909]
                    List.zip eigs npeigs
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Eigenvalues are not equal")
                }

                test "np.linalg.eig reference test 2" {
                    let A =
                        [[|3.0; 1.0; 1.0|] |> Vector
                         [|1.0; 4.0; 0.0|] |> Vector
                         [|1.0; 0.0; 4.0|] |> Vector]
                        |> Matrix.FromRowVectors
                    
                    let D, _ = A.Hessenberg // To tridiagonal form
                    let eigs = ShiftedQR D None // To diagonal form

                    let npeigs = [2.; 5.; 4.]
                    List.zip eigs npeigs
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Eigenvalues are not equal")                    
                }

                test "np.linalg.eig reference test 3" {
                    let A = 
                        [[|4.0; 1.0; 0.0|] |> Vector
                         [|1.0; 3.0; 1.0|] |> Vector
                         [|0.0; 1.0; 4.0|] |> Vector]
                        |> Matrix.FromRowVectors
                       
                    let D, _ = A.Hessenberg // To tridiagonal form
                    let eigs = ShiftedQR D None // To diagonal form

                    let npeigs = [2.; 5.; 4.]
                    List.zip eigs npeigs
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Eigenvalues are not equal")
                }

                testProp "Trace of A = Sum of eigenvalues" <| fun (As: SymmetricPositiveDefiniteMatrix) ->
                    let (SymmetricPositiveDefiniteMatrix A) = As
                    let D, _ = A.Hessenberg // To tridiagonal form
                    let eigs = ShiftedQR D None // To diagonal form
                    let sumeig = eigs |> List.sum
                    let trace = A.Trace

                    Expect.floatClose Accuracy.medium sumeig trace "Trace did not equal the sum of eigenvalues"
            ]

            //test
        ]