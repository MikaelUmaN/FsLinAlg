namespace FsLinAlg.Test

open Expecto

open FsLinAlg

module Eigenvalue =

    let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let eigenvalueTests =
        testList "Eigenvalue" [
            testList "Svd" [

                test "Golub-Kahan 3rd ed. p. 456 reference test" {
                    let A =
                        [[|1.0; 1.0; 0.0; 0.0|] |> Vector
                         [|0.0; 2.0; 1.0; 0.0|] |> Vector
                         [|0.0; 0.0; 3.0; 1.0|] |> Vector
                         [|0.0; 0.0; 0.0; 4.0|] |> Vector]
                        |> Matrix.FromRowVectors
                    
                    let D = SvdSteps A // To tridiagonal form
                    let s = 
                        D.D.Data 
                        |> Array.toList
                        |> List.sortDescending

                    let nps = [4.26000668; 3.10734857; 2.11178459; 0.85854166]
                    List.zip s nps
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Singular values are not equal")                    
                }

                test "np.linalg.svd reference test" {
                    let A =
                        [[|1.0; 0.0; 0.0|] |> Vector
                         [|0.0; 1.0; -0.81|] |> Vector
                         [|0.0; 0.0; 1.0|] |> Vector]
                        |> Matrix.FromRowVectors
                    
                    let D = SvdSteps A // To tridiagonal form
                    let s = 
                        D.D.Data 
                        |> Array.toList
                        |> List.sortDescending

                    let nps = [1.4838999; 1.; 0.6738999]
                    List.zip s nps
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Singular values are not equal")                    
                }

                testProp "Svd step yields orthogonal U and V" <| fun (Bs: BidiagonalMatrix) ->
                    let (BidiagonalMatrix B) = Bs
                    let _, U, V = SvdStep B

                    Expect.isTrue U.IsOrthogonal "Matrix U is not orthogonal"
                    Expect.isTrue V.IsOrthogonal "Matrix V is not orthogonal"

                testProp "Svd steps yield diagonal matrix" <| fun (Bs: BidiagonalMatrix) ->
                    let (BidiagonalMatrix B) = Bs
                    let D = SvdSteps B

                    Expect.isTrue D.IsDiagonal "Matrix D is not diagonal"
            ]

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