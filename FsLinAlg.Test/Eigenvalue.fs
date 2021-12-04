namespace FsLinAlg.Test

open Expecto

open FsLinAlg

module Eigenvalue =

    let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    // TODO: Some instabilities with these tests.
    let shiftedQrTests =
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

    [<Tests>]
    let eigenvalueTests =
        testList "Eigenvalue" [
            testList "Svd" [

                test "Golub, Van Loan 3rd ed. p. 456 reference test" {
                    let A =
                        [[|1.0; 1.0; 0.0; 0.0|] |> Vector
                         [|0.0; 2.0; 1.0; 0.0|] |> Vector
                         [|0.0; 0.0; 3.0; 1.0|] |> Vector
                         [|0.0; 0.0; 0.0; 4.0|] |> Vector]
                        |> Matrix.FromRowVectors
                    
                    let D, Ut, V = SvdSteps A // To tridiagonal form
                    let s = 
                        D.D.Data 
                        |> Array.toList
                        |> List.sortDescending

                    let nps = [4.26000668; 3.10734857; 2.11178459; 0.85854166]
                    List.zip s nps
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Singular values are not equal")                    
                }

                test "np.lingalg.svd reference test with diagonal zero element" {
                    let A =
                        [[|1.; 1.; 0.; 0.|] |> Vector
                         [|0.; 2.; 1.; 0.|] |> Vector
                         [|0.; 0.; 0.; 1.|] |> Vector
                         [|0.; 0.; -0.; -4.|] |> Vector]
                        |> Matrix.FromRowVectors

                    let D, Ut, V = SvdSteps A
                    let s = 
                        D.D.Data 
                        |> Array.toList
                        |> List.sortDescending

                    let nps = [4.12310563; 2.44948974; 1.; 0.]
                    List.zip s nps
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Singular values are not equal")  
                }

                test "np.linalg.svd eigenvalue reference test" {
                    let A =
                        [[|1.0; 0.0; 0.0|] |> Vector
                         [|0.0; 1.0; -0.81|] |> Vector
                         [|0.0; 0.0; 1.0|] |> Vector]
                        |> Matrix.FromRowVectors
                    
                    let D, Ut, V = SvdSteps A
                    let s = 
                        D.D.Data 
                        |> Array.toList
                        |> List.sortDescending

                    let nps = [1.4838999; 1.; 0.6738999]
                    List.zip s nps
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Singular values are not equal")                    
                }

                test "np.lingalg.svd reference test" {
                    let A =
                        [[|1.; 2.; 3.; 12.|] |> Vector
                         [|4.; -5.; 6.; 11. |] |> Vector
                         [|7.; 8.; -9.; 10.|] |> Vector
                         [|10.; -11.; 12.; 9.|] |> Vector
                         [|-2.; 5.; 7.; -3.|] |> Vector]
                        |> Matrix.FromRowVectors
                    
                    let D, Ut, V = Svd A
                    
                    let s = 
                        D.D.Data 
                        |> Array.toList
                        |> List.sortDescending
                    let nps = [26.25064582; 19.15201094; 9.62822768;  5.86526235]
                    List.zip s nps
                    |> List.iter (fun (x, y) -> Expect.floatClose Accuracy.medium x y <| "Singular values are not equal") 
                
                    let A2 = Ut * D * V
                    Expect.equal A2 A "Reconstituted matrix U' * D * V does not equal A"
                }

                testProp "Svd on bidiagonal matrix yields diagonal D and orthogonal U and V with D + E =  U' * B * V" <| fun (Bs: BidiagonalMatrix) ->
                    let (BidiagonalMatrix B) = Bs
                    let D, Ut, V = SvdSteps B

                    Expect.isTrue Ut.IsOrthogonal "Matrix U is not orthogonal"
                    Expect.isTrue V.IsOrthogonal "Matrix V is not orthogonal"
                    Expect.isTrue D.IsDiagonal "Matrix D is not diagonal"

                    let D2 = Ut * B * V
                    let diff = D - D2 //TODO: Tolerance based on frobenius norm as stated on p.455 Golub, Van Loan 3rd Ed.
                    diff |> Matrix.iter (fun f -> Expect.isTrue (f < 1e-3) "D-D2 element expected to be zero")
            ]
        ]