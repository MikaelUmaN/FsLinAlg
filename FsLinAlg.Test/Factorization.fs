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
                test "scip.lingalg.hessenberg reference test" {
                    (*
                        From Python. Some difference in precision.

                        B = np.array([
                            [4.0, 0.9999960883, 5.741236271e-19, 1.0],
                            [0.9999960883, 5.0, 0.5, 0.5],
                            [5.741236271e-19, 0.5, 4.0, 0.497615855],
                            [1.0, 0.5, 0.497615855, 5.0]
                        ])

                        H = sc.hessenberg(B)
                        array([[ 4.00000000e+00, -1.41421080e+00, -1.03997563e-23,
                                -3.75089678e-18],
                               [-1.41421080e+00,  5.50000000e+00,  7.05420933e-01,
                                 1.92714487e-16],
                               [ 0.00000000e+00,  7.05420933e-01,  4.00000001e+00,
                                 1.68861110e-03],
                               [ 0.00000000e+00,  0.00000000e+00,  1.68861110e-03,
                                 4.49999999e+00]])

                    *)
                    let data =
                        [[4.0; 0.9999960883; 5.741236271e-19; 1.0]
                         [0.9999960883; 5.0; 0.5; 0.5]
                         [5.741236271e-19; 0.5; 4.0; 0.497615855]
                         [1.0; 0.5; 0.497615855; 5.0]]
                        |> List.map (List.toArray >> fun d -> Vector(d))
                    let A = Matrix.FromRowVectors(data)
                    let H, _ = A.Hessenberg
                    Expect.isTrue H.IsTridiagonal "Matrix is not tridiagonal"
                }

                test "scipy.linalg.hessenberg reference test 2" {

                    (*
                        From Python.

                        A = np.array([
                                [7.0, 0.0, -1.0, -0.5, -0.5, -1.0, -0.5],
                                [0.0, 8.0, 1.115662563, 2.775557562e-17, 2.775557562e-17, 2.220446049e-16, 2.775557562e-17],
                                [-1.0, 1.115662563, 9.52809281, 0.7395713062, -1.221245327e-15, -2.164934898e-15, 1.33226763e-15],
                                [-0.5, 2.775557562e-17, 0.7395713062, 6.471699392, 0.2352372022, 4.163336342e-17, 1.387778781e-16],
                                [-0.5, 2.775557562e-17, 0.0, 0.2352372022, 7.177027183, 0.4032641895, -8.326672685e-16],
                                [-1.0, 2.220446049e-16, 5.551115123e-17, -2.775557562e-17, 0.4032641895, 7.216453431, -0.05823811227],
                                [-0.5, 2.775557562e-17, 0.0, 5.551115123e-17, 0.0, -0.05823811227, 5.24387933]])

                        H = scipy.linalg.hessenberg(A)
                    *)
                    let data =
                        [[7.0; 0.0; -1.0; -0.5; -0.5; -1.0; -0.5]
                         [0.0; 8.0; 1.115662563; 2.775557562e-17; 2.775557562e-17; 2.220446049e-16;
                          2.775557562e-17]
                         [-1.0; 1.115662563; 9.52809281; 0.7395713062; -1.221245327e-15;
                          -2.164934898e-15; 1.33226763e-15]
                         [-0.5; 2.775557562e-17; 0.7395713062; 6.471699392; 0.2352372022;
                          4.163336342e-17; 1.387778781e-16]
                         [-0.5; 2.775557562e-17; 0.0; 0.2352372022; 7.177027183; 0.4032641895;
                          -8.326672685e-16]
                         [-1.0; 2.220446049e-16; 5.551115123e-17; -2.775557562e-17; 0.4032641895;
                          7.216453431; -0.05823811227]
                         [-0.5; 2.775557562e-17; 0.0; 5.551115123e-17; 0.0; -0.05823811227; 5.24387933]]
                        |> List.map (List.toArray >> fun d -> Vector(d))
                    let A = Matrix.FromRowVectors(data)
                    let H, _ = A.Hessenberg
                    Expect.isTrue H.IsTridiagonal "Matrix is not tridiagonal"
                }

                testProp "Upper Hessenberg result" <| fun (As: SquareMatrix) ->
                    let (SquareMatrix A) = As
                    let H, _ = A.Hessenberg
                    Expect.isTrue (H.IsHessenberg()) "Matrix is not of hessenberg form"

                testProp "Tridiagonal result" <| fun (As: SymmetricPositiveDefiniteMatrix) ->
                    let (SymmetricPositiveDefiniteMatrix A) = As
                    let H, _ = A.Hessenberg
                    Expect.isTrue H.IsTridiagonal "Matrix is not tridiagonal"
            ]

            testList "Bidiagonalization" [
                testProp "Bidiagonal result" <| fun (Ar: TallThinMatrix) ->
                        let (TallThinMatrix A) = Ar
                        let B, _, _ = A.Bidiagonalize

                        Expect.isTrue (B.IsBidiagonal()) "Matrix is not bidiagonal"
            ]
        ]