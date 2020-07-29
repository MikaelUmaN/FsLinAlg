namespace FsLinAlg.Test

open Expecto

open FsLinAlg.DataStructure
open FsLinAlg.Factorization

module Factorization =

    //let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let factorizationTests =
        testList "A=PLU" [
            test "Matrix Multiplication Equivalency" {
                
                let A =
                    Matrix.FromRowVectors([
                        V [| 3.; 2.; 10. |]
                        V [| 1.; 5.; 1. |]
                        V [| 2.; 4.; -2. |]
                    ])
                let P, L, U = A.LU

                let ``A'`` = P * L * U

                Expect.equal ``A'`` A "Multiplication PLU should have equaled A" // No float close array check?
            }
        ]