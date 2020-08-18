namespace FsLinAlg.Test

open Expecto

open FsLinAlg
open FsLinAlg.DataStructure
open FsLinAlg.Factorization
open DataStructure

module Factorization =

    let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let factorizationTests =
        testList "A=PLU" [
            testProp "Matrix Multiplication Equivalency" <| fun (As: SquareMatrix) ->
                let (SquareMatrix A) = As
                
                try
                    let P, L, U = A.LU
                    let ``PA*`` = L * U
                    let PA = P * A
    
                    let d = PA.Data |> Seq.cast<float>
                    let ``d*`` =  ``PA*``.Data |> Seq.cast<float>
                    let data = Seq.zip d ``d*`` |> Seq.toList
                    data
                    |> Seq.iter (fun (``pa*``, pa) -> Expect.floatClose Accuracy.high ``pa*`` pa "Elements do not match")
                with
                | :? LinearDependenceException as ex -> () // Expected not to work on linearly dependent matrices.
        ]