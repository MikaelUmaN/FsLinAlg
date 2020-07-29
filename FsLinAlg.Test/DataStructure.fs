namespace FsLinAlg.Test

open Expecto

open FsLinAlg.DataStructure

module DataStructure =

    //let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let properties =
        testList "Matrix creation consistency tests" [
            testProperty "Array2D equivalence" <| fun (m: float[,]) ->
                let mat = Matrix m
                let mm = mat.Data
                m = mm
        ]