namespace FsLinAlg.Test

open Expecto
open FsCheck

open FsLinAlg

module DataStructure =

    //let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let properties =
        testList "Matrix creation consistency tests" [
            testProp "Array2D equivalence" <| fun (m: NormalFloat[,]) ->
                let fm = m |> Array2D.map (fun v -> v.Get)
                let mat = Matrix fm
                let mm = mat.Data
                fm = mm
        ]