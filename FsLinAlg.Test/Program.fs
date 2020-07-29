namespace FsLinAlg.Test

module Program =

    open System

    open Expecto
    open FsCheck        

    let float2d = Arb.from<float[,]>

    type CustomGenerators =
        static member Array2D() =
            { new Arbitrary<float[,]>() with
                override x.Generator = 
                    float2d.Generator
                    |> Gen.filter (fun a ->
                        (Array2D.length1 a < 2 || Array2D.length1 a > 10) && (Array2D.length2 a < 2 || Array2D.length2 a > 10)
                    )
                    |> Gen.filter (fun a ->
                        let m = Array2D.length1 a
                        [|0..m-1|]
                        |> Array.forall (fun i -> a.[i, *] |> Array.forall (Double.IsNaN >> not))
                    )
                override x.Shrinker t = float2d.Shrinker t }

    [<EntryPoint>]
    let main args =
        Arb.register<CustomGenerators>() |> ignore
        runTestsInAssemblyWithCLIArgs [] args