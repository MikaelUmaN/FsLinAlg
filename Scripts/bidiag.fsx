#r "../FsLinAlg/bin/Debug/net5.0/FsLinAlg.dll"

open System
open FsLinAlg

// Example Golub Van Loan 3rd Ed. p. 252
let A = [
            Vec [|1.;2.;3.|]
            Vec [|4.;5.;6.|]
            Vec [|7.;8.;9.|]
            Vec [|10.;11.;12.|]
        ]
        |> Matrix.FromRowVectors

A.Bidiagonalize 