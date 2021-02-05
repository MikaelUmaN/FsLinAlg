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

let B = A.Bidiagonalize

// Example Golub Van Loan 3rd Ed. p. 457
let D = [
            Vec [|1.;1.;0.;0.|]
            Vec [|0.;2.;1.;0.|]
            Vec [|0.;0.;3.;1.|]
            Vec [|0.;0.;0.;4.|]
        ]
        |> Matrix.FromRowVectors

let BD, Ut, V = SvdStep D
//let BDD, Utss, Vss = SvdStep BD
//let BDDD, Utsss, Vsss = SvdStep BDD
//let BDDDD, Utssss, Vssss = SvdStep BDDD

let Utt = 
    array2D [[0.443645989; 0.07739799259; -0.003066782529; 0.8928484653]
             [0.8962021181; -0.03831424661; 0.00151814612; -0.4419858338]
             [0.0; -0.996263805; -0.0002966377327; 0.08636169728]
             [0.0; 0.0; 0.999994101; 0.003434809553]]
    |> Matrix