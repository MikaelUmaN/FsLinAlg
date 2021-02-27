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



let A =
    [[|1.0; 1.0; 0.0; 0.0|] |> Vector
     [|0.0; 2.0; 1.0; 0.0|] |> Vector
     [|0.0; 0.0; 3.0; 1.0|] |> Vector
     [|0.0; 0.0; 0.0; 4.0|] |> Vector]
    |> Matrix.FromRowVectors

// Example 1, 1 is zero. zero out the entire row by a series of givens rotations.
let B =
    [[|1.0; 1.0; 0.0; 0.0; 0.0|] |> Vector
     [|0.0; 0.0; 1.0; 0.0; 0.0|] |> Vector
     [|0.0; 0.0; 3.0; 1.0; 0.0|] |> Vector
     [|0.0; 0.0; 0.0; 4.0; 1.0|] |> Vector
     [|0.0; 0.0; 0.0; 0.0; 3.0|] |> Vector]
     |> Matrix.FromRowVectors


// This works for zeroing out a row.
let G2 = givensMatrix B.N 1 2 B.[2, 2] B.[1, 2]
let gB = G2.T * B

let G3 = givensMatrix B.N 1 3 gB.[3, 3] gB.[1, 3]
let ggB = G3.T * gB

let G4 = givensMatrix B.N 1 4 ggB.[4, 4] ggB.[1, 4]
let gggB = G4.T * ggB

// Zero out the last column...
// works.
let B =
    [[|1.0; 1.0; 0.0; 0.0; 0.0|] |> Vector
     [|0.0; 5.0; 1.0; 0.0; 0.0|] |> Vector
     [|0.0; 0.0; 3.0; 1.0; 0.0|] |> Vector
     [|0.0; 0.0; 0.0; 4.0; 1.0|] |> Vector
     [|0.0; 0.0; 0.0; 0.0; 0.0|] |> Vector]
     |> Matrix.FromRowVectors

let V = givensMatrix B.N (B.N-2) (B.N-1) B.[B.N-2, B.N-2] B.[B.N-2, B.N-1]
let Bv = B * V.T

let V2 = givensMatrix B.N (B.N-3) (B.N-1) Bv.[B.N-3, B.N-3] Bv.[B.N-3, B.N-1]
let Bvv = Bv * V2.T

let V3 = givensMatrix B.N (B.N-4) (B.N-1) Bvv.[B.N-4, B.N-4] Bvv.[B.N-4, B.N-1]
let Bvvv = Bvv * V3.T

let V4 = givensMatrix B.N (B.N-5) (B.N-1) Bvvv.[B.N-5, B.N-5] Bvvv.[B.N-5, B.N-1]
let Bvvvv = Bvvv * V4.T