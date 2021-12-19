#r "/home/jovyan/work/FsLinAlg/FsLinAlg/bin/Debug/net6.0/FsLinAlg.dll"

open System
open FsLinAlg

// Matrix that fails...???
let A = [ 
        Vec [|1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0|]
        Vec [|0.0; 1.0; -0.4146722985; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0|]
        Vec [|0.0; 0.0; 1.0; -0.9820034309; 0.0; 0.0; 0.0; 0.0; 0.0|]
        Vec [|0.0; 0.0; 0.0; 1.0; 0.5544461998; 0.0; 0.0; 0.0; 0.0|]
        Vec [|0.0; 0.0; 0.0; 0.0; 1.0; -0.8537384044; 0.0; 0.0; 0.0|]
        Vec [|0.0; 0.0; 0.0; 0.0; 0.0; 1.0; -0.5603501902; 0.0; 0.0|]
        Vec [|0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.270416723; 0.0|]
        Vec [|0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; -0.8415624747|]
        Vec [|0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0|]
        ] |> Matrix.FromRowVectors

let D, Ut, V = SvdSteps A













// Example Golub Van Loan 3rd Ed. p. 252
let A = [
            Vec [|1.;2.;3.|]
            Vec [|4.;5.;6.|]
            Vec [|7.;8.;9.|]
            Vec [|10.;11.;12.|]
        ]
        |> Matrix.FromRowVectors

let B, ukst, vkst = A.Bidiagonalize

// TODO: rev them in function?
let uks = List.rev ukst
let vks = List.rev vkst

// For the U' matrix, we get M by M
let m = B.M
let Qu = Matrix.I B.M
// n equals the dimension to create, for U' it will be m
let n = m
// r equals the number of vectors we have to work with.
let r = uks.Length
// starting at the end, the vector with the least number of updates.
let j = r-1

// backaccumulate, first householder vector -> matrix update.
let u1 = uks.[j]

// Because we norm the vector, we don't follow golub-kahan formula where they append back the 1.
//let vu1 = Vector.concat [Vector.e 0 1; u1.[j..]]
// First update, lower right corner.
Qu.[j.., j..] <- (Matrix.I (n-j) - 2.*u1*u1.T) * Qu.[j.., j..]

// Next update, is a bigger corner.
let j = r-2
let u2 = uks.[j]
Qu.[j.., j..] <- (Matrix.I (n-j) - 2.*u2*u2.T) * Qu.[j.., j..]

// final update, last corner.
let j = r-3
let u3 = uks.[j]
Qu.[j.., j..] <- (Matrix.I (n-j) - 2.*u3*u3.T) * Qu.[j.., j..]



// Create V using same method.
// For the U' matrix, we get M by M
let n = B.N
let Qv = Matrix.I n
// n equals the dimension to create, for U' it will be m
//let n = n

// vks is special... it needs an identity vector at the start...
let vks2 = List.concat [[Vector.e 0 n]; vks]

// r equals the number of vectors we have to work with.
let r = vks2.Length
// starting at the end, the vector with the least number of updates.
let j = r-1

// backaccumulate, first householder vector -> matrix update.
let v1 = vks2.[j]

// First update, lower right corner.
Qv.[j.., j..] <- (Matrix.I (n-j) - 2.*v1*v1.T) * Qv.[j.., j..]

// Next update is just a do nothing.
let j = r-2
let v2 = vks2.[j]
Qv.[j.., j..]  <- (Matrix.I (n-j) - 2.*v2*v2.T) * Qv.[j.., j..]

// Does it work...
Qu.T * A * Qv

// Yes! it works and they are orthogonal. Great.

/// Data structure for accumulation of Householder vectors.
/// Accepts a list of vectors, each vector is used once in the accumulation.
/// Optionally accepts the size of the Householder matrix, else it is inferred from the largest vector.
/// The vectors are supplied in the order from smallest to largest.
type HouseholderAccumulation(vs: List<Vector>, ?N: int) =
    do
        if vs.IsEmpty then
            raise <| ArgumentException($"vs must contain at least one vector but was empty")

    let r = vs.Length
    let n = if N.IsNone then vs.[^0].Length else N.Value

    /// Accumulates the vectors into the full matrix representation.
    member _.Accumulate =
        let Q = Matrix.I n
        for j in 0..r-1 do
            let v = vs.[j]
            Q.[^j+1.., ^j+1..]  <- (Matrix.I (j+2) - 2.*v*v.T) * Q.[^j+1.., ^j+1..]
        Q

let h = HouseholderAccumulation(ukst)
let Qu = h.Accumulate

let hv = HouseholderAccumulation(vkst, vkst.[^0].Length + 1)
let Qv = hv.Accumulate

// Does it work...
let B2 = Qu.T * A * Qv
B2 = B

let vv = Vector.e 0 3
let vvv = -1. * Vector.e 0 3