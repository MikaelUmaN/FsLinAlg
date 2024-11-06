#r "../FsLinAlg/bin/Debug/net8.0/FsLinAlg.dll"

open System
open FsLinAlg


// Example from Golub, Van Loand 3rd Ed. p. 216
let x = Vec [| 1.; 2.; 3.; 4.; |]
//let x = [Vec [| 1.; 2.; 3.; 4.; |]] |> Matrix.FromColumnVectors

// radian to degrees...
let radianToDegree r =
    180./Math.PI * r

let degreeToRadian d =
    d * Math.PI/180.

// From cos theta to theta in radians and then degrees.
let toDegrees c s =
    let thetaCc = acos c |> radianToDegree //63
    let thetaCs = asin c |> radianToDegree
    let thetaCcn = acos -c |> radianToDegree
    let thetaCsn = asin -c |> radianToDegree
    
    let thetaSc = acos s |> radianToDegree
    let thetaSs = asin s |> radianToDegree //63
    let thetaScn = acos -s |> radianToDegree
    let thetaSsn = asin -s |> radianToDegree

    (thetaCc, thetaCs, thetaCcn, thetaCsn, thetaSc, thetaSs, thetaScn, thetaSsn)

// gives back cos theta, sin theta
let (c, s) = givensNumbers 2. 4.
let (c1, s1) = givensNumbers -4. 2.
let (c2, s2) = givensNumbers -1. -2.

let (thetaCc, thetaCs, thetaCcn, thetaCsn, thetaSc, thetaSs, thetaScn, thetaSsn) = toDegrees c s

let G = givens 4 1 3 c s
G * x

x.T*G.T

// Example from wikipedia.
let A = array2D [
                [|6.;5.;0.|]
                [|5.;1.;4.|]
                [|0.;4.;3.|]
            ] |> Matrix

// from left: 0, 1 <- zero 1
let (c, s) = givensNumbers A[0, 0] A[1, 0]
let G1 = givens A.N 0 1 c s
let A2 = G1 * A

// from left: 1, 2
let (c1, s1) = givensNumbers A2[1, 1] A2[2, 1]
let G2 = givens A.N 1 2 c1 s1
let A3 = G2 * A2