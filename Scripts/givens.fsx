#r "../FsLinAlg/bin/Debug/net6.0/FsLinAlg.dll"

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
let A = [
            Vec [|6.;5.;0.|]
            Vec [|5.;1.;4.|]
            Vec [|0.;4.;3.|]
        ]
        |> Matrix.FromColumnVectors

// from left: 0, 1 <- zero 1
let (c, s) = givensNumbers 6. 5.
let G1 = givens A.N 0 1 c s
let A2 = G1 * A

// from left: 1, 2
let (c1, s1) = givensNumbers -2.4327 4.
let G2 = givens A.N 1 2 c1 s1
let A3 = G2 * A2



// try from right side..
// Example from wikipedia.
let A = [
            Vec [|6.;5.;0.|]
            Vec [|5.;1.;4.|]
            Vec [|0.;4.;3.|]
        ]
        |> Matrix.FromColumnVectors

// from right: 1, 2 <- zero 4
let (c, s) = givensNumbers 1. 4.
let G1 = givens A.N 1 2 c s
let A2 = A * G1.T

// from right: 2, 3 <- zero -3.152963125
let (c1, s1) = givensNumbers 0. -3.152963125
let G2 = givens A.N 0 2 c1 s1
let A3 = A2 * G2.T
