#r "../FsLinAlg/bin/Debug/net8.0/FsLinAlg.dll"

open FsLinAlg

// Matrix to zero out elements on (same as wikipedia).
let data: float[,] = array2D [ 
                                [6.; 5.; 0.]
                                [5.; 4.; 1.]
                                [0.; 4.; 3.] ]
let A = Matrix data

// zeroing out the 5, (1, 0), by rotation in the (i, k) = (0, 1) plane, affecting rows/columns 0 and 1.
let (c, s) = givensNumbers A[0, 0] A[1, 0]
let G = givens 3 0 1 c s
let GA = G * A

// zeroing out the 4, (2, 1), by rotation in the (i, k) = (0, 2) plane, affecting rows/columns 0 and 2.
// Note: this won't triangularize the matrix as we now affect the (2, 0) element.
let (c1, s1) = givensNumbers GA[0, 1] GA[2, 1]
let G2 = givens 3 0 2 c1 s1

let GGA = G2 * GA

// Now rotating in the (i, k) = (1, 2) plane.
let (c1_2, s1_2) = givensNumbers GA[1, 1] GA[2, 1]
let G2_2 = givens 3 1 2 c1_2 s1_2

let GGA_2 = G2_2 * GA

// Now, create a type to represent a givens rotation but not materialize it as a matrix.
type GivensRotation(i: int, k: int, c: float, s: float) =
    member _.I = i
    member _.K = k
    member _.C = c
    member _.S = s


    static member (*) (rotation: GivensRotation, A: Matrix) =
        let G = givens A.N rotation.I rotation.K rotation.C rotation.S
        G * A

    static member (*) (A: Matrix, rotation: GivensRotation) =
        let G = givens A.N rotation.I rotation.K rotation.C rotation.S
        A * G

let rot1 = GivensRotation(0, 1, c, s)
let GA_2 = rot1 * A

// Produces the same result.
let GA_diff = GA_2 - GA
