namespace FsLinAlg

open System

/// Contains datastructures that have specialized usages.
[<AutoOpen>]
module SpecializedStructure =


    /// TODO: We want an interface here that enables us to do cost-efficient multiplications
    /// with givens matrices, as if they were any matrices at all.

    /// The product G(i, k, θ)x represents a counterclockwise rotation of the vector x in the (i, k) 
    /// plane of θ radians, hence the name Givens rotation.
    /// 
    /// When a Givens rotation matrix, G(i, j, θ), multiplies another matrix, A, from the left, 
    /// G*A, only rows i and j of A are affected.
    /// 
    /// The main use of Givens rotations in numerical linear algebra is to introduce zeros in vectors or matrices. 
    /// This effect can, for example, be employed for computing the QR decomposition of a matrix. 
    /// One advantage over Householder transformations is that they can easily be parallelised, 
    /// and another is that often for very sparse matrices they have a lower operation count.
    let givens n i k c s =
        if k <= i then
            raise <| ArgumentException($"k must be strictly larger than i, but {i} >= {k}")
        else
            let G = Matrix.I n
            G.[i, i] <- c
            G.[k, k] <- c
            G.[i, k] <- s
            G.[k, i] <- -s
            G

    /// Returns numbers (c, s) to be used in a Givens rotation matrix.
    /// asin s = theta
    /// acos c = theta
    /// Note: acos and asin use a different reference for the angle so the value will only be the same
    /// when they are in the first quadrant.
    /// Note: Does not follow Golub Van Loan 3rd Ed. Eq. 5.1.7, instead follows Wikipedia.
    let givensNumbers a b =
        if b = 0. then
            signv a, 0.
        elif a = 0. then
            0., signv b
        elif abs a > abs b then
            let t = b / a
            let u = signv a * sqrt(1. + t**2.)
            let c = 1./u
            c, c*t
        else
            let t = a / b
            let u = signv b * sqrt(1. + t**2.)
            let s = 1./u
            s*t, s
            
    /// Returns a given rotation matrix.
    /// Note: Does not follow Golub Van Loan 3rd Ed. Eq. 5.1.7, instead follows Wikipedia.
    /// They are equivalent up to a transpose operation.
    let givensMatrix n i k a b =
        let (c, s) = givensNumbers a b
        givens n i k c s

    /// Zeros the entire row, because zeroing one element affects the nearby element
    /// on the same row.
    /// <returns>
    /// Tuple of matrix B* (with row i zeroed) and matrix U representing the
    /// operation of givens matrices pre-multiplying B to give B*.
    /// </returns>
    let zeroRow (B: Matrix) (i: int) =
        let GB, Us =
            List.fold (fun (GB: Matrix, Us: List<Matrix>) k -> 
                let U = givensMatrix B.N i k GB.[k, k] GB.[i, k] 
                (U.T * GB, U::Us)) (B, []) [i+1..B.N-1]

        // Accumulate givens matrices (G1 * G2 * G3) and transpose to represent the operation:
        // G3' * G2' * G1' * B22
        let Utd = (List.fold (fun (U: Matrix) (u: Matrix) -> u * U) (Matrix.I B.M) Us).T

        GB, Utd


    /// Zeros the entire column, because zeroing one element offsets the nearby element
    /// on the same column.
    /// <returns>
    /// Tuple of matrix B* (with column k zeroed) and matrix V representing the 
    /// operation of givens matrices post-multiplying B to give B*.
    /// </returns>
    let zeroColumn (B: Matrix) (k: int) =
        let BG, Vs =
            List.fold (fun (BG: Matrix, Vs: List<Matrix>) i -> 
                let V = givensMatrix B.N i k BG.[i, i] BG.[i, k] 
                (BG * V.T, V.T::Vs)) (B, []) [k-1..-1..0]


        // Accumulate givens matrices (G1 * G2 * G3) to represent the operation:
        // B * G1 * G2 * G3
        let Vd = List.fold (fun (V: Matrix) (v: Matrix) -> V * v) (Matrix.I B.N) (Vs |> List.rev)

        BG, Vd