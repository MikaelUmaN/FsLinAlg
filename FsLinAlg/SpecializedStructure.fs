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
            (sign a |> float, 0.)
        elif a = 0. then
            (0., sign b |> float)
        elif abs a > abs b then
            let t = b / a
            let u = (sign a |> float) * sqrt(1. + t**2.)
            let c = 1./u
            (c, c * t)
        else
            let t = a / b
            let u = (sign b |> float) * sqrt(1. + t**2.)
            let s = 1./u
            (s*t, s)

    /// Returns a given rotation matrix.
    /// Note: Does not follow Golub Van Loan 3rd Ed. Eq. 5.1.7, instead follows Wikipedia.
    /// They are equivalent up to a transpose operation.
    let givensMatrix n i k a b =
        let (c, s) = givensNumbers a b
        givens n i k c s