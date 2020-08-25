namespace FsLinAlg

module Equation =

    /// Performs backsubstitution on a upper-triangular system.
    /// Solves Ax=b for x.
    let backSubstitute (A: Matrix) (b: Vector) =
        let m = A.M
        let x = Vector.zero m

        x.[m-1] <- b.[m-1] / A.[m-1, m-1]
        if m = 0 then
            x
        else
            let idx = [m-2..-1..0]
            List.fold (fun (x: Vector) (i: int) -> x.[i] <- (b.[i] - A.[i, i+1..] *+ x.[i+1..]) / A.[i, i]; x) x idx
