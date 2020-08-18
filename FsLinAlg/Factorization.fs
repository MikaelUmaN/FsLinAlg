namespace FsLinAlg

open DataStructure

module Factorization =

    let PLU (A: Matrix) =
        let L = Matrix.I A.M
        let P = Matrix.I A.M
        let U = A.Clone()
        let m = A.M

        for k in 0..m-2 do

            // Pivoting
            let i = (U.[k.., k] |> Vector.map abs |> Vector.findMaxIndex) + k
            let up = U.[i, *]
            U.[i, *] <- U.[k, *]
            U.[k, *] <- up

            let lp = L.[i, ..k-1]
            L.[i, ..k-1] <- L.[k, ..k-1]
            L.[k, ..k-1] <- lp

            let pp = P.[i, *]
            P.[i, *] <- P.[k, *]
            P.[k, *] <- pp

            for j in k+1..m-1 do
                if U.[k, k] = 0. then
                     raise linDep
                else
                    L.[j, k] <- U.[j, k] / U.[k, k]
                    U.[j, k..m-1] <- U.[j, k..m-1] - L.[j, k]*U.[k, k..m-1]
        P, L, U

    type Matrix with
    
        /// Gaussian elimination with partial pivoting.
        member this.LU = PLU this