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

    /// Reduced QR, Q m by n
    let QR (A: Matrix) =
        let R = A.Clone()
        let m = A.M
        let n = A.N

        let rec triangulate k Qp =
            let x = R.[k.., k]
            let sn = float(sign x.[0])
            let vk = sn * x.Norm * Vector.e 0 x.Length + x

            // Forming R
            let qk =
                if vk.Norm <> 0. then
                    let u = vk / vk.Norm
                    let q = 2.*u*u.T
                    R.[k.., k..] <- R.[k.., k..] - q * R.[k.., k..]
                    q
                else
                    Matrix.I (n-k)

            // Forming Q
            let qkm = Matrix.I m
            qkm.[k.., k..] <- qk
            let Qn = Matrix.I m - qkm
            let Qk = Qp * Qn.T

            if k <> n-1 then
                // Because of machine precision rounding errors.
                R.[k+1.., k] <- m-1-k |> Vector.zero
                triangulate (k+1) Qk
            else
                Qk

        let Q0 = Matrix.I m
        let Q = triangulate 0 Q0

        Q, R

    type Matrix with
    
        /// Gaussian elimination with partial pivoting.
        member this.LU = PLU this

        /// Reduced QR using Householder reflections for triangularization.
        member this.QR = QR this