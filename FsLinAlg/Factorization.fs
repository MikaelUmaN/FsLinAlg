namespace FsLinAlg

[<AutoOpen>]
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
                if isZeroStrict U.[k, k] then
                     raise linDep
                else
                    L.[j, k] <- U.[j, k] / U.[k, k]
                    U.[j, k..m-1] <- U.[j, k..m-1] - L.[j, k]*U.[k, k..m-1]
        P, L, U

    let private QRinner(A: Matrix) =
        let R = A.Clone()
        let m = A.M
        let n = A.N

        let rec triangulate k uks =
            let x = R.[k.., k]
            let sn = float(sign x.[0])
            let vk = sn * x.Norm * Vector.e 0 x.Length + x

            // Forming R
            let uk =
                if not <| isZeroStrict vk.Norm then
                    let u = vk / vk.Norm
                    let q = 2.*u*u.T
                    R.[k.., k..] <- R.[k.., k..] - q * R.[k.., k..]
                    u
                else
                    Vector.zero vk.Length

            let uksn = uk::uks
            if k <> n-1 then triangulate (k+1) uksn else uksn

        let uks = triangulate 0 [] |> List.zip [n-1..-1..0]
        R, uks

    /// Returns Product Q.T * b and reduced R from QR factorization.
    let QRb (A: Matrix) (bs: Vector) =
        let b = bs.Clone()
        let n = A.N

        let R, uks = QRinner A
        let ruks = uks |> List.rev
        let qb = List.fold (fun (x: Vector) (k, v: Vector) -> x.[k..] <- x.[k..] - 2.*v*(v.T *+ x.[k..]); x)
        
        let Qtb = qb b ruks
        let Rr = R.[..n-1, ..n-1]
        Qtb, Rr

    /// Reduced QR, Q m by n
    let QR (A: Matrix) =
        let m = A.M
        let n = A.N

        let R, uks = QRinner A
        let qx = List.fold (fun (x: Vector) (k, v: Vector) -> x.[k..] <- x.[k..] - 2.*v*(v.T *+ x.[k..]); x)
        let qks = 
            [for i in 0..n-1 -> 
                let ei = Vector.e i m
                qx ei uks]
        let Q = qks |> Matrix.FromColumnVectors

        let Rr = R.[..n-1, ..n-1]
        Q, Rr

    /// Reduction to Hessenberg form using two-sided Householder reflections.
    /// Tridiagonal for symmetric matrices.
    /// Returns the reduced matrix as well as vectors that can be used
    /// to reconstruct Q.
    let Hessenberg (A: Matrix) =
        let m = A.M
        if not A.IsSquare then raise <| invDimMsg "Expected square matrix"
        let issym = A.IsSymmetric

        let H = A.Clone()
        if m <= 2 then 
            H, []
        else
            let rec triangulate k uks =
                let x = H.[k+1.., k]
                let sn = float(sign x.[0])
                let vt = x.Norm * Vector.e 0 x.Length + x
                let vk = if sn < 0. then -1. * vt else vt

                let uk = 
                    if not <| isZeroStrict vk.Norm then
                        let u = vk / vk.Norm
                        H.[k+1.., k..] <- H.[k+1.., k..] - 2.*u*(u.T*H.[k+1.., k..]) // Householder refl applied on the left.
                        let rs = if issym then k else 0
                        H.[rs.., k+1..] <- H.[rs.., k+1..] - 2.*(H.[rs.., k+1..]*u)*u.T // Householder refl applied to the right.
                        u
                    else 
                        Vector.zero vk.Length
                if k = m-3 then uk::uks else triangulate (k+1) (uk::uks)

            let uks = triangulate 0 []
            H, uks

    /// Reduction to bidiagonal form using Householder reflections for rectangular matrix A: m by n,
    /// where m >= n. Returns the reduced matrix as well as vectors that can be used
    /// to reconstruct Q.
    let Bidiagonalize (A: Matrix) =
        if A.M < A.N then raise <| invDimMsg $"Expected M >= N, but {A.M} <= {A.N}"

        let B = A.Clone()
        let n = B.N

        let rec bidiagonalize k uks vks =
            // Introduce zeros in column k.
            let x = B.[k.., k]
            let sn = float(sign x.[0])
            let ut = sn * x.Norm * Vector.e 0 x.Length + x
            let uk = if sn < 0. then -1. * ut else ut
            let ukn = uk / uk.Norm

            // U' * A
            B.[k.., k..] <- B.[k.., k..] - 2.*ukn*(ukn.T*B.[k.., k..])

            if k < n-2 then
                // Introduce zeros in row k.
                let y = B.[k, k+1..]
                let sny = float(sign y.[0])
                let vt = sny * y.Norm * Vector.e 0 y.Length + y
                let vk = if sny < 0. then -1. * vt else vt
                let vkn = vk / vk.Norm

                // A * V
                B.[k.., k+1..] <- B.[k.., k+1..] - (B.[k.., k+1..]*2.*vkn)*vkn.T

                bidiagonalize (k+1) (uk::uks) (vk::vks)
            elif k < n-1 then
                bidiagonalize (k+1) (uk::uks) vks
            else
                (uk::uks), vks

        let uks, vks = bidiagonalize 0 [] []
        B, uks, vks

    /// A = R.T R
    /// Elimination on a symmetric matrix.
    /// Returns upper-triangular matrix R.
    let Cholesky (A: Matrix) =
        let R = A.UpperTriangular
        let m = R.M

        for k in 0..m-1 do
            for j in k+1..m-1 do
                if isZeroStrict R.[k, k] then
                    raise notPosDef
                R.[j, j..] <- R.[j, j..] - (R.[k, j]*R.[k, j..])/R.[k, k]
            if R.[k, k] <= 0. || isZeroStrict R.[k, k] then
                raise notPosDef
            R.[k, k..] <- R.[k, k..]/sqrt(R.[k, k])   
        R

    type Matrix with
    
        /// Gaussian elimination with partial pivoting.
        member this.LU = PLU this

        /// Reduced QR using Householder reflections for triangularization.
        member this.QR = QR this

        member this.QRb b = QRb this b

        /// Returns upper-triangular matrix R for which A=R.T * R
        member this.Cholesky = Cholesky this

        member this.Hessenberg = Hessenberg this

        member this.Bidiagonalize = Bidiagonalize this