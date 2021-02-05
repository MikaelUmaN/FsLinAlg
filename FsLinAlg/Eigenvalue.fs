namespace FsLinAlg

[<AutoOpen>]
module Eigenvalue =

    /// Numerically stable computation of eigenvalue closer to T.[1, 1] from Trefethen, Bau p. 222
    /// [a(m-1) b(m-1)
    ///  b(m-1) a(m)]
    let wilkinsonShift (T: Matrix) =
        let d = (T.[0, 0] - T.[1, 1]) / 2.
        // If zero, arbitrarily set to one
        let ds = if d < 0. then -1. else 1.
        T.[1, 1] - ds*T.[0, 1]**2. / (abs ds + sqrt(d**2. + T.[0, 1]**2.))
        
    let SvdStep (A: Matrix) =
        let B = A.Clone()
        let m = B.M
        let n = B.N

        let rec inner (B: Matrix) y z k (Uts: Matrix list) (Vs: Matrix list) =
            // Right-side givens rotation.
            let V = givensMatrix n k (k+1) y z 
            let BV = B * V.T

            // Left-side givens rotation
            let y = BV.[k, k]
            let z = BV.[k+1, k]
            let Ut = givensMatrix m k (k+1) y z
            let UtBV = Ut * BV

            if k < n-2 then
                let y = UtBV.[k, k+1]
                let z = UtBV.[k, k+2]
                inner UtBV y z (k+1) (Ut::Uts) (V::Vs)
            else
                UtBV, (Ut::Uts), (V::Vs)

        // Eigenvalue estimate
        let T = 
            let a = B.[m-1, m-1]**2. + B.[m-2, m-1]
            let b = B.[m-1, m-1]*B.[m-2, m-1] 
            array2D [[ a; b ]; [b; a]]
            |> Matrix
        let mu = wilkinsonShift T
        let y = T.[0, 0] - mu
        let z = T.[0, 1]

        let UtBV, Uts, Vs = inner B y z 0 [] []

        let I = Matrix.I m
        let Ut = List.fold (fun (U: Matrix) (ut: Matrix) -> ut * U) I Uts
        let V = List.fold (fun (V: Matrix) (v: Matrix) -> V * v) I Vs

        UtBV, Ut, V

    /// QR Algorithm with shifts for revealing eigenvalues in
    /// the diagonal of the tridiagonal matrix S using
    /// similarity transforms.
    let ShiftedQR (S: Matrix) maxIterOpt =
        if S.M = 2 then // Special case where Tridiagonal = Diagonal
            S.D
            |> Vector.toList
        else
            let maxItr = defaultArg maxIterOpt 10000
            let A = S.Clone()

            let rec qr (Ap: Matrix) k =
                if k >= maxItr then
                    raise <| maxIter k
                else
                    let mu = wilkinsonShift Ap.[Ap.M-2.., Ap.M-2..]
                    let shift = mu * (Matrix.I Ap.M)
                    let Q, R = (Ap - shift).QR
                    let An = R*Q + shift

                    let shouldDeflate = isZeroStrict An.[An.M-2, An.M-1]
                    if shouldDeflate then
                        if An.M = 2 then
                            [An.[0, 0]; An.[1, 1]]
                        else
                            let A2e = [An.[An.M-1, An.M-1]]
                            let A1 = An.[..An.M-2, ..An.M-2]
                            let A1e = qr A1 0
                            [A1e; A2e] |> List.collect (id)
                    else
                        qr An (k+1)

            qr A 0

    /// Arnoldi iteration for eigenvalue computation of non-hermitian matrices.
    /// Iterative projection onto successive Krylov subspaces.
    /// The eigenvalues of Hessenberg matrix Hn are the Arnoldi eigenvalue estimates of A.
    let Arnoldi (A: Matrix) maxIterOpt =
        let m = A.M
        let maxItr = defaultArg maxIterOpt m

        let b = (Vector.randn m).AsUnitVector

        let subv (v: Vector, hs: list<float>) (qj: Vector) =
            let hjn = (qj.T*v).AsScalar
            let s = hjn*qj
            v - s, hjn::hs

        let rec inner (qs: list<Vector>) (qn: Vector) n (hvs: list<Vector>) =
            let v = (A * qn).AsVector
            let (hnn, hs) = List.fold subv (v, []) qs
            let qnn = v / hnn.AsUnitVector
            let hv = Vector(hs |> List.rev |> List.toArray)
            let hvsn = hvs@[hv]

            if n < m-1 then
                inner (qs@[qnn]) qnn (n+1) hvsn
            else
                let Hnn = Matrix.FromColumnVectors hvsn
                // TODO: non-tridiagonal solver for eigenvalues
                Hnn

        inner [b] b 0

    /// Lancoz iteration for eigenvalue computation of hermitian matrices.
    /// Iterative projection onto successive Krylov subspaces.
    let Lancoz (A: Matrix) maxIterOpt =
        let m = A.M
        let maxItr = defaultArg maxIterOpt m

        let q1 = (Vector.randn m).AsUnitVector

        let betas = [0.]
        let q0 = Vector.zero m

        let rec inner  (qn: Vector) (qnp: Vector) (betas: list<float>) (alphas: list<float>) n =
            let v = (A * qn).AsVector
            let an = (qn.T*v).AsScalar
            let v = v - betas.[0]*qnp - an*qn
            let bn = v.Norm
            let qnn = v / bn

            if n < m-1 then
                inner qnn qn (bn::betas) (an::alphas) (n+1) 
            else
                let b = betas |> List.rev |> List.skip 1
                let a = alphas |> List.rev
                let H = Matrix.CreateTridiagonal b a b 
                ShiftedQR H

        inner q1 q0 betas [] 0