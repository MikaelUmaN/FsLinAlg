namespace FsLinAlg

[<AutoOpen>]
module Eigenvalue =

    /// QR Algorithm with shifts for revealing eigenvalues in
    /// the diagonal of the tridiagonal matrix S using
    /// similarity transforms.
    let ShiftedQR (S: Matrix) maxIterOpt =
        if S.M = 2 then // Special case where Tridiagonal = Diagonal
            S.D
            |> Vector.toList
        else
            let wilkinsonShift (B: Matrix) =
                let d = (B.[0, 0] - B.[1, 1]) / 2.
                let ds = if d < 0. then -1. else 1.
                B.[1, 1] - ds*B.[0, 1]**2. / (abs ds + sqrt(d**2. + B.[0, 1]**2.))

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