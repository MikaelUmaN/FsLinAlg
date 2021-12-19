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
        T.[1, 1] - ds*T.[1, 0]**2. / (abs d + sqrt(d**2. + T.[1, 0]**2.))
        
    /// Finds eigenvalue of symmetric 2 by 2 matrix that is closer to T.[1, 1].
    /// Uses characteristic polynomial
    let private symmetricEigenvalue (S: Matrix) =
        let r1 = S.[1, 1] + S.[1, 0]
        let r2 = S.[1, 1] - S.[1, 0]
        if abs(r1 - S.[1, 1]) < abs(r2 - S.[1, 1]) then
            r1
        else
            r2

    /// Takes one step in the Svd algorithm whereby A is being diagonlized.
    /// We are implicitly operating on the tridiagonal matrix A'A.
    /// <param name="A">Bidiagonal Matrix m by n with m >= n.</param>
    /// <returns>Triple of matrices U' * A * V, U', V</returns>
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
                inner UtBV y z (k+1) (Ut::Uts) (V.T::Vs)
            else
                UtBV, (Ut::Uts), (V.T::Vs)

        // Eigenvalue estimate
        let T = 
            if m = 2 && n = 2 then
                B.T * B
            else
                let t11 = B.[n-2, n-2]**2. + B.[n-3, n-2]**2.
                let t = B.[n-2, n-2] * B.[n-2, n-1]
                let t22 = B.[n-1, n-1]**2. + B.[n-2, n-1]**2.
                array2D [[ t11; t ]; [t; t22]]
                |> Matrix
        let mu = wilkinsonShift T // Wilkinson shift is defined as the eigenvalue closer to t22.
        let y = T.[0, 0] - mu
        let z = T.[0, 1]

        let UtBV, Uts, Vs = inner B y z 0 [] []

        let I = Matrix.I m
        let Ut = List.fold (fun (U: Matrix) (ut: Matrix) -> ut * U) I (List.rev Uts)
        let V = List.fold (fun (V: Matrix) (v: Matrix) -> V * v) I (List.rev Vs)

        UtBV, Ut, V

    let rec SvdSteps (A: Matrix) =
        if A.M <> A.N then raise <| invDimMsg $"Square matrix expected but received: M = {A.M} <> N = {A.N}"
    
        let B = A.Clone()
        let n = B.N
    
        let rec inner (B: Matrix) (Uts: Matrix list) (Vs: Matrix list) =
    
            // 1. Check for convergence.
            // Zero any superdiagonal element that is negligible compared to its closest diagonal elements.
            for i in 0..n-2 do
                if insignificantComparedTo (abs(B.[i, i+1])) (abs(B.[i, i]) + abs(B.[i+1, i+1])) then
                    B.[i, i+1] <- 0.
    
            // Find the largest diagonal from the bottom-right. Base case is a complete diagonal matrix,
            // in which case we are done.
            let qOpt = ([n-2..-1..0] |> List.tryFind (fun i -> not (isZero B.[i, i+1])))
            let q = (if qOpt.IsSome then qOpt.Value + 1 else 0)
            if q = 0 then
                B, Uts, Vs
            else
    
                // Find the smallest diagonal from the top-left. Base case is a size zero matrix.
                let pOpt = [0..q] |> List.tryFind (fun i -> not (isZero B.[i, i+1]))
                let p = if pOpt.IsSome then pOpt.Value else 0
    
                // B22 matrix. Must have non-zero superdiagonal.
                // This is the matrix to perform work on.
                let Bd = B.[p..q, p..q]
    
                // 2. Detect decoupling.
                // Zero any superdiagonal element if the diagonal element is zero.
                let zeroDiagOpt = [0..Bd.N-1] |> List.tryFind (fun i -> isZero Bd.[i, i])
                match zeroDiagOpt with
                | Some(i) ->
                    // This addresses into B rather than Bd.
                    let k = i + p
    
                    let Utd, Vd =
                        if i = Bd.N-1 then
                            let BG, Vd = zeroColumn B k
                            B.[*, *] <- BG
                            None, Some(Vd)
                        else
                            let GB, Utd = zeroRow B k
                            B.[*, *] <- GB
                            Some(Utd), None
                    let Utsn = if Utd.IsSome then (Utd.Value::Uts) else Uts
                    let Vsn = if Vd.IsSome then (Vd.Value::Vs) else Vs
    
                    if k = B.N-1 || k = 0 then
                        // Upper left diagonal element is zero, or lower right diagonal element is zero, 
                        // no need to decouple into sub-problems.
                        // The next iteration, this part will be in the (diagonal) B11 or B33 matrix.
    
                        inner B Utsn Vsn
                    else
                        // Diagonal element in the middle is zero
                        // Decouple into two separate problems.
                        let B1, Uts1, Vs1 = SvdSteps B.[..k, ..k]
                        let B2, Uts2, Vs2 = SvdSteps B.[k+1.., k+1..]
    
                        // Overwrite B with solutions to subproblems.
                        B.[..k, ..k] <- B1
                        B.[k+1.., k+1..] <- B2
    
                        // Reshape matrices U and V to full dimension.
                        let Uts1n = 
                            Uts1 
                            |> List.map (fun ut ->
                                let Ut = Matrix.I n
                                Ut.[..k, ..k] <- ut
                                Ut)
                        let Vs1n = 
                            Vs1 
                            |> List.map (fun v ->
                                let V = Matrix.I n
                                V.[..k, ..k] <- v
                                V)
                        let Uts2n = 
                            Uts2 
                            |> List.map (fun ut ->
                                let Ut = Matrix.I n
                                Ut.[..k, ..k] <- ut
                                Ut)
                        let Vs2n = 
                            Vs2 
                            |> List.map (fun v ->
                                let V = Matrix.I n
                                V.[..k, ..k] <- v
                                V)
    
                        // By convention, apply B1 first, then B2.
                        let Utsnn = Uts2n @ Uts1n @ Utsn
                        let Vsnn = Vs2n @ Vs1n @ Vsn
    
                        // And we are done.
                        B, Utsnn, Vsnn
                | None ->
                    let UtBV, Utd, Vd = SvdStep Bd
                    B.[p..q, p..q] <- UtBV
    
                    // Reshape matrices U and V to full dimension.
                    let Ut = Matrix.I n
                    Ut.[p..q, p..q] <- Utd
                    let V = Matrix.I n
                    V.[p..q, p..q] <- Vd
    
                    // Recurse until convergence.
                    inner B (Ut::Uts) (V::Vs)
    
        inner B [] []

    /// Svd (reduced) algorithm following Golub-Kahan.
    /// Diagonalizes A to D = U' * A * V. This is the Singular Value Decomposition.
    /// <param name="A">An m by n matrix with m >= n.</param>
    /// <returns>A triple of matrices D (n by n), U (m by n) and V (n by n) such that A = U * D * V'</returns> 
    let Svd (A: Matrix) =
        if A.M < A.N then raise <| invDimMsg $"M must be >= N but was: {A.M} < {A.N}"

        let B = A.Clone()
        let m = B.M
        let n = B.N

        // Bidiagonalize, A = U * [B; 0] * V'
        let Ba, Ua, Va = Bidiagonalize B
        let Ub = Ua.Accumulate
        let Vb = Va.Accumulate

        // Skip superflous zero rows, these do not contribute to singular values.
        let Bb = Ba.[..n-1, ..n-1]

        // Find singular values from bidiagonal matrix.
        let Db, Uts, Vs = SvdSteps Bb

        // Sign-adjust if necessary.
        // Singular values are square root of eigenvalues of A'A, 
        // by definition positive.
        let eyes = 
            (Db.D.Data |> Array.map signv) 
            |> Array.toList 
            |> Matrix.CreateDiagonal
        let Vss = eyes::Vs
        let D = Db * eyes

        // Accumulate U' and V.
        let I = Matrix.I n
        let Utsn = List.fold (fun (U: Matrix) (ut: Matrix) -> U * ut) I Uts
        let Vsn = List.fold (fun (V: Matrix) (v: Matrix) -> V * v) I (List.rev Vss)

        // Form full left and right singular vectors.
        let Utsnn = Matrix.I m
        Utsnn.[..n-1, ..n-1] <- Utsn
        let Ut = Utsnn * Ub.T
        let U = Ut.T.[*, ..n-1]
        let V = Vb * Vsn

        D, U, V

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