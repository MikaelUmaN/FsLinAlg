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
                let b = ds * B.[0, 1]**2.
                if isZeroStrict b then B.[1, 1] else B.[1, 1] - b/(abs d + sqrt (d**2. + B.[0, 1]**2.))

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

                    let shouldDeflate =
                        [for j in 0..An.M-2 -> if isZeroStrict An.[j, j+1] then Some(j) else None]
                        |> List.tryFind Option.isSome

                    match shouldDeflate with
                    | Some(oj) ->
                        let j = oj |> Option.get
                        if An.M = 2 then
                            [An.[0, 0]; An.[1, 1]]
                        else
                            let A1 = An.[..j, ..j]
                            let A2 = An.[j+1.., j+1..]
                            let A1e = if A1.IsScalar then [A1.AsScalar] else qr A1 0
                            let A2e = if A2.IsScalar then [A2.AsScalar] else qr A2 0
                            [A1e; A2e] |> List.collect (id)
                    | None -> qr An (k+1)

            qr A 0

