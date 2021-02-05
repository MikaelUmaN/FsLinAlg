namespace FsLinAlg
open System
open MathNet.Numerics.Distributions

/// Contains the fundamental datastructures.
[<AutoOpen>]
module DataStructure =
   
    let private randnf = Normal(0., 1.)

    /// Vector that is either a column- or a row-vector.
    /// Provides a transpose operation to change the orientation.
    type Vector(data: float[], ?isColumnVector: bool) =
        let isColumnVector = defaultArg isColumnVector true

        let m = if isColumnVector then data.Length else 1
        let n = if isColumnVector then 1 else data.Length

        member _.M = m
        member _.N = n

        member _.IsColumnVector = isColumnVector

        member _.Data = data

        member _.Length = data.Length

        member this.T = Vector(this.Data, not this.IsColumnVector)

        member this.Norm =
            this |> Vector.map (fun x -> x**2.) |> Vector.sum |> sqrt

        member this.AsMatrix =
            if not this.IsColumnVector then 
                Matrix.FromRowVectors [this] 
            else 
                Matrix.FromColumnVectors [this]

        member this.AsUnitVector =
            let n = this.Norm
            this |> Vector.map (fun x -> x/n)

        /// Get/Set value by row.
        member _.Item
            with get(row) =
                data.[row]
            and set(row) x =
                data.[row] <- x

        /// Column vector slices.
        member this.GetSlice(sr, er) =
            let sr, er =
                defaultArg sr 0,
                defaultArg er (this.Length-1)
            data.[sr..er] |> Vector
            
        member this.SetSlice(sr, er, x: Vector) =
            let sr, er =
                defaultArg sr 0,
                defaultArg er (this.Length-1)
            data.[sr..er] <- x.Data

        override this.Equals other =
            match other with
            | :? Vector as v ->
                if v.Length <> this.Length then
                    false
                else
                    let xy = List.zip (this |> Vector.toList) (v |> Vector.toList)
                    xy
                    |> List.forall (fun (x, y) -> if isZero y then isZero (x-y) else relEq x y)
            | _ -> false

        /// Unary ops
        static member (~-) (v: Vector) = v.Data |> Array.map (~-) |> Vector

        /// Vector-scalar ops
        static member (*) (v: Vector, x: float) =
            v.Data |> Array.map (fun y -> y * x) |> Vector
        static member (*) (x: float, v: Vector) = v * x
        static member (/) (v: Vector, x: float) = v * (1./x)       
        static member (+) (v: Vector, x: float) =
            v.Data |> Array.map (fun y -> y + x) |> Vector
        static member (+) (x: float, v: Vector) = v + x
        static member (-) (v: Vector, x: float) = v + (-x)
        static member (-) (x: float, v: Vector) = x + (-v)

        /// Elementwise vector ops
        static member (+) (x: Vector, y: Vector) =
            Array.zip x.Data y.Data |> Array.map (fun (a, b) -> a + b) |> Vector
        static member (-) (x: Vector, y: Vector) = x + (-y)
        static member elementWiseProduct (x: Vector, y: Vector) =
            Array.zip x.Data y.Data |> Array.map (fun (a, b) -> a * b) |> Vector
        static member (/) (x: Vector, y: Vector) =
            Array.zip x.Data y.Data |> Array.map (fun (a, b) -> a / b) |> Vector

        /// Vector dot-product
        static member ( *+) (x: Vector, y: Vector) =
            Array.zip x.Data y.Data |> Array.sumBy (fun (a, b) -> a * b)

        /// Inner (dot-product) or outer product.
        static member (*) (x: Vector, y: Vector) = x.AsMatrix * y.AsMatrix

        static member zero n = Vector(Array.zeroCreate n)
        static member e i n =
            let d = Array.zeroCreate n
            d.[i] <- 1.
            Vector(d)

        static member exists pred (v: Vector) = v.Data |> Array.exists pred
        static member toList (v: Vector) = v.Data |> Array.toList
        static member toArray (v: Vector) = v.Data
        static member sum (v: Vector) = v.Data |> Array.sum
        static member norm (v: Vector) = v.Norm
        static member map f (v: Vector) = v.Data |> Array.map f |> Vector
        static member length (v: Vector) = v.Data |> Array.length
        static member max (v: Vector) = v.Data |> Array.max
        static member min (v: Vector) = v.Data |> Array.min
        static member findIndex pred (v: Vector) = v.Data |> Array.findIndex pred
        static member findMaxIndex (v: Vector) = v |> Vector.findIndex (fun x -> x = (v |> Vector.max))
        static member findMinIndex (v: Vector) = v |> Vector.findIndex (fun x -> x = (v |> Vector.min))
        static member randn m = randnf.Samples() |> Seq.take m |> Seq.toArray |> Vector

        interface ICloneable with
            member this.Clone() =
                this.Data |> Array.copy |> Vector :> obj

        member this.Clone() = (this :> ICloneable).Clone() :?> Vector


    /// Row-major matrix of size m (rows) by n (columns).
    and Matrix(data: float[,]) =
        let m = Array2D.length1 data 
        let n = Array2D.length2 data
        
        let (|Tall|Wide|Square|) (matrix: Matrix) = 
            if matrix.IsTall then 
                Tall
            elif matrix.IsWide then
                Wide
            else
                Square

        member _.M = Array2D.length1 data
        member _.N = Array2D.length2 data
        member _.Dimensions = (m, n)

        static member FromColumnVectors (columns: list<Vector>) =
            if columns.IsEmpty then
                failwith "Column vector list cannot be empty"
            else
                let lens = columns |> List.map (Vector.length)
                if List.zip lens.[..^1] lens.[1..] |> List.sumBy (fun (i, j) -> i - j) <> 0 then failwith "Column vectors have inconsistent length" else
                    let n = columns.Length
                    let m = columns.[0].Length
                    let data = Array2D.init m n (fun i j -> columns.[j].[i]) 
                    Matrix(data)

        static member FromRowVectors (rows: list<Vector>) =
            if rows.IsEmpty then
                failwith "Row vector list cannot be empty"
            else
                let lens = rows |> List.map (Vector.length)
                if List.zip lens.[..^1] lens.[1..] |> List.sumBy (fun (i, j) -> i - j) <> 0 then failwith "Row vectors have inconsistent length" else
                    let m = rows.Length
                    let n = rows.[0].Length
                    let data = Array2D.init m n (fun i j -> rows.[i].[j]) 
                    Matrix(data)

        static member CreateTridiagonal (subdiagonal: list<float>) (diagonal: list<float>) (superdiagonal: list<float>) =
            let n = diagonal.Length
            let d = Array2D.zeroCreate<float> n n
            for i in 0..n-2 do
                d.[i+1, i] <- subdiagonal.[i]  
                d.[i, i] <- diagonal.[i]  
                d.[i, i+1] <- superdiagonal.[i]
            d.[n-1, n-1] <- diagonal.[n-1]
            Matrix(d)

        static member I(dim) =
            let eye = Array2D.zeroCreate dim dim
            for i in 0..dim-1 do
                eye.[i, i] <- 1.
            eye |> Matrix

        /// Returns the length of the maximum dimension.
        static member maxDim (A: Matrix) = max A.M A.N

        member _.Data = data

        member _.Item
            with get(row, column) =
                data.[row, column]
            and set(row, column) x =
                data.[row, column] <- x

        /// Sub-matrices.
        member _.GetSlice(sr, er, sc, ec) =
            let sr, er, sc, ec =
                defaultArg sr 0,
                defaultArg er (m-1),
                defaultArg sc 0,
                defaultArg ec (n-1)
            data.[sr..er, sc..ec] |> Matrix
            
        member _.SetSlice(sr, er, sc, ec, x: Matrix) =
            let sr, er, sc, ec =
                defaultArg sr 0,
                defaultArg er (m-1),
                defaultArg sc 0,
                defaultArg ec (n-1)
            data.[sr..er, sc..ec] <- x.Data

        /// Row vector slices.
        member _.GetSlice(r, sc, ec) =
            let sc, ec =
                defaultArg sc 0,
                defaultArg ec (n-1)
            data.[r, sc..ec] |> Vector

        member _.SetSlice(r, sc, ec, x: Vector) =
            let sc, ec =
                defaultArg sc 0,
                defaultArg ec (n-1)
            data.[r, sc..ec] <- x.Data

        /// Column vector slices.
        member _.GetSlice(sr, er, c) =
            let sr, er =
                defaultArg sr 0,
                defaultArg er (m-1)
            data.[sr..er, c] |> Vector
        
        member _.SetSlice(sr, er, c, x: Vector) =
            let sr, er =
                defaultArg sr 0,
                defaultArg er (m-1)
            data.[sr..er, c] <- x.Data

        member this.Rows = seq { for i in 0..m-1 -> this.[i, *].T }

        member this.Columns = seq { for i in 0..n-1 -> this.[*, i] }

        /// Transpose of the matrix.
        member this.T = this.Columns |> Seq.toList |> Matrix.FromRowVectors

        /// Sum of the diagonal.
        member this.Trace = this.D |> Vector.sum
            
        /// L(p, q) norm.
        member this.Lpq (p: int) (q: int) =
            let pf, qf = float(p), float(q)
            (Seq.sumBy (fun c -> 
                    c 
                    |> Vector.map (pwn pf)
                    |> Vector.sum
                    |> pwn (qf/pf)) this.Columns)
            |> pwn (1./qf)

        /// Frobenius norm = L(2, 2) norm.
        member this.FrobeniusNorm = this.Lpq 2 2

        /// The diagonal of the matrix.
        member _.D =
            let d = min m n |> int
            [| for i in 0..d-1 -> data.[i, i] |]
            |> Vector

        /// Returns a clone matrix that is upper triangular (diagonal included).
        member this.UpperTriangular =
            let (A: Matrix) = this.Clone()
            A.Rows
            |> Seq.mapi (fun i r -> r.[..i-1] <- Vector.zero i; r)
            |> Seq.toList
            |> Matrix.FromRowVectors

        member _.AsVector =
            if m = 1 then 
                Vector(data.[0, *], isColumnVector=false)
            elif n = 1 then
                Vector(data.[*, 0])
            else
                raise invDim

        member this.IsScalar = this.IsSquare && m = 1

        member this.AsScalar = if this.IsScalar then data.[0, 0] else raise invDim

        member _.IsSquare = m = n

        member this.IsSymmetric = this.IsSquare && this = this.T

        /// If true, then M > N
        member _.IsTall = m > n
            
        /// If true, then N > M
        member _.IsWide = n > m

        /// Generalized check for if this matrix is in hessenberg form.
        /// Allows tall and wide matrices, surplus elements must be zero.
        member this.IsHessenberg(?u: bool) =
            let cu = defaultArg u true
            let d = if this.IsTall then n else m

            let idxs = 
                [for i in 2..d-1 -> 
                    [for j in 0..i-2 -> 
                        if cu then (i, j) else (j, i)]]
                |> List.collect (id)
            
            let allIdxs =
                if cu then 
                    idxs 
                else
                    match this with
                    | Square -> idxs
                    | Tall ->
                        let addIdx =
                            [for i in n..m-1 ->
                                [for j in 0..n-1 ->
                                    (i, j)]]
                            |> List.collect (id)
                        List.concat [idxs; addIdx]
                    | Wide ->
                        let addIdx =
                            [for i in 0..m-1 ->
                                [for j in m..n-1 ->
                                    (i, j)]]
                            |> List.collect (id)
                        List.concat [idxs; addIdx]
            allIdxs
            |> List.forall (fun (i, j) -> isZero data.[i, j])

        /// Generalized check for tridiagonal matrices.
        /// Allows tall and wide matrices, surplus elements must be zero.
        member this.IsTridiagonal = this.IsHessenberg() && this.IsHessenberg(false)

        /// Generalized check for bidiagonal matrices.
        /// Allows tall and wide matrices, surplus elements must be zero.
        /// <param name="u">If true, then the check is for upper-bidiagonality (default=true).</param>
        member this.IsBidiagonal(?u: bool) =
            let cu = defaultArg u true
            let isTri = this.IsTridiagonal
            let (rShift, cShift) = if cu then (0, -1) else (-1, 0)
            let idxs = [for i in 1..n-1 -> (i + rShift, i + cShift)]
            isTri && idxs |> List.forall (fun (i, j) -> isZero data.[i, j])

        member this.IsOrthogonal =
            if m <> n then
                false
            else
                seq { 
                    for i in 0..m-2 do
                        for j in i+1..n-1 do
                            yield (this.[i, *] *+ this.[j, *])
                }
                |> Seq.toList
                |> List.forall isZero

        override this.Equals other =
            match other with
            | :? Matrix as A ->
                if A.Dimensions <> this.Dimensions then
                    false
                else
                    let xy = List.zip (this.Data |> Seq.cast<float> |> Seq.toList) (A.Data |> Seq.cast<float> |> Seq.toList)
                    xy
                    |> List.forall (fun (x, y) -> if isZero y then isZero (x-y) else relEq x y)
            | _ -> false

        override _.ToString() = sprintf "%A" data

        /// Unary ops
        static member (~-) (M: Matrix) = M.Data |> Array2D.map (~-) |> Matrix

        /// Matrix-scalar ops
        static member (*) (M: Matrix, x: float) =
            M.Data |> Array2D.map (fun y -> y * x) |> Matrix
        static member (*) (x: float, M: Matrix) = M * x
        static member (/) (M: Matrix, x: float) =
            M.Data |> Array2D.map (fun y -> y / x) |> Matrix
        static member (/) (x: float, M: Matrix) = M / x

        /// Matrix-matrix ops
        static member (*) (A: Matrix, B: Matrix) =
            if A.N <> B.M then
                failwithf "Inconsistent dimensions for matrix multiplication: %d and %d" A.N B.M
            else
                Array2D.init A.M B.N (fun i j -> A.[i, *] *+ B.[*, j])
                |> Matrix

        static member (+) (A: Matrix, B: Matrix) =
            Array2D.init A.M B.N (fun i j -> A.[i, j] + B.[i, j])
            |> Matrix

        static member (-) (A: Matrix, B: Matrix) =
            -B + A

        /// Matrix-vector ops
        static member (*) (v: Vector, M: Matrix) = v.AsMatrix * M

        static member (*) (M: Matrix, v: Vector) = M * v.AsMatrix

        static member map f (M: Matrix) = M.Data |> Array2D.map f |> Matrix
        static member toVector (M: Matrix) = M.AsVector
        static member toScalar (M: Matrix) = M.AsScalar
        static member lpq p q (M: Matrix) = M.Lpq p q
        static member frobeniusNorm (M: Matrix) = M.FrobeniusNorm

        interface ICloneable with
            member this.Clone() =
                this.Data |> Array2D.copy |> Matrix :> obj

        member this.Clone() = (this :> ICloneable).Clone() :?> Matrix

    type Vec = Vector

    /// Constructor for a vector.
    let V (d: float[]) = Vector(d)

    type Mat = Matrix

