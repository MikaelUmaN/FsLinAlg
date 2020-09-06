namespace FsLinAlg
open System

/// Contains the fundamental datastructures.
[<AutoOpen>]
module DataStructure =
    
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
                    let vn = v.Norm
                    let dn = this - v |> Vector.norm
                    if isZero vn then isZero dn else dn / vn |> isZero
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
        static member sum (v: Vector) = v.Data |> Array.sum
        static member norm (v: Vector) = v.Norm
        static member map f (v: Vector) = v.Data |> Array.map f |> Vector
        static member length (v: Vector) = v.Data |> Array.length
        static member max (v: Vector) = v.Data |> Array.max
        static member min (v: Vector) = v.Data |> Array.min
        static member findIndex pred (v: Vector) = v.Data |> Array.findIndex pred
        static member findMaxIndex (v: Vector) = v |> Vector.findIndex (fun x -> x = (v |> Vector.max))
        static member findMinIndex (v: Vector) = v |> Vector.findIndex (fun x -> x = (v |> Vector.min))

        interface ICloneable with
            member this.Clone() =
                this.Data |> Array.copy |> Vector :> obj

        member this.Clone() = (this :> ICloneable).Clone() :?> Vector


    /// Row-major matrix of size m (rows) by n (columns).
    and Matrix(data: float[,]) =
        
        let m = Array2D.length1 data 
        let n = Array2D.length2 data
        
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

        member this.IsHessenberg(?u: bool) =
            if this.IsSquare then

                let cu = defaultArg u true
                
                let idxs = 
                    [for i in 2..m-1 -> 
                        [for j in 0..i-2 -> 
                            if cu then (i, j) else (j, i)]]
                    |> List.collect (id)
                
                idxs
                |> List.forall (fun (i, j) -> isZero data.[i, j])
            else
                false

        member this.IsTridiagonal = this.IsHessenberg() && this.IsHessenberg(false)

        override this.Equals other =
            match other with
            | :? Matrix as A ->
                if A.Dimensions <> this.Dimensions then
                    false
                else
                    let fn = A.FrobeniusNorm
                    let dn = this - A |> Matrix.frobeniusNorm
                    if isZero fn then isZero dn else dn / fn |> isZero
            | _ -> false

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

