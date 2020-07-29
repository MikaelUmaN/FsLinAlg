namespace FsLinAlg
open System

module DataStructure =
    
    /// Column vector.
    type Vector(data: float[]) =
        member _.Length = data.Length

        member _.Data = data

        /// Get/Set value by row.
        member _.Item
            with get(row) =
                data.[row]
            and set(row) x =
                data.[row] <- x

        /// Column vector slices.
        member this.GetSlice
            with get(sr, er) =
                let sr, er =
                    defaultArg sr 0,
                    defaultArg er this.Length-1
                data.[sr..er] |> Vector
            and set(sr, er) (x: Vector) =
                let sr, er =
                    defaultArg sr 0,
                    defaultArg er this.Length-1
                data.[sr..er] <- x.Data

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
        static member (*) (x: Vector, y: Vector) =
            Array.zip x.Data y.Data |> Array.map (fun (a, b) -> a * b) |> Vector
        static member (/) (x: Vector, y: Vector) =
            Array.zip x.Data y.Data |> Array.map (fun (a, b) -> a / b) |> Vector

        /// Vector dot-product
        static member ( *+) (x: Vector, y: Vector) =
            Array.zip x.Data y.Data |> Array.sumBy (fun (a, b) -> a * b)

        static member length (v: Vector) = v.Data |> Array.length
        static member max (v: Vector) = v.Data |> Array.max
        static member min (v: Vector) = v.Data |> Array.min
        static member findIndex pred (v: Vector) = v.Data |> Array.findIndex pred
        static member findMaxIndex (v: Vector) = v |> Vector.findIndex (fun x -> x = (v |> Vector.max))
        static member findMinIndex (v: Vector) = v |> Vector.findIndex (fun x -> x = (v |> Vector.min))

    type Vec = Vector

    /// Constructor for a vector.
    let V (d: float[]) = Vector(d)

    /// Row-major matrix of size m (rows) by n (columns).
    type Matrix(data: float[,]) =
        
        let m = Array2D.length1 data 
        let n = Array2D.length2 data
        
        member _.M = Array2D.length1 data
        member _.N = Array2D.length2 data
        member this.Dimensions = (this.M, this.N)

        static member FromColumnVectors (columns: list<Vec>) =
            if columns.IsEmpty then
                failwith "Column vector list cannot be empty"
            else
                let lens = columns |> List.map (Vector.length)
                if List.zip lens.[..^1] lens.[1..] |> List.sumBy (fun (i, j) -> i - j) <> 0 then failwith "Column vectors have inconsistent length" else
                    let n = columns.Length
                    let m = columns.[0].Length
                    let data = Array2D.init m n (fun i j -> columns.[j].[i]) 
                    Matrix(data)

        static member FromRowVectors (rows: list<Vec>) =
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

        member this.Data = data
        
        member _.Item
            with get(row, column) =
                data.[row, column]
            and set(row, column) x =
                data.[row, column] <- x

        /// Sub-matrices.
        member this.GetSlice(sr, er, sc, ec) =
            let sr, er, sc, ec =
                defaultArg sr 0,
                defaultArg er this.M-1,
                defaultArg sc 0,
                defaultArg ec this.N-1
            data.[sr..er, sc..ec] |> Matrix
            
        member this.SetSlice(sr, er, sc, ec, x: Matrix) =
            let sr, er, sc, ec =
                defaultArg sr 0,
                defaultArg er this.M-1,
                defaultArg sc 0,
                defaultArg ec this.N-1
            data.[sr..er, sc..ec] <- x.Data

        /// Row vector slices.
        member this.GetSlice(r, sc, ec) =
            let sc, ec =
                defaultArg sc 0,
                defaultArg ec this.N-1
            data.[r, sc..ec] |> Vector

        member this.SetSlice(r, sc, ec, x: Vector) =
            let sc, ec =
                defaultArg sc 0,
                defaultArg ec this.N-1
            data.[r, sc..ec] <- x.Data

        /// Column vector slices.
        member this.GetSlice(sr, er, c) =
            let sr, er =
                defaultArg sr 0,
                defaultArg er this.M-1
            data.[sr..er, c] |> Vector
        
        member this.SetSlice(sr, er, c, x: Vector) =
            let sr, er =
                defaultArg sr 0,
                defaultArg er this.M-1
            data.[sr..er, c] <- x.Data

        /// Unary ops
        static member (~-) (M: Matrix) = M.Data |> Array2D.map (~-) |> Matrix

        /// Matrix-scalar ops
        static member (*) (M: Matrix, x: float) =
            M.Data |> Array2D.map (fun y -> y * x) |> Matrix
        static member (*) (x: float, M: Matrix) = M * x

        /// Matrix-vector ops
        //static member (*) (M: Matrix, v: Vector) =
        //static member (*) (x: float, v: Matrix) = v * x

        /// Matrix-matrix ops
        static member (*) (A: Matrix, B: Matrix) =
            if A.N <> B.M then
                failwithf "Inconsistent dimensions for matrix multiplication: %d and %d" A.N B.M
            else
                Array2D.init A.M B.N (fun i j -> A.[i, *] *+ B.[*, j])
                |> Matrix

        interface ICloneable with
            member this.Clone() =
                this.Data |> Array2D.copy |> Matrix :> obj

        member this.Clone() = (this :> ICloneable).Clone() :?> Matrix

    type Mat = Matrix

