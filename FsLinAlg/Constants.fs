namespace FsLinAlg

open System

[<AutoOpen>]
module Constants =

    /// Float epsilon.
    /// The smallest representable change from zero.
    let eps = Double.Epsilon

    /// Calculates machine epsilon valid on the executing machine,
    /// representing the upper bound of the relative error due to rounding in floating-point arithmetic.
    let private calcMachEps =
        let rec inner e = if e + 1.0 <> 1.0 then inner (e/2.) else e
        inner 1e-4

    /// Calculated machine epsilon
    let machEps = calcMachEps

    /// Relaxed tolerance threshold.
    /// Tolerance should be based on the condition number of the calculation.
    let tol = 1e-5
    let tolStrict = machEps * 2.**8.

    let isZero x = abs x < tol
    let isZeroStrict x = abs x < tolStrict

    let insignificantComparedTo x y = x <= tol * y
    let insignificantComparedToStrict x y = x <= tolStrict * y

    let relEq x y = if isZero y then isZero x else (x-y) / y |> isZero
    let relEqStrict x y = if isZeroStrict y then isZeroStrict x else (x-y) / y |> isZeroStrict

[<AutoOpen>]
module Utils =

    let pwn p x = x**p

    let sqr x = pwn 2. x