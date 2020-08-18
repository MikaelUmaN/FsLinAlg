namespace FsLinAlg

[<AutoOpen>]
module Exceptions =
    exception LinearDependenceException of string

    let linDep = LinearDependenceException("Columns/Rows are linearly dependent")
