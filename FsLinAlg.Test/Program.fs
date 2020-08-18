namespace FsLinAlg.Test
open System
open Expecto
open FsCheck

open FsLinAlg.DataStructure
open DataStructure

module Program =

    [<EntryPoint>]
    let main args =
        runTestsInAssemblyWithCLIArgs [] args