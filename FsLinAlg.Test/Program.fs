namespace FsLinAlg.Test
open System
open Expecto
open FsCheck

open FsLinAlg.DataStructure
open DataStructure

module Program =

    [<EntryPoint>]
    let main args =
        runTestsInAssemblyWithCLIArgs [Sequenced] args

        // Run specific tests.
        //runTestsWithCLIArgs [] args Factorization.factorizationTests