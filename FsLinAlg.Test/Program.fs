namespace FsLinAlg.Test
open System
open Expecto
open FsCheck

open FsLinAlg.DataStructure
open DataStructure

module Program =

    [<EntryPoint>]
    let main args =
        
        // Normal run
        //runTestsInAssemblyWithCLIArgs [] args

        // For easier debugging.
        //runTestsInAssemblyWithCLIArgs [Sequenced] args

        // Run specific tests.
        //runTestsWithCLIArgs [] args Factorization.factorizationTests

        //runTestsWithCLIArgs [Sequenced] args Factorization.factorizationTests

        //runTestsWithCLIArgs [Sequenced] args SpecializedStructure.specializedStructureTests

        runTestsWithCLIArgs [Sequenced] args Eigenvalue.eigenvalueTests