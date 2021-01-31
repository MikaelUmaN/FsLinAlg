namespace FsLinAlg.Test

open System

open Expecto

open FsLinAlg

module SpecializedStructure =

    let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let specializedStructureTests =
        testList "Givens" [
            testProp "acos c = asin s" <| fun (a: NonNanNonZeroFloat, b: NonNanNonZeroFloat) ->
                let (c, s) = givensNumbers (a.Get) (b.Get)
                let thetaC = acos c
                let thetaS = asin s

                // why not work...
                if isZero (abs(thetaC - thetaS)) then
                    Expect.equal thetaC thetaS "Theta angles do not agree"
                else
                    Expect.equal thetaC (Math.PI - thetaS) "Theta angles do not agree (shifted)"
        ]