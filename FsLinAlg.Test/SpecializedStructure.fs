namespace FsLinAlg.Test

open System

open Expecto

open FsLinAlg

module SpecializedStructure =

    let config = { FsCheckConfig.defaultConfig with maxTest = 10000 }

    [<Tests>]
    let specializedStructureTests =
        testList "Givens" [
            testProp "Angles agree after reflection" <| fun (a: NonNanNonZeroFloat, b: NonNanNonZeroFloat) ->
                let radianToDegree r =
                    180./Math.PI * r

                let degreeToRadian d =
                    d * Math.PI/180.

                let (c, s) = givensNumbers (a.Get) (b.Get)
                let thetaC = acos c
                let thetaS = asin s

                let thetaCc = radianToDegree thetaC
                let thetaSs =
                    if c >= 0. && s >= 0. then
                        // First quadrant.
                        radianToDegree thetaS
                    elif c < 0. && s >= 0. then
                        // c in second quadrant, s in first.
                        radianToDegree (Math.PI - thetaS)
                    elif c < 0. && s < 0. then
                        // c in second, s in third.
                        radianToDegree (Math.PI + thetaS)
                    else
                        // c in first, s in third
                        radianToDegree (-thetaS)
                Expect.floatClose Accuracy.medium thetaCc thetaSs "Theta angles do not agree"
        ]