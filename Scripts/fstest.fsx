// Just tests for basic F# functionality.

let arr = [| 3; 4; 2; 10 |]

let d = 2

// From end indexing needs parenthesis if arithmetic.
let fromEnd = arr[^(d-1)..]
let fromEnd2 = arr[^(2-1)..]
