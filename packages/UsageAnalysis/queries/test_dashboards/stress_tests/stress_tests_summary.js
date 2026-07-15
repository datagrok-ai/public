//language: javascript
//input: dataframe result
//output: dataframe out

// Passthrough: project datasync of plain (post-process-free) queries loses the
// result dataframe client-side; routing through a JS step makes sync reliable.
let out = result;
