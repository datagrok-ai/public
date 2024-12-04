//tags: Func
//help-url: https://datagrok.ai/help/datagrok/functions/function
// Executing functions

// Calling "RDup" function that belongs to the "RSrcipts" package
// It is implemented in R and returns a duplicated string
// https://github.com/datagrok-ai/public/blob/master/packages/Samples/scripts/r/r_dup.R
let res1 = await grok.functions
  .call("Samples:RDup", {s: "Foo"});
grok.shell.info(res1);
  

// This is a core function, package name could be omitted
let res2 = await grok.functions
  .call("Sin", {x: 1.5});
grok.shell.info(res2);