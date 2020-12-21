// Executing functions

// Calling "RDup" function that belongs to the "RSrcipts" package
// It is implemented in R and returns a duplicated string
// https://github.com/datagrok-ai/public/blob/master/packages/RScripts/scripts/r_dup.R
grok.functions
  .call("RScripts:RDup", {s: "Foo"})
  .then((res) => {
    grok.shell.info(res);
  });

// This is a core function, package name could be omitted
grok.functions
  .call("Sin", {x: 1.5})
  .then((res) => {
    grok.shell.info(res);
  });