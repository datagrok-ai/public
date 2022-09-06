//tags: Func
//help-url: https://datagrok.ai/help/datagrok/functions/function
// Executing functions

// Calling "RDup" function that belongs to the "RSrcipts" package
// It is implemented in R and returns a duplicated string
// https://github.com/datagrok-ai/public/blob/master/packages/DemoScripts/scripts/r/r_dup.R
grok.functions
  .call("DemoScripts:RDup", {s: "Foo"})
  .then((res) => {
    grok.shell.info(res);
  });

// This is a core function, package name could be omitted
grok.functions
  .call("Sin", {x: 1.5})
  .then((res) => {
    grok.shell.info(res);
  });