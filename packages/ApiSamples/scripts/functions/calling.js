//help-url: https://datagrok.ai/help/datagrok/functions/function
// Executing functions

// Calling "RDup" function that is implemented in R and returns a duplicated string
let res1 = await grok.functions
  .call('ApiSamples:RDup', {s: 'Foo'});
grok.shell.info(res1);
  

// This is a core function, package name could be omitted
let res2 = await grok.functions
  .call('Sin', {x: 1.5});
grok.shell.info(res2);