// Demonstrates function result caching.
// The following code uses dynamic function registration, for package functions you can do that:
//     //input: int x
//     //output: int bar
//     //meta.cache: true

//grok.functions.clientCache.clear();
let counts = {};

grok.functions.register({
  signature: 'int bar(int x)',
  run: (x) => {
    grok.shell.info('Evaluating bar(' + x + ')');
    return x * x;
  },
  options: { cache: 'true',  'cache.invalidateOn': '* * * * *'}  // invalidate every minute
});

grok.functions.register({
  signature: 'int foo(int x)',
  run: (x) => { grok.shell.info('Evaluating foo(' + x + ')'); return x; },
  options: { cache: 'true'}
});

// results are cached
grok.functions.call('bar', {x: 7}).then((result) => grok.shell.info(result));
grok.functions.call('bar', {x: 7}).then((result) => grok.shell.info(result));
grok.functions.call('bar', {x: 7}).then((result) => grok.shell.info(result));

grok.functions.call('bar', {x: 8}).then((result) => grok.shell.info(result));
grok.functions.call('bar', {x: 8}).then((result) => grok.shell.info(result));

grok.functions.call('foo', {x: 12}).then((result) => grok.shell.info(result));
grok.functions.call('foo', {x: 12}).then((result) => grok.shell.info(result));

grok.shell.info(counts);