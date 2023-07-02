// Demonstrates function result caching.
// The following code uses dynamic function registration, for package functions you can do that:
//     //input: int x
//     //output: int bar
//     //meta.cache: true

let counts = {};

grok.functions.register({
  signature: 'int bar(int x)',
  run: (x) => {
    grok.shell.info('Evaluating bar(' + x + ')');

    // for internal diagnostics
    let key = '_' + x;
    if (!counts[key])
      counts[key] = 0;
    counts[key]++;

    return x * x;
  },
  options: { cache: 'true' }
});

// results are cached
grok.functions.call('bar', {x: 7}).then((result) => grok.shell.info(result));
grok.functions.call('bar', {x: 7}).then((result) => grok.shell.info(result));
grok.functions.call('bar', {x: 7}).then((result) => grok.shell.info(result));

grok.functions.call('bar', {x: 8}).then((result) => grok.shell.info(result));
grok.functions.call('bar', {x: 8}).then((result) => grok.shell.info(result));

grok.shell.info(counts);