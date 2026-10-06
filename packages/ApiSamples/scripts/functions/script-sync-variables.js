// Synchronous GrokScript evaluation. With a variables map the expression
// is evaluated over a fresh context built from that map only; without one
// it is evaluated over the shell context.

grok.shell.info(grok.functions.scriptSync('2 + 3'));

grok.shell.info(grok.functions.scriptSync('a + b', {a: 2, b: 3}));
