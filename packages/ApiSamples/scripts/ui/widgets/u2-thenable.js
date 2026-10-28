//api: DG.U2.Scope.thenable
// Did a handler answer a promise? DG.U2.Scope.thenable answers the promise-like it is given, else undefined, and
// never throws: an object whose `then` cannot be read (a proxy that answers every name) is just a value.
const row = grok.data.demo.demog(10).row(0);
const answers = {'a promise': Promise.resolve(42), 'a thenable': {then: (ok) => ok(1)}, 'a row': row, 'a text': 'done'};
for (const [what, answer] of Object.entries(answers))
  grok.shell.info(`${what}: ${DG.U2.Scope.thenable(answer) === undefined ? 'a value' : 'waited for'}`);
