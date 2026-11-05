//api: DG.U2.debug, DG.U2.Scope.bind
// Leak check with DG.U2.debug: snapshot, build, dispose, diff. A scope's bound handler stops once
// the scope is disposed.
const debug = DG.U2.debug;
const wasEnabled = debug.enabled;
if (!wasEnabled)
  debug.enable();
const before = debug.snapshot();

const control = new DG.U2.Control();
const onClick = control.scope.bind(() => grok.shell.info('clicked'));
control.scope.listen(document, 'click', onClick);
grok.shell.info(`live after build: ${debug.diff(before).total}`);

control.dispose();
grok.shell.info(`live after dispose: ${debug.diff(before).total}, handler answers ${onClick()}`);
if (!wasEnabled)
  debug.disable();
