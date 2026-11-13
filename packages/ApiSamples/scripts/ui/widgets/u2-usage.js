//api: DG.U2.Usage, DG.U2.U2UsageError
// How u2 answers a misuse it detects. In `report` mode (the default for now) it logs it once per call site and keeps
// it in Usage.reports(); in `throw` mode it throws a U2UsageError whose message names the fix.
const usage = DG.U2.Usage;
const mode = usage.mode;
const session = DG.U2.Scope.root('sample session');
session.dispose();

usage.clear();
session.effect(() => {});   // an effect on a disposed owner does nothing: reported, with the fix
for (const r of usage.reports())
  grok.shell.info(`${r.code} x${r.count}: ${r.message}`);

usage.mode = 'throw';
try {
  session.setTimeout(() => {}, 100);
}
catch (e) {
  grok.shell.warning(`${e instanceof DG.U2.U2UsageError ? e.code : 'unexpected'}: ${e.message}`);
}
finally {
  usage.mode = mode;
}
