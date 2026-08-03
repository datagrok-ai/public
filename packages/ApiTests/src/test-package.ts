import type * as _DG from 'datagrok-api/dg';

// Shared, mutable reference to the ApiTests package, populated by whichever test entry
// runs: package-test.ts (browser) or package-test-node.ts (Node). Test files import
// `_package` from here instead of from a specific entry, so they don't transitively pull
// a whole suite graph — and its import-time browser/Dart side effects (e.g. cache.ts's
// top-level grok.functions.register) — into the Node runner.
export let _package: _DG.Package;

export function setTestPackage(p: _DG.Package): void {
  _package = p;
}
