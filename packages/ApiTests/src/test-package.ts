import type * as _DG from 'datagrok-api/dg';

// Test files take `_package` from here rather than from package-test.ts, so importing one
// suite doesn't transitively pull in the whole suite graph and its import-time side effects.
export let _package: _DG.Package;

export function setTestPackage(p: _DG.Package): void {
  _package = p;
}
