// UI-FACADE §1.1 — the WHOLE CRUD app, in one body line.
//
// In a package it is an app function, and that is all of it:
//
//   //name: Issues
//   //meta.role: app
//   //meta.browsePath: Grit
//   //input: string path {meta.url: true; optional: true}
//   //output: view result
//   export async function issuesApp(): Promise<DG.ViewBase> {
//     return (await domains.table('grit.issue')).app();
//   }
//
// Free: the list page (cards / brief / editable grid, search, New, the row-cap
// banner), a row page per entity (form, FK-inverted child tabs as editable grids
// with the parent FK prefilled, history, permission-gated actions), URL round trip
// and cold deep links, ONE unsaved-changes gate on every navigation path,
// beforeunload — and every page reports itself through getWidgetStatus() /
// getFunctions() with zero extra code.
//
// ApiTests ships exactly this app over its own fixture table: Browse | Apps |
// Domain Items. Here is the same expression, run from a script.
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

// The example, verbatim (modulo the fixture table):
grok.shell.addView((await domains.table('apitests.item')).app());
