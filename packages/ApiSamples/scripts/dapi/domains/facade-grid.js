// UI-FACADE §1.5 — an editable grid over a filtered subset, in TWO lines.
//
// Free: platform grid decoration (the table's own handler customizes it), batch
// editing with amber / red / conflict markers, the save bar, one-transaction save
// with the 409 flows, permission gating down to per-column writability, `defaults`
// prefilling added rows, and service columns that never reach a saved project or
// an export.
//
// The factory is SYNCHRONOUS on the handle — the widget exists at once and its
// rows load inside it; `grid.ready` is the one place that boundary shows.
//
// In a package: `import {domains} from '@datagrok-libraries/domain-ui';`
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

// The example, verbatim (modulo the fixture table and its columns):
const items = await domains.table('apitests.item');
grok.shell.newView('My items',
  [items.grid({query: {filter: 'quantity > 0'}, defaults: {quantity: 1}}).root]);

// The machine surface of a grid is its status, its functions and its dataFrame —
// no per-cell inputs (the frame is a first-class platform object already).
const grid = items.grid();
await grid.ready;
grok.shell.info(`${grid.getWidgetStatus().description}\n` +
  `Actions: ${grid.getFunctions().map((f) => f.name).join(', ')}`);
