// UI-FACADE §1.8 / §1.9 — a composed page: a form and a grid with ONE dirty
// state, one status bar, one Save/Discard and one unsaved-changes gate; then a
// domain grid next to a regular viewer, live on the SAME DataFrame.
//
// A View is a container widget: what its children edit aggregates up to the page,
// their functions merge onto its ribbon, and a DG.Viewer needs no adapter at all
// (it is already a widget — a raw element enters via DG.Widget.fromRoot).
//
// In a package: `import {domains} from '@datagrok-libraries/domain-ui';`
await grok.functions.call('ApiTests:getDT', {rows: 1});
const {domains} = apitests.domainUi;

// §1.8 — master-detail. The doc's example pairs a project form with the issues of
// that project; the fixture's parent/child pair is item / item_event.
const [items, events] = await Promise.all(
  [domains.table('apitests.item'), domains.table('apitests.item_event')]);
const row = (await items.query({limit: 1}))[0];
if (row == null) {
  grok.shell.info('apitests.item has no rows yet — run transaction.js first.');
  return;
}

grok.shell.addView(domains.view([
  items.form({id: row.id}),
  events.grid({query: {filter: `item_id = "${row.id}"`}, defaults: {item_id: row.id}}),
], {name: row.name}));

// §1.9 — a domain grid next to a regular viewer, both on one frame.
const grid = await events.grid({query: {filter: `item_id = "${row.id}"`}}).ready;
const chart = DG.Viewer.barChart(grid.dataFrame, {split: 'kind'});
grok.shell.addView(domains.view([grid, chart], {name: 'Triage'}));
