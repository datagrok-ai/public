// DomainObjectHandler — the per-table handler every domain table gets for free.
// It is reflective (columns, labels, choices and capabilities come from the
// runtime registry) and delegating: whatever you do NOT override renders exactly
// like the platform, so subclassing is never a regression.

// A complete handler: override only what is genuinely custom.
class IssueHandler extends DG.DomainObjectHandler {
  constructor() { super('grit.issue'); }
  renderCard(x) {
    const row = this.rowOf(x);
    return ui.divV([ui.divText(row.semValue, {style: {fontWeight: 'bold'}}), ui.divText(row.values.title ?? '')],
      'd4-gallery-item');
  }
}
// Registration is the standard mechanism — no new API. NB: no undo, the handler
// stays for the session; platform CRUD commands (Edit/Delete/Share) stay too.
DG.ObjectHandler.register(new IssueHandler());

// Everything else works without a single line of custom code:
const handler = new DG.DomainObjectHandler('grit.issue');
const issue = (await grok.dapi.domains.table('grit.issue').query({limit: 1}))[0];
if (issue == null)
  throw new Error('grit.issue has no rows — create one first (/domains/grit/issue).');
const row = await handler.getById(issue.id);                 // DG.DomainRow, one round trip
// ...or, for rows a query already returned (a list rendering a thousand cards):
const local = handler.rowFrom(issue);                        // same DomainRow, no round trip

const props = await handler.getProperties();                 // registry Property metadata
const caps = await handler.capabilities();                   // server-truth capabilities
const perms = await row.permissions();                       // per-row edit/delete/share
const tabs = await handler.getDetailTabs(row);               // FK-inverted child tables
const actions = await handler.getRibbonActions(row);         // only what this user may do

grok.shell.newView('grit.issue handler', [
  ui.h2(handler.getCaption(local)),
  handler.renderCard(row),
  ui.tableFromMap({
    'Columns': props.map((p) => p.name).join(', '),
    'Writable': caps.writableColumns.join(', ') || '(read-only)',
    'Row permissions': Object.keys(perms).filter((k) => perms[k]).join(', ') || '(none)',
    'Detail tables': tabs.map((t) => t.table).join(', ') || '(none)',
    'Actions': actions.map((a) => a.name).join(', '),
  }),
  await handler.renderEditor(row),                           // reflective form, writable columns only
]);
