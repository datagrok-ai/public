// DomainQuery — the one serializable representation of what a user is looking at: the
// platform function's parameters, a URL deep link, and a REST spec, all the same object.

const table = 'grit.issue'; // any '<schema>.<table>' registered on this server

const q = new DG.DomainQuery({
  schema: 'grit', table: 'issue',
  filters: ['status = "open"'],   // smart-filter grammar, or a JSON condition node
  orderBy: ['!created_on'],
  limit: 100,
});

// A deep link: one URL key per list element, so a shared link (or a data-synced
// dashboard) can re-parameterize `filters[0]` alone. 'view'/'entity' are reserved for
// the app's own UI state and are ignored on the way back — they never enter the query.
const url = new URLSearchParams({...q.toUrlParams(), view: 'grid'}).toString();
const restored = DG.DomainQuery.fromUrlParams('grit', 'issue',
  Object.fromEntries(new URLSearchParams(url)));
grok.shell.info(`${url}\n${JSON.stringify(restored.toParams())}`);

// Silent read: the same rows over REST — no history, no workspace entry.
const rows = await grok.dapi.domains.table(table).queryDf(q.toSpec());
grok.shell.info(`${rows.rowCount} rows`);

// Reproducible run: through the DomainQuery function, so the frame carries a creation
// script (Source-pane refresh, data sync, URL parameters) and the platform handles the
// result — this is what "Open in Table View" does in the Domain View.
const recorded = await q.run();
grok.shell.info(recorded.tags['.script']);

// The fluent builder exports the same state: each top-level AND conjunct becomes one
// filter element, as a condition node whose value is bound server-side.
const builder = grok.dapi.domains.table(table).query().where('status', '=', 'open').top(10);
grok.shell.info(JSON.stringify(DG.DomainQuery.fromBuilder(builder).toParams()));

// And what a Domain View is showing right now IS a DomainQuery.
const view = DG.DomainView.create({schema: 'grit', table: 'issue'});
grok.shell.addView(view);
await DG.delay(2000);
grok.shell.info(new URLSearchParams(DG.DomainQuery.fromParams(view.query).toUrlParams()).toString());
