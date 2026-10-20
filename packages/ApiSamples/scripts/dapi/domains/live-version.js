// A live list without polling the rows: `table.version()` answers the table's change token,
// `{seq, at}`. `seq` moves by exactly ONE per write transaction that touched the table's rows —
// a 50-row batch and a 3-op transaction each move it once, a rolled-back or validate-only write
// not at all — and `at` is when it last moved (null until the first write).
//
// So "did anything change?" costs one indexed read per open list instead of a count + max
// (updated_on) aggregate over the whole filtered set. What it does NOT tell you is what changed:
// on a move, re-read the rows (or ask auditLog for the diff).
//
// Needs View on the securing entity: a caller that reaches rows only through per-row grants is
// refused, and has to fall back to counting.

const items = grok.dapi.domains.table('apitests.item');

const df = await items.queryDf({sort: '!created_on', limit: 50});
const view = grok.shell.addTableView(df);
view.name = 'apitests.item (live)';

let held = (await items.version()).seq;
grok.shell.info(`watching apitests.item from seq ${held}`);

const timer = setInterval(async () => {
  const now = await items.version();
  if (now.seq === held)
    return;
  held = now.seq;
  // Something was written — re-read, and rebind the view to the fresh frame.
  const fresh = await items.queryDf({sort: '!created_on', limit: 50});
  view.dataFrame = fresh;
  grok.shell.info(`apitests.item changed (seq ${now.seq}, at ${now.at}) — ${fresh.rowCount} rows`);
}, 3000);

// One subscription, one teardown: the poll must not outlive the view it feeds.
view.subs.push(grok.events.onViewRemoved.subscribe((v) => {
  if (v === view)
    clearInterval(timer);
}));

// Write something, and the next tick picks it up.
const [ins] = await items.insert({sku: `LV-${Date.now()}`, name: 'Live sample'});
setTimeout(() => items.delete(ins.id), 10000);
