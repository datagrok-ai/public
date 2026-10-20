// Caption projection: `captions: ['<ref_column>']` adds ONE service column per ref column,
// `~caption_<column>` (DG.domainCaptionColumn), carrying the TARGET ROW'S DISPLAY NAME — what
// the target declares as its name column, else its business key. A list renders `item_id` as
// 'SKU-42' without a second round trip and without the caller knowing which column of the
// target is its name.
//
// A caption is not an expand. `expand: ['item_id']` brings the target's declared COLUMNS across
// as `item_id.<name>` and costs a join plus a projection of each; a caption brings the name and
// nothing else. Ask for fields with expand, ask for a label with captions.
//
// Never on by default: each caption is one more LEFT JOIN per page.

const events = grok.dapi.domains.table('apitests.item_event');
const items = grok.dapi.domains.table('apitests.item');
const caption = DG.domainCaptionColumn('item_id'); // '~caption_item_id'

const [parent] = await items.insert({sku: `CAP-${Date.now()}`, name: 'Caption sample'});
await events.insert({item_id: parent.id, kind: 'created', amount: 1});

const rows = await events.query({filter: `item_id = "${parent.id}"`, captions: ['item_id']});
// apitests.item declares no isName column, so its display name is its business key (sku).
grok.shell.info(`${rows.length} event(s) of "${rows[0][caption]}"`);

// Independent of `columns`: a caption may be asked for a ref column that is not projected, and
// asking for one never projects the ref column itself.
const [narrow] = await events.query(
  {filter: `item_id = "${parent.id}"`, columns: ['kind'], captions: ['item_id']});
grok.shell.info(`${narrow.kind} on ${narrow[caption]} — no item_id in the row: ${!('item_id' in narrow)}`);

// In a frame the column is a plain nullable string, tagged out of both exports like every
// other `~` service column — toCsv() never carries it.
const df = await events.queryDf({filter: `item_id = "${parent.id}"`, captions: ['item_id']});
grok.shell.info(`frame columns: ${df.columns.names().join(', ')}; csv has the caption: ` +
  `${df.toCsv().includes('~caption_')}`);

// Null means "no name to show", never "the row is missing": a target the caller may not View
// reads null exactly like a target that does not exist. Nothing about a caption tells a caller
// whether a row they cannot see is there — and a caption is never editable and never a field,
// so it is absent from access().fields.
// Refusals are equally silent about the schema: an unknown column, a column that is not a ref,
// a nested 'a.b' name and a repeated name all reject with a DomainFilterError.
try {
  await events.query({filter: `item_id = "${parent.id}"`, captions: ['kind']});
} catch (e) {
  grok.shell.info(`a non-ref caption is refused like an unknown one: ${e.message}`);
}

await items.delete(parent.id); // cascades item_event
