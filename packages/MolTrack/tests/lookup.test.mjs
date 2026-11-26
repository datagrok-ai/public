// Run: pnpm --dir packages/MolTrack run test:unit
import assert from 'node:assert/strict';
import {test} from 'node:test';
import {lookupFirstByHash} from '../src/ems/lookup.ts';

// a server that answers 100 rows where a read names no limit, filtered by `=` over a list and by `>`
const server = (rows) => async ({filter: [{property, value}, after], limit = 100}) => rows
  .filter((r) => value.includes(r[property]) && r.molregno > after.value).slice(0, limit);

// under a stereo- or tautomer-insensitive rule several compounds share a hash: 1500 rows over 300 hashes run past
// the server's default 100 and past one page, and the first registered (lowest molregno) is the match
test('every registered structure is recognised, not only the first 100', async () => {
  const rows = Array.from({length: 1500}, (_, i) => ({molregno: i + 1, hash_all_layers: `h${Math.floor(i / 5)}`}));
  const found = await lookupFirstByHash(Array.from({length: 300}, (_, i) => `h${i}`), 'ALL_LAYERS', server(rows));
  assert.equal(found.size, 300, `recognised ${found.size} of 300 registered structures`);
  assert.equal(found.get('h7').molregno, 36, 'the lowest molregno wins');
});
