import assert from 'node:assert/strict';
import {test} from 'node:test';
import {lookupFirstByHash} from '../src/ems/lookup.ts';

const server = (rows) => async ({filter: [{property, value}, after], limit = 100}) => rows
  .filter((r) => value.includes(r[property]) && r.molregno > after.value).slice(0, limit);

test('every registered structure is recognised, not only the first 100', async () => {
  const rows = Array.from({length: 1500}, (_, i) => ({molregno: i + 1, hash_all_layers: `h${Math.floor(i / 5)}`}));
  const found = await lookupFirstByHash(Array.from({length: 300}, (_, i) => `h${i}`), 'ALL_LAYERS', server(rows));
  assert.equal(found.size, 300, `recognised ${found.size} of 300 registered structures`);
  assert.equal(found.get('h7').molregno, 36, 'the lowest molregno wins');
});
