/* Headless regression test for the registered-structure lookup of MolTrack's import (EQ-065).
   Run: pnpm --dir packages/MolTrack run test:unit
   registration.ts is loaded as is; the platform modules it imports are replaced by a stub whose compounds table
   answers like the server: a filter by `=` over a list, rows in no order unless `sort` names one, and a default
   limit of 100 when no limit is named. */
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {test} from 'node:test';
import {pathToFileURL} from 'node:url';

const SERVER_DEFAULT_LIMIT = 100;
const real = pathToFileURL(new URL('../src/ems/', import.meta.url).pathname).href;

const stub = `
export const grok = {}, cond = (property, operator, value) => ({property, operator, value});
export const isEmpty = (v) => v == null || v === '', span = () => null;
export const HASH_COLUMNS = {ALL_LAYERS: 'hash_mol'}, IDENTITY_COLUMNS = [];
export const Chemistry = {}, MolTrackTables = {};
export const moltrackDb = {compounds: {async query(spec) {
  const server = globalThis.__moltrackServer;
  server.calls.push(spec);
  const values = new Set(spec.filter.value);
  const hit = server.rows.filter((r) => values.has(r[spec.filter.property]));
  if (spec.sort === 'molregno') hit.sort((a, b) => a.molregno - b.molregno); else hit.reverse();
  const at = spec.offset ?? 0;
  return hit.slice(at, at + (spec.limit ?? ${SERVER_DEFAULT_LIMIT}));
}}};`;
const hook = `
const STUB = ${JSON.stringify('data:text/javascript,' + encodeURIComponent(stub))};
export async function resolve(specifier, context, next) {
  if (context.parentURL?.endsWith('/registration.ts')) {
    if (specifier === './lookup')
      return next(${JSON.stringify(real)} + 'lookup.ts', context);
    return {url: STUB, shortCircuit: true};
  }
  return next(specifier, context);
}`;
register('data:text/javascript,' + encodeURIComponent(hook));

const {BulkRegistration} = await import(real + 'registration.ts');

// under a stereo- or tautomer-insensitive rule several compounds share a hash: 1500 rows over 300 hashes run past
// the server's default 100 and past one page, and the first registered (lowest molregno) is the match
test('every registered structure is recognised, not only the first 100', async () => {
  const rows = Array.from({length: 1500}, (_, i) => ({molregno: i + 1, hash_mol: `h${Math.floor(i / 5)}`}));
  globalThis.__moltrackServer = {rows, calls: []};
  const found = await BulkRegistration._registered(Array.from({length: 300}, (_, i) => `h${i}`), 'ALL_LAYERS');
  assert.equal(found.size, 300, `recognised ${found.size} of 300 registered structures`);
  assert.equal(found.get('h7').molregno, 36, 'the lowest molregno wins');
});
