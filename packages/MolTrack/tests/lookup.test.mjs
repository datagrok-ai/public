/* Headless regression test for the registered-structure lookup of MolTrack's import (EQ-065).
   Run from the repo root: node --experimental-strip-types --test packages/MolTrack/tests/lookup.test.mjs
   registration.ts is loaded as is; the platform modules it imports are replaced by a stub whose compounds table
   answers like the server: a filter by `=` over a list, `sort`, and a default limit of 100 when no limit is named. */
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
  const hit = server.rows.filter((r) => values.has(r[spec.filter.property])).sort((a, b) => a.molregno - b.molregno);
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

async function lookup(rows, hashes) {
  const server = globalThis.__moltrackServer = {rows, calls: []};
  const found = await BulkRegistration._registered(hashes, 'ALL_LAYERS');
  return {found, calls: server.calls};
}

test('every registered structure is recognised, not only the first 100', async () => {
  const rows = Array.from({length: 150}, (_, i) => ({molregno: i + 1, hash_mol: `h${i + 1}`}));
  const {found} = await lookup(rows, rows.map((r) => r.hash_mol));
  assert.equal(found.size, 150, `recognised ${found.size} of 150 registered structures`);
});

test('hashes shared by several compounds can outnumber the hashes asked for', async () => {
  const rows = Array.from({length: 1500}, (_, i) => ({molregno: i + 1, hash_mol: `h${i % 300}`}));
  const hashes = Array.from({length: 300}, (_, i) => `h${i}`);
  const {found} = await lookup(rows, hashes);
  assert.equal(found.size, 300, `recognised ${found.size} of 300 shared hashes`);
  assert.equal(found.get('h7').molregno, 8, 'the lowest molregno wins');
});

test('a long list is read in bounded chunks', async () => {
  const rows = Array.from({length: 1200}, (_, i) => ({molregno: i + 1, hash_mol: `h${i + 1}`}));
  const {found, calls} = await lookup(rows, rows.map((r) => r.hash_mol));
  assert.equal(found.size, 1200, `recognised ${found.size} of 1200 registered structures`);
  assert.ok(calls.every((c) => c.filter.value.length <= 500), 'a filter holds at most 500 values');
});
