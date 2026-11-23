import {test} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const {ImportLinks} = await import('../src/dg/import-links.js');

// a compound table keyed by a whole number the app shows as DG-000000, as MolTrack's Grok ID is
const KEY = {name: 'molregno', propertyType: 'int', type: 'int', format: 'DG-000000', friendlyName: 'Grok ID'};
const link = {
  prop: {name: 'compound_id', semType: 'moltrack.compound'},
  target: {address: 'moltrack.compound', info: {singularName: 'Compound', pluralName: 'Compounds'}, properties: [KEY]},
  fields: [{prop: KEY, path: 'molregno'}],
};
const STORED = [{id: 'c72', molregno: 72}, {id: 'c73', molregno: 73}];

/** The server's `=` over a number key: what it keeps of the stored rows. */
function serve(asked) {
  grok.dapi.domains.table = () => ({query: async ({filter}) => {
    asked.push(filter);
    const wanted = new Set(filter[0].value);
    return STORED.filter((r) => wanted.has(r.molregno));
  }});
}

test('a key typed as its format shows it finds its record, as the bare number does', async () => {
  const asked = [];
  serve(asked);
  const found = await ImportLinks.resolve(link, link.fields, [['DG-000072'], ['dg-000073'], ['72'], ['nonesuch']]);
  assert.deepEqual(found.map((m) => m.id), ['c72', 'c73', 'c72', null],
    'a Grok ID typed as DG-000072 finds the compound stored as 72');
  assert.equal(found[3].error, 'No compound with Grok ID nonesuch', 'a value that parses to nothing is a row error');
  assert.deepEqual(asked[0][0].value.sort((a, b) => a - b), [72, 73], 'the server is asked for the numbers');
});
