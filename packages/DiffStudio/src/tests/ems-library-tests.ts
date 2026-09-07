// Tests for the EMS `diffstudio.library_model` domain table and library seeding.
// These are real integration tests that require an EMS-enabled server (dev); on a server
// without the domain schema the `before` probe flips `emsOk` and each test returns early.

import * as grok from 'datagrok-api/grok';

import {category, expect, test, before} from '@datagrok-libraries/test/src/test';

category('EMS Library', () => {
  const lib = () => grok.dapi.domains.table('diffstudio.library_model');
  const SRC = '#name: t\n#equations:\n  dx/dt = -x\n#argument: t\n  t0 = 0\n  t1 = 1\n  h = 0.1\n#inits:\n  x = 1';
  let emsOk = true;

  before(async () => {
    try {
      await lib().query({limit: 1});
    } catch {
      emsOk = false; // no EMS / diffstudio schema on this server — the checks below are skipped
    }
  });

  test('library_model CRUD round-trip', async () => {
    if (!emsOk)
      return;
    const name = `test-${Date.now()}`;
    const [ins] = await lib().insert({name, source: SRC, category: 'test'});
    expect(ins.id != null, true, 'insert returned no id');
    try {
      const got = await lib().get(ins.id);
      expect(got?.name, name);
      expect(got?.source, SRC);
      expect((await lib().query({filter: `name = "${name}"`})).length, 1);
      await lib().update(ins.id, {description: 'updated'}, {version: got!.version});
      expect((await lib().get(ins.id))?.description, 'updated');
    } finally {
      await lib().delete(ins.id);
    }
    expect(await lib().get(ins.id), null, 'row must be gone after delete');
  });

  test('library_model rejects a missing required column', async () => {
    if (!emsOk)
      return;
    let rejected = false;
    try {
      await lib().insert({name: `test-${Date.now()}`} as any); // no `source` (required)
    } catch {
      rejected = true;
    }
    expect(rejected, true, 'insert without the required "source" column must be rejected');
  });

  test('seedLibraryModels populates the curated library', async () => {
    if (!emsOk)
      return;
    const report = await grok.functions.call('DiffStudio:seedLibraryModels');
    expect(typeof report === 'string', true, 'seed did not return a report');
    const names = new Set((await lib().query({})).map((m) => m.name));
    for (const expected of ['Pollution', 'Bioreactor', 'Nimotuzumab'])
      expect(names.has(expected), true, `library missing "${expected}" after seeding`);
  });
});
