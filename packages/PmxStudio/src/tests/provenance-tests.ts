import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {exampleStudy} from './example-study';
import {PLAN_TAG, atlasNote, parsePlan, provenance, sha256Hex} from '../plan-file';
import {parseStudy} from '../study-file';

category('Provenance', () => {
  test('sha256 of a known string', async () => {
    expect(await sha256Hex('abc'), 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');
  });

  test('the current provenance is a hash and a version, cached', async () => {
    const p = await provenance();
    expect(/^[0-9a-f]{64}$/.test(p.atlas), true);
    expect(p.app.length > 0, true);
    expect((await provenance()).atlas, p.atlas);
  });

  test('a different atlas hash gives the sentence, the same or none says nothing', async () => {
    const now = await provenance();
    expect((await atlasNote({atlas: 'f'.repeat(64), app: '0'}))?.includes('different version of the atlas'), true);
    expect(await atlasNote(now), null);
    expect(await atlasNote(undefined), null);
  });

  const study = (extra: object) => JSON.stringify({format: 'pmx-study', version: 1, name: 's', plan: null, looks: '',
    humane: null, rows: [], log: [], ...extra});

  test('a study keeps its provenance and an old one has none', async () => {
    const p = {atlas: 'a'.repeat(64), app: '1.2.3'};
    expect(parseStudy(study({provenance: p})).provenance, p);
    expect(parseStudy(study({})).provenance, undefined);
    expect('provenance' in parseStudy(study({provenance: {atlas: 1}})), false);
  });

  test('a plan round-trips its provenance; an old or damaged one reads unchanged', async () => {
    const plan = parsePlan((await exampleStudy()).getTag(PLAN_TAG)!);
    const p = {atlas: 'b'.repeat(64), app: '1.2.3'};
    expect(parsePlan(JSON.stringify({...plan, provenance: p})).provenance, p);
    const {provenance: _, ...old} = plan as any;
    expect(JSON.stringify(parsePlan(JSON.stringify(old))), JSON.stringify(old));
    expect('provenance' in parsePlan(JSON.stringify({...plan, provenance: 'x'})), false);
    expect('provenance' in parsePlan(JSON.stringify({...plan, provenance: {atlas: 1}})), false);
  });
});
