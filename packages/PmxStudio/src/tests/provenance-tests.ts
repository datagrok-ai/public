import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {armsFrame} from '../arms-input';
import {exampleStudy} from './example-study';
import {FROZEN_ATLAS_NOTE, PLAN_ATLAS_NOTE, PLAN_TAG, atlasNote, parsePlan, planOf, provenance, sha256Hex}
  from '../plan-file';
import {PlanFrames} from '../frames';
import {parseStudy, stamped, studyOf} from '../study-file';

const same = (a: unknown, b: unknown, msg?: string) => expect(JSON.stringify(a), JSON.stringify(b), msg);

category('Provenance', () => {
  test('sha256 of a known string', async () => {
    expect(await sha256Hex('abc'), 'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');
  });

  test('a CRLF text, once its line ends are LF, has the LF hash', async () => {
    expect(await sha256Hex('a\r\nb\r\n'.replace(/\r\n/g, '\n')), await sha256Hex('a\nb\n'));
  });

  test('the current provenance is a hash and a version, cached', async () => {
    const p = await provenance();
    expect(/^[0-9a-f]{64}$/.test(p.atlas), true);
    expect(p.app.length > 0, true);
    expect((await provenance()).atlas, p.atlas);
  });

  test('a different atlas hash gives the sentence, the same or none says nothing', async () => {
    const now = await provenance();
    expect((await atlasNote({atlas: 'f'.repeat(64), app: '0'}, PLAN_ATLAS_NOTE)), PLAN_ATLAS_NOTE);
    expect((await atlasNote({atlas: 'f'.repeat(64), app: '0'}, FROZEN_ATLAS_NOTE)), FROZEN_ATLAS_NOTE);
    expect(await atlasNote(now, PLAN_ATLAS_NOTE), null);
    expect(await atlasNote(undefined, PLAN_ATLAS_NOTE), null);
  });

  const study = (extra: object) => JSON.stringify({format: 'pmx-study', version: 1, name: 's', plan: null, looks: '',
    humane: null, rows: [], log: [], ...extra});

  test('a study keeps its provenance and an old one has none', async () => {
    const p = {atlas: 'a'.repeat(64), app: '1.2.3'};
    same(parseStudy(study({provenance: p})).provenance, p);
    expect(parseStudy(study({})).provenance, undefined);
    expect('provenance' in parseStudy(study({provenance: {atlas: 1}})), false);
  });

  test('the study export writes the current provenance', async () => {
    const df = await exampleStudy();
    same((await stamped(studyOf(df, 2000))).provenance, await provenance());
  });

  test('a plan round-trips its provenance; an old or damaged one reads unchanged', async () => {
    const plan = parsePlan((await exampleStudy()).getTag(PLAN_TAG)!);
    const p = {atlas: 'b'.repeat(64), app: '1.2.3'};
    same(parsePlan(JSON.stringify({...plan, provenance: p})).provenance, p);
    const {provenance: _, ...old} = plan as any;
    same(parsePlan(JSON.stringify(old)), old);
    expect('provenance' in parsePlan(JSON.stringify({...plan, provenance: 'x'})), false);
    expect('provenance' in parsePlan(JSON.stringify({...plan, provenance: {atlas: 1}})), false);
  });

  test('planOf writes the current provenance', async () => {
    const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxPlanStudy'})[0].prepare({line: 'MC38', endpoint: 1000,
      power: 80, correction: 'holm', compareSchedules: false,
      arms: armsFrame([{name: 'Cpd-A', therapy: 'small_molecule', metric: 'TGI', effect: 60}])});
    await fc.call();
    const f = Object.fromEntries(['plan', 'armResults', 'powerCurve'].map((o) => [o, fc.outputs[o]])) as
      unknown as PlanFrames;
    same((await planOf('P', fc, f)).provenance, await provenance());
  }, {timeout: 120000});
});
