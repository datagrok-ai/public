import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {exampleStudy} from './example-study';
import {PLAN_TAG} from '../plan-file';
import {cruxCalls} from '../crux';
import {hubView} from '../hub-view';
import {bring, readFile, suggestedStage} from '../open-study';
import {gpu} from '../gpu';

const CSV = 'Mouse ID,Arm,Day,Volume (mm3),Body weight\n' +
  'm1,Vehicle,0,100,20\nm1,Vehicle,3,180,20.4\nm1,Vehicle,12,900,21\n' +
  'm2,Cpd-1,0,105,20.1\nm2,Cpd-1,3,140,19.8\nm2,Cpd-1,12,300,19.5\n';

category('Hub', () => {
  test('WebGPU: one device opened once (or none, without an error, where the browser has none)', async () => {
    const d = await gpu();
    expect(await gpu() === d, true);
    expect(d === null || typeof d.createBuffer === 'function', true);
  });

  test('a table of measurements: its columns by name, animals, groups and last day', async () => {
    const r = readFile('study.csv', CSV);
    expect(r.kind, 'table');
    if (r.kind !== 'table')
      return;
    expect(r.animals, 2);
    expect(r.groups, 2);
    expect(r.lastDay, 12);
    expect(r.missing.length, 0);
    expect(r.weight, true);
    expect(r.found.join('|'), 'Animal ← Mouse ID|Group ← Arm|Day ← Day|Volume ← Volume (mm3)');
  });

  test('a table without a volume column says so', async () => {
    const r = readFile('study.tsv', 'animal\tgroup\tday\nm1\tA\t0\n');
    expect(r.kind === 'table' && r.missing.join(), 'Volume');
  });

  test('a saved plan, and the stage its study opens in from the dates alone', async () => {
    const text = (await exampleStudy()).getTag(PLAN_TAG)!;
    const plan = readFile('MC38.pmxplan', text);
    expect(plan.kind, 'plan');
    if (plan.kind !== 'plan')
      return;
    expect(plan.end !== null && plan.end > 5, true);
    const finished = bring({}, [plan, readFile('study.csv', CSV.replace(/,12,/g, `,${plan.end},`))]).brought;
    const running = bring({}, [plan, readFile('study.csv', CSV)]).brought;
    expect(suggestedStage(finished), 'analyze');
    // the CSV ends on day 12: still running when the plan's last day is past day 13
    expect(plan.end! > 13, true, `the example plan ends on day ${plan.end}`);
    expect(suggestedStage(running), 'monitor');
    expect(suggestedStage({plan}), 'plan');
  });

  test('what is not a study file is set aside, by name', async () => {
    const {brought, bad} = bring({}, [readFile('notes.docx', 'x'), readFile('broken.pmxplan', '{'),
      readFile('empty.csv', ''), readFile('cut.pmxstudy', '{"format": "pmx-study", "version": 1, "rows": [')]);
    expect(bad.join(), 'notes.docx,broken.pmxplan,empty.csv,cut.pmxstudy');
    expect(Object.keys(brought).length, 0);
  });

  test('the home: three stages, each open', async () => {
    const root = (await hubView()).root;
    const cards = [...root.querySelectorAll<HTMLButtonElement>('.pmx-card')];
    expect(cards.map((c) => c.dataset.stage).join(), 'plan,monitor,analyze');
    expect(cards.some((c) => c.disabled), false);
    expect(root.querySelectorAll('.pmx-link, .pmx-arc').length, 3);
  });

  test('the home opens without the engine: no calculation, no workers', async () => {
    const before = cruxCalls();
    const root = (await hubView()).root;
    const studies = root.querySelector<HTMLElement>('.pmx-studies')!;
    for (let i = 0; i < 100 && studies.hidden; i++) // the user's studies read
      await delay(50);
    await delay(500);
    expect(studies.hidden, false);
    expect(cruxCalls() - before, 0);
  });
});
