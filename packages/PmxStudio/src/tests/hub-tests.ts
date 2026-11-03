import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {exampleStudy} from '../analyze';
import {PLAN_TAG} from '../plan-file';
import {MONITOR_READY} from '../ui';
import {hubView} from '../hub-view';
import {bring, readFile, suggestedStage} from '../open-study';

const CSV = 'Mouse ID,Arm,Day,Volume (mm3),Body weight\n' +
  'm1,Vehicle,0,100,20\nm1,Vehicle,3,180,20.4\nm1,Vehicle,12,900,21\n' +
  'm2,Cpd-1,0,105,20.1\nm2,Cpd-1,3,140,19.8\nm2,Cpd-1,12,300,19.5\n';

category('Hub', () => {
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
    expect(suggestedStage(finished, true), 'analyze');
    expect(suggestedStage(running, true), plan.end! - 1 > 12 ? 'monitor' : 'analyze');
    // without Monitor, a running study still opens: in Analyze
    expect(suggestedStage(running, false), 'analyze');
    expect(suggestedStage({plan}, true), 'plan');
  });

  test('what is not a study file is set aside, by name', async () => {
    const {brought, bad} = bring({}, [readFile('notes.docx', 'x'), readFile('broken.pmxplan', '{'),
      readFile('empty.csv', '')]);
    expect(bad.join(), 'notes.docx,broken.pmxplan,empty.csv');
    expect(Object.keys(brought).length, 0);
  });

  test('the home: three stages, Monitor open only once it is built', async () => {
    const root = (await hubView()).root;
    const cards = [...root.querySelectorAll<HTMLButtonElement>('.pmx-card')];
    expect(cards.map((c) => c.dataset.stage).join(), 'plan,monitor,analyze');
    expect(cards[1].disabled, !MONITOR_READY);
    expect(root.querySelectorAll('.pmx-link, .pmx-arc').length, 3);
  });
});
