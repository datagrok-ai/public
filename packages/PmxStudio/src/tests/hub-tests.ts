import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {exampleStudy} from './example-study';
import {PLAN_TAG} from '../plan-file';
import {cruxCalls} from '../crux';
import {PackageFunctions} from '../package';
import {hubView} from '../hub-view';
import {LibraryFiles, indexPath, indexSaved, myLibrary} from '../library';
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

  test('Your studies from the index: a study is read only when its row is missing or behind', async () => {
    // three plans of no study ID (each its own row: rows of one ID are one study)
    const plan = JSON.parse((await exampleStudy()).getTag(PLAN_TAG)!);
    const files: {[name: string]: {text: string, updated: number}} = {
      'a.pmxplan': {text: JSON.stringify({...plan, name: 'Plan A', study: undefined}), updated: 1},
      'b.pmxplan': {text: JSON.stringify({...plan, name: 'Plan B', study: undefined}), updated: 1}};
    const store = new Map<string, string>();
    let reads = 0;
    const fs: LibraryFiles = {
      list: async (_folder, _recursive, ext) => ext !== 'pmxplan' ? [] : Object.entries(files)
        .map(([fileName, f]) => ({fileName, fullPath: fileName, updatedOn: {valueOf: () => f.updated}}) as any),
      readAsText: async (f) => {
        if (typeof f === 'string') {
          if (!store.has(f))
            throw new Error('no file');
          return store.get(f)!;
        }
        reads++;
        return files[f.fileName].text;
      },
      writeAsText: async (path, text) => void store.set(path, text),
    };
    const rows = (es: {name: string, line: string, end: number | null}[]) =>
      es.map((e) => `${e.name}|${e.line}|${e.end}`).sort().join();
    const first = await myLibrary(fs);
    expect(reads, 2);
    expect(rows(first), `Plan A|MC38|${plan.predicted.studyEnd},Plan B|MC38|${plan.predicted.studyEnd}`);
    reads = 0;
    expect(rows(await myLibrary(fs)), rows(first), 'the index gives other rows than the files');
    expect(reads, 0, 'a study was read with its row in the index');
    files['b.pmxplan'].updated = 2; // changed outside the app
    await myLibrary(fs);
    expect(reads, 1);
    reads = 0;
    // saved by the app: its row written as it is saved, the file not read to list it
    files['c.pmxplan'] = {text: JSON.stringify({...plan, name: 'Plan C', study: undefined}), updated: 3};
    await indexSaved('home:Home/PMX Studio/plans/c.pmxplan', files['c.pmxplan'].text, fs);
    expect((await myLibrary(fs)).length, 3);
    expect(reads, 0, 'a plan saved by the app was read to list it');
    // a damaged index is rebuilt from the files
    store.set(indexPath(), '{"version": 1, "rows": [');
    expect(rows(await myLibrary(fs)).split(',').length, 3);
    expect(reads, 3);
    expect(JSON.parse(store.get(indexPath())!).version, 1);
  });

  test('the app opens on its home without the engine: no calculation, no workers', async () => {
    // what earlier tests left running settles first (their calls counted before the app opens)
    let seen = -1;
    for (let i = 0; i < 60 && cruxCalls() !== seen; i++) {
      seen = cruxCalls();
      await delay(500);
    }
    const before = cruxCalls();
    // through the app's own function: what it starts beside the home counts too
    const root = (await PackageFunctions.pmxStudioApp()).root;
    const studies = root.querySelector<HTMLElement>('.pmx-studies')!;
    for (let i = 0; i < 100 && studies.hidden; i++) // the user's studies listed
      await delay(50);
    await delay(500);
    expect(studies.hidden, false);
    expect(cruxCalls() - before, 0);
  });
});
