import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {exampleStudy} from './example-study';
import {PLAN_TAG} from '../plan-file';
import {cruxCalls} from '../crux';
import {PackageFunctions} from '../package';
import {hubView} from '../hub-view';
import {LibraryFiles, indexPath, libraryReader, myLibrary} from '../library';
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
    const files: {[name: string]: {text: string, updated: number | null}} = {
      'a.pmxplan': {text: JSON.stringify({...plan, name: 'Plan A', study: undefined}), updated: 1},
      'b.pmxplan': {text: JSON.stringify({...plan, name: 'Plan B', study: undefined}), updated: 1}};
    const store = new Map<string, string>();
    let reads = 0;
    let [writable, listable] = [true, {pmxplan: true, pmxstudy: true} as {[ext: string]: boolean}];
    const fs: LibraryFiles = {
      list: async (_folder, _recursive, ext) => {
        if (!listable[ext])
          throw new Error('not now');
        return ext !== 'pmxplan' ? [] : Object.entries(files).map(([fileName, f]) => ({fileName, fullPath: fileName,
          updatedOn: f.updated === null ? null : {valueOf: () => f.updated}}) as any);
      },
      // the index from the store; any study's file, by its info or by its path, counts as a read
      readAsText: async (f) => {
        if (f === indexPath()) {
          if (!store.has(f))
            throw new Error('no file');
          return store.get(f)!;
        }
        reads++;
        return files[typeof f === 'string' ? f.split('/').pop()! : f.fileName].text;
      },
      writeAsText: async (path, text) => {
        if (!writable)
          throw new Error('read only');
        store.set(path, text);
      },
    };
    const rows = (es: {name: string, line: string, end: number | null}[]) =>
      es.map((e) => `${e.name}|${e.line}|${e.end}`).sort().join();
    const indexed = () => JSON.parse(store.get(indexPath())!).rows;
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
    // a new save: the next listing reads it once, dated as its file, and the one after none
    files['c.pmxplan'] = {text: JSON.stringify({...plan, name: 'Plan C', study: undefined}), updated: 3};
    expect((await myLibrary(fs)).length, 3);
    expect(reads, 1);
    expect(indexed()['c.pmxplan'].updated, 3, 'the row is not dated as its file');
    reads = 0;
    expect((await myLibrary(fs)).length, 3);
    expect(reads, 0, 'a listed plan was read again');
    // saved again: its row is older than the file, which is read
    files['c.pmxplan'] = {text: JSON.stringify({...plan, name: 'Plan C2', study: undefined}), updated: 4};
    expect(rows(await myLibrary(fs)).includes('Plan C2|'), true, 'an old row taken for a newer file');
    expect(reads, 1);
    reads = 0;
    // a folder that fails to list fails the reading and drops no row, one or both
    for (const failing of [['pmxstudy'], ['pmxplan', 'pmxstudy']]) {
      failing.forEach((x) => listable[x] = false);
      const failed = await myLibrary(fs).then(() => false, () => true);
      failing.forEach((x) => listable[x] = true);
      expect(failed, true, `read as if listed when ${failing} failed to list`);
      expect(Object.keys(indexed()).length, 3, `rows dropped when ${failing} failed to list`);
    }
    // a file with no date is read every time, and leaves no row
    files['d.pmxplan'] = {text: JSON.stringify({...plan, name: 'Plan D', study: undefined}), updated: null};
    reads = 0;
    await myLibrary(fs);
    await myLibrary(fs);
    expect(reads, 2);
    expect('d.pmxplan' in indexed(), false);
    delete files['d.pmxplan'];
    // a deleted file's row leaves the index
    delete files['c.pmxplan'];
    expect((await myLibrary(fs)).length, 2);
    expect(Object.keys(indexed()).join(), 'a.pmxplan,b.pmxplan');
    // a damaged index, or one of another version, is rebuilt from the files
    for (const damaged of ['{"version": 1, "rows": [', '{"version": 1, "rows": []}', '{"version": 2, "rows": {}}',
      JSON.stringify({version: 1, rows: {'a.pmxplan': {name: 'Plan A', updated: 1}, 'b.pmxplan': {name: 'Plan B',
        updated: 2}}})]) {
      store.set(indexPath(), damaged);
      reads = 0;
      expect(rows(await myLibrary(fs)), rows(first), damaged);
      expect(reads, 2, damaged);
      expect(JSON.parse(store.get(indexPath())!).version, 1);
      // and rebuilt for good: the next load reads nothing
      reads = 0;
      await myLibrary(fs);
      expect(reads, 0, `not rebuilt: ${damaged}`);
    }
    // an index that cannot be written: the list is read as before
    store.clear();
    writable = false;
    reads = 0;
    expect(rows(await myLibrary(fs)), rows(first));
    expect(reads, 2);
  });

  test('the app opens on its home without the engine: no calculation, no workers', async () => {
    // what earlier tests left running settles first (their calls counted before the app opens)
    let seen = -1;
    for (let i = 0; i < 60 && cruxCalls() !== seen; i++) {
      seen = cruxCalls();
      await delay(500);
    }
    expect(cruxCalls() === seen, true, 'earlier engine calls never settled');
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

  test('Your studies are read once when asked twice, and again when asked after the folders were listed', async () => {
    let [listings, shown, warned, listable] = [0, 0, 0, true];
    let open = () => {}; // holds the one file's read: the reading has listed the folders and not ended
    const fs: LibraryFiles = {
      list: async (_folder, _recursive, ext) => {
        listings++;
        if (!listable)
          throw new Error('not now');
        return ext !== 'pmxplan' ? [] : [{fileName: 'a.pmxplan', fullPath: 'a.pmxplan', updatedOn: null} as any];
      },
      readAsText: (f) => typeof f === 'string' ? Promise.reject(new Error('no index')) :
        new Promise((r) => open = () => r('{}')),
      writeAsText: async () => {},
    };
    const read = libraryReader(() => shown++, fs, () => warned++);
    read();
    read(); // as the Overview opens: before the folders are listed
    await delay(20);
    open();
    await delay(20);
    expect(shown, 1);
    expect(listings, 2, 'read twice for two asks before the listing'); // (two folders a reading)
    read();
    await delay(20);
    read(); // a save meanwhile: the folders are listed already
    open();
    await delay(20);
    expect(shown, 2);
    expect(listings, 6, 'an ask after the listing was lost');
    open();
    await delay(20);
    expect(shown, 3);
    expect(listings, 6);
    // a reading that fails to list keeps the list shown, says so, and the next ask reads again
    listable = false;
    read();
    await delay(20);
    expect(shown, 3, 'a failed reading replaced the list');
    expect(warned, 1);
    listable = true;
    read();
    await delay(20);
    open();
    await delay(20);
    expect(shown, 4);
    expect(warned, 1);
  });
});
