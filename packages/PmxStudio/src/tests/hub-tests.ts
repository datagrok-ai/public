import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';

import {exampleStudy} from './example-study';
import {PLAN_TAG} from '../plan-file';
import {cruxCalls} from '../crux';
import {placeOf} from '../stages';
import {PackageFunctions} from '../package';
import {PICK_WIDTH, ctdOfStudy, demoDated, demoStart, hubView, libraryWidths} from '../hub-view';
import {INDEX_VERSION, LibraryFiles, entryOf, indexPath, libraryReader, myLibrary} from '../library';
import * as DG from 'datagrok-api/dg';
import {bring, fileText, readFile, suggestedStage} from '../open-study';
import {mergedInWorkbook, mergedNames, mergedWarning, sheetCells} from '../xlsx-raw';
import {ANIMALS_007_AND_7, GROUPS_01_AND_1, GROUPS_1_2_3} from './xlsx-fixtures';
import {VOLUME_TAG, parseStudy, studyOf, tableOf} from '../study-file';
import {studySummary} from '../study-reader';
import {gpu} from '../gpu';

const CSV = 'Mouse ID,Arm,Day,Volume (mm3),Body weight\n' +
  'm1,Vehicle,0,100,20\nm1,Vehicle,3,180,20.4\nm1,Vehicle,12,900,21\n' +
  'm2,Cpd-1,0,105,20.1\nm2,Cpd-1,3,140,19.8\nm2,Cpd-1,12,300,19.5\n';

category('Hub', () => {
  test('the demo is dated for today: day 0 a Monday, its plan\'s Mon / Wed / Fri, nothing measured after today',
    async () => {
      // Wednesday 25 November 2026 is day 9 (day 0 Monday the 16th); Friday 27th is day 11 (Wednesday's values the
      // last); Monday 30th is day 7, its day 9 not measured yet
      expect(demoStart(7, new Date(2026, 10, 25)), '2026-11-16');
      expect(demoStart(7, new Date(2026, 10, 27)), '2026-11-16');
      expect(demoStart(7, new Date(2026, 10, 30)), '2026-11-23');
      expect(demoStart(7, new Date(2026, 11, 6)), '2026-11-23');
      // Analyze's, ended at day 21: Monday 30th is day 21, Friday 27th day 25
      expect(demoStart(21, new Date(2026, 10, 30)), '2026-11-09');
      expect(demoStart(21, new Date(2026, 10, 27)), '2026-11-02');
      const plan = '{"start":"2026-11-30","startWeekday":"Monday","measurementDays":"Mon, Wed, Fri"}';
      expect(demoDated(plan, 7, new Date(2026, 10, 27)),
        '{"start":"2026-11-16","startWeekday":"Monday","measurementDays":"Mon, Wed, Fri"}');
      const study = '{"start":"2026-11-30","rows":[{"day":7},{"day":9}]}';
      expect(demoDated(study, 7, new Date(2026, 10, 27)), '{"start":"2026-11-16","rows":[{"day":7},{"day":9}]}');
      expect(demoDated(study, 7, new Date(2026, 10, 30)), '{"start":"2026-11-23","rows":[{"day":7}]}');
    });

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
    expect(r.weight, true);
    expect(studySummary(r.read).join('|'),
      'Read 2 groups: Vehicle (1 animal), Cpd-1 (1 animal), days 0–12.|Volumes as ' +
      'given.|No unit in the column names: assumed g.');
    expect(r.df.columns.names().join(), 'animal,group,day,volume,body_weight');
    expect(r.df.getTag(VOLUME_TAG), JSON.stringify({formula: null, source: 'given'}));
  });

  test('a table without a volume column says so', async () => {
    const r = readFile('study.tsv', 'animal\tgroup\tday\nm1\tA\t0\n');
    expect(r.kind === 'table' && r.read.error, 'Not read: no column is the tumor volume, nor its length and width.');
  });

  test('a table read with an answer: its volumes from length and width, kept in the study it opens as', async () => {
    const csv = 'Animal,Group,Day,Length,Width\nm1,V,0,10,8\nm2,C,0,9,7\n';
    const r = readFile('caliper.csv', csv, {formula: 'pi/6*LW2'});
    expect(r.kind, 'table');
    if (r.kind !== 'table')
      return;
    expect(Math.abs(r.df.col('volume')!.get(0) - Math.PI / 6 * 640) < 1e-9, true);
    const saved = parseStudy(JSON.stringify(studyOf(r.df, null)));
    expect(JSON.stringify(saved.provenance), JSON.stringify({volume: {formula: 'pi/6*LW2', source: 'computed'}}));
    expect(tableOf(saved).getTag(VOLUME_TAG), JSON.stringify({formula: 'pi/6*LW2', source: 'computed'}));
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
    // each with what is wrong, in a sentence
    const why = (name: string, text: string) => {
      const r = readFile(name, text);
      return r.kind === 'bad' ? r.why : '';
    };
    expect(why('notes.docx', 'x').includes('.pmxstudy'), true, why('notes.docx', 'x'));
    expect(why('broken.pmxplan', '{').startsWith('Not a PMX Studio plan'), true, why('broken.pmxplan', '{'));
  });

  test('Your studies: the columns fill the width exactly, with room for Delete / Restore', async () => {
    const sum = (w: {[c: string]: number | undefined}) => Object.values(w).reduce((a: number, b) => a + (b ?? 0), 0);
    expect([sum(libraryWidths(true)), sum(libraryWidths(false))].join(), '100,100');
    expect(libraryWidths(true).action! >= 7, true, 'Restore read "Rest…" at 5');
    // the tick's column a fixed width outside the shares: as one it fell to 9 px at 900 px (QA 0.1.3)
    expect(['pick' in libraryWidths(true), PICK_WIDTH].join(), 'false,32px');
  });

  test('Your studies: what each group receives in one format, the groups apart by "; ", no control', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M2', 'M3']),
      DG.Column.fromStrings('group', ['Vehicle', 'Cpd-A 10 mg/kg', 'Cpd-A 10 mg/kg + Cpd-B 5 mg/kg']),
      DG.Column.fromList('double', 'day', [0, 0, 0]), DG.Column.fromList('double', 'volume', [100, 90, 95])]);
    expect(entryOf(JSON.stringify(studyOf(df, null)), 1)!.treatment, 'Cpd-A 10 mg/kg; Cpd-A 10 mg/kg + Cpd-B 5 mg/kg');
  });

  test('Your studies: a drug listed once however written, a name with no letter no drug (QA 0.1.3 #6, #8)', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M2', 'M3']),
      DG.Column.fromStrings('group', ['Trametinib 1 mg/kg', 'trametinib 3 mg/kg', '1']),
      DG.Column.fromList('double', 'day', [0, 0, 0]), DG.Column.fromList('double', 'volume', [100, 90, 95])]);
    expect(entryOf(JSON.stringify(studyOf(df, null)), 1)!.drugs, 'Trametinib');
  });

  test('the home: three stages, each open', async () => {
    const root = (await hubView()).root;
    const cards = [...root.querySelectorAll<HTMLElement>('.pmx-card')];
    expect(cards.map((c) => c.dataset.stage).join(), 'plan,monitor,analyze');
    // each card: the stage's button and its demo's, apart (no button inside a button)
    expect(cards.every((c) => c.querySelectorAll(':scope > button.pmx-card-main, .pmx-card-foot > button.pmx-demo')
      .length === 2 && !c.querySelector('button button')), true);
    expect(cards.some((c) => c.querySelector<HTMLButtonElement>('.pmx-card-main')!.disabled), false);
    expect(root.querySelectorAll('.pmx-link, .pmx-arc').length, 3);
  });

  test('Your studies from the index: a study is read only when its row is missing or behind', async () => {
    // three plans of no study ID (each its own row: rows of one ID are one study)
    const plan = JSON.parse((await exampleStudy()).getTag(PLAN_TAG)!);
    const files: {[name: string]: {text: string | null, updated: number | null}} = { // (text null: fails to read)
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
        const text = files[typeof f === 'string' ? f.split('/').pop()! : f.fileName].text;
        if (text === null)
          throw new Error('not now');
        return text;
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
    // a file that is no plan is read once and not again until it changes; one that fails to read, every time
    files['e.pmxplan'] = {text: 'not a plan', updated: 5};
    reads = 0;
    expect((await myLibrary(fs)).length, 3);
    await myLibrary(fs);
    expect(reads, 1, 'a file that is no plan was read again');
    files['e.pmxplan'].updated = 6;
    await myLibrary(fs);
    expect(reads, 2, 'a file that is no plan, changed, was not read again');
    files['e.pmxplan'] = {text: null, updated: 7};
    await myLibrary(fs);
    await myLibrary(fs);
    expect(reads, 4, 'a file that failed to read was not read again');
    // one that was no plan and is now a valid one, at a new date, is listed
    files['e.pmxplan'] = {text: 'not a plan', updated: 8};
    await myLibrary(fs);
    files['e.pmxplan'] = {text: JSON.stringify({...plan, name: 'Plan E', study: undefined}), updated: 9};
    expect(rows(await myLibrary(fs)).includes('Plan E|'), true, 'a file become a plan was not listed');
    // a plan of a future version (ours, unreadable here) leaves no row: the next version reads it, at the same date
    files['e.pmxplan'] = {text: JSON.stringify({...plan, version: 99}), updated: 10};
    reads = 0;
    await myLibrary(fs);
    await myLibrary(fs);
    expect(reads, 2, 'a plan of a future version was not read again');
    expect(indexed()['e.pmxplan']?.unreadable === true, false, 'a plan of a future version marked unreadable');
    delete files['e.pmxplan'];
    // a deleted file's row leaves the index
    delete files['c.pmxplan'];
    expect((await myLibrary(fs)).length, 2);
    expect(Object.keys(indexed()).join(), 'a.pmxplan,b.pmxplan');
    // a damaged index, or one of another version, is rebuilt from the files
    for (const damaged of ['{"version": 1, "rows": [', `{"version": ${INDEX_VERSION}, "rows": []}`,
      '{"version": 99, "rows": {}}', JSON.stringify({version: 1, rows: {'a.pmxplan': {name: 'Plan A', updated: 1},
        'b.pmxplan': {name: 'Plan B', updated: 2}}})]) {
      store.set(indexPath(), damaged);
      reads = 0;
      expect(rows(await myLibrary(fs)), rows(first), damaged);
      expect(reads, 2, damaged);
      expect(JSON.parse(store.get(indexPath())!).version, INDEX_VERSION);
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

  // (its workers may be made ready once the home is idle, cruxWarm: no calculation runs on them)
  test('the app opens on its home without the engine: no calculation', async () => {
    // what earlier tests left running settles first (their calls counted before the app opens)
    let seen = -1;
    for (let i = 0; i < 60 && cruxCalls() !== seen; i++) {
      seen = cruxCalls();
      await delay(500);
    }
    expect(cruxCalls() === seen, true, 'earlier engine calls never settled');
    const before = cruxCalls();
    // through the app's own function: what it starts beside the home counts too (with no page of the app open: one
    // open would take it, and no view is answered)
    [...grok.shell.views].filter((v) => placeOf(v)).forEach((v) => v.close());
    const root = (await PackageFunctions.pmxStudioApp())!.root;
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

  const xlsx = (b64: string, name: string) =>
    new File([Uint8Array.from(atob(b64), (c) => c.charCodeAt(0))], name);

  test('names a sheet writes differently but reads as one number are found; others are not', async () => {
    const m = mergedNames([['Animal', 'Group', 'Day'], ['m1', '01', '0'], ['m2', '1', '0'], ['m3', '2', '0']]);
    expect(JSON.stringify(m), JSON.stringify([{role: 'group', texts: ['01', '1'], as: '1'}]));
    expect(mergedNames([['Animal', 'Group'], ['m1', '01'], ['m2', 'B']]).length, 0, 'a text column keeps both');
    expect(mergedNames([['Animal', 'Day'], ['m1', '01'], ['m1', '1']]).length, 0, 'only group and animal columns');
    expect(mergedWarning(m, false).includes('"01" and "1" (read as 1)'), true, mergedWarning(m, false));
    expect(mergedWarning(m, true).includes('01'), false, 'blinded: no group named');
    const many = ['1', '2', '3', '4'].map((as) => ({role: 'animal' as const, texts: [`0${as}`, as], as}));
    expect(mergedWarning(many, false).includes('and 1 more'), true, mergedWarning(many, false));
    expect(mergedWarning([], false), '');
    expect(mergedNames([['Group'], ['+1'], ['1']]).length, 0, 'a sign: no finding rather than a false one');
  });

  test('a sheet with prefixed elements and cells without references reads in order', async () => {
    const ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main';
    const sheet = `<x:worksheet xmlns:x="${ns}"><x:sheetData>` +
      '<x:row><x:c t="s"><x:v>0</x:v></x:c><x:c t="s"><x:v>1</x:v></x:c></x:row>' +
      '<x:row r="2"><x:c><x:v>1</x:v></x:c><x:c r="C2" t="inlineStr"><x:is><x:t>01</x:t></x:is></x:c></x:row>' +
      '<x:row><x:c t="s"><x:v>2</x:v></x:c></x:row></x:sheetData></x:worksheet>';
    expect(JSON.stringify(sheetCells(sheet, ['Animal', 'Group', 'm2'])),
      JSON.stringify([['Animal', 'Group'], ['1', '', '01'], ['m2']]));
  });

  test('a workbook that does not read as a zip gives no finding', async () => {
    expect((await mergedInWorkbook(new Uint8Array([1, 2, 3]), 'Sheet')).length, 0);
    const cut = Uint8Array.from(atob(GROUPS_01_AND_1), (c) => c.charCodeAt(0)).slice(0, 2000);
    expect((await mergedInWorkbook(cut, 'Sheet')).length, 0);
  });

  const bytesOf = (b64: string) => Uint8Array.from(atob(b64), (c) => c.charCodeAt(0));

  test('Excel workbooks as written: 01 and 1 found in shared texts, 007 and 7 inline, 1, 2, 3 nothing', async () => {
    expect(JSON.stringify(await mergedInWorkbook(bytesOf(GROUPS_01_AND_1), 'Sheet')),
      JSON.stringify([{role: 'group', texts: ['01', '1'], as: '1'}]));
    expect(JSON.stringify(await mergedInWorkbook(bytesOf(ANIMALS_007_AND_7), 'Sheet')),
      JSON.stringify([{role: 'animal', texts: ['007', '7'], as: '7'}]));
    expect((await mergedInWorkbook(bytesOf(GROUPS_1_2_3), 'Sheet')).length, 0);
    expect((await mergedInWorkbook(bytesOf(GROUPS_01_AND_1), 'No such sheet')).length, 0);
  });

  test('Excel: groups 01 and 1 are read as one, and the summary says so; groups 1, 2, 3 say nothing', async () => {
    // the reading itself needs PowerPack; without it, the test above still checks the workbooks
    if (!DG.Func.find({package: 'PowerPack', name: 'XlsxFileHandler'}).length)
      return;
    const typo = await fileText(xlsx(GROUPS_01_AND_1, 'typo.xlsx'));
    const r = readFile('typo.xlsx', typo.text, {}, typo.merged);
    // PowerPack reads the column as numbers (pmx-app-notes/xlsx-numeric-group-names.md): one group, said
    expect(r.kind === 'table' ? JSON.stringify(r.read.groups) : r.kind, JSON.stringify([['1', 2]]));
    expect(JSON.stringify(typo.merged), JSON.stringify([{role: 'group', texts: ['01', '1'], as: '1'}]));
    const plain = await fileText(xlsx(GROUPS_1_2_3, 'plain.xlsx'));
    expect(plain.merged.length, 0);
    const animal = await fileText(xlsx(ANIMALS_007_AND_7, 'animal.xlsx'));
    expect(JSON.stringify(animal.merged), JSON.stringify([{role: 'animal', texts: ['007', '7'], as: '7'}]));
  });

  // QA 0.1.3 r5: a study whose kept analysis predates the growth estimate read "p = 0.088" alone, BW nadir "—"
  test('CTD summary: each study\'s rows computed from its data as Analyze does, not an old kept analysis',
    async () => {
      const s = studyOf(await exampleStudy(), null);
      const old = {id: '', name: s.name, compounds: [], animals: '', line: 'MC38', endpointDay: 14, correction: '',
        groups: [{name: 'Vehicle', control: true, given: '', route: '', n: 10, nadir: null, lost: null},
          {name: 'old arm', control: false, given: '', route: '', n: 10, p: 0.088, nadir: null, lost: null}]};
      const ctd = await ctdOfStudy({...s, exposure: {effect: {day: 14, tgi: {}, ctd: old}}});
      const treated = ctd!.groups.filter((g) => !g.control);
      expect(treated.length > 0 && treated.every((g) => !!g.rate && typeof g.p === 'number'), true,
        'each treated group: its estimate with its interval, and its p');
      expect(treated.every((g) => g.nadir !== null), true, 'body weight as weighed');
      expect(treated.some((g) => g.name === 'old arm'), false, 'not the kept rows');
    });
});
