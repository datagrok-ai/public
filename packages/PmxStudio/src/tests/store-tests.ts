import * as grok from 'datagrok-api/grok';
import {category, delay, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import * as ui from 'datagrok-api/ui';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import {armedGates, saveGate} from '../native';

import {exampleStudy} from './example-study';
import {PLAN_TAG, SavedPlan, parsePlan, planKey} from '../plan-file';
import {SavedStudy, parseStudy, studyOf} from '../study-file';
import {detailsForm, emptyMeta} from '../study-meta';
import * as DG from 'datagrok-api/dg';
import {entryOf, libraryOf} from '../library';
import {storageOnFrom} from '../storage-setting';
import {ROW_TAG, domainStorageOn, entryText, firstOpen, firstOpenNotice, importFile, importFolders,
  isReadOnly, keptRows, keptRowsLine, newRowNote, rowsOf, saveCopy,
  linkedRow, list, memoryTables, nextCode, open, remove,
  restore, rowOfPlan, rowOfTable, save, saveError, setRowOn, withCode} from '../store';

/** The example's plan and its study, both under one study ID (as Plan and Monitor save them). */
async function example(): Promise<{plan: ReturnType<typeof parsePlan>, study: SavedStudy}> {
  const df = await exampleStudy();
  const plan = parsePlan(df.getTag(PLAN_TAG)!);
  // (the user's own: a file of another's is brought in as a copy)
  plan.study = {...plan.study ?? emptyMeta(), id: 'PMX-2026-0007', project: 'IO', owner: grok.shell.user?.login ?? ''};
  const study = studyOf(df, 1500);
  study.meta = {...plan.study};
  study.plan = plan;
  return {plan, study};
}

/** A large study with text in many scripts: `animals` animals measured on 20 days, a note on every row. */
function bigStudy(study: SavedStudy, animals: number): SavedStudy {
  const rows = Array.from({length: animals * 20}, (_, i) => ({animal: `Ratón-${Math.floor(i / 20)} µ`,
    group: i % 2 ? 'Véhicule — 対照' : 'Cpd-Ω 10 mg/kg', day: i % 20, volume: 100 + (i % 997) / 7, weight: 20.125}));
  return {...study, name: 'Étude α 😀 中文', rows, extra: {comment: rows.map((_, i) => `nota «${i}» ñé "q"\n\t\\`)},
    log: [{when: '2026-11-24T10:00:00.000Z', who: 'demo', what: 'Révélé — 群'}]};
}

category('Store', () => {
  test('a plan saved as a row opens as it was saved: the same text, so the same plan key', async () => {
    const t = memoryTables();
    const {plan} = await example();
    const saved = await save('plan', plan, null, {}, t);
    expect(saved.version, 1);
    const back = await open('plan', saved.id, t);
    expect(JSON.stringify(back.doc) === JSON.stringify(plan), true);
    expect(planKey(back.doc) === planKey(plan), true);
    expect(back.saved.id === saved.id && back.saved.version === 1, true);
    const row = (await t.plan.get(saved.id))!;
    expect([row.code, row.line, row.project, row.end_day === plan.predicted.studyEnd].join('|'),
      `PMX-2026-0007|${plan.predicted.line}|IO|true`);
  });

  test('an update carries the version it was opened at; a save made elsewhere since is refused', async () => {
    const t = memoryTables();
    const {plan} = await example();
    const first = await save('plan', plan, null, {}, t);
    const second = await save('plan', {...plan, name: 'Renamed'}, first, {}, t);
    expect(second.version, 2);
    expect((await open('plan', first.id, t)).doc.name, 'Renamed');
    let code = '';
    await save('plan', plan, first, {}, t).catch((e) => code = e.code);
    expect(code, 'version-conflict');
  });

  test('the Overview lists a plan and the study followed from it once, at the study, newest first', async () => {
    const t = memoryTables();
    const {plan, study} = await example();
    const p = await save('plan', plan, null, {}, t);
    const s = await save('study', study, null, {planId: p.id}, t);
    const other = await save('plan', {...plan, name: 'Other', study: {...plan.study!, id: 'PMX-2026-0008'}}, null, {},
      t);
    const all = await list('', 1000, t);
    expect(all.map((e) => `${e.id}:${e.kind}:${e.stage}:${e.row.id === (e.kind === 'study' ? s.id : other.id)}`)
      .join(), 'PMX-2026-0008:plan:1:true,PMX-2026-0007:study:2:true');
    expect((await t.study.get(s.id))!.plan_id === p.id, true);
    expect(all[1].text === '' && all[1].file === null && all[1].updated > 0, true);
    expect((await list('other', 1000, t)).map((e) => e.name).join(), 'Other');
  });

  test('a blinded study\'s row names no compound until its groups are revealed', async () => {
    const t = memoryTables();
    const {study} = await example();
    const open1 = await save('study', study, null, {}, t);
    const shown = (await t.study.get(open1.id))!;
    expect(shown.drugs!.length > 0 && shown.compounds !== '', true);
    const blinded = {...study, meta: {...study.meta!, blinded: true}};
    const b = await save('study', blinded, null, {}, t);
    const row = (await t.study.get(b.id))!;
    expect([row.compounds, row.drugs!.length, row.treatment].join('|'), '|0|');
    expect(row.document.includes('"blinded":true'), true);
    // nor does a search by a compound find it, in the list or in the trash
    const drug = shown.drugs![0];
    const hits = async (deleted: boolean) => (await list(drug, 1000, t, deleted)).map((e) => e.row.id);
    expect((await hits(false)).includes(open1.id), true, 'the revealed study is found by its drug');
    await remove('study', open1.id, t);
    expect([(await hits(false)).includes(b.id), (await hits(true)).includes(b.id)].join(), 'false,false');
    await remove('study', b.id, t);
    expect((await hits(true)).includes(b.id), false);
  });

  test('study IDs come from the lab-wide counter, never twice', async () => {
    const t = memoryTables();
    const year = new Date().getFullYear();
    expect([await nextCode(t), await nextCode(t)].join(), `PMX-${year}-0001,PMX-${year}-0002`);
  });

  test('a removed plan leaves the list and no longer opens', async () => {
    const t = memoryTables();
    const {plan} = await example();
    const p = await save('plan', plan, null, {}, t);
    await remove('plan', p.id, t);
    expect((await list('', 1000, t)).length, 0);
    await expectExceptionAsync(async () => {
      await open('plan', p.id, t);
    });
    // in the trash until brought back
    expect((await list('', 1000, t, true)).map((e) => e.row.id).join(), p.id);
    await restore('plan', p.id, t);
    expect((await list('', 1000, t)).length, 1);
    expect((await list('', 1000, t, true)).length, 0);
  });

  test('a view still open on a deleted plan does not save into it; once restored it does', async () => {
    const t = memoryTables();
    const {plan} = await example();
    const p = await save('plan', plan, null, {}, t);
    await remove('plan', p.id, t);
    let said = '';
    await save('plan', plan, p, {}, t).catch((e) => said = saveError(e, 'plan'));
    expect(said.includes('was deleted'), true, said);
    await restore('plan', p.id, t);
    expect((await save('plan', plan, p, {}, t)).version, 2);
  });

  test('a document that is not of the table\'s kind is refused before anything is written', async () => {
    const t = memoryTables();
    const {plan} = await example();
    await expectExceptionAsync(async () => {
      await save('study', plan as unknown as SavedStudy, null, {}, t);
    });
    expect((await t.study.query({})).length, 0);
  });

  test('on the server: a plan row saves, opens as saved, updates at its version and goes to the trash', async () => {
    // the package's own schema (databases/pmx), deployed when the package is published; no study ID is taken
    const {plan} = await example();
    // (no study ID: one that a row on the server has already would be saved under a new one)
    const named = {...plan, name: `Store test ${crypto.randomUUID()}`, study: {...plan.study!, id: ''}};
    const saved = await save('plan', named);
    try {
      expect(JSON.stringify((await open('plan', saved.id)).doc) === JSON.stringify(named), true);
      expect((await list(named.name)).some((e) => e.row.id === saved.id), true);
      const next = await save('plan', {...named, name: named.name + ' (2)'}, saved);
      expect(next.version > saved.version, true);
      let code = '';
      await save('plan', named, saved).catch((e) => code = e.code ?? e.constructor?.name);
      expect(code !== '', true, 'a stale version was not refused');
    } finally {
      await remove('plan', saved.id);
    }
    expect((await list(named.name)).length, 0);
    // in the trash (Show deleted), and back
    expect((await list(named.name, 10, undefined, true)).some((e) => e.row.id === saved.id), true);
    await restore('plan', saved.id);
    expect((await list(named.name)).some((e) => e.row.id === saved.id), true);
    await remove('plan', saved.id);
  });

  test('a file brought in twice is one row; another file of a taken ID is a new row; a colleague’s is a copy',
    async () => {
      const t = memoryTables('ana');
      const {plan, study} = await example();
      const text = JSON.stringify(study, null, 2);
      const first = await importFile('study', text, t);
      const again = await importFile('study', text, t);
      expect([first.imported, again.imported, again.saved.id === first.saved.id].join(), 'true,false,true');
      expect((await t.study.query({})).length, 1);
      await importFile('plan', JSON.stringify(plan), t);
      await importFile('plan', JSON.stringify(plan), t);
      expect((await t.plan.query({})).length, 1);
      // another file of a study ID already there (other text): a new row under a new ID, the first untouched
      const before = (await t.study.get(first.saved.id))!.document;
      const second = await importFile('study', JSON.stringify(study), t);
      const row2 = (await t.study.get(second.saved.id))!;
      expect([second.imported, second.renumbered, second.saved.id !== first.saved.id,
        row2.code !== study.meta!.id].join(), 'true,true,true,true');
      expect((await t.study.get(first.saved.id))!.document, before);
      // a colleague's file (another owner): the user's copy, with an ID of its own, marked as Ana's
      const anas = {...study, name: 'Ana study', meta: {...study.meta!, id: 'PMX-2026-0050', owner: 'ana.lopez'}};
      const copy = await importFile('study', JSON.stringify(anas), t, 'ben.ruiz');
      const row3 = (await t.study.get(copy.saved.id))!;
      expect([copy.copyOf, row3.owner, row3.copy_owner, row3.code !== 'PMX-2026-0050'].join('|'),
        'ana.lopez|ben.ruiz|ana.lopez|true');
      expect((await open('study', copy.saved.id, t)).doc.meta?.copyOf, 'ana.lopez');
      const other = JSON.parse(text);
      other.meta = {...other.meta, id: 'PMX-2026-0099'};
      expect((await importFile('study', JSON.stringify(other), t)).imported, true);
      await expectExceptionAsync(async () => {
        await importFile('study', text.slice(0, 100), t);
      });
    });

  test('the model lists by the name given to a model outside the atlas; its line stays the key', async () => {
    const t = memoryTables();
    const {plan} = await example();
    const own = {...plan, history: {name: 'PDX-CRC-117 (vendor)'} as SavedPlan['history'],
      predicted: {...plan.predicted, line: 'new:pdx:ab12cd'}};
    const s = await save('plan', own, null, {}, t);
    const row = (await t.plan.get(s.id))!;
    expect([row.model, row.line].join('|'), 'PDX-CRC-117 (vendor)|new:pdx:ab12cd');
    expect((await list('crc-117', 1000, t)).map((e) => e.model).join(), 'PDX-CRC-117 (vendor)');
    const atlas = await save('plan', plan, null, {}, t);
    expect((await t.plan.get(atlas.id))!.model, plan.predicted.line);
  });

  test('a large study in many scripts comes back byte for byte (memory, then the server)', async () => {
    const {study} = await example();
    const big = bigStudy({...study, meta: {...study.meta!, id: ''}}, 400);
    const text = JSON.stringify(big);
    expect(text.length > 500000, true, `${text.length} characters`);
    const t = memoryTables();
    const m = await save('study', big, null, {}, t);
    expect((await t.study.get(m.id))!.document === text, true);
    // (opened, it is what the file's own reading gives)
    expect(JSON.stringify((await open('study', m.id, t)).doc) === JSON.stringify(parseStudy(text)), true);
    const saved = await save('study', big);
    try {
      expect((await domainRow(saved.id)) === text, true, 'the server changed the document text');
      expect(JSON.stringify((await open('study', saved.id)).doc) === JSON.stringify(parseStudy(text)), true);
    } finally {
      await remove('study', saved.id);
    }
  });
});

const domainRow = async (id: string) => (await grok.dapi.domains.table('pmx.study').get(id))?.document as string;

category('Store: bringing the files in', () => {
  test('a Home with plans and studies: the first open lists them with no click, moves and deletes nothing',
    async () => {
      const t = memoryTables();
      const {plan, study} = await example();
      const other = {...plan, name: 'Another plan', study: {...plan.study!, id: 'PMX-2026-0101'}};
      const files: {[path: string]: string} = {'plans/a.pmxplan': JSON.stringify(plan),
        'plans/b.pmxplan': JSON.stringify(other), 'studies/s.pmxstudy': JSON.stringify(study, null, 2)};
      const before = JSON.stringify(files);
      const info = (path: string) => ({fileName: path.split('/')[1], path} as unknown as DG.FileInfo);
      const called = new Set<string>();
      const api = {
        list: async (folder: string) => Object.keys(files)
          .filter((k) => folder.endsWith(k.split('/')[0] + '/')).map(info),
        readAsText: async (f: DG.FileInfo) => files[(f as unknown as {path: string}).path],
      };
      // every call the import makes is recorded: only listing and reading are expected
      const fake = new Proxy(api, {get: (o, k: string) => {
        called.add(k);
        return (o as any)[k];
      }});
      // the Overview's first opening: what it listed from the files, it lists from the rows at once
      const {entries, report} = await firstOpen(t, fake);
      const fromFiles = libraryOf(Object.values(files).map((x) => entryOf(x, 0))).map((e) => e.name).sort();
      expect(entries.map((e) => e.name).sort().join('|'), fromFiles.join('|'));
      expect(firstOpenNotice(report).info,
        'Your 3 plans and studies are now saved in Datagrok (the files are left as ' +
        'they were).');
      expect(firstOpenNotice(report).warning === null, true);
      // opened again (another session before the setting was kept): nothing doubles
      expect((await firstOpen(t, fake)).entries.length, entries.length);
      expect((await t.plan.query({})).length + (await t.study.query({})).length, 3);
      expect([...called].sort().join(), 'list,readAsText');
      expect(JSON.stringify(files), before, 'a file was changed');
    });

  test('every plan and study file once, damaged ones named, the files only read', async () => {
    const t = memoryTables();
    const {plan, study} = await example();
    const files: {[path: string]: string} = {
      'plans/a.pmxplan': JSON.stringify(plan), 'plans/bad.pmxplan': '{"format": "pmx-plan"',
      'studies/s.pmxstudy': JSON.stringify(study), 'studies/s copy.pmxstudy': JSON.stringify(study, null, 1)};
    const info = (path: string) => ({fileName: path.split('/')[1], path} as unknown as DG.FileInfo);
    const fake = {
      list: async (folder: string) => Object.keys(files)
        .filter((k) => folder.endsWith(k.split('/')[0] + '/')).map(info),
      readAsText: async (f: DG.FileInfo) => files[(f as unknown as {path: string}).path],
    };
    const first = await importFolders(t, fake);
    expect(JSON.stringify(first), JSON.stringify({imported: 3, already: 0, renumbered: ['s copy.pmxstudy'],
      damaged: ['bad.pmxplan']}));
    const again = await importFolders(t, fake);
    expect([again.imported, again.already].join(), '0,3');
    expect((await list('', 1000, t)).length, 2, 'the plan and its study: one study; the second file: another');
  });
});

category('Store: the views\' side', () => {
  test('a save gate lets go when its view closes, however the views were saved and closed before', async () => {
    const before = armedGates.size;
    const asks = () => {
      const e = new Event('beforeunload', {cancelable: true});
      window.dispatchEvent(e);
      return e.defaultPrevented;
    };
    // a study saved as new (Analyze), then a plan saved with its button (Plan): both changed again, then saved or
    // discarded (while changed, closing asks first), both closed
    let changed = true;
    const views = ['study', 'plan'].map((what) => {
      const host = ViewHost.create({name: `PMX gate test ${what}`, content: () => ui.div()});
      grok.shell.addView(host.view);
      saveGate(host, what, async () => {}, () => changed);
      return host.view;
    });
    expect([armedGates.size - before, asks()].join('|'), '2|true');
    changed = false;
    views.forEach((v) => v.close());
    await delay(100);
    expect([armedGates.size - before, asks()].join('|'), '0|false');
    // a gate whose view never opened, or went without its removal being seen: never asks
    const host = ViewHost.create({name: 'PMX gate test unseen', content: () => ui.div()});
    saveGate(host, 'plan', async () => {}, () => true);
    expect(asks(), false);
    host.dispose();
    expect(armedGates.size - before, 0);
  });

  test('a save never replaces another row: a finished file of a running study’s ID is a new row', async () => {
    const t = memoryTables();
    const {study} = await example();
    const running = await save('study', {...study, rows: study.rows.filter((r) => (r.day ?? 0) <= 9)}, null, {}, t);
    const before = (await t.study.get(running.id))!.document;
    // the finished table, opened from a file (no row): its plan's ID is the running study's
    const finished = await save('study', study, null, {}, t);
    expect(finished.id !== running.id && !!finished.code && finished.code !== study.meta!.id, true);
    expect((await t.study.get(running.id))!.document, before, 'the running study was changed');
    expect((await t.study.get(finished.id))!.code, finished.code);
    expect((await list('', 1000, t)).length, 2);
    // the view takes the new details and row: its next save updates that row, under the new ID, not a third
    expect(finished.meta?.id, finished.code);
    const df = DG.DataFrame.create(1);
    setRowOn(df, finished);
    const again = await save('study', {...study, meta: finished.meta}, rowOfTable(df), {}, t);
    expect([again.id === finished.id, again.version, again.code === undefined].join(), 'true,2,true');
    expect([(await t.study.query({})).length, (await t.study.get(finished.id))!.code === finished.code].join(),
      '2,true');
  });

  test('a plan made in the app is the user\'s, whatever its owner field says; from a colleague\'s file, a copy',
    async () => {
      const t = memoryTables();
      const {plan} = await example();
      // a lab member types their PI as the study's owner: still their own new plan, not a copy of the PI's
      const typed = {...plan, study: {...plan.study!, id: 'PMX-2026-0030', owner: 'Dr Smith'}};
      const own = await save('plan', typed, null, {}, t, 'ben.ruiz');
      const row = (await t.plan.get(own.id))!;
      expect([own.copyOf === undefined, row.copy_owner ?? '', row.code, row.owner].join('|'),
        'true||PMX-2026-0030|Dr Smith');
      // the same plan opened from Dr Smith's file (owner as read at load): the user's copy
      const copy = await save('plan', {...typed, name: 'From the file'}, null, {fileOwner: 'Dr Smith'}, t, 'ben.ruiz');
      expect([copy.copyOf, copy.meta?.owner, (await t.plan.get(copy.id))!.copy_owner].join('|'),
        'Dr Smith|ben.ruiz|Dr Smith');
      // a row's next save is never a copy, whatever is passed
      const again = await save('plan', typed, own, {fileOwner: 'Dr Smith'}, t, 'ben.ruiz');
      expect([again.id === own.id, again.copyOf === undefined].join('|'), 'true|true');
    });

  test('a colleague\'s plan or study opened from a file is first saved as the user\'s copy', async () => {
    const t = memoryTables();
    const {plan, study} = await example();
    // (from the Open dialog, a drop or the Overview: no row yet, and someone else's file)
    const anas = {...plan, study: {...plan.study!, id: 'PMX-2026-0010', owner: 'ana.lopez'}};
    const p = await save('plan', anas, null, {fileOwner: 'ana.lopez'}, t, 'ben.ruiz');
    const row = (await t.plan.get(p.id))!;
    expect([p.copyOf, p.meta?.owner, p.meta?.copyOf, row.owner, row.copy_owner, row.code !== 'PMX-2026-0010',
      row.code === p.code].join('|'), 'ana.lopez|ben.ruiz|ana.lopez|ben.ruiz|ana.lopez|true|true');
    expect(newRowNote('plan', p, 'PMX-2026-0010'), `Saved as your copy of ana.lopez's plan, ${p.code}.`);
    // the user's own keeps its ID (free here); saved again at its row, nothing changes hands
    const own = {...study, meta: {...study.meta!, id: 'PMX-2026-0020', owner: 'ben.ruiz'}};
    const s = await save('study', own, null, {fileOwner: 'ben.ruiz'}, t, 'ben.ruiz');
    expect([s.code === undefined, s.copyOf === undefined, (await t.study.get(s.id))!.code].join('|'),
      'true|true|PMX-2026-0020');
  });

  test('a deep link by the study ID: its one row, or every row with it to choose from; by row id as it is',
    async () => {
      const t = memoryTables();
      const {plan} = await example();
      const a = await save('plan', plan, null, {}, t);
      const b = await save('plan', {...plan, name: 'Again'}, a, {}, t);
      expect((await rowsOf('plan', plan.study!.id, t)).map((r) => r.id).join(','), b.id);
      // (a colleague's row shared with the user under the same ID: both offered, with whose they are)
      await t.plan.insert({name: 'Theirs', code: plan.study!.id, owner: 'ana.lopez', document: '{}'});
      const two = await rowsOf('plan', plan.study!.id, t);
      expect([two.length, two.some((r) => r.owner === 'ana.lopez' && r.name === 'Theirs')].join('|'), '2|true');
      expect((await rowsOf('plan', '0b6c1a4e-5d2f-4c8a-9e1b-2f3a4b5c6d7e', t))[0].id,
        '0b6c1a4e-5d2f-4c8a-9e1b-2f3a4b5c6d7e');
      await expectExceptionAsync(async () => {
        await rowsOf('plan', 'PMX-1999-0001', t);
      });
    });

  test('storage off: the Overview says how many rows are kept, and nothing when there are none', async () => {
    const t = memoryTables();
    expect(keptRowsLine(await keptRows(t)) === null, true);
    const {plan, study} = await example();
    await save('plan', plan, null, {}, t);
    await save('study', study, null, {}, t);
    // a plan and its study: one study, as the list shows it
    expect(keptRowsLine(await keptRows(t)), '1 study is kept in Datagrok storage: switch ' +
      '"Save in Datagrok" back on to see or export it.');
    expect(keptRowsLine(1002)!.startsWith('More than 1000 studies are kept'), true, 'past the count');
  });

  test('the blinding note shows with blinding and storage on, not with blinding off', async () => {
    const shown = (blinded: boolean) => [...detailsForm({...emptyMeta(), blinded}, []).root
      .querySelectorAll<HTMLElement>('.pmx-details-hint')].some((e) =>
      e.textContent!.includes('can still see the groups') && e.style.display !== 'none');
    expect([domainStorageOn(), shown(true), shown(false)].join('|'), 'true|true|false');
  });

  test('a study table keeps the row it came from; a plan opened from a row saves back to it', async () => {
    const t = memoryTables();
    const {plan, study} = await example();
    const s = await save('study', study, null, {}, t);
    const df = DG.DataFrame.create(1);
    expect(rowOfTable(df) === null, true);
    setRowOn(df, s);
    expect(JSON.stringify(rowOfTable(df)), JSON.stringify(s));
    df.setTag(ROW_TAG, '{"id": 3}');
    expect(rowOfTable(df) === null, true, 'a damaged tag names no row');
    const p = await save('plan', plan, null, {}, t);
    const opened = (await open('plan', p.id, t)).doc;
    expect(rowOfPlan(opened)?.id === p.id, true);
    // read as the Overview reads an entry: its row's text, and the row remembered
    const [entry] = await list('', 1000, t);
    expect((await entryText(entry, t)).length > 0, true);
    const next = await save('plan', opened, rowOfPlan(opened), {}, t);
    expect(rowOfPlan(opened)?.version, next.version);
  });

  test('a refused save says why in the user\'s words', async () => {
    expect(saveError({code: 'version-conflict'}, 'study').startsWith('This study was saved elsewhere'), true);
    expect(saveError({code: 'forbidden'}, 'plan').includes('shared with you to read'), true);
    expect(saveError(new Error('offline'), 'plan'), 'The plan could not be saved: offline');
  });

  test('domain storage is on unless the setting says false; a deep link names its row; an ID is kept', async () => {
    expect([storageOnFrom({}), storageOnFrom(null), storageOnFrom({domainStorage: true}),
      storageOnFrom({domainStorage: false}), storageOnFrom({domainStorage: 'false'})].join(),
    'true,true,true,false,false');
    // (this server's setting may have been switched off by hand: the link follows it)
    const link = linkedRow('?study=abc');
    expect(domainStorageOn() ? link?.kind === 'study' && link.id === 'abc' : link === null, true);
    const m = await withCode({...emptyMeta(), id: 'PMX-2026-0042'});
    expect(m.id, 'PMX-2026-0042');
  });

  test('a study shared to read: the save is refused as read-only, Save as my copy forks it, marked', async () => {
    const t = memoryTables();
    const {study} = await example();
    const source = await save('study', study, null, {}, t);
    await t.study.update(source.id, {owner: 'Ana'}, {version: source.version});
    const at = {id: source.id, version: source.version + 1};
    const before = (await t.study.get(source.id))!.document;
    // the server refuses the grantee's update (fake: every update forbidden)
    const update = t.study.update;
    t.study.update = async () => {
      throw Object.assign(new Error('Forbidden'), {code: 'forbidden'});
    };
    let refused: unknown = null;
    await save('study', {...study, rows: study.rows.slice(1)}, at, {}, t).catch((e) => refused = e);
    expect(isReadOnly(refused), true);
    const {saved, doc} = await saveCopy('study', {...study, rows: study.rows.slice(1)}, at, t, 'Ben');
    t.study.update = update;
    const copy = (await t.study.get(saved.id))!;
    expect([copy.name, copy.copy_of === source.id, copy.copy_owner, copy.owner].join('|'),
      `${study.name} (copy)|true|Ana|Ben`, 'the copy is Ben’s, of Ana’s study');
    expect(doc.meta!.id !== study.meta!.id && copy.code === doc.meta!.id, true, 'the copy has a study ID of its own');
    expect((await t.study.get(source.id))!.document, before, 'the original is not touched');
    // listed apart: both rows, the copy marked
    const listed = await list('', 1000, t);
    expect(listed.length, 2);
    expect(listed.find((e) => e.row.id === saved.id)?.copyOf?.owner, 'Ana');
    expect(listed.find((e) => e.row.id === source.id)?.copyOf === undefined, true);
    // a second copy of the same study is numbered, and so is a copy of the copy
    const second = await saveCopy('study', study, at, t, 'Ben');
    const third = await saveCopy('study', doc, saved, t, 'Ben');
    expect([second.doc.name, third.doc.name].join('|'), `${study.name} (copy 2)|${study.name} (copy 3)`);
  });

  test('a copy\'s table keeps the copy\'s name when saved again, though its plan has the original\'s', async () => {
    const df = await exampleStudy();
    const plan = JSON.parse(df.getTag(PLAN_TAG)!);
    df.setTag(PLAN_TAG, JSON.stringify({...plan, name: 'MC38 full study'}));
    df.name = 'MC38 full study (copy 2)'; // as Save as my copy names it
    expect(studyOf(df, null).name, 'MC38 full study (copy 2)');
    df.name = 'measurements';
    expect(studyOf(df, null).name, 'MC38 full study', 'any other table still goes by its plan\'s name');
  });

  test('the way back: a study saved as a row exports to a file and comes back the same document', async () => {
    const [t1, t2] = [memoryTables(), memoryTables()];
    const {plan, study} = await example();
    const s = await save('study', study, null, {}, t1);
    const p = await save('plan', plan, null, {}, t1);
    const doc = (await open('study', s.id, t1)).doc;
    const planDoc = (await open('plan', p.id, t1)).doc;
    // as Download (.pmxstudy) writes it, and as Save writes a file (pretty-printed)
    for (const text of [JSON.stringify(doc), JSON.stringify(doc, null, 2)]) {
      const t = memoryTables();
      const back = await importFile('study', text, t);
      expect(JSON.stringify((await open('study', back.saved.id, t)).doc), JSON.stringify(doc));
    }
    const planBack = await importFile('plan', JSON.stringify(planDoc), t2);
    expect(JSON.stringify((await open('plan', planBack.saved.id, t2)).doc), JSON.stringify(planDoc));
  });
});
