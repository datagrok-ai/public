import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  contentHash, evidenceCsv, evidenceFileName, evidenceHash, EvidenceInputs, evidenceMarkdown, evidenceSources,
  evidenceTableRows, MAX_SNAPSHOT_ROWS, nextGeneratedVersion, supersededOnRegenerate, compareEvidence, describeChanges,
  shouldNotify,
} from '../model/evidence-tables';
import {mergeFeeds, programFeedEntries} from '../model/feed';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programContext} from '../portal';
import {buildPortfolioWorkspace, subspaces} from '../workspace';
import {deleteRunTable, saveModelRun, selectRunForCriterion} from '../model-runs';
import {snapshotsOf, takeSnapshot} from '../snapshots';
import {documentCardsOf, documentsOf, folderFiles} from '../documents';
import {checkOf, evidenceChangesOf, evidenceTableOf, generatedOf, generateEvidenceTableFor} from '../evidence-tables';
import {programFeedOf} from '../feed';
import {provenanceOf} from '../provenance';
import {PortalView} from '../views/portal-view';
import {_package} from '../package-test';

const inputs = (over: Partial<EvidenceInputs> = {}): EvidenceInputs => ({
  program: {code: 'P-1', name: 'Program one'},
  criterion: {id: 'c-pd', code: 'pd', name: 'Process description', stage: 'Process locked'},
  run: {runId: 'r-4', title: 'Feed 4', model: 'Pkg:FedBatchTiter', selectedOn: '2026-10-02T10:00:00.000Z', selectedBy: 'Ann',
    rationale: 'higher titer', inputs: [{name: 'feedRate', caption: 'Feed rate', units: '% v/v per day', value: 4},
      {name: 'data', value: 'saved table'}], outputs: [{name: 'titer', caption: 'Final titer', units: 'g/L', value: 3.14159}]},
  snapshot: {id: 's-1', label: 'P-1 2026-10-02 09:00', takenOn: '2026-10-02T09:00:00.000Z', sha256: 'abc', rows: 2,
    criterionState: 'in_progress', values: [{lot: 'L-1', purpose: 'tox', released: true}, {lot: 'L-2', purpose: 'gmp', released: false}]},
  conditions: [{id: 'c-tox', code: 'tox', name: 'Tox lot', state: 'not_met', summary: '0 of 1 rows'}],
  decisions: [{id: 'd-1', decision: 'hold', on: '2026-10-01T08:00:00.000Z', from: 'D4', to: 'D4', rationale: 'wait, "for" tox'}],
  ...over,
});

category('Portal: evidence table rules', () => {
  test('sources: one per run, snapshot, data condition and decision, hashed; the overall hash ignores order', async () => {
    const s = evidenceSources(inputs());
    expect(s.map((x) => x.key).join(), 'run,snapshot,condition:tox,decision:d-1');
    expect(s.every((x) => /^[0-9a-f]{14}$/.test(x.hash)), true);
    expect(contentHash('a'), contentHash('a'));
    expect(contentHash('a') === contentHash('b'), false);
    expect(evidenceHash(s), evidenceHash([...s].reverse()));
    // another run selected changes the run's hash and nothing else
    const other = evidenceSources(inputs({run: {...inputs().run!, runId: 'r-2', title: 'Feed 2'}}));
    expect(s.filter((x, i) => x.hash !== other[i].hash).map((x) => x.key).join(), 'run');
    // a verdict that changes changes the condition's hash
    const verdict = evidenceSources(inputs({conditions: [{...inputs().conditions[0], state: 'met'}]}));
    expect(verdict[2].hash === s[2].hash, false);
    expect(evidenceSources(inputs({run: null, snapshot: null, conditions: [], decisions: []})).length, 0);
  });

  test('rows: the run with units, the snapshot as frozen, the data conditions, the decisions', async () => {
    const rows = evidenceTableRows(inputs());
    const line = (section: string) => rows.filter((r) => r.section === section)
      .map((r) => [r.item, r.value, r.units, r.state].filter((x) => x).join(' ')).join(' | ');
    expect(line('run'), 'selected FedBatchTiter Ann | rationale higher titer');
    expect(line('run input'), 'Feed rate 4 % v/v per day | data saved table');
    expect(line('run output'), 'Final titer 3.142 g/L');
    expect(line('snapshot'), 'pd 2 rows in_progress');
    expect(line('snapshot value'), 'L-1 purpose: tox; released: true | L-2 purpose: gmp; released: false');
    expect(line('data condition'), 'tox 0 of 1 rows not met');
    expect(line('decision'), 'hold D4 → D4 — wait, "for" tox');
    const empty = evidenceTableRows(inputs({run: null, snapshot: null, conditions: [], decisions: []}));
    expect(empty.map((r) => r.item).join(' | '), 'no run selected for the criterion | no snapshot of the program | ' +
      'the Process locked gate has no data condition | no decision cites it yet');
    const many = evidenceTableRows(inputs({snapshot: {...inputs().snapshot!, rows: 40,
      values: Array.from({length: 40}, (_, i) => ({lot: `L-${i}`}))}}));
    expect(many.filter((r) => r.section === 'snapshot value').length, MAX_SNAPSHOT_ROWS + 1);
    expect(many[many.length - 3].item, `${40 - MAX_SNAPSHOT_ROWS} more rows in the snapshot`);
  });

  test('the CSV quotes what needs it; the Markdown has a table per section and says what it was generated from', async () => {
    const rows = evidenceTableRows(inputs());
    const csv = evidenceCsv(rows);
    expect(csv.split('\n')[0], 'section,source,ref,item,value,units,state,on');
    expect(csv.includes(',hold,"D4 → D4 — wait, ""for"" tox",'), true, csv);
    expect(DG.DataFrame.fromCsv(csv).rowCount, rows.length);
    const md = evidenceMarkdown(inputs(), rows, {version: '2', generatedOn: '2026-10-02T10:31:05.000Z', generatedBy: 'Ann', hash: 'h1'});
    expect(md.startsWith('# Evidence table, Process description\n'), true);
    expect(md.includes('Generated 2026-10-02 10:31 UTC by Ann from evidence h1'), true);
    expect(md.includes('## Run input\n\n| Source | Item | Value | Units | State | On |'), true);
    expect(md.includes('| Feed 4 | Feed rate | 4 | % v/v per day |  |  |'), true);
  });

  test('versions count up per criterion; a new one supersedes older generated drafts, not approved ones', async () => {
    expect(nextGeneratedVersion([]), '1');
    expect(nextGeneratedVersion(['1', '3', null, 'x']), '4');
    expect(evidenceFileName('process_description', 'EXB 201', '2', '2026-10-02T10:31:05.123Z', 'md'),
      'evidence-process_description-EXB_201-v2-20261002T103105Z.md');
    const docs = [{id: 'a', status: 'draft', generated: true, criterionId: 'c'}, {id: 'b', status: 'approved', generated: true, criterionId: 'c'},
      {id: 'c', status: 'draft', generated: false, criterionId: 'c'}, {id: 'd', status: 'rejected', generated: true, criterionId: 'c'},
      {id: 'e', status: 'draft', generated: true, criterionId: 'x'}, {id: 'f', status: 'submitted', generated: true, criterionId: 'c'}];
    expect(supersededOnRegenerate(docs, 'c').join(), 'a,d');
  });

  test('what changed since written: another run, a newer snapshot, a verdict, a decision added or gone', async () => {
    const then = evidenceSources(inputs());
    expect(compareEvidence(then, evidenceSources(inputs())).length, 0);
    expect(describeChanges([]), 'Evidence unchanged since written');
    const now = evidenceSources(inputs({
      run: {...inputs().run!, runId: 'r-2', title: 'Feed 2', selectedOn: '2026-10-03T00:00:00.000Z'},
      snapshot: {...inputs().snapshot!, id: 's-2', label: 'P-1 2026-10-03', takenOn: '2026-10-03T00:00:00.000Z', sha256: 'def'},
      conditions: [{...inputs().conditions[0], state: 'met'}],
      decisions: [{id: 'd-2', decision: 'advance', on: '2026-10-03T00:00:00.000Z'}],
    }));
    expect(describeChanges(compareEvidence(then, now)), 'Evidence changed since written: another run selected: Feed 2 (was Feed 4); ' +
      'a newer snapshot: P-1 2026-10-03 (was P-1 2026-10-02 09:00); tox is now met (was not met); ' +
      'decision hold 2026-10-01 no longer cites it; decision added: advance 2026-10-03');
    expect(compareEvidence(then, evidenceSources(inputs({run: null}))).map((c) => `${c.change} ${c.text}`).join(),
      'removed the run selection was removed (was Feed 4)');
    // told once per state of the evidence
    expect(shouldNotify({state: 'stale', hash: 'h2'}, {notified: 'h1'}), true);
    expect(shouldNotify({state: 'stale', hash: 'h2'}, {notified: 'h2'}), false);
    expect(shouldNotify({state: 'current', hash: 'h3'}, null), false);
  });

  test('a stale document is an entry in the activity feed', async () => {
    const e = mergeFeeds(programFeedEntries({program: {id: 'p', code: 'P-1'}, stages: [], decisions: [], evidence: [], runs: [],
      selections: [], documents: [{id: 'd', title: 'Evidence table, PD', version: '1', status: 'draft', linkedBy: 'u',
        linkedOn: '2026-10-02T08:00:00Z', staleOn: '2026-10-02T09:00:00Z', stale: 'another run selected: Feed 2 (was Feed 4)'}],
      snapshots: [], roleChanges: [], syncs: [], news: []}));
    expect(e.map((x) => x.text).join(' | '), 'Evidence changed since written: Evidence table, PD v1: another run selected: Feed 2 ' +
      '(was Feed 4) | Document linked: Evidence table, PD v1');
  });
});

/** Everything the platform tests made, removed in reverse in `after`. */
const undo: (() => Promise<unknown>)[] = [];
let program: ProgramRow;

async function removeFolder(id: string): Promise<void> {
  for (const c of await subspaces(id))
    await removeFolder(c.id);
  const space = await grok.dapi.spaces.find(id).catch(() => null);
  if (space)
    await grok.dapi.spaces.delete(space);
}

async function saveRun(feedRate: number, title: string) {
  const f = DG.Func.find({package: 'ProjectPortal', name: 'FedBatchTiter'})[0];
  const call = await f.prepare({feedRate, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
  const r = await saveModelRun(call, 'EX-BIOC', program.code, {unitOperation: 'EX-PROD', runType: 'EX-SIM', title, link: false});
  for (const id of r.tableIds)
    undo.push(() => deleteRunTable(id));
  undo.push(async () => grok.dapi.functions.calls.delete(await grok.dapi.functions.calls.find(r.runId)));
  undo.push(() => portalDb.modelRuns.delete(r.modelRun.id));
  return r;
}

category('Portal: evidence tables', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    const bioc = await portfolioByCode('EX-BIOC');
    const d4 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D4'}))!;
    const [ins] = await portalDb.programs.insert({portfolio_id: bioc.id, code: `TEST-ET-${Date.now() % 1e8}`,
      name: 'Evidence table test', stage_model_id: bioc.stage_model_id, stage_id: d4.id, source: 'manual', modality: 'mAb'});
    undo.push(() => portalDb.programs.delete(ins.id));
    await buildPortfolioWorkspace('EX-BIOC', false, ins.id);
    program = await portalDb.programs.get(ins.id);
    undo.push(async () => {
      for (const d of await documentsOf(ins.id))
        await portalDb.documents.delete(d.id).catch(() => undefined);
      for (const s of await snapshotsOf(ins.id))
        await portalDb.programSnapshots.delete(s.id).catch(() => undefined);
      if (program.space_id)
        await removeFolder(program.space_id);
    });
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
  });

  test('the process description\'s table: the selected run with units, a snapshot, linked as draft v1 in the folder', async () => {
    expect(!!program.space_id, true, 'the test program has a folder');
    const feed4 = await saveRun(4, 'Feed 4');
    const ctx = await programContext(program.id);
    const pd = ctx.criteria.find((c) => c.code === 'process_description')!;
    await selectRunForCriterion(program.id, pd.id, feed4.modelRun.id, 'higher titer');
    await takeSnapshot(program, {label: 'before the table'});
    const df: DG.DataFrame = await grok.functions.call('ProjectPortal:generateEvidenceTable',
      {portfolio: 'EX-BIOC', program: program.code, criterion: 'process_description'});
    const at = (section: string, item: string) => {
      for (let i = 0; i < df.rowCount; i++) {
        if (df.get('section', i) === section && df.get('item', i) === item)
          return `${df.get('value', i)} ${df.get('units', i)}`.trim();
      }
      return 'none';
    };
    expect(at('run input', 'Feed rate'), '4 % v/v per day');
    expect(at('run', 'selected'), 'FedBatchTiter');
    expect(at('snapshot', 'process_description'), '0 rows');
    expect(at('data condition', 'the Process Lock gate has no data condition'), '');
    expect(at('decision', 'no decision cites it yet'), '');
    const docs = await documentsOf(program.id);
    expect(docs.length, 1);
    const d = docs[0];
    const g = generatedOf(d)!;
    expect(`${d.title} · v${d.doc_version} · ${d.status} · ${d.kind}`, 'Evidence table, Locked process description · v1 · draft · other');
    expect(`${g.criterionCode} ${g.sources.map((s) => s.kind).join('+')} ${g.rows === df.rowCount}`, 'process_description run+snapshot true');
    expect(d.gate_criterion_id, pd.id);
    const files = (await folderFiles(program)).map((f) => f.name);
    expect(files.filter((f) => f.startsWith('evidence-process_description-') && f.includes('-v1-')).length, 2, files.join());
    expect((await grok.dapi.files.readAsText(d.path)).includes('| Feed 4 | Feed rate | 4 | % v/v per day |'), true);
    expect((await evidenceTableOf(d)).rowCount, df.rowCount);
    // provenance runs through the generated document: up to the run and the snapshot it was written from
    const prov = await provenanceOf('document', d.id, 1);
    expect(prov.nodes.filter((n) => n.direction === 'up' && n.relation === 'generated into').map((n) => n.kind).sort().join(),
      'model_run,snapshot');
    // the card shows it, with its table
    const card = (await documentCardsOf(await programContext(program.id))).find((c) => c.code === 'process_description')!;
    expect(card.documents.map((x) => x.title).join(), 'Evidence table, Locked process description');
    const v = new PortalView();
    await v.load('EX-BIOC');
    await v.showProgram('EX-BIOC', program.id);
    const q = (n: string) => v.detailRoot.querySelector(`[data-u2-name="${n}"]`) as HTMLElement | null;
    expect(!!q('generate-evidence-process-description'), true, 'the card offers Generate evidence table');
    expect(!!q('table-document-process-description-0'), true, 'the generated document opens as a table');
  });

  test('generating again: v2, the v1 draft superseded; the IND quality table carries its gate\'s data condition', async () => {
    const r = await generateEvidenceTableFor(program, 'process_description');
    expect(`${r.version} ${r.superseded}`, '2 1');
    const docs = await documentsOf(program.id);
    expect(docs.map((d) => `v${d.doc_version} ${d.status}`).sort().join(), 'v1 superseded,v2 draft');
    const ind = await generateEvidenceTableFor(program, 'ind_quality');
    const rows = ind.df.col('section')!.toList().map((s, i) => `${s}:${ind.df.get('item', i)}:${ind.df.get('state', i)}`);
    expect(rows.includes('data condition:tox_lot:not met'), true, rows.join(' | '));
    expect(rows.includes('run:no run selected for the criterion:'), true, rows.join(' | '));
    let refused = '';
    await generateEvidenceTableFor(program, 'tox_lot').catch((e) => refused = e.message);
    expect(refused, 'criterion tox_lot is a data condition; evidence tables are generated for document criteria');
  });

  test('another run selected: the table is stale, its author told once, in the feed; regenerating makes it current', async () => {
    const feed2 = await saveRun(2, 'Feed 2');
    const ctx = await programContext(program.id);
    const pd = ctx.criteria.find((c) => c.code === 'process_description')!;
    const v2 = (await documentsOf(program.id)).find((d) => generatedOf(d)?.criterion === pd.id && d.status === 'draft')!;
    expect(v2.doc_version, '2');
    expect((await evidenceChangesOf(ctx)).get(v2.id)?.changes.length, 0, 'current before anything changes');
    await selectRunForCriterion(program.id, pd.id, feed2.modelRun.id, 'lower feed holds the titer');
    const found = (await evidenceChangesOf(await programContext(program.id))).get(v2.id)!;
    expect(describeChanges(found.changes), 'Evidence changed since written: another run selected: Feed 2 (was Feed 4)');
    const report: string = await grok.functions.call('ProjectPortal:checkEvidenceChanges', {code: 'EX-BIOC'});
    expect(report.includes(`${program.code}: Evidence table, Locked process description v2: Evidence changed since written`), true, report);
    expect(/told: \S+: notified/.test(report), true, report);
    const check = checkOf(await portalDb.documents.get(v2.id))!;
    expect(`${check.state} ${check.notified === check.hash}`, 'stale true');
    // told once: the same state of the evidence is not told again
    const again: string = await grok.functions.call('ProjectPortal:checkEvidenceChanges', {code: 'EX-BIOC'});
    expect(again.includes('told:'), false, again);
    const bioc = await portfolioByCode('EX-BIOC');
    const feed = await programFeedOf(bioc, program);
    expect(feed.some((x) => x.text.startsWith('Evidence changed since written: Evidence table, Locked process description v2')), true);
    // the card says so and offers a new version
    const v = new PortalView();
    await v.load('EX-BIOC');
    await v.showProgram('EX-BIOC', program.id);
    const q = (n: string) => v.detailRoot.querySelector(`[data-u2-name="${n}"]`) as HTMLElement | null;
    const card = q('document-card-process-description')?.textContent ?? '';
    expect(card.includes('stale') && card.includes('another run selected: Feed 2 (was Feed 4)'), true, card);
    expect(!!q('regenerate-document-process-description-0'), true);
    const r = await generateEvidenceTableFor(program, pd.id);
    expect(`${r.version} ${r.superseded}`, '3 1');
    const after = await evidenceChangesOf(await programContext(program.id));
    expect(after.get(r.document.id)?.changes.length, 0);
    expect(after.has(v2.id), false, 'a superseded version is not checked');
    expect((await grok.functions.call('ProjectPortal:checkEvidenceChanges', {code: 'EX-BIOC'}) as string)
      .includes(`${program.code}: Evidence table, Locked process description`), false);
  });
});
