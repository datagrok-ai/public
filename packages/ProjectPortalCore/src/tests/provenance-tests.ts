import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  nodeKey, provenanceGraph, ProvenanceRecords, provenanceTree, startKey, summarizeProvenance, traverse,
} from '../model/provenance';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, recordDecision, setCriterion, stagesOf} from '../portal';
import {deleteRunTable, saveModelRun, selectRunForCriterion} from '../model-runs';
import {linkDocumentTo, PACKAGE_FILES} from '../documents';
import {provenanceOf} from '../provenance';
import {openProvenanceDialog} from '../views/provenance-dialog';
import {PortalView} from '../views/portal-view';
import {_package} from '../package-test';

/** One program: a run of model M reads T1 and writes T2, is selected for the process description,
 * which document D answers; the advance cites the run, the document, a snapshot of the data, the
 * criterion and a link; a run cited and deleted since. */
const records = (): ProvenanceRecords => ({
  programs: [{id: 'p', code: 'P-1', name: 'Alpha'}],
  decisions: [{id: 'x', programId: 'p', decision: 'advance', on: '2026-10-02T10:00:00Z', from: 'D4', to: 'D5'},
    {id: 'y', programId: 'p', decision: 'hold', on: '2026-10-03T10:00:00Z'}],
  evidence: [
    {id: 'e1', decisionId: 'x', kind: 'model_run', target: 'call-1', name: 'Feed 4'},
    {id: 'e2', decisionId: 'x', kind: 'file', target: 'X:pd.md', name: 'PD', context: {document: 'd'}},
    {id: 'e3', decisionId: 'x', kind: 'snapshot', target: 's', name: 'Snap'},
    {id: 'e4', decisionId: 'x', kind: 'criterion', target: 'c-pd', name: 'Process description'},
    {id: 'e5', decisionId: 'x', kind: 'url', target: 'https://example.org/r', name: 'Review'},
    {id: 'e6', decisionId: 'y', kind: 'model_run', target: 'call-gone', name: 'Old run'},
  ],
  runs: [{id: 'r', programId: 'p', runId: 'call-1', title: 'Feed 4', model: {code: 'Pkg:M', name: 'Titer model'},
    inputs: [{id: 't1', name: 'Feed plan'}], outputs: [{id: 't2', name: 'Profile'}]}],
  selections: [{programId: 'p', criterionId: 'c-pd', runId: 'r'}],
  documents: [{id: 'd', programId: 'p', title: 'PD', version: '1.0', status: 'approved', criterionId: 'c-pd', path: 'X:pd.md'}],
  snapshots: [{id: 's', programId: 'p', label: 'Snap', takenOn: '2026-10-02T09:00:00Z'}],
  criteria: [{id: 'c-pd', code: 'process_description', name: 'Process description', stage: 'D5'},
    {id: 'c-tox', code: 'tox', name: 'Tox lot', stage: 'D3'}],
  datasets: [{programId: 'p', name: 'Lots'}],
  compounds: [{id: 'k', code: 'EXC-1', programId: 'p', role: 'lead'}],
});

category('Portal: provenance rules', () => {
  test('the graph: data, runs, documents, criteria, decisions and the program', async () => {
    const g = provenanceGraph(records());
    const rel = (from: string, to: string) => g.edges.find((e) => e.from === from && e.to === to)?.relation ?? 'none';
    expect(rel('compound:k', 'dataset:p:data'), 'is in');
    expect(rel('dataset:p:data', 'snapshot:s'), 'frozen in');
    expect(rel('table:t1', 'model_run:r'), 'input of');
    expect(rel('model_run:r', 'table:t2'), 'output');
    expect(rel('model:Pkg:M', 'model_run:r'), 'model of');
    expect(rel('model_run:r', 'document:d'), 'selected for its criterion');
    expect(rel('model_run:r', 'criterion:c-pd'), 'none', 'a document answers it, so the run feeds the document');
    expect(rel('document:d', 'criterion:c-pd'), 'answers');
    expect(rel('document:d', 'decision:x'), 'evidence for');
    expect(rel('url:https://example.org/r', 'decision:x'), 'evidence for');
    expect(rel('decision:x', 'program:p'), 'moved');
    // a cited run deleted since is a node named as the evidence froze it
    expect(g.nodes.get('model_run:call-gone')?.detail, 'as recorded');
    expect(g.nodes.has('criterion:c-tox'), true);
  });

  test('a generated document rests on what it was written from: its run, snapshot, data conditions, decisions', async () => {
    const r = records();
    r.criteria.push({id: 'c-dc', code: 'titer', name: 'Titer above 2 g/L', stage: 'D5'});
    r.documents.push({id: 'g', programId: 'p', title: 'Evidence table', version: '1', status: 'draft', criterionId: 'c-tox',
      path: 'X:evidence.md', generatedFrom: [{kind: 'run', ref: 'call-1'}, {kind: 'snapshot', ref: 's'},
        {kind: 'data_condition', ref: 'c-dc'}, {kind: 'decision', ref: 'y'}, {kind: 'run', ref: 'call-gone'}]});
    const g = provenanceGraph(r);
    const into = g.edges.filter((e) => e.to === 'document:g' && e.relation === 'generated into').map((e) => e.from).sort();
    expect(into.join(), 'criterion:c-dc,decision:y,model_run:r,snapshot:s', 'a run no longer recorded is not linked');
    const up = traverse(g, nodeKey('document', 'g'), 2).nodes.filter((n) => n.direction === 'up').map((n) => n.key);
    expect(up.includes('model:Pkg:M') && up.includes('table:t1'), true, 'through the run to its model and input');
    const down = traverse(g, nodeKey('model_run', 'r'), 1).nodes.filter((n) => n.direction === 'down').map((n) => n.key);
    expect(down.includes('document:g'), true, 'the run leads to the table generated from it');
  });

  test('a walk from a decision: what it rests on, up to the data, and the program it moved', async () => {
    const r = records();
    const g = provenanceGraph(r);
    const p = traverse(g, nodeKey('decision', 'x'), 4);
    const up = p.nodes.filter((n) => n.direction === 'up').map((n) => n.key);
    for (const k of ['model_run:r', 'document:d', 'snapshot:s', 'criterion:c-pd', 'model:Pkg:M', 'table:t1', 'dataset:p:data',
      'compound:k'])
      expect(up.includes(k), true, k);
    expect(p.nodes.filter((n) => n.direction === 'down').map((n) => n.key).join(), 'program:p');
    expect(p.nodes.some((n) => n.key === 'decision:y'), false, 'another decision is not on the path');
    const tree = provenanceTree(p, 'up');
    expect(tree.filter((t) => t.level === 0).map((t) => t.node.key).sort().join(),
      'criterion:c-pd,document:d,model_run:r,snapshot:s,url:https://example.org/r');
    expect(tree.find((t) => t.node.key === 'table:t1')!.level, 1);
    expect(traverse(g, nodeKey('decision', 'x'), 1).nodes.length, 7);
    expect(traverse(g, 'program:none').nodes.length, 0);
  });

  test('from a run down to its document, the decision and the program; starts from evidence and run ids', async () => {
    const r = records();
    const g = provenanceGraph(r);
    const p = traverse(g, nodeKey('model_run', 'r'), 4);
    const down = provenanceTree(p, 'down').map((t) => `${t.level}:${t.node.key}`);
    expect(down.join(' '), '0:table:t2 0:document:d 1:criterion:c-pd 0:decision:x 1:program:p');
    expect(startKey(g, r, 'evidence', 'e2'), 'document:d');
    expect(startKey(g, r, 'evidence', 'e1'), 'model_run:r');
    expect(startKey(g, r, 'model_run', 'call-1'), 'model_run:r');
    expect(startKey(g, r, 'compound', 'EXC-1'), 'compound:k');
    expect(startKey(g, r, 'snapshot', 'nope'), null);
    const text = summarizeProvenance(traverse(g, nodeKey('document', 'd'), 2));
    expect(text.split('\n')[0], 'Document: PD (v1.0 approved)');
    expect(text.includes('    Model: Titer model · model of'), true, text);
    expect(text.includes('Used by:\n  Criterion: Process description (D5 · process_description) · answers'), true, text);
  });
});

/** Everything the platform test made, removed in reverse in `after`. */
const undo: (() => Promise<unknown>)[] = [];

category('Portal: provenance', () => {
  let program: ProgramRow;
  let documentId: string;
  let runRowId: string;

  before(async () => {
    await loadAllPortfolios(_package);
    const bioc = await portfolioByCode('EX-BIOC');
    const st = await stagesOf(bioc.stage_model_id);
    const [ins] = await portalDb.programs.insert({portfolio_id: bioc.id, code: `TEST-PV-${Date.now() % 1e9}`, name: 'Provenance test',
      stage_model_id: bioc.stage_model_id, stage_id: st.find((s) => s.code === 'D4')!.id, source: 'manual', modality: 'mAb'});
    undo.push(() => portalDb.programs.delete(ins.id));
    program = await portalDb.programs.get(ins.id);
    const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
    const call = await f.prepare({feedRate: 4, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
    const r = await saveModelRun(call, 'EX-BIOC', program.code, {unitOperation: 'EX-PROD', runType: 'EX-SIM', title: 'PV feed 4', link: false});
    for (const id of r.tableIds)
      undo.push(() => deleteRunTable(id));
    undo.push(async () => grok.dapi.functions.calls.delete(await grok.dapi.functions.calls.find(r.runId)));
    undo.push(() => portalDb.modelRuns.delete(r.modelRun.id));
    runRowId = r.modelRun.id;
    const pd = (await portalDb.gateCriterions.getByKey({stage_id: st.find((s) => s.code === 'D5')!.id, code: 'process_description'}))!;
    await selectRunForCriterion(program.id, pd.id, runRowId, 'the higher titer');
    const doc = await linkDocumentTo(program, {path: `${PACKAGE_FILES}documents/ex-bioc/exb-201-process-description-v0.3.md`,
      criterion: pd.id, title: 'PV process description', version: '1.0', status: 'approved'});
    documentId = doc.document.id;
    undo.push(async () => {
      const by = {filter: `program_id = "${program.id}"`};
      for (const d of await portalDb.documents.query(by)) await portalDb.documents.delete(d.id);
      for (const a of await portalDb.criterionArtifacts.query(by)) await portalDb.criterionArtifacts.delete(a.id);
      for (const s of await portalDb.criterionStatuses.query(by)) await portalDb.criterionStatuses.delete(s.id);
    });
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
    document.querySelectorAll('[data-u2-name="provenance-dialog"]').forEach((d) => d.closest('.d4-dialog')?.remove());
  });

  test('a document rests on the run selected for its criterion and its model; the run wrote its tables', async () => {
    const p = await provenanceOf('document', documentId);
    const up = p.nodes.filter((n) => n.direction === 'up');
    expect(up.find((n) => n.kind === 'model_run')?.name, 'PV feed 4');
    expect(up.some((n) => n.kind === 'model' && n.id === 'ProjectPortalCore:FedBatchTiter'), true, up.map((n) => n.key).join());
    expect(p.nodes.filter((n) => n.direction === 'down').map((n) => `${n.kind} ${n.relation}`).join(), 'criterion answers');
    const fromRun = await provenanceOf('model_run', (await portalDb.modelRuns.get(runRowId)).run_id!);
    expect(fromRun.nodes.some((n) => n.direction === 'down' && n.key === `document:${documentId}`), true);
    // the saved call's output table (the time course) is what the run wrote
    expect(fromRun.nodes.some((n) => n.direction === 'down' && n.kind === 'table' && n.relation === 'output'), true,
      fromRun.nodes.map((n) => n.key).join());
  });

  test('from a decision on EX-SMD: its readiness and snapshot up to the program data, and the program it moved', async () => {
    const smd = await portfolioByCode('EX-SMD');
    const st = await stagesOf(smd.stage_model_id);
    const code = `TEST-PV-${Date.now() % 1e9}`;
    const [row] = await portalDb.programs.insert({portfolio_id: smd.id, code, name: code, stage_model_id: smd.stage_model_id,
      stage_id: st.find((s) => s.code === 'TV')!.id, source: 'manual'});
    try {
      const c = (await portalDb.gateCriterions.query({filter: `stage_id = "${st.find((s) => s.code === 'HF')!.id}" and code = "assay_ready"`}))[0];
      await setCriterion(row.id, c.id, 'met');
      const done = await recordDecision(row.id, 'advance', null, 'go', {source: 'approver', snapshot: true});
      const p = await provenanceOf('decision', done.transitionId!);
      const up = p.nodes.filter((n) => n.direction === 'up').map((n) => n.kind).sort();
      expect(up.join(), 'criterion,dataset,snapshot');
      expect(p.nodes.filter((n) => n.direction === 'down').map((n) => n.key).join(), `program:${row.id}`);
      const evidence = (await portalDb.evidences.query({filter: `program_id = "${row.id}" and kind = "snapshot"`}))[0];
      expect((await provenanceOf('evidence', evidence.id)).start, `snapshot:${evidence.target}`);
      const result: any = await grok.functions.call('ProjectPortalCore:provenance', {kind: 'program', id: row.id, depth: 2});
      expect(`${result.start} ${result.nodes.length} ${result.edges.length}`, `program:${row.id} 4 3`);
      const dlg = await openProvenanceDialog('decision', done.transitionId!);
      const q = (n: string) => dlg.root.querySelector(`[data-u2-name="${n}"]`) as HTMLElement | null;
      expect(q('provenance-start')?.textContent, `Decision: advance ${new Date().toISOString().slice(0, 10)} (TV → HF)`);
      expect(q('provenance-up')?.textContent?.includes('Snapshot: '), true, q('provenance-up')?.textContent ?? 'no tree');
      expect(q('provenance-down')?.textContent?.includes(`Program: ${code}`), true);
      dlg.close();
      const v = new PortalView();
      await v.load('EX-SMD');
      await v.showProgram('EX-SMD', row.id);
      expect(!!v.detailRoot.querySelector('[data-u2-name="provenance-decision-0"]'), true);
      expect(!!v.detailRoot.querySelector('[data-u2-name="provenance-evidence-0-0"]'), true);
      expect(!!v.detailRoot.querySelector('[data-u2-name="program-provenance-button"]'), true);
    } finally {
      const by = {filter: `program_id = "${row.id}"`};
      for (const e of await portalDb.evidences.query(by)) await portalDb.evidences.delete(e.id).catch(() => undefined);
      for (const s of await portalDb.programSnapshots.query(by)) await portalDb.programSnapshots.delete(s.id).catch(() => undefined);
      for (const t of await portalDb.stageTransitions.query(by)) await portalDb.stageTransitions.delete(t.id).catch(() => undefined);
      for (const s of await portalDb.criterionStatuses.query(by)) await portalDb.criterionStatuses.delete(s.id).catch(() => undefined);
      await portalDb.programs.delete(row.id);
    }
  });
});
