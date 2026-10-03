import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {buildManifest, ManifestInput, manifestCsv, manifestFileName, manifestStages, summarizeManifest} from '../model/data-package';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programContext} from '../portal';
import {loadExampleDocuments} from '../documents';
import {saveModelRun} from '../model-runs';
import {dataPackageOf} from '../data-package';
import {openDataPackageDialog} from '../views/data-package-dialog';
import {filterValue, ROW_LIMIT} from '../rows';
import {_package} from '../package-test';

const stages = [{id: 's1', code: 'A', ordinal: 1}, {id: 's2', code: 'B', ordinal: 2}, {id: 's3', code: 'C', ordinal: 3}];
const input = (over: Partial<ManifestInput> = {}): ManifestInput => ({
  stages, stageId: 's1',
  criteria: [
    {id: 'c3', stageId: 's3', code: 'late', name: 'Late', kind: 'document', required: true},
    {id: 'c2', stageId: 's2', code: 'assay', name: 'Assay valid', kind: 'data_condition', required: true},
    {id: 'c1', stageId: 's1', code: 'start', name: 'Started', kind: 'decision', required: false},
  ],
  states: {c1: 'met', c2: 'not_met'},
  evaluated: new Set(['c2']),
  documents: [
    {id: 'd1', title: 'Report, "final"', kind: 'development_report', version: '1.0', status: 'approved', path: 'S:/r.md',
      linkedBy: 'u1', linkedOn: '2026-10-01T10:00:00Z', approvedBy: 'u2', approvedOn: '2026-10-02T09:00:00Z'},
    {id: 'd2', title: 'Evidence table: Assay valid', kind: 'other', version: '2', status: 'draft', path: 'S:/e.md',
      criterionId: 'c2', linkedBy: 'u1', linkedOn: '2026-10-02T11:00:00Z',
      generated: {sources: ['run', 'snapshot', 'run'], state: 'stale'}},
  ],
  runs: [{id: 'r1', runId: 'call-1', title: 'Feed 4', model: 'Titer', modelVersion: '1.2', stageId: 's2',
    unitOperation: 'Production', runType: 'Simulated', purpose: 'characterization', by: 'u1', on: '2026-09-30T08:00:00Z',
    selectedFor: ['c2']},
  {id: 'r2', title: 'Feed 2', model: 'Titer', selectedFor: [], isReference: true}],
  datasets: [{code: 'program', name: 'Process data', connection: 'System:Datagrok'}],
  snapshots: [{id: 'sn1', label: 'At gate B', by: 'u2', on: '2026-10-02T12:00:00Z'}],
  names: new Map([['u1', 'Ann'], ['u2', 'Bo']]),
  ...over,
});

category('Portal: data package rules', () => {
  test('criteria to the next gate, in stage order, with state and basis', async () => {
    expect([...manifestStages(stages, 's1')].join(), 's1,s2');
    expect([...manifestStages(stages, null)].join(), 's1', 'an unstaged program: the first stage');
    expect([...manifestStages(stages, 's3')].join(), 's1,s2,s3');
    const lines = buildManifest(input()).filter((l) => l.section === 'criterion');
    expect(lines.map((l) => `${l.stage}/${l.criterion}:${l.status}`).join(), 'A/start:met,B/assay:not met');
    expect(lines[0].source, 'decision, optional · recorded');
    expect(lines[1].source, 'data condition · evaluated from data now');
  });

  test('documents, evidence tables, runs, datasets and snapshots, with provenance', async () => {
    const lines = buildManifest(input());
    expect(lines.map((l) => l.section).join(), 'criterion,criterion,document,evidence_table,run,run,dataset,snapshot');
    const doc = lines.find((l) => l.section === 'document')!;
    expect([doc.status, doc.version, doc.by, doc.on].join(), 'approved,1.0,Bo,2026-10-02', 'approved: by the approver');
    const table = lines.find((l) => l.section === 'evidence_table')!;
    expect([table.criterion, table.status, table.source].join('|'), 'assay|draft, evidence changed|generated from run, snapshot');
    const [selected, reference] = lines.filter((l) => l.section === 'run');
    expect([selected.criterion, selected.status, selected.source, selected.by, selected.ref, selected.stage].join('|'),
      'assay|selected|Titer · Production · Simulated · characterization|Ann|call-1|B');
    expect([reference.status, reference.ref].join(), 'reference,r2', 'no call id: the row id');
    expect(lines.find((l) => l.section === 'dataset')!.source, 'System:Datagrok');
    expect(summarizeManifest(lines), '2 criteria · 1 document · 1 evidence table · 2 runs · 1 dataset · 1 snapshot');
    expect(summarizeManifest([]), 'Nothing recorded yet');
  });

  test('CSV: a header, a line each, quoted where needed; the file name', async () => {
    const csv = manifestCsv(buildManifest(input()));
    const rows = csv.trimEnd().split('\n');
    expect(rows[0], 'Section,Item,Stage,Criterion,Status,Version,Source,By,On,Reference');
    expect(rows.length, 9);
    expect(rows.some((r) => r.startsWith('Document,"Report, ""final""",')), true, rows.join('\n'));
    expect(rows.some((r) => r.startsWith('Evidence table,Evidence table: Assay valid,,assay,"draft, evidence changed"')), true);
    expect(manifestFileName('EX B/1', new Date('2026-10-03T01:00:00Z')), 'EX_B_1-data-package-2026-10-03.csv');
  });
});

category('Portal: data package', () => {
  let code = '';
  let programId = '';
  const undo: (() => Promise<any>)[] = [];

  before(async () => {
    await loadAllPortfolios(_package);
    await loadExampleDocuments(_package);
    const bioc = await portfolioByCode('EX-BIOC');
    const d5 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D5'}))!;
    code = `DP-${Date.now() % 1e9}`;
    const [p] = await portalDb.programs.insert({portfolio_id: bioc.id, code, name: 'Data package test',
      stage_model_id: bioc.stage_model_id, stage_id: d5.id, source: 'manual', modality: 'mAb'});
    programId = p.id;
    undo.push(() => portalDb.programs.delete(p.id));
    const f = DG.Func.find({package: 'ProjectPortalCore', name: 'FedBatchTiter'})[0];
    const call = await f.prepare({feedRate: 3, temperature: 36.5, days: 14, seedingDensity: 0.5}).call();
    const saved = await saveModelRun(call, 'EX-BIOC', code, {unitOperation: 'EX-PROD', runType: 'EX-SIM',
      purpose: 'characterization', title: 'Package run'});
    undo.push(async () => {
      const back = await grok.dapi.functions.calls.find(saved.runId).catch(() => null);
      for (const id of (back?.options?.['portalTables'] ?? []) as string[]) {
        const t = await grok.dapi.tables.find(id).catch(() => null);
        if (t)
          await grok.dapi.tables.delete(t);
      }
      for (const r of await portalDb.modelRuns.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT}))
        await portalDb.modelRuns.delete(r.id);
      if (back)
        await grok.dapi.functions.calls.delete(back);
    });
  });

  after(async () => {
    for (const u of undo.splice(0).reverse())
      await u().catch(() => {});
  });

  test('EXB-201: its criteria to D5, its example documents, its process dataset', async () => {
    const bioc = await portfolioByCode('EX-BIOC');
    const exb201 = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(bioc.id)} and code = "EXB-201"`, limit: ROW_LIMIT}))[0];
    const p = await dataPackageOf(exb201.id);
    const ctx = await programContext(exb201.id);
    const upTo = new Set(ctx.stages.filter((s) => s.ordinal <= 5).map((s) => s.id));
    expect(p.lines.filter((l) => l.section === 'criterion').length, ctx.criteria.filter((c) => upTo.has(c.stageId)).length,
      'every criterion of D1 to D5');
    const docs = p.lines.filter((l) => l.section === 'document').map((l) => `${l.name} ${l.version} ${l.status}`);
    expect(docs.includes('Process description, Antibody Alpha 0.3 draft'), true, docs.join('; '));
    expect(docs.includes('Upstream development report, Antibody Alpha 1.0 approved'), true, docs.join('; '));
    expect(p.lines.some((l) => l.section === 'dataset' && l.source === 'System:Datagrok'), true);
  });

  test('a saved run is a line with its model, unit operation, run type and who ran it', async () => {
    const p = await dataPackageOf(programId);
    const run = p.lines.find((l) => l.section === 'run')!;
    expect(run.name, 'Package run');
    expect(run.status, 'saved');
    expect(run.stage, 'D5');
    expect(run.source, 'Fed-batch titer (example) · Production bioreactor · Simulated · characterization');
    expect(run.by, DG.User.current().friendlyName);
    expect(p.lines.filter((l) => l.section === 'criterion').every((l) => ['D1', 'D2', 'D3', 'D4', 'D5', 'D6'].includes(l.stage!)), true);
  });

  test('the functions: a table and CSV of the same lines', async () => {
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:dataPackage', {portfolio: 'EX-BIOC', program: code});
    const csv: string = await grok.functions.call('ProjectPortalCore:dataPackageCsv', {portfolio: 'EX-BIOC', program: code});
    expect(df.name, `${code} data package`);
    expect(df.columns.names().join(), 'Section,Item,Stage,Criterion,Status,Version,Source,By,On,Reference');
    expect(csv.trimEnd().split('\n').length, df.rowCount + 1);
    expect(df.col('Section')!.toList().includes('Run'), true);
  });

  test('the dialog lists the lines and offers the CSV', async () => {
    const dlg = await openDataPackageDialog(programId);
    try {
      const root = dlg.root;
      expect(root.querySelector('[data-u2-name="data-package-summary"]')!.textContent!.includes('1 run'), true);
      expect(root.querySelectorAll('[data-u2-name^="data-package-line-"]').length > 1, true);
      expect(!!root.querySelector('[data-u2-name="data-package-csv"]'), true);
    } finally {
      dlg.close();
    }
  });
});
