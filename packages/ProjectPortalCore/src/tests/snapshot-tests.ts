import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  checkTable, describeSnapshot, snapshotCriteria, snapshotEvidence, snapshotFileName,
} from '../model/snapshots';
import {computeReadiness, ReadinessCriterion} from '../model/readiness';
import {portalDb, ProgramRow, ProgramSnapshotRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, recordDecision, setCriterion, stagesOf} from '../portal';
import {buildPortfolioWorkspace} from '../workspace';
import {frozenColumns, frozenCriteria, sha256, snapshotsOf, snapshotTable, takeSnapshot} from '../snapshots';
import {evidenceOf} from '../evidence';
import {openSnapshotDialog} from '../views/snapshot-dialog';
import {PortalView} from '../views/portal-view';
import {_package} from '../package-test';

const stages = [{id: 's1', code: 'A', ordinal: 1}, {id: 's2', code: 'B', ordinal: 2}];
const criteria = [
  {id: 'k2', stageId: 's2', code: 'potent', name: 'Potent', kind: 'data_condition', required: true},
  {id: 'k1', stageId: 's1', code: 'doc', name: 'Doc', kind: 'document', required: true},
  {id: 'k3', stageId: 's2', code: 'extra', name: 'Extra', kind: 'document', required: false},
];
const readinessCriteria: ReadinessCriterion[] = [
  {id: 'k2', code: 'potent', name: 'Potent', kind: 'data_condition', required: true,
    condition: {where: [{column: 'ic50', op: '<', value: 1}]}},
  {id: 'k3', code: 'extra', name: 'Extra', kind: 'document', required: false},
];

category('Portal: snapshot rules', () => {
  test('every criterion in stage order, recorded or as evaluated for the next gate', async () => {
    const r = computeReadiness('s2', readinessCriteria, {}, [{compound: 'X1', ic50: 0.2}]);
    const frozen = snapshotCriteria(criteria, stages, {k1: 'met', k2: 'not_met'}, r);
    expect(frozen.map((c) => `${c.stage}:${c.code}:${c.state}:${c.basis}`).join(' '),
      'A:doc:met:recorded B:potent:met:evaluated B:extra:missing:recorded');
    expect(!!frozen[1].summary, true);
    expect(snapshotCriteria(criteria, stages, {}, null)[1].state, 'not_started');
  });

  test('a file name per program and second; the table check; the summary line', async () => {
    expect(snapshotFileName('EXS-KIN1', '2026-10-02T10:31:05.123Z'), 'snapshot-EXS-KIN1-2026-10-02T103105Z.csv');
    expect(snapshotFileName('A B/C', '2026-10-02T10:31:05Z'), 'snapshot-A_B_C-2026-10-02T103105Z.csv');
    expect(checkTable('x', 'x', true), 'intact');
    expect(checkTable('x', 'y', true), 'changed');
    expect(checkTable('x', null, true), 'missing');
    expect(checkTable(null, null, false), 'not kept');
    const s = {takenOn: '2026-10-02T10:31:05Z', label: 'L', stageName: 'Alpha',
      criteria: snapshotCriteria(criteria, stages, {k1: 'met'}, null), tableRows: 12};
    expect(describeSnapshot(s), '2026-10-02 10:31 · Alpha · 1 of 2 required met or waived · 12 rows');
    expect(describeSnapshot({...s, stageName: null, criteria: [], tableRows: null}), '2026-10-02 10:31 · unstaged · no required criteria · no table');
    const e = snapshotEvidence('id-1', s);
    expect(`${e.kind} ${e.target} ${e.name} ${e.context!.stage} ${e.context!.rows}`, 'snapshot id-1 L Alpha 12');
  });
});

/** Removes [s]: its file, then its row. Tests only; the app never deletes a snapshot. */
async function removeSnapshot(s: ProgramSnapshotRow): Promise<void> {
  if (s.table_file && s.table_folder)
    await grok.dapi.spaces.id(s.table_folder).files.delete(s.table_file).catch(() => undefined);
  await portalDb.programSnapshots.delete(s.id).catch(() => undefined);
}

category('Portal: program snapshots', () => {
  let kin1: ProgramRow;

  before(async () => {
    await loadAllPortfolios(_package);
    const smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`}))[0];
    await buildPortfolioWorkspace('EX-SMD', false, kin1.id);
    kin1 = await portalDb.programs.get(kin1.id);
  });

  test('a snapshot of EXS-KIN1 keeps its SAR table in the program folder, fingerprinted; a changed file is noticed', async () => {
    const {snapshot: s, notes} = await takeSnapshot(kin1, {note: 'package test'});
    try {
      expect(notes.length, 0, notes.join('; '));
      // unstaged on a fresh stand, at Hit-to-Lead once the example history is seeded
      const stage = (await stagesOf(kin1.stage_model_id)).find((x) => x.id === kin1.stage_id)?.code ?? 'unstaged';
      expect(`${s.status} ${s.stage_code ?? 'unstaged'} ${s.dataset}`, `active ${stage} SAR table`);
      expect(s.table_file!.startsWith('snapshot-EXS-KIN1-'), true, s.table_file);
      expect(s.table_folder, kin1.space_id);
      expect(frozenCriteria(s).some((c) => c.code === 'assay_ready'), true);
      const {df, check} = await snapshotTable(s);
      expect(check, 'intact');
      expect(df!.rowCount, s.table_rows);
      expect(df!.rowCount > 0, true);
      expect(df!.columns.names().join(), frozenColumns(s).join());
      expect((await snapshotsOf(kin1.id))[0].id, s.id);
      await grok.dapi.spaces.id(s.table_folder!).files.writeString(s.table_file!, 'compound\nEXC-9999\n');
      expect((await snapshotTable(s)).check, 'changed');
      expect(s.table_sha256 === await sha256('compound\nEXC-9999\n'), false);
    } finally {
      await removeSnapshot(s);
    }
  });

  test('the functions: snapshotProgram, then programSnapshots lists it', async () => {
    const report: string = await grok.functions.call('ProjectPortalCore:snapshotProgram',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', label: 'TEST snapshot'});
    const s = (await snapshotsOf(kin1.id)).find((x) => x.label === 'TEST snapshot')!;
    try {
      expect(report.startsWith('EXS-KIN1: snapshot "TEST snapshot" taken:'), true, report);
      const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:programSnapshots', {portfolio: 'EX-SMD', program: 'EXS-KIN1'});
      expect(df.col('label')!.toList().includes('TEST snapshot'), true);
      const dlg = await openSnapshotDialog(s.id);
      const q = (n: string) => dlg.root.querySelector(`[data-u2-name="${n}"]`) as HTMLElement | null;
      expect(q('snapshot-table')?.textContent?.includes('the file holds what was frozen'), true, q('snapshot-table')?.textContent ?? '');
      expect(q('snapshot-criteria')?.textContent?.includes('Primary assay validated'), true);
      dlg.close();
    } finally {
      if (s)
        await removeSnapshot(s);
    }
  });

  test('a decision takes a snapshot now as evidence; without a workspace the table is fingerprinted, not kept', async () => {
    const smd = await portfolioByCode('EX-SMD');
    const st = await stagesOf(smd.stage_model_id);
    const code = `TEST-SN-${Date.now() % 1e9}`;
    const [row] = await portalDb.programs.insert({portfolio_id: smd.id, code, name: code, stage_model_id: smd.stage_model_id,
      stage_id: st.find((s) => s.code === 'TV')!.id, source: 'manual'});
    try {
      const c = (await portalDb.gateCriterions.query({filter: `stage_id = "${st.find((s) => s.code === 'HF')!.id}" and code = "assay_ready"`}))[0];
      await setCriterion(row.id, c.id, 'met');
      const done = await recordDecision(row.id, 'advance', null, 'go', {source: 'approver', items: [], snapshot: true});
      expect(done.evidence?.written, 1, JSON.stringify(done.evidence));
      const [s] = await snapshotsOf(row.id);
      expect(`${s.stage_code} ${s.table_rows} ${s.table_file ?? 'no file'}`, 'TV 0 no file');
      expect(s.note?.includes('has no workspace folder'), true, s.note);
      const lines = [...(await evidenceOf(row.id)).values()].flat();
      expect(lines.map((l) => `${l.row.kind} ${l.row.target === s.id} ${l.state}`).join(), 'snapshot true live');
      const v = new PortalView();
      await v.load('EX-SMD');
      await v.showProgram('EX-SMD', row.id);
      const q = (n: string) => v.detailRoot.querySelector(`[data-u2-name="${n}"]`) as HTMLElement | null;
      expect(q('snapshots-heading')?.textContent, 'Snapshots (1)');
      expect(q('snapshot-0')?.textContent?.includes('Target Validation'), true, q('snapshot-0')?.textContent ?? '');
      expect(q('decision-evidence-0')?.textContent?.includes('Snapshot: '), true, q('decision-evidence-0')?.textContent ?? '');
    } finally {
      const by = {filter: `program_id = "${row.id}"`};
      for (const e of await portalDb.evidences.query(by)) await portalDb.evidences.delete(e.id).catch(() => undefined);
      for (const s of await portalDb.programSnapshots.query(by)) await removeSnapshot(s);
      for (const t of await portalDb.stageTransitions.query(by)) await portalDb.stageTransitions.delete(t.id).catch(() => undefined);
      for (const s of await portalDb.criterionStatuses.query(by)) await portalDb.criterionStatuses.delete(s.id).catch(() => undefined);
      await portalDb.programs.delete(row.id);
    }
  });
});
