import * as grok from 'datagrok-api/grok';
import {before, category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, stagesOf} from '../portal';
import {checkRecord, DecisionView, describeDecision, describeRecord, RECORD_NOTE, searchDecisions} from '../model/decision-records';
import {decisionsOf, recordGateDecision} from '../decision-records';
import {linkEvidence} from '../model/evidence';
import {withRestrictedUser} from './restricted-user';
import {_package} from '../package-test';

const ok = {stageAuthority: 'feed', isApprover: true, programStatus: 'active', outcome: 'go', rationale: 'Review passed',
  gateStageId: 's2', stages: [{id: 's1', ordinal: 1}, {id: 's2', ordinal: 2}], currentStageId: 's1'};

const view = (id: string, on: string, extra: Partial<DecisionView> = {}): DecisionView => ({id, portfolio: 'EX', programId: 'p1',
  programCode: 'P-1', programName: 'Program One', decision: 'record', outcome: 'go', from: {code: 'D1', name: 'Discovery'},
  to: {code: 'D2', name: 'Development'}, on, by: null, source: 'approver', rationale: null, evidence: [], ...extra});

category('Portal: decision records, rules', () => {
  test('a record needs a registry-owned stage, an approver, an outcome, a rationale and a gate', async () => {
    expect(checkRecord(ok).length, 0);
    expect(checkRecord({...ok, stageAuthority: 'portal'})[0], 'Decisions here move the stage: record them as gate decisions');
    expect(checkRecord({...ok, isApprover: false}).some((r) => r.includes('approvers')), true);
    expect(checkRecord({...ok, outcome: 'maybe'}).some((r) => r.includes('Go, Hold, No-go')), true);
    expect(checkRecord({...ok, rationale: '  '}).some((r) => r.includes('rationale')), true);
    expect(checkRecord({...ok, gateStageId: null}).some((r) => r.includes('no gate')), true);
    expect(checkRecord({...ok, gateStageId: 'other'}).some((r) => r.includes('different stage model')), true);
    expect(checkRecord({...ok, programStatus: 'terminated'}).some((r) => r.includes('terminated')), true);
  });

  test('a record says the registry sets the stage', async () => {
    expect(describeRecord('no_go', 'Development', 'Yield too low'),
      `No-go at the gate to Development — Yield too low (${RECORD_NOTE})`);
    expect(describeDecision(view('d1', '2026-10-01T00:00:00Z')).includes(RECORD_NOTE), true);
  });

  test('search: words in rationale or evidence, program, gate, portfolio; newest first', async () => {
    const all = [
      view('d1', '2026-09-01T00:00:00Z', {rationale: 'Titer met the target', evidence: [{kind: 'url', name: 'Minutes', target: 'https://x/m'}]}),
      view('d2', '2026-09-05T00:00:00Z', {decision: 'sync', outcome: null, rationale: 'Registry feed'}),
      view('d3', '2026-09-03T00:00:00Z', {programCode: 'P-2', programName: 'Second', to: {code: 'D3', name: 'Launch'}, outcome: 'hold'}),
      view('d4', '2026-09-04T00:00:00Z', {portfolio: 'OTHER'}),
    ];
    const ids = (q: Parameters<typeof searchDecisions>[1]) => searchDecisions(all, q).map((d) => d.id).join();
    expect(ids({}), 'd2,d4,d3,d1');
    expect(ids({text: 'titer minutes'}), 'd1');
    expect(ids({text: 'hold'}), 'd3');
    expect(ids({program: 'p-2'}), 'd3');
    expect(ids({program: 'second'}), 'd3');
    expect(ids({gate: 'launch'}), 'd3');
    expect(ids({gate: 'D2', withoutSync: true, portfolio: 'EX'}), 'd1');
    expect(ids({withEvidence: true}), 'd1');
  });
});

/** A manual program of EX-BIOC at D2 for [fn]; it and everything written for it are removed after. */
async function scratchAtD2(fn: (p: ProgramRow, bioc: PortfolioRow, stage: (code: string) => string) => Promise<void>): Promise<void> {
  const bioc = await portfolioByCode('EX-BIOC');
  const stages = await stagesOf(bioc.stage_model_id);
  const stage = (code: string) => stages.find((s) => s.code === code)!.id;
  const code = `TEST-REC-${Date.now() % 1e9}`;
  const [row] = await portalDb.programs.insert({portfolio_id: bioc.id, code, name: code, stage_model_id: bioc.stage_model_id,
    stage_id: stage('D2'), source: 'manual'});
  try {
    await fn(await portalDb.programs.get(row.id), bioc, stage);
  } finally {
    const by = {filter: `program_id = "${row.id}"`};
    for (const e of await portalDb.evidences.query(by)) await portalDb.evidences.delete(e.id).catch(() => undefined);
    for (const t of await portalDb.stageTransitions.query(by)) await portalDb.stageTransitions.delete(t.id).catch(() => undefined);
    await portalDb.programs.delete(row.id);
  }
}

category('Portal: decision records where the registry owns the stage', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  test('a record keeps outcome, rationale and evidence, and the program stays where the registry put it', async () => {
    await scratchAtD2(async (p, bioc, stage) => {
      const stamp = `rec${Date.now() % 1e6}`;
      const r = await recordGateDecision(p.id, {outcome: 'go', rationale: `Tox package complete ${stamp}`,
        items: [linkEvidence(`https://example.org/minutes/${stamp}`)]});
      expect(r.evidence.written, 1, r.evidence.errors.join('; '));
      expect(r.line.includes(RECORD_NOTE), true, r.line);
      const after = await portalDb.programs.get(p.id);
      expect(after.stage_id, stage('D2'), 'the record moved the program');
      expect(after.status ?? 'active', p.status ?? 'active', 'the record changed the status');
      const t = await portalDb.stageTransitions.get(r.transitionId);
      expect(`${t.decision} ${t.outcome} ${t.from_stage_id === stage('D2')} ${t.to_stage_id === stage('D3')}`, 'record go true true');

      const found = searchDecisions(await decisionsOf('EX-BIOC'), {text: stamp});
      expect(found.map((d) => d.programCode).join(), p.code);
      expect(found[0].evidence.length, 1);
      expect(searchDecisions(await decisionsOf(), {program: p.code, gate: 'D3', withEvidence: true}).length, 1);
      const df = await grok.functions.call('ProjectPortal:findDecisions', {text: stamp});
      expect(df.rowCount, 1);
    });
  });

  test('where the portal owns the stage, a record is refused: decide the gate instead', async () => {
    const smd = await portfolioByCode('EX-SMD');
    const [p] = await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`});
    await expectExceptionAsync(async () => {
      await recordGateDecision(p.id, {outcome: 'go', rationale: 'x'});
    });
  });

  test('a viewer can neither record in the app nor write the gate log on the server', async () => {
    await scratchAtD2(async (p, bioc, stage) => {
      const ran = await withRestrictedUser('pprec', async (u) => {
        const viewers = await grok.dapi.groups.find(bioc.viewers!);
        const personal = await grok.dapi.groups.find(u.group);
        await grok.dapi.groups.includeTo(personal, viewers);
        try {
          let refused = false;
          await u.asUser(() => portalDb.stageTransitions.insert({program_id: p.id, from_stage_id: stage('D2'),
            to_stage_id: stage('D3'), decision: 'record', outcome: 'go', rationale: 'by a viewer', source: 'approver'}))
            .catch(() => refused = true);
          expect(refused, true, 'a viewer wrote the gate log');
        } finally {
          await grok.dapi.groups.excludeFrom(personal, viewers);
        }
      });
      expect(ran !== null, true, 'no restricted session could be created on this stand');
      expect(checkRecord({...ok, isApprover: false}).length > 0, true);
    });
  });
});
