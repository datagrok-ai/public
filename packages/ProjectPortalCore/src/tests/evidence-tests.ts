import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  describeEvidence, EvidenceItem, linkEvidence, orderEvidence, readinessEvidence, suggestEvidence, targetState,
  uniqueEvidence, validateEvidence,
} from '../model/evidence';
import {computeReadiness, ReadinessCriterion} from '../model/readiness';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, recordDecision, setCriterion, stagesOf} from '../portal';
import {autoAdvance} from '../auto-advance';
import {evidenceOf, selectedRunsOf} from '../evidence';
import {PortalView} from '../views/portal-view';
import {gpc3AtHitFinding} from './auto-advance-tests';
import {_package} from '../package-test';

const criteria: ReadinessCriterion[] = [
  {id: 'k1', code: 'potent', name: 'Potent', kind: 'data_condition', required: true,
    condition: {where: [{column: 'ic50', op: '<', value: 1}]}},
  {id: 'k2', code: 'report', name: 'Report', kind: 'document', required: true},
  {id: 'k3', code: 'extra', name: 'Extra', kind: 'data_condition', required: false, condition: {measured: ['herg']}},
];
const readiness = computeReadiness('s2', criteria, {k2: 'met'}, [{compound: 'X1', ic50: 0.3}]);
const run = {criterionId: 'k2', criterionCode: 'report', criterionName: 'Report', runId: 'call-1', title: 'Titer at 3%',
  model: 'Pkg:Titer', selectedOn: '2026-10-02T10:00:00.000Z', selectedBy: 'Ann', rationale: 'closest'};

category('Portal: evidence rules', () => {
  test('an item needs a known kind, a target and a name; a URL is http(s)', async () => {
    expect(validateEvidence({kind: 'url', target: 'https://example.org/a', name: 'a', role: 'supporting'}).length, 0);
    expect(validateEvidence({kind: 'url', target: 'ftp://x', name: 'x'})[0], '"ftp://x" is not an http(s) address');
    expect(validateEvidence({kind: 'note' as any, target: 'x', name: 'x'})[0].startsWith('evidence kind "note" is not one of'), true);
    expect(validateEvidence({kind: 'file', target: ' ', name: ''}).join('; '), 'evidence needs a target; evidence needs a name');
  });

  test('a typed link is a URL or a file path, named by its last segment', async () => {
    const u = linkEvidence(' https://example.org/reviews/hf-gate ');
    expect(`${u.kind} ${u.target} ${u.name}`, 'url https://example.org/reviews/hf-gate example.org/…/hf-gate');
    expect(linkEvidence('https://example.org').name, 'example.org');
    const f = linkEvidence('System:AppData/Reports/gate.pdf');
    expect(`${f.kind} ${f.name} ${f.role}`, 'file gate.pdf supporting');
  });

  test('the gate\'s readiness is one item per required criterion, frozen with its state and what the data said', async () => {
    const items = readinessEvidence(readiness, 'Beta');
    expect(items.map((e) => `${e.kind}:${e.target}:${e.role}`).join(' '), 'criterion:k1:readiness criterion:k2:readiness');
    expect(items[0].context!.state, 'met');
    expect(items[0].context!.basis, 'evaluated');
    expect(String(items[0].context!.summary).length > 0, true);
    expect(items[1].context!.basis, 'recorded');
    expect(readinessEvidence(null, 'x').length, 0);
  });

  test('the suggestion adds the run selected for a criterion of the gate, not of another stage', async () => {
    const other = {...run, criterionId: 'elsewhere', runId: 'call-2'};
    const s = suggestEvidence(readiness, 'Beta', [run, other]);
    expect(s.map((e) => `${e.kind}:${e.target}`).join(' '), 'criterion:k1 criterion:k2 model_run:call-1');
    expect(s[2].name, 'Titer at 3%');
    expect(s[2].context!.criterion, 'report');
    expect(uniqueEvidence([...s, s[2]]).length, 3);
  });

  test('a target renamed or gone keeps its frozen name; a URL is not checked', async () => {
    expect(targetState('model_run', 'A', 'A'), 'live');
    expect(targetState('model_run', 'A', 'B'), 'renamed');
    expect(targetState('model_run', 'A', null), 'gone');
    expect(targetState('file', 'A', undefined), 'unchecked');
    expect(targetState('url', 'A', null), 'unchecked');
    const e = {kind: 'model_run' as const, name: 'Titer at 3%', context: {criterion: 'report', model: 'Pkg:Titer'}};
    expect(describeEvidence(e), 'Run: Titer at 3% · for report, Pkg:Titer');
    expect(describeEvidence(e, 'renamed', 'Titer v2'), 'Run: Titer at 3% · for report, Pkg:Titer (now Titer v2)');
    expect(describeEvidence(e, 'gone'), 'Run: Titer at 3% · for report, Pkg:Titer (gone; as recorded)');
    expect(describeEvidence({kind: 'criterion', name: 'Potent', context: {state: 'not_met', summary: '0 of 2'}}),
      'Criterion: Potent · not met: 0 of 2');
    const ordered = orderEvidence([{role: 'readiness', kind: 'criterion', name: 'B'}, {role: 'supporting', kind: 'url', name: 'A'}]);
    expect(ordered.map((x) => x.name).join(), 'A,B');
  });
});

/** A manual program of EX-SMD at Target Validation for [fn]; it and everything written for it are removed after. */
async function scratchAtTv(fn: (p: ProgramRow, smd: PortfolioRow, stage: (code: string) => string) => Promise<void>): Promise<void> {
  const smd = await portfolioByCode('EX-SMD');
  const stages = await stagesOf(smd.stage_model_id);
  const stage = (code: string) => stages.find((s) => s.code === code)!.id;
  const code = `TEST-EV-${Date.now() % 1e9}`;
  const [row] = await portalDb.programs.insert({portfolio_id: smd.id, code, name: code, stage_model_id: smd.stage_model_id,
    stage_id: stage('TV'), source: 'manual'});
  try {
    await fn(await portalDb.programs.get(row.id), smd, stage);
  } finally {
    const by = {filter: `program_id = "${row.id}"`};
    for (const e of await portalDb.evidences.query(by)) await portalDb.evidences.delete(e.id).catch(() => undefined);
    for (const t of await portalDb.stageTransitions.query(by)) await portalDb.stageTransitions.delete(t.id).catch(() => undefined);
    for (const s of await portalDb.criterionStatuses.query(by)) await portalDb.criterionStatuses.delete(s.id).catch(() => undefined);
    for (const r of await portalDb.modelRuns.query(by)) await portalDb.modelRuns.delete(r.id).catch(() => undefined);
    await portalDb.programs.delete(row.id);
  }
}

const assayReady = async (smd: PortfolioRow, stage: (code: string) => string) =>
  (await portalDb.gateCriterions.query({filter: `stage_id = "${stage('HF')}" and code = "assay_ready"`}))[0];

category('Portal: decisions with evidence', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  test('an advance without a list rests on the gate\'s readiness, written beside the gate-log entry', async () => {
    await scratchAtTv(async (p, smd, stage) => {
      const c = await assayReady(smd, stage);
      await setCriterion(p.id, c.id, 'met');
      const done = await recordDecision(p.id, 'advance', null, 'assay validated');
      expect(done.evidence?.written, 1, JSON.stringify(done.evidence));
      const rows = await portalDb.evidences.query({filter: `program_id = "${p.id}"`});
      expect(rows.map((r) => `${r.kind} ${r.role} ${r.target === c.id} ${r.stage_transition_id === done.transitionId}`).join(),
        'criterion readiness true true');
      expect(rows[0].name, c.name);
      expect((rows[0].context as any).state, 'met');
      expect((rows[0].context as any).stage, 'Hit Finding');
    });
  });

  test('an explicit list replaces the suggestion; a bad item is refused, the decision stands', async () => {
    await scratchAtTv(async (p, smd, stage) => {
      await setCriterion(p.id, (await assayReady(smd, stage)).id, 'met');
      const items: EvidenceItem[] = [linkEvidence('https://example.org/gate-review'),
        {kind: 'url', target: 'not a link', name: 'bad', role: 'supporting'}];
      const done = await recordDecision(p.id, 'advance', null, 'reviewed', {source: 'approver', items});
      expect(done.evidence?.written, 1);
      expect(done.evidence?.errors[0], 'bad: "not a link" is not an http(s) address');
      expect((await portalDb.programs.get(p.id)).stage_id, stage('HF'));
      const none = await recordDecision(p.id, 'hold', null, 'budget', {source: 'approver', items: []});
      expect(none.evidence?.written, 0, 'an empty list writes none');
    });
  });

  test('addEvidence and decisionEvidence: a run renamed and then deleted keeps its frozen name', async () => {
    await scratchAtTv(async (p, smd, stage) => {
      await setCriterion(p.id, (await assayReady(smd, stage)).id, 'met');
      await recordDecision(p.id, 'advance', null, 'go', {source: 'approver', items: []});
      const [r] = await portalDb.modelRuns.insert({program_id: p.id, run_id: `${p.code}-call`, model_name: 'ProjectPortalCore:FedBatchTiter',
        title: 'Titer at 3%'});
      const added: string = await grok.functions.call('ProjectPortalCore:addEvidence',
        {portfolio: 'EX-SMD', program: p.code, kind: 'model_run', target: `${p.code}-call`});
      expect(added.startsWith(`${p.code}: 1 piece(s) of evidence added to the advance of`), true, added);
      await grok.functions.call('ProjectPortalCore:addEvidence',
        {portfolio: 'EX-SMD', program: p.code, kind: 'url', target: 'https://example.org/minutes'});
      const state = async () => {
        const lines = [...(await evidenceOf(p.id)).values()].flat();
        return lines.filter((l) => l.row.kind === 'model_run').map((l) => `${l.row.name} ${l.state}`).join();
      };
      expect(await state(), 'Titer at 3% live');
      await portalDb.modelRuns.update(r.id, {title: 'Titer v2'});
      expect(await state(), 'Titer at 3% renamed');
      await portalDb.modelRuns.delete(r.id);
      expect(await state(), 'Titer at 3% gone');
      const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:decisionEvidence', {portfolio: 'EX-SMD', program: p.code});
      const rows = Array.from({length: df.rowCount}, (_, i) => `${df.get('decision', i)} ${df.get('kind', i)} ${df.get('state', i)}`);
      expect(rows.join(' | '), 'advance model_run gone | advance url unchecked');
      let refused = '';
      await grok.functions.call('ProjectPortalCore:addEvidence', {portfolio: 'EX-SMD', program: p.code, kind: 'url', target: 'nope'})
        .catch((e: any) => refused = String(e?.message ?? e));
      expect(refused.includes('is not an http(s) address'), true, refused);
    });
  });

  test('the gate panel suggests readiness, attaches a link, and the history lists each decision\'s evidence', async () => {
    await scratchAtTv(async (p, smd, stage) => {
      await setCriterion(p.id, (await assayReady(smd, stage)).id, 'met');
      const v = new PortalView();
      await v.load('EX-SMD');
      await v.showProgram('EX-SMD', p.id);
      const q = (name: string) => v.detailRoot.querySelector(`[data-u2-name="${name}"]`) as HTMLElement | null;
      // on a u2 input the name is on its root, and its value follows the editor's input events
      const type = (name: string, value: string) => {
        const editor = q(name)!.querySelector('input, textarea') as HTMLInputElement | HTMLTextAreaElement;
        editor.value = value;
        editor.dispatchEvent(new Event('input'));
      };
      expect(q('evidence-readiness')?.textContent, 'Readiness of Hit Finding: 1 criterion');
      expect((q('evidence-readiness')!.querySelector('input') as HTMLInputElement).checked, true);
      type('evidence-link-input', 'https://example.org/hf-review');
      q('add-evidence-link-button')!.click();
      expect(q('evidence-link-0')?.textContent, 'Link: example.org/…/hf-review');
      type('rationale-input', 'assay validated');
      q('advance-button')!.click();
      await DG.delay(100);
      for (let i = 0; i < 100 && !q('decision-evidence-0'); i++)
        await DG.delay(100);
      expect(q('decision-evidence-0')?.textContent?.includes('Link: example.org/…/hf-review'), true, q('transition-history')?.textContent ?? '');
      expect(q('decision-evidence-0')?.textContent?.includes('Criterion: '), true);
      expect(!!v.detailRoot.querySelector('[data-u2-name^="evidence-0-url-"]'), true);
    });
  });

  test('an advance by data records its readiness evidence the same way', async () => {
    await gpc3AtHitFinding(async (p, stage) => {
      await autoAdvance('EX-SMD', false);
      const t = (await portalDb.stageTransitions.query({filter: `program_id = "${p.id}"`}))
        .find((x) => x.to_stage_id === stage('H2L') && x.source === 'data')!;
      const rows = (await portalDb.evidences.query({filter: `stage_transition_id = "${t.id}"`}))
        .map((r) => `${r.kind}:${(r.context as any).code}:${r.role}:${(r.context as any).basis}`).sort();
      expect(rows.join(' '), 'criterion:potency:readiness:evaluated criterion:validated_series:readiness:evaluated');
    });
    expect((await selectedRunsOf('00000000-0000-0000-0000-000000000000')).length, 0);
  });
});
