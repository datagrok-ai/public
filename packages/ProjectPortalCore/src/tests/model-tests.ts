import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {PortfolioConfig, StageModelConfig, validatePortfolioSet, validateStageModel} from '../model/config';
import {planSync, PlanInput, ProgramRef, resolveStage, StageRef} from '../model/sync-plan';
import {checkDecision, GateInput, GateStage} from '../model/gates';
import {findProgramFolder, programFolderName} from '../model/workspace-plan';
import {readPortfolioFiles} from '../portal';
import {_package} from '../package-test';

const stages = (codes: string[]): StageRef[] =>
  codes.map((c, i) => ({id: `s-${c}`, code: c, name: `Stage ${c}`, ordinal: i + 1, aliases: [`SG${c}`]}));

const model = (code: string, stageCodes: string[]): StageModelConfig => ({code, name: code, setting: 'discovery',
  stages: stageCodes.map((c, i) => ({code: c, name: `${code} ${c}`, ordinal: i + 1}))});

const portfolio = (code: string, m: StageModelConfig, extra: Partial<PortfolioConfig> = {}): PortfolioConfig => ({
  code, name: code, setting: 'discovery', stageAuthority: 'portal', stageModel: m, feeds: [], ...extra});

const plan = (over: Partial<PlanInput>) => planSync({portfolioCode: 'P', system: 'sys', stageAuthority: 'feed',
  stages: stages(['1', '2', '6a']), existing: [], elsewhere: [], rows: [], ...over});

category('Portal: configuration', () => {
  test('the shipped portfolios load as one set', async () => {
    const configs = await readPortfolioFiles(_package);
    expect(configs.map((c) => c.code).sort().join(','), 'EX-BIOC,EX-BIOD,EX-SMD');
    expect(validatePortfolioSet(configs).join('\n'), '');
  });

  test('stage codes are scoped to their model', async () => {
    // 'LO' is lead optimization in both discovery models; that must not collide
    const errors = validatePortfolioSet([
      portfolio('A', model('M-A', ['TV', 'LO'])), portfolio('B', model('M-B', ['TI', 'LO']))]);
    expect(errors.length, 0, errors.join('; '));
  });

  test('one spelling may not name two stages of a model', async () => {
    const m = model('M', ['H2L', 'LO']);
    m.stages[1].aliases = ['h2l'];
    expect(validateStageModel(m).some((e) => e.includes('would resolve to both')), true);
  });

  test('a shared stage model must be declared identically', async () => {
    const a = model('SHARED', ['A', 'B']);
    const b = model('SHARED', ['A', 'B', 'C']);
    const errors = validatePortfolioSet([portfolio('P1', a), portfolio('P2', b)]);
    expect(errors.some((e) => e.includes('declared differently')), true, errors.join('; '));
    expect(validatePortfolioSet([portfolio('P1', a), portfolio('P2', model('SHARED', ['A', 'B']))]).length, 0);
  });

  test('handoff must leave from a terminal stage to a configured portfolio', async () => {
    const m = model('M', ['A', 'B']);
    const errors = validatePortfolioSet([portfolio('P', m, {handoff: {to: 'NOPE', fromStage: 'B'}})]);
    expect(errors.some((e) => e.includes('must be a terminal stage')), true);
    expect(errors.some((e) => e.includes('not configured')), true);
  });

  test('feed authority needs a feed', async () => {
    const errors = validatePortfolioSet([portfolio('P', model('M', ['A']), {stageAuthority: 'feed'})]);
    expect(errors.some((e) => e.includes('needs at least one feed')), true);
  });
});

category('Portal: sync plan', () => {
  test('feed stage resolves by code, name or alias; case and spacing ignored', async () => {
    const s = stages(['1', '6a']);
    expect(resolveStage(s, '6A')?.code, '6a');
    expect(resolveStage(s, 'SG 6a')?.code, '6a');
    expect(resolveStage(s, 'stage 1')?.code, '1');
    expect(resolveStage(s, 'Full Development'), null);
  });

  test('a new program lands in its stage and logs a sync transition', async () => {
    const p = plan({rows: [{code: 'X', name: 'X', stage: 'SG2', status: 'Active'}],
      statusMap: {Active: 'active'}});
    expect(p.inserts.length, 1);
    expect(p.inserts[0].stage_id, 's-2');
    expect(p.inserts[0].external_ref, 'sys:X');
    expect(p.stageChanges.length, 1);
  });

  test('portal authority keeps the registry stage as text and does not place the program', async () => {
    // a discovery portfolio reading a development registry: its gate code means nothing in the discovery model
    const p = plan({stageAuthority: 'portal', rows: [{code: 'X', name: 'X', source_stage: '3'}]});
    expect(p.inserts[0].stage_id, null);
    expect(p.inserts[0].source_stage, '3');
    expect(p.stageChanges.length, 0);
    expect(p.unresolvedStages.length, 0);
  });

  test('an unresolvable spelling is reported and never clears a known stage', async () => {
    const existing: ProgramRef[] = [{id: 'p1', code: 'X', portfolioCode: 'P', stageId: 's-2', status: 'active',
      source: 'feed', externalRef: 'sys:X', name: 'X', sourceStage: 'Full Development'}];
    const p = plan({existing, rows: [{code: 'X', name: 'X', stage: 'Full Development'}]});
    expect(p.unresolvedStages.map((u) => u.spelling).join(), 'Full Development');
    expect(p.updates.length, 0, JSON.stringify(p.updates));
    expect(p.unchanged.join(), 'X');
  });

  test('a stage move in the registry becomes an update and a transition', async () => {
    const existing: ProgramRef[] = [{id: 'p1', code: 'X', portfolioCode: 'P', stageId: 's-1', status: 'active',
      source: 'feed', externalRef: 'sys:X', name: 'X', sourceStage: '1'}];
    const p = plan({existing, rows: [{code: 'X', name: 'X', stage: '2'}]});
    expect(p.updates[0].changes.stage_id, 's-2');
    expect(p.stageChanges[0].fromStageId, 's-1');
  });

  test('a program arriving from another portfolio links back as its successor', async () => {
    const elsewhere: ProgramRef[] = [{id: 'disc-1', code: 'X', portfolioCode: 'EX-BIOD', status: 'completed',
      source: 'feed', externalRef: 'sys:X'}];
    const p = plan({elsewhere, rows: [{code: 'X', name: 'X', stage: '1'}]});
    expect(p.inserts[0].predecessor_id, 'disc-1');
    expect(p.lineage.length, 1);
    expect(p.overlaps.length, 0);
  });

  test('two feeds delivering the same key at once is an overlap', async () => {
    const elsewhere: ProgramRef[] = [{id: 'disc-1', code: 'X', portfolioCode: 'EX-BIOD', status: 'active',
      source: 'feed', externalRef: 'sys:X'}];
    const p = plan({elsewhere, rows: [{code: 'X', name: 'X', stage: '1'}]});
    expect(p.overlaps.map((o) => o.otherPortfolio).join(), 'EX-BIOD');
  });

  test('a program gone from its feed is transitioned when another portfolio has it, else flagged', async () => {
    const existing: ProgramRef[] = [
      {id: 'a', code: 'MOVED', portfolioCode: 'P', status: 'active', source: 'feed', externalRef: 'sys:MOVED'},
      {id: 'b', code: 'GONE', portfolioCode: 'P', status: 'active', source: 'feed', externalRef: 'sys:GONE'},
      {id: 'c', code: 'MINE', portfolioCode: 'P', status: 'active', source: 'manual'}];
    const elsewhere: ProgramRef[] = [{id: 'z', code: 'MOVED', portfolioCode: 'EX-BIOC', status: 'active',
      source: 'feed', externalRef: 'sys:MOVED', predecessorId: 'a'}];
    const p = plan({existing, elsewhere, rows: []});
    const byCode = Object.fromEntries(p.missing.map((m) => [m.code, m.status]));
    expect(byCode['MOVED'], 'transitioned');
    expect(byCode['GONE'], 'missing_from_feed');
    expect('MINE' in byCode, false, 'manually created programs are not the feed\'s to flag');
  });

  test('dirty rows: blank and NULL codes rejected, duplicates rejected, unknown status warned', async () => {
    const p = plan({rows: [{code: ' '}, {code: 'NULL'}, {code: 'A', status: 'Weird'}, {code: 'A'}],
      statusMap: {Active: 'active'}});
    expect(p.rejected.length, 3);
    expect(p.inserts.length, 1);
    expect(p.inserts[0].status, 'active');
    expect(p.warnings.length, 1);
  });

  test('unmapped columns are kept as details', async () => {
    const p = plan({rows: [{code: 'A', molecule_code: 'MOL-1', registries: 'x,y'}]});
    expect(p.inserts[0].details?.molecule_code, 'MOL-1');
  });
});

category('Portal: gate rules', () => {
  const gs: GateStage[] = [
    {id: 'tv', code: 'TV', name: 'Target Validation', ordinal: 1},
    {id: 'hf', code: 'HF', name: 'Hit Finding', ordinal: 2},
    {id: 'h2l', code: 'H2L', name: 'Hit-to-Lead', ordinal: 3, isTerminal: true}];
  const base = (over: Partial<GateInput>): GateInput => ({decision: 'advance', stageAuthority: 'portal',
    programStatus: 'active', stages: gs, currentStageId: 'tv', criteria: [
      {id: 'c1', stageId: 'hf', code: 'assay', name: 'Assay', required: true},
      {id: 'c2', stageId: 'hf', code: 'nice', name: 'Nice to have', required: false}],
    criterionStates: {}, isApprover: true, ...over});

  test('advance needs the target stage\'s required criteria', async () => {
    const v = checkDecision(base({}));
    expect(v.allowed, false);
    expect(v.outstanding.join(), 'assay');
    expect(checkDecision(base({criterionStates: {c1: 'met'}})).allowed, true);
    expect(checkDecision(base({criterionStates: {c1: 'waived'}})).toStageId, 'hf');
    expect(checkDecision(base({criterionStates: {c1: 'not_met'}})).allowed, false);
  });

  test('advance moves one stage at a time', async () => {
    const v = checkDecision(base({targetStageId: 'h2l', criterionStates: {c1: 'met'}}));
    expect(v.allowed, false);
    expect(v.reasons.some((r) => r.includes('one stage at a time')), true);
  });

  test('an unstaged program advances into the first stage', async () => {
    const v = checkDecision(base({currentStageId: null}));
    expect(v.toStageId, 'tv');
    expect(v.allowed, true);
  });

  test('a stage of another model is unreachable', async () => {
    const v = checkDecision(base({targetStageId: 'ex-5', criterionStates: {c1: 'met'}}));
    expect(v.allowed, false);
    expect(v.reasons.some((r) => r.includes('different stage model')), true);
  });

  test('where the registry owns the stage, the portal records nothing', async () => {
    const v = checkDecision(base({stageAuthority: 'feed', criterionStates: {c1: 'met'}}));
    expect(v.allowed, false);
    expect(v.reasons.some((r) => r.includes('registry owns')), true);
  });

  test('only approvers decide', async () => {
    expect(checkDecision(base({isApprover: false, criterionStates: {c1: 'met'}})).allowed, false);
  });

  test('the terminal stage hands off rather than advances', async () => {
    const v = checkDecision(base({currentStageId: 'h2l'}));
    expect(v.reasons.some((r) => r.includes('hand it off')), true);
  });

  test('recycle goes back with a rationale; hold and terminate need one too', async () => {
    expect(checkDecision(base({decision: 'recycle', currentStageId: 'hf', targetStageId: 'tv'})).allowed, false);
    expect(checkDecision(base({decision: 'recycle', currentStageId: 'hf', targetStageId: 'tv',
      rationale: 'assay failed'})).allowed, true);
    expect(checkDecision(base({decision: 'recycle', currentStageId: 'tv', targetStageId: 'hf',
      rationale: 'x'})).allowed, false);
    expect(checkDecision(base({decision: 'hold'})).allowed, false);
    expect(checkDecision(base({decision: 'terminate', rationale: 'deprioritized'})).allowed, true);
  });

  test('a terminated program takes no decisions', async () => {
    expect(checkDecision(base({programStatus: 'terminated', criterionStates: {c1: 'met'}})).allowed, false);
  });
});

category('Portal: workspace plan', () => {
  test('a program folder leads with its code', async () => {
    expect(programFolderName({code: 'EXS-KIN1', name: 'Kinase inhibitor'}), 'EXS-KIN1 · Kinase inhibitor');
    expect(programFolderName({code: 'EXS-KIN1', name: 'EXS-KIN1'}), 'EXS-KIN1');
    expect(programFolderName({code: 'EXS-KIN1', name: ' '}), 'EXS-KIN1');
  });

  test('the recorded folder wins over the name; the name is matched as the server stores it', async () => {
    const program = {id: 'p', code: 'EXS-KIN1', name: 'Kinase inhibitor', spaceId: 's1'};
    const renamed = {id: 's1', name: 'EXSKIN1Oldname', friendlyName: 'EXS-KIN1 · Old name'};
    const byName = {id: 's2', name: 'EXSKIN1Kinaseinhibitor', friendlyName: null};
    expect(findProgramFolder(program, [byName, renamed])?.space.id, 's1');
    expect(findProgramFolder({...program, spaceId: null}, [byName, renamed])?.space.id, 's2');
    expect(findProgramFolder({...program, spaceId: 'gone'}, [renamed]), null);
  });

  test('workstreams must name distinct spaces under a root', async () => {
    const m = model('M', ['A']);
    const errs = (extra: Partial<PortfolioConfig>) => validatePortfolioSet([portfolio('P', m, extra)]).join('; ');
    expect(errs({workspaceRoot: 'R', workstreams: ['Chemistry', 'DMPK']}), '');
    expect(errs({workspaceRoot: 'R', workstreams: ['Comp Chem', 'CompChem']}).includes('declared twice'), true);
    expect(errs({workspaceRoot: 'R', workstreams: ['--']}).includes('no letters or digits'), true);
    expect(errs({workstreams: ['Chemistry']}).includes('need a workspaceRoot'), true);
  });
});
