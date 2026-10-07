import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {FlowFoldersConfig, flowFolderPaths, flowFoldersDue, flowsForProgram, TableFlow, validateFlowFolders,
  validateFlowSet} from '../model/flow-folders';
import {noNewFoldersReason, programFolderName} from '../model/workspace-plan';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, syncPortfolio} from '../portal';
import {buildPortfolioWorkspace} from '../workspace';
import {filterValue, ROW_LIMIT} from '../rows';
import {_package} from '../package-test';

const config: FlowFoldersConfig = {
  placement: 'process_type',
  flows: [
    {code: 'DS1', name: 'DS one', processType: 'Drug Substance', steps: [
      {step: 2, unitOperation: 'PROD', folderName: '02 - Production'},
      {step: 1, unitOperation: 'SEED', folderName: '01 - Seed'}]},
    {code: 'DS2', name: 'DS two', processType: 'Drug Substance', steps: [{step: 1, unitOperation: 'SEED', folderName: '01 - Seed'}]},
    {code: 'DP1', name: 'DP one', processType: 'Drug Product', steps: [{step: 1, unitOperation: 'FILL', folderName: '01 - Fill'}]},
  ],
  rules: [
    {modality: 'mAb', flow: 'DS1', priority: 1},
    {flow: 'DS2'},
    {modality: 'Fusion protein', flow: 'DS2', priority: 5},
    {flow: 'DP1'},
  ],
};

const paths = (modality: string | null, table: TableFlow[] = [], placement?: 'process_type' | 'program') =>
  flowFolderPaths(flowsForProgram(modality, table, config), placement).map((p) => p.join('/')).join(',');

category('Portal: flow folder rules', () => {
  test('a declared flow is picked by modality, highest priority first, one per process type', async () => {
    expect(paths('MAB'), 'Drug Product/01 - Fill,Drug Substance/01 - Seed,Drug Substance/02 - Production',
      'mAb: DS1 over the any-modality DS2, steps in order');
    expect(flowsForProgram('fusion-protein', [], config).map((f) => f.code).join(), 'DP1,DS2', 'modality spelled loosely');
    expect(flowsForProgram(null, [], config).map((f) => f.code).join(), 'DP1,DS2', 'no modality: any-modality rules only');
    expect(flowsForProgram('mAb', [], {rules: [], flows: []}).length, 0);
    expect(flowsForProgram('mAb', [], null).length, 0);
  });

  test('a program_flow row wins its process type; a manual row wins over a rule row', async () => {
    const table: TableFlow[] = [
      {code: 'LOADED', processType: 'drug substance', rowSource: 'rule', steps: [{step: 1, folderName: '01 - Thaw'}]},
      {code: 'MINE', processType: 'Drug Substance', rowSource: 'manual', steps: [{step: 1, folderName: '01 - Vial'}]},
    ];
    const flows = flowsForProgram('mAb', table, config);
    expect(flows.map((f) => `${f.code}:${f.source}`).join(), 'DP1:portfolio,MINE:table');
    expect(paths('mAb', table, 'program'), '01 - Fill,01 - Vial', 'placement program: straight under the program');
  });

  test('flows from a stage on; folders shared by two flows listed once', async () => {
    expect(flowFoldersDue(undefined, null), true);
    expect(flowFoldersDue(2, 1), false);
    expect(flowFoldersDue(2, 2), true);
    expect(flowFoldersDue(2, null), false, 'an unstaged program has not reached it');
    const twice = flowFolderPaths([
      {code: 'A', processType: 'X', source: 'portfolio', steps: [{step: 1, folderName: '01 - Seed'}]},
      {code: 'B', processType: 'X', source: 'portfolio', steps: [{step: 1, folderName: '01 seed'}]},
    ]);
    expect(twice.length, 1);
  });

  test('an ended program gets no new folders; one on hold or active does', async () => {
    expect(noNewFoldersReason({code: 'P1', status: 'completed'}), 'P1 is completed: no new folders');
    expect(noNewFoldersReason({code: 'P2', status: 'terminated'}), 'P2 is terminated: no new folders');
    for (const status of ['active', 'on_hold', null])
      expect(noNewFoldersReason({code: 'P3', status}), null, String(status));
  });

  test('flowFolders are checked', async () => {
    expect(validateFlowFolders(config, ['D1'], 'Root').join('; '), '');
    expect(validateFlowFolders(undefined, [], undefined).length, 0);
    const bad: FlowFoldersConfig = {placement: 'deep' as any, fromStage: 'D9',
      flows: [{code: 'F 1', name: '', processType: 'A/B', steps: [
        {step: 1, unitOperation: '', folderName: 'x/y'}, {step: 1, unitOperation: 'U', folderName: '--'}]}],
      rules: [{flow: 'NOPE', priority: 1.5}]};
    const errors = validateFlowFolders(bad, ['D1'], undefined);
    for (const part of ['workspaceRoot', 'placement', 'fromStage', 'flow code', 'name is required', 'processType',
      'used once', 'unitOperation is required', "may not contain '/'", 'letter or digit', 'not a declared flow', 'whole number'])
      expect(errors.some((e) => e.includes(part)), true, `${part}: ${errors.join('; ')}`);
    expect(validateFlowFolders({flows: config.flows}, [], 'Root').some((e) => e.includes('no rule')), true);
    expect(validateFlowSet([{code: 'P', flowFolders: config, unitOperations: ['SEED', 'PROD']}]).join('; '),
      'P: flow DP1 step 1: unit operation FILL is not in the catalogue of any portfolio of the set');
  });
});

let ready: Promise<void> | undefined;
const setUp = () => ready ??= (async () => {
  await loadAllPortfolios(_package);
  await syncPortfolio('EX-BIOC');
})();

const ROOT = 'Example Biologics Development';
const DS = ['01 - Seed train', '02 - Production bioreactor', '03 - Harvest', '04 - Protein A capture', '05 - UF-DF'];

const subspaces = async (id: string) =>
  await grok.dapi.spaces.id(id).children.ofTypes('Project', false).list({pageSize: 1000}) as DG.Project[];
/** Deletes a folder this test made, its children first. */
async function removeFolder(id: string): Promise<void> {
  for (const c of await subspaces(id))
    await removeFolder(c.id);
  const space = await grok.dapi.spaces.find(id).catch(() => null);
  if (space)
    await grok.dapi.spaces.delete(space);
}

const named = (spaces: DG.Project[], name: string) =>
  spaces.find((s) => (s.friendlyName ?? s.name).replace(/[^A-Za-z0-9]/g, '').toLowerCase() ===
    name.replace(/[^A-Za-z0-9]/g, '').toLowerCase());

category('Portal: flow folders', () => {
  before(async () => { await setUp(); });

  test('a preview lists the unit-operation folders of the antibodies past D1, and of no other program', async () => {
    const preview = await buildPortfolioWorkspace('EX-BIOC', true);
    const flowPaths = preview.entries.map((e) => e.path).filter((p) => /\/\d\d - /.test(p));
    const base = `${ROOT}/${programFolderName({code: 'EXB-201', name: 'Antibody Alpha'})}/Drug Substance/`;
    for (const f of DS)
      expect(flowPaths.includes(base + f), true, `${base + f} in ${flowPaths.join(', ')}`);
    // scratch programs other tests may leave behind are not the example's
    const programsWith = [...new Set(flowPaths.map((p) => p.split('/')[1].split(' · ')[0]))].filter((c) => c.startsWith('EXB-')).sort();
    expect(programsWith.join(), 'EXB-201,EXB-204', 'D1 antibodies, the bispecific, the ADC and the fusion protein get none');
    // the note lists scratch programs too, possibly ahead of the example's
    expect(preview.notes.some((n) => n.startsWith('flow folders: ') && n.includes('EXB-201 (EX-MAB-DS)')), true,
      preview.notes.join('; '));
  });

  test('a build makes them under the Drug Substance workstream; a second build creates nothing', async () => {
    const first = await buildPortfolioWorkspace('EX-BIOC');
    expect(first.failed, 0, first.entries.filter((e) => e.status !== 'ok').map((e) => `${e.path}: ${e.status}`).join('; '));
    const bioc = await portfolioByCode('EX-BIOC');
    const alpha = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(bioc.id)} and code = "EXB-201"`, limit: ROW_LIMIT}))[0];
    const ds = named(await subspaces(alpha.space_id!), 'Drug Substance')!;
    const folders = (await subspaces(ds.id)).map((s) => s.friendlyName ?? s.name);
    for (const f of DS)
      expect(!!named(await subspaces(ds.id), f), true, `${f} in ${folders.join(', ')}`);
    const again = await buildPortfolioWorkspace('EX-BIOC');
    expect(again.created, 0, again.entries.filter((e) => e.action === 'create').map((e) => e.path).join());
  });

  test('a portfolio without flowFolders gets none, and reads no flow', async () => {
    const preview = await buildPortfolioWorkspace('EX-SMD', true);
    expect(preview.notes.some((n) => n.startsWith('flow folders')), false, preview.notes.join('; '));
    expect(preview.entries.some((e) => /\/\d\d - /.test(e.path)), false);
  });

  test('a program_flow row a loader wrote wins its process type', async () => {
    const bioc = await portfolioByCode('EX-BIOC');
    const d4 = (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code: 'D4'}))!;
    const op = (await portalDb.unitOperations.getByKey({code: 'EX-PROD'}))!;
    const stamp = Date.now() % 1e9;
    const code = `WSF-${stamp}`;
    const [flow] = await portalDb.processFlows.insert({code: `EX-TEST-FLOW-${stamp}`, name: 'Test flow',
      process_type: 'Drug Substance'});
    const [program] = await portalDb.programs.insert({portfolio_id: bioc.id, code, name: code,
      stage_model_id: bioc.stage_model_id, stage_id: d4.id, source: 'manual', modality: 'mAb'});
    try {
      await portalDb.processSteps.insert({process_flow_id: flow.id, unit_operation_id: op.id, step_number: 1,
        folder_name: '01 - Loaded step'});
      await portalDb.programFlows.insert({program_id: program.id, process_flow_id: flow.id, source: 'rule'});
      const preview = await buildPortfolioWorkspace('EX-BIOC', true, program.id);
      const base = `${ROOT}/${code}/`;
      const mine = preview.entries.filter((e) => e.path.startsWith(base)).map((e) => e.path.slice(base.length));
      expect(mine.includes('Drug Substance/01 - Loaded step'), true, mine.join());
      expect(mine.some((p) => p.includes('Seed train')), false, 'the declared flow was not overridden');
    } finally {
      await portalDb.programs.delete(program.id);
      await portalDb.processFlows.delete(flow.id);
    }
  });

  test('an ended program: what it has is reused, nothing new is made, the preview lists the rest as skipped', async () => {
    const bioc = await portfolioByCode('EX-BIOC');
    const stage = async (code: string) => (await portalDb.stages.getByKey({stage_model_id: bioc.stage_model_id, code}))!.id;
    const code = `WSE-${Date.now() % 1e9}`;
    const [program] = await portalDb.programs.insert({portfolio_id: bioc.id, code, name: code,
      stage_model_id: bioc.stage_model_id, stage_id: await stage('D1'), source: 'manual', modality: 'mAb', status: 'completed'});
    let spaceId: string | null = null;
    try {
      const none = await buildPortfolioWorkspace('EX-BIOC', false, program.id);
      expect(none.entries.map((e) => `${e.action} ${e.path}${e.reason ? ` (${e.reason})` : ''}`).join(),
        `reuse ${ROOT},skip ${ROOT}/${code} (${code} is completed: no new folders)`, 'the root only, which is there');
      expect((await portalDb.programs.get(program.id)).space_id ?? null, null, 'no folder, none recorded');

      // its folder from when it was active, at D1; then it ends at D4, past the flow folders' stage
      await portalDb.programs.update(program.id, {status: 'active'});
      const built = await buildPortfolioWorkspace('EX-BIOC', false, program.id);
      expect(built.failed, 0, built.entries.filter((e) => e.status !== 'ok').map((e) => `${e.path}: ${e.status}`).join('; '));
      spaceId = (await portalDb.programs.get(program.id)).space_id!;
      await portalDb.programs.update(program.id, {status: 'completed', stage_id: await stage('D4')});

      const preview = await buildPortfolioWorkspace('EX-BIOC', true, program.id);
      const skipped = preview.entries.filter((e) => e.action === 'skip');
      expect(preview.created, 0, preview.entries.filter((e) => e.action === 'create').map((e) => e.path).join());
      expect(preview.reused, built.created + built.reused, 'what it has is listed as there');
      expect(skipped.some((e) => /\/Drug Substance\/01 - /.test(e.path)), true, skipped.map((e) => e.path).join());
      expect(skipped.every((e) => e.reason === `${code} is completed: no new folders`), true);
      expect(preview.skipped, skipped.length);
      const again = await buildPortfolioWorkspace('EX-BIOC', false, program.id);
      expect(again.created, 0, again.entries.filter((e) => e.action === 'create').map((e) => e.path).join());
      expect(again.skipped, preview.skipped);
      const folders = await subspaces(spaceId);
      expect(!!named(await subspaces(named(folders, 'Drug Substance')!.id), DS[0]), false, 'no unit-operation folder made');
    } finally {
      if (spaceId)
        await removeFolder(spaceId).catch(() => {});
      await portalDb.programs.delete(program.id);
    }
  });
});
