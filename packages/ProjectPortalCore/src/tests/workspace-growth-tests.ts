import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {findOrphans, programsWithoutFolder, stageFolderPaths, StageFolders, validateStageFolders} from '../model/workspace-plan';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programContext, recordDecision, setCriterion, syncPortfolio} from '../portal';
import {buildMissingFolders, buildPortfolioWorkspace, findRoot, summarizeFolderGrowth, workspaceOrphans} from '../workspace';
import {recordedSync} from '../sync-runs';
import {programFolderName} from '../model/workspace-plan';
import {_package} from '../package-test';

const stages: StageFolders[] = [
  {code: 'A', ordinal: 1, folders: []},
  {code: 'B', ordinal: 2, folders: ['Biology/Assay']},
  {code: 'C', ordinal: 3, folders: ['Package', 'biology / assay']},
  {code: 'D', ordinal: 4, folders: ['Package/Handover']},
];

category('Portal: stage folders and orphans', () => {
  test('a program earns the folders of every stage it has reached', async () => {
    const at = (o: number | null) => stageFolderPaths(stages, o).map((p) => p.join('/')).join(',');
    expect(at(null), '', 'an unstaged program has none');
    expect(at(1), '');
    expect(at(2), 'Biology/Assay');
    expect(at(3), 'Biology/Assay,Package', 'the same folder spelled again is one folder');
    expect(at(4), 'Biology/Assay,Package,Package/Handover');
  });

  test('a stage folder must name spaces, shallowly, once', async () => {
    expect(validateStageFolders('X', ['A/B', 'C']).join(), '');
    expect(validateStageFolders('X', ['']).length, 1);
    expect(validateStageFolders('X', ['a/b/c/d']).length, 1);
    expect(validateStageFolders('X', ['A|B']).length, 1);
    expect(validateStageFolders('X', ['A/--']).length, 1);
    expect(validateStageFolders('X', ['Dmpk', 'DMPK']).length, 1);
  });

  test('orphans: a folder no program claims, and an ended program\'s folder; nothing else', async () => {
    const children = [{id: 's1', name: 'P1Alpha', friendlyName: 'P1 · Alpha'}, {id: 's2', name: 'Stray'},
      {id: 's3', name: 'P2', friendlyName: 'P2'}];
    const programs = [{id: 'p1', code: 'P1', name: 'Alpha', spaceId: 's1', status: 'active'},
      {id: 'p2', code: 'P2', name: 'P2', status: 'terminated'}];
    const o = findOrphans('Root', children, programs);
    expect(o.map((x) => `${x.kind}:${x.path}`).join(), 'park:Root/P2,unclaimed:Root/Stray');
  });

  test('a sync builds folders only for programs that have none and have not ended or been parked', async () => {
    const children = [{id: 's1', name: 'P1', friendlyName: 'P1'}, {id: 's4', name: 'P4Delta', friendlyName: 'P4 · Delta'}];
    const programs = [{id: 'p1', code: 'P1', name: 'P1', spaceId: 's1', status: 'active'},
      {id: 'p2', code: 'P2', name: 'P2', status: 'active'},
      {id: 'p3', code: 'P3', name: 'P3', status: 'completed'},
      {id: 'p4', code: 'P4', name: 'Delta', status: 'active'},
      {id: 'p5', code: 'P5', name: 'P5', spaceId: 'gone', status: 'active'},
      {id: 'p6', code: 'P6', name: 'P6', spaceId: 'sp', status: 'active'}];
    const w = programsWithoutFolder(programs, children, new Set(['sp']));
    expect(w.build.map((p) => p.code).join(), 'P2,P5', 'P1 by id and P4 by name have one; P5\'s is gone');
    expect(w.ended.map((p) => p.code).join(), 'P3');
    expect(w.parked.map((p) => p.code).join(), 'P6');
  });

  test('the sync report says what was built, and that an unbuilt workspace was left alone', async () => {
    const g = {portfolio: 'EX-SMD', dryRun: false, built: ['A', 'B'], folders: 9, failed: [], ended: [], parked: [], newPrograms: []};
    expect(summarizeFolderGrowth(g), 'EX-SMD folders: built for A, B (9 folder(s))');
    expect(summarizeFolderGrowth({...g, built: []}), '', 'nothing missing, nothing said');
    expect(summarizeFolderGrowth({...g, dryRun: true, built: [], newPrograms: ['C']}), 'EX-SMD folders: would build for the new program(s) C');
    expect(summarizeFolderGrowth({...g, notBuilt: 'workspace "R" is not built yet'}), 'EX-SMD folders: workspace "R" is not built yet');
  });
});

let ready: Promise<void> | undefined;
const setUp = () => ready ??= (async () => {
  await loadAllPortfolios(_package);
  for (const code of ['EX-SMD', 'EX-BIOD', 'EX-BIOC'])
    await syncPortfolio(code);
})();

const subspaces = async (id: string) =>
  await grok.dapi.spaces.id(id).children.ofTypes('Project', false).list({pageSize: 1000}) as DG.Project[];
const named = (spaces: DG.Project[], name: string) =>
  spaces.find((s) => (s.friendlyName ?? s.name).replace(/[^A-Za-z0-9]/g, '').toLowerCase() ===
    name.replace(/[^A-Za-z0-9]/g, '').toLowerCase());

/** Deletes a folder this test made, its children first. */
async function removeFolder(id: string): Promise<void> {
  for (const c of await subspaces(id))
    await removeFolder(c.id);
  const space = await grok.dapi.spaces.find(id).catch(() => null);
  if (space)
    await grok.dapi.spaces.delete(space);
}

category('Portal: workspace grows with the stage', () => {
  before(async () => { await setUp(); });

  test('a preview of a program at lead optimization lists the folders of the stages it passed', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const lo = (await portalDb.stages.getByKey({stage_model_id: sm.stage_model_id, code: 'LO'}))!;
    const code = `WSG-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: code,
      stage_model_id: sm.stage_model_id, stage_id: lo.id, source: 'manual'});
    try {
      const preview = await buildPortfolioWorkspace('EX-SMD', true);
      const base = `Example Small Molecule Discovery/${programFolderName({code, name: code})}/`;
      const mine = preview.entries.filter((e) => e.path.startsWith(base)).map((e) => e.path.slice(base.length)).sort();
      expect(mine.join(), 'Biology,Biology/Assay Development,Chemistry,CompChem,DMPK,DMPK/In vivo PK',
        'workstreams once, stage folders beneath them, nothing from candidate selection');
    } finally {
      await portalDb.programs.delete(ins.id);
    }
  });

  test('a registry-placed program gets its stage folders from the builder; a second build creates nothing', async () => {
    const first = await buildPortfolioWorkspace('EX-BIOC');
    expect(first.failed, 0, first.entries.filter((e) => e.status !== 'ok').map((e) => `${e.path}: ${e.status}`).join('; '));
    const bioc = await portfolioByCode('EX-BIOC');
    const alpha = (await portalDb.programs.query({filter: `portfolio_id = "${bioc.id}" and code = "EXB-201"`}))[0];
    const folder = await subspaces(alpha.space_id!);
    const regulatory = named(folder, 'Regulatory')!;
    expect(!!named(await subspaces(regulatory.id), 'IND'), true, 'D3 passed: Regulatory/IND');
    expect(!!named(await subspaces(named(folder, 'Drug Substance')!.id), 'Process Lock'), false, 'D5 not reached');
    const again = await buildPortfolioWorkspace('EX-BIOC');
    expect(again.created, 0, again.entries.filter((e) => e.action === 'create').map((e) => e.path).join());
  });

  test('advancing a program grows its workspace by the stage it entered', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const code = `WSA-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: code,
      stage_model_id: sm.stage_model_id, source: 'manual'});
    let spaceId: string | undefined;
    try {
      await buildPortfolioWorkspace('EX-SMD', false, ins.id);
      spaceId = (await portalDb.programs.get(ins.id)).space_id;
      expect(!!spaceId, true, 'the program has a folder');
      const kickoff = await recordDecision(ins.id, 'advance', null, 'kick-off');
      expect(kickoff.workspace === undefined, true, 'target validation adds no folder');
      const ctx = await programContext(ins.id);
      const hf = ctx.stages.find((s) => s.code === 'HF')!;
      await setCriterion(ins.id, ctx.criteria.find((c) => c.stageId === hf.id)!.id, 'met');
      const done = await recordDecision(ins.id, 'advance', null, 'assay validated');
      expect(done.workspace?.includes('created 1'), true, done.workspace);
      const biology = named(await subspaces(spaceId!), 'Biology')!;
      expect(!!named(await subspaces(biology.id), 'Assay Development'), true);
    } finally {
      if (spaceId)
        await removeFolder(spaceId);
      await portalDb.programs.delete(ins.id);
    }
  });

  test('a sync builds the folder of a program new since the last build, and of no other', async () => {
    await buildPortfolioWorkspace('EX-SMD');
    const sm = await portfolioByCode('EX-SMD');
    const root = (await findRoot(sm.workspace_root!))!;
    const lo = (await portalDb.stages.getByKey({stage_model_id: sm.stage_model_id, code: 'LO'}))!;
    const code = `SYF-${Date.now() % 1e9}`;
    const [made] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: code,
      stage_model_id: sm.stage_model_id, stage_id: lo.id, source: 'manual'});
    const [ended] = await portalDb.programs.insert({portfolio_id: sm.id, code: `${code}-T`, name: `${code}-T`,
      stage_model_id: sm.stage_model_id, status: 'terminated', source: 'manual'});
    let spaceId: string | undefined;
    let runId: string | null = null;
    try {
      const before = (await subspaces(root.id)).map((c) => c.id).sort();
      const preview = await buildMissingFolders('EX-SMD', true);
      expect(preview.built.join(), code, 'only the new program; the ended one gets none');
      expect(preview.ended.includes(`${code}-T`), true);
      expect((await subspaces(root.id)).length, before.length, 'a preview makes nothing');

      const r = await recordedSync('EX-SMD', 'manual');
      runId = r.runId;
      expect(r.error, null);
      expect(r.folders?.join(), code, r.report);
      expect(r.report.includes(`EX-SMD folders: built for ${code} (`), true, r.report);
      spaceId = (await portalDb.programs.get(made.id)).space_id;
      expect(!!spaceId, true, 'the folder is recorded on the program');
      const after = (await subspaces(root.id)).map((c) => c.id).filter((id) => id !== spaceId).sort();
      expect(after.join(), before.join(), 'no other folder made or removed under the root');
      expect(!!named(await subspaces(named(await subspaces(spaceId!), 'DMPK')!.id), 'In vivo PK'), true, 'stage folders came with it');
      expect((await portalDb.programs.get(ended.id)).space_id ?? null, null);
      if (runId)
        expect(((await portalDb.syncLogs.get(runId)).changes as any)?.folders?.join(), code, 'the sync history names it');
      expect((await buildMissingFolders('EX-SMD', true)).built.length, 0, 'a second sync builds nothing');
    } finally {
      if (spaceId)
        await removeFolder(spaceId);
      if (runId)
        await portalDb.syncLogs.delete(runId).catch(() => undefined);
      await portalDb.programs.delete(made.id);
      await portalDb.programs.delete(ended.id);
    }
  });

  test('the orphan report names what prune-or-park would look at, and touches nothing', async () => {
    await buildPortfolioWorkspace('EX-BIOC');
    const before = (await buildPortfolioWorkspace('EX-BIOC', true)).reused;
    const r = await workspaceOrphans('EX-BIOC');
    const bioc = await portfolioByCode('EX-BIOC');
    const ended = (await portalDb.programs.query({filter: `portfolio_id = "${bioc.id}"`}))
      .filter((p) => ['terminated', 'transitioned', 'missing_from_feed', 'completed'].includes(p.status ?? ''));
    expect(r.orphans.filter((o) => o.kind === 'park').length, ended.length, r.orphans.map((o) => o.path).join());
    expect((await buildPortfolioWorkspace('EX-BIOC', true)).reused, before, 'the report changed the workspace');
  });
});
