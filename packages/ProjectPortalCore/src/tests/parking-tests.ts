import * as grok from 'datagrok-api/grok';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {checkParking, checkUnpark, isParkingFolder, parkDay, planParking, summarizeParking} from '../model/parking';
import {findOrphans, Orphan} from '../model/workspace-plan';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {buildPortfolioWorkspace, findRoot, subspaces, summarizeWorkspace, workspaceOrphans} from '../workspace';
import {parkedFoldersOf, parkOrphans, summarizeParkResult, unparkFolder} from '../parking';
import {withRestrictedUser} from './restricted-user';
import {_package} from '../package-test';

const orphans: Orphan[] = [
  {id: 'a', path: 'Root/Old notes', kind: 'unclaimed', reason: 'no program of the portfolio claims it'},
  {id: 'b', path: 'Root/EXS-ION4 · Ion channel', kind: 'park', reason: 'EXS-ION4 is terminated'},
];

category('Portal: parking rules', () => {
  test('parking moves each orphan into _Parked/<day> under the root; one already parked is skipped', async () => {
    const plan = planParking('Root', orphans, '2026-10-02', ['a']);
    expect(plan.moves.map((m) => `${m.from} → ${m.to}`).join('; '),
      'Root/EXS-ION4 · Ion channel → Root/_Parked/2026-10-02/EXS-ION4 · Ion channel');
    expect(JSON.stringify(plan.skipped), '[{"path":"Root/Old notes","reason":"already parked"}]');
    expect(planParking('Root', orphans, '2026-10-02', [], ['Root/Old notes']).moves.map((m) => m.id).join(), 'a', 'only by path');
    expect(planParking('Root', orphans, '2026-10-02', [], ['b']).moves.map((m) => m.id).join(), 'b', 'only by id');
    expect(parkDay(new Date('2026-10-02T23:30:00Z')), '2026-10-02');
  });

  test('a dry run is anyone\'s; parking and unparking are the approvers\'', async () => {
    expect(checkParking({isApprover: false, dryRun: true}).length, 0);
    expect(checkParking({isApprover: false, dryRun: false})[0], 'Only the portfolio\'s approvers park or unpark folders; preview with dryRun');
    expect(checkParking({isApprover: true, dryRun: false}).length, 0);
  });

  test('a folder goes back only while parked, and not over a folder of the same name', async () => {
    const p = {id: 'r', name: 'Old notes', status: 'parked' as const, fromPath: 'Root/Old notes', parkedPath: 'Root/_Parked/2026-10-02/Old notes'};
    expect(checkUnpark(p, []).length, 0);
    expect(checkUnpark(p, [{id: 'x', name: 'Oldnotes', friendlyName: 'Old notes'}])[0],
      'a folder named Old notes is back under the root; move or rename it first');
    expect(checkUnpark({...p, status: 'restored'}, [])[0], 'Old notes is not parked (restored)');
  });

  test('the parking folder is never an orphan itself, under either spelling', async () => {
    expect(isParkingFolder({id: '1', name: 'Parked'}), true);
    expect(isParkingFolder({id: '1', name: 'x', friendlyName: '_Parked'}), true);
    const found = findOrphans('Root', [{id: '1', name: 'Parked', friendlyName: '_Parked'}, {id: '2', name: 'Oldnotes', friendlyName: 'Old notes'}], []);
    expect(found.map((o) => o.path).join(), 'Root/Old notes');
    expect(found[0].id, '2', 'the report carries the folder\'s id');
  });

  test('the summary says where each folder goes and that nothing is deleted', async () => {
    const plan = planParking('Root', orphans, '2026-10-02');
    expect(summarizeParking('EX-SMD', plan, true).split('\n')[0],
      'EX-SMD workspace "Root": 2 folder(s) would be parked in _Parked/2026-10-02; nothing deleted');
    expect(summarizeParking('EX-SMD', plan, false, ['x: not parked']).split('\n').slice(-1)[0], '  ! x: not parked');
  });
});

/** Every folder under [id] as indented names, to compare a tree before and after. */
async function tree(id: string, depth = 2, indent = ''): Promise<string[]> {
  const out: string[] = [];
  for (const c of (await subspaces(id)).sort((a, b) => (a.friendlyName || a.name).localeCompare(b.friendlyName || b.name))) {
    out.push(`${indent}${c.friendlyName || c.name}`);
    if (depth > 1)
      out.push(...await tree(c.id, depth - 1, `${indent}  `));
  }
  return out;
}

/** Deletes a folder this test made, its children first. */
async function removeFolder(id: string): Promise<void> {
  for (const c of await subspaces(id))
    await removeFolder(c.id);
  const space = await grok.dapi.spaces.find(id).catch(() => null);
  if (space)
    await grok.dapi.spaces.delete(space);
}

category('Portal: park orphans', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  test('an unclaimed folder on EX-SMD is previewed, parked, left out of the report and the build, unparked; the tree is as before', async () => {
    const smd = await portfolioByCode('EX-SMD');
    const root = (await findRoot(smd.workspace_root!))!;
    const before = (await tree(root.id)).join('\n');
    const name = `TEST-PARK-${Date.now() % 1e8}`;
    const made = await grok.dapi.spaces.id(root.id).addSubspace(name);
    const parkedRows: string[] = [];
    try {
      const path = `${smd.workspace_root}/${made.friendlyName || made.name}`;
      const preview = await parkOrphans('EX-SMD', true, [made.id]);
      expect(preview.plan.moves.length, 1, summarizeParkResult(preview));
      expect(preview.plan.moves[0].to, `${smd.workspace_root}/_Parked/${parkDay(new Date())}/${made.friendlyName || made.name}`);
      expect((await subspaces(root.id)).some((c) => c.id === made.id), true, 'a preview moves nothing');

      const done = await parkOrphans('EX-SMD', false, [made.id]);
      expect(done.errors.join('; '), '');
      expect(done.parked.length, 1);
      parkedRows.push(done.parked[0].id);
      expect(`${done.parked[0].from_path} ${done.parked[0].status}`, `${path} parked`);
      expect((await subspaces(root.id)).some((c) => c.id === made.id), false, 'moved off the root');
      expect(!!(await grok.dapi.spaces.find(made.id).catch(() => null)), true, 'moved, not deleted');
      const report = await workspaceOrphans('EX-SMD');
      expect(report.orphans.some((o) => o.id === made.id || o.path.includes('Parked')), false,
        report.orphans.map((o) => o.path).join());
      const build = summarizeWorkspace(await buildPortfolioWorkspace('EX-SMD', true));
      expect(build.includes('Parked'), false, build);
      expect((await parkOrphans('EX-SMD', true, [made.id])).plan.moves.length, 0, 'not offered again');

      const back = await unparkFolder(done.parked[0].id);
      expect(back.line.endsWith(`${path}: unparked`), true, back.line);
      expect(back.notes.join('; ').includes('_Parked was empty: removed'), true, back.notes.join('; '));
      expect((await subspaces(root.id)).some((c) => c.id === made.id), true, 'back under the root');
      expect((await portalDb.parkedFolders.get(done.parked[0].id)).status, 'restored');
    } finally {
      await removeFolder(made.id);
      for (const id of parkedRows)
        await portalDb.parkedFolders.delete(id).catch(() => undefined);
    }
    expect((await tree(root.id)).join('\n'), before, 'EX-SMD\'s tree is as it was');
  });

  test('the folder of an ended program is parked, not rebuilt by the builder, and unparked', async () => {
    const smd = await portfolioByCode('EX-SMD');
    const code = `TEST-PK-${Date.now() % 1e8}`;
    const [r] = await portalDb.programs.insert({portfolio_id: smd.id, code, name: code, stage_model_id: smd.stage_model_id,
      source: 'manual'});
    let spaceId: string | null = null;
    const parkedRows: string[] = [];
    try {
      await buildPortfolioWorkspace('EX-SMD', false, r.id);
      spaceId = (await portalDb.programs.get(r.id)).space_id ?? null;
      expect(!!spaceId, true, 'the program has a folder');
      await portalDb.programs.update(r.id, {status: 'terminated'});
      const done = await parkOrphans('EX-SMD', false, [spaceId!]);
      expect(`${done.parked.length} ${done.parked[0]?.kind} ${done.parked[0]?.reason}`, `1 park ${code} is terminated`);
      parkedRows.push(done.parked[0].id);
      const build = summarizeWorkspace(await buildPortfolioWorkspace('EX-SMD', true));
      expect(build.includes(`+ ${smd.workspace_root}/${code}`), false, build);
      expect(build.includes(`${code}: folder parked`), true, build);
      await unparkFolder(spaceId!);
      expect((await workspaceOrphans('EX-SMD')).orphans.some((o) => o.id === spaceId), true, 'an orphan again, until someone acts');
    } finally {
      if (spaceId)
        await removeFolder(spaceId);
      for (const id of parkedRows)
        await portalDb.parkedFolders.delete(id).catch(() => undefined);
      await portalDb.programs.delete(r.id);
    }
  });

  test('someone who is not an approver may preview but not park', async () => {
    const ran = await withRestrictedUser('pppark', async (user) => {
      let refused = '';
      try {
        await user.asUser(() => parkOrphans('EX-SMD', false));
      } catch (e: any) {
        refused = e?.message ?? String(e);
      }
      expect(refused, 'Only the portfolio\'s approvers park or unpark folders; preview with dryRun');
      return true;
    });
    if (ran === null)
      console.log('skipped: no self-signup on this stand');
    const smd = await portfolioByCode('EX-SMD');
    expect((await parkedFoldersOf(smd)).length, 0, 'nothing is left parked on EX-SMD');
  });
});

