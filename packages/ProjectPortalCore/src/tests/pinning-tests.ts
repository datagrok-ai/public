import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {assembleTiles, checkGroupPin, describeTile, PinHolders, pinsOf, programAddress, TileInput, tileTitle, unseenNews,
  whyNotPinnable} from '../model/pinning';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {myPins, pinByCode} from '../pinning';
import {pinnedFrame, pinnedTiles} from '../pinned-programs';
import {COL, portfolioTables} from '../portfolio-table';
import {_package} from '../package-test';

const stages = [{code: 'H2L', name: 'Hit to lead'}, {code: 'LO', name: 'Lead optimization'}, {code: 'CS', name: 'Candidate'}];
const input = (over: Partial<TileInput> = {}): TileInput => ({
  programId: 'p1', code: 'EXS-KIN1', name: 'Kinase', spaceId: 's1', portfolioCode: 'EX-SMD',
  portfolioName: 'Example Small Molecule Discovery', status: 'active', stages, stage: 'LO',
  next: {code: 'CS', name: 'Candidate'}, readiness: {score: 0.75, requiredMet: 3, requiredTotal: 4, gaps: ['hERG missing']},
  attention: 2, follows: [], ...over,
});
const mine = (...ids: string[]) => new Map<string, PinHolders>(ids.map((id) => [id, {mine: true, groups: []}]));

category('Portal: pinning rules', () => {
  test('only a program with a live folder can be pinned, and it says why not', async () => {
    expect(whyNotPinnable({spaceId: 's', spaceLive: true}), null);
    expect(whyNotPinnable({spaceId: null, spaceLive: false})!.startsWith('It has no workspace folder yet'), true);
    expect(whyNotPinnable({spaceId: 's', spaceLive: false})!.startsWith('Its workspace folder is gone'), true);
  });

  test('pinning for a group is the approvers\'', async () => {
    expect(checkGroupPin({isApprover: true}).length, 0);
    expect(checkGroupPin({isApprover: false})[0], 'Only the portfolio\'s approvers pin a program for a group');
  });

  test('pins are the folders among the user\'s and their groups\' favorites', async () => {
    const pins = pinsOf([{id: 's1', isSpace: true}, {id: 'q', isSpace: false}],
      [{group: 'Chemists', entities: [{id: 's1', isSpace: true}, {id: 's2', isSpace: true}]},
        {group: 'Biology', entities: [{id: 's2', isSpace: true}]}]);
    expect([...pins.keys()].sort().join(), 's1,s2', 'a favorite that is not a folder is no pin');
    expect(JSON.stringify(pins.get('s1')), '{"mine":true,"groups":["Chemists"]}');
    expect(JSON.stringify(pins.get('s2')), '{"mine":false,"groups":["Chemists","Biology"]}');
  });

  test('a tile: place on the track, readiness for the next gate, attention, its address', async () => {
    const [t] = assembleTiles([input()], mine('s1'));
    expect(t.stage, 'LO · Lead optimization');
    expect(`${t.step}/${t.steps}`, '2/3');
    expect(t.track.map((s) => s.state).join(), 'done,current,ahead');
    expect(t.next, 'CS · Candidate');
    expect(t.readinessLine, '3 of 4 required met');
    expect(t.ready, false);
    expect(t.address, '/apps/ProjectPortalCore/Portal/EX-SMD/EXS-KIN1');
    expect(describeTile(t), 'EXS-KIN1 · Kinase (EX-SMD) · LO · Lead optimization [2/3] · active · CS · Candidate: 75% · 2 need attention');
    expect(tileTitle({code: 'X-1', name: 'X-1'}), 'X-1', 'a name that is the code is said once');
    expect(programAddress('EX SMD', 'A/B'), '/apps/ProjectPortalCore/Portal/EX%20SMD/A%2FB');
  });

  test('tiles follow the pins exactly, sorted by portfolio then code', async () => {
    const a = input({programId: 'a', code: 'B-2', spaceId: 'sa', portfolioName: 'Beta'});
    const b = input({programId: 'b', code: 'A-1', spaceId: 'sb', portfolioName: 'Beta'});
    const c = input({programId: 'c', code: 'Z-9', spaceId: 'sc', portfolioName: 'Alpha'});
    const none = input({programId: 'd', code: 'N-0', spaceId: null});
    expect(assembleTiles([a, b, c, none], mine('sa', 'sb', 'sc')).map((t) => t.code).join(), 'Z-9,A-1,B-2');
    expect(assembleTiles([a, b, c], mine('sb')).map((t) => t.code).join(), 'A-1', 'unpinned programs have no tile');
    const group = assembleTiles([a], new Map([['sa', {mine: false, groups: ['Chemists']}]]))[0];
    expect(`${group.pinnedByMe} ${group.pinnedFor.join()}`, 'false Chemists', 'a group\'s pin gives a tile too');
  });

  test('unstaged, at the last stage, ready', async () => {
    const [u] = assembleTiles([input({stage: null, next: stages[0], readiness: null})], mine('s1'));
    expect(`${u.stage} ${u.step} ${u.track.every((s) => s.state === 'ahead')}`, 'Unstaged 0 true');
    const [last] = assembleTiles([input({stage: 'CS', next: null, readiness: null})], mine('s1'));
    expect(`${last.readinessLine} ${last.readiness}`, 'No further gate null');
    const [ready] = assembleTiles([input({readiness: {score: 1, requiredMet: 4, requiredTotal: 4, gaps: []}})], mine('s1'));
    expect(ready.ready, true);
  });

  test('new data is news that landed after the user\'s last visit', async () => {
    expect(unseenNews([]), null);
    expect(unseenNews([{newsOn: '2026-10-01T10:00:00Z', seenOn: '2026-10-02T09:00:00Z', lastNews: 'old'}]), null, 'seen since');
    expect(unseenNews([{newsOn: '2026-10-02T10:00:00Z', seenOn: '2026-10-02T09:00:00Z', lastNews: 'EXS-KIN1: 3 new rows'},
      {newsOn: '2026-10-02T08:00:00Z', seenOn: null, lastNews: 'earlier'}]), 'EXS-KIN1: 3 new rows', 'the latest unseen');
    const [t] = assembleTiles([input({follows: [{newsOn: '2026-10-02T10:00:00Z', seenOn: null, lastNews: 'x'}]})], mine('s1'));
    expect(`${t.followed} ${t.newData}`, 'true x');
  });
});

/** The program used here; EXS-KIN1 is the BDD scenario's. */
const CODE = 'EXS-PRO2';
const SCRATCH = 'EXS-PINX';

category('Portal: pinned programs', () => {
  let program: ProgramRow;
  let space: DG.Project;
  let viewers: DG.Group;
  /** How the favorites stood, to put them back. */
  let had = {mine: false, group: false};

  const favorite = async (group?: DG.Group) =>
    (await grok.dapi.entities.getFavorites(group ?? DG.User.current().group)).some((e) => e.id === space.id);

  before(async () => {
    await loadAllPortfolios(_package);
    const pf = await portfolioByCode('EX-SMD');
    program = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "${CODE}"`}))[0];
    space = await grok.dapi.spaces.find(program.space_id!);
    viewers = await grok.dapi.groups.find(pf.viewers!);
    had = {mine: await favorite(), group: await favorite(viewers)};
  });

  after(async () => {
    for (const [on, group] of [[had.mine, undefined], [had.group, viewers]] as [boolean, DG.Group | undefined][]) {
      if (on !== await favorite(group))
        await (on ? DG.Favorites.add(space, group) : DG.Favorites.remove(space, group));
    }
    for (const r of await portalDb.programs.query({filter: `code = "${SCRATCH}"`}))
      await portalDb.programs.delete(r.id);
  });

  test('pin and unpin: a favorite on the folder, a tile, the table\'s pinned column', async () => {
    await pinByCode('EX-SMD', CODE, false);
    expect(await favorite(), false);
    expect(await pinByCode('EX-SMD', CODE, true), `${CODE}: pinned; its folder is in the Workspace and its tile on Home`);
    expect(await favorite(), true, 'the platform favorite is the pin');
    expect(await pinByCode('EX-SMD', CODE, true), `${CODE}: already pinned`);
    const tile = (await pinnedTiles()).find((t) => t.code === CODE);
    expect(tile?.portfolioCode, 'EX-SMD');
    expect(tile?.pinnedByMe, true);
    const df = pinnedFrame(await pinnedTiles());
    expect(df.col('program')!.toList().includes(CODE), true);
    const table = (await portfolioTables('EX-SMD')).programs;
    const row = table.col(COL.program)!.toList().indexOf(CODE);
    expect(table.get(COL.pinned, row), true, 'pinned column');
    expect(await pinByCode('EX-SMD', CODE, false), `${CODE}: unpinned`);
    expect(await favorite(), false);
    expect((await pinnedTiles()).some((t) => t.code === CODE), false);
  });

  test('unpinning in the platform takes the tile away', async () => {
    await pinByCode('EX-SMD', CODE, true);
    expect((await pinnedTiles()).some((t) => t.code === CODE), true);
    await DG.Favorites.remove(space);   // what the Workspace's Unpin and the star do
    expect((await myPins()).get(space.id)?.mine ?? false, false);
    expect((await pinnedTiles()).some((t) => t.code === CODE), false);
  });

  test('an approver pins for a group: the group\'s favorite, a tile marked for it', async () => {
    const name = viewers.friendlyName || viewers.name;
    expect(await pinByCode('EX-SMD', CODE, true, name), `${CODE} for ${name}: pinned; its folder is in the Workspace and its tile on Home`);
    expect(await favorite(viewers), true);
    expect((await pinnedTiles()).find((t) => t.code === CODE)?.pinnedFor.includes(name), true);
    await pinByCode('EX-SMD', CODE, false, name);
    expect(await favorite(viewers), false);
  });

  test('a program without a folder cannot be pinned', async () => {
    const pf = await portfolioByCode('EX-SMD');
    await portalDb.programs.insert({portfolio_id: pf.id, code: SCRATCH, name: SCRATCH, stage_model_id: pf.stage_model_id,
      source: 'manual'});
    let error = '';
    try {
      await pinByCode('EX-SMD', SCRATCH, true);
    } catch (e: any) {
      error = e?.message ?? String(e);
    }
    expect(error.startsWith(`${SCRATCH} cannot be pinned: It has no workspace folder yet`), true, error);
    expect(await pinByCode('EX-SMD', SCRATCH, false), `${SCRATCH}: not pinned (no folder)`);
  });

  test('a program\'s address opens its page', async () => {
    const pf = await portfolioByCode('EX-SMD');
    const kin = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "EXS-KIN1"`}))[0];
    const view: DG.View = await grok.functions.call('ProjectPortalCore:ProjectPortalApp', {path: 'EX-SMD/EXS-KIN1'});
    try {
      expect(view.path.endsWith('/Portal/EX-SMD/EXS-KIN1'), true, view.path);
      // the program page is in the context panel, not the view
      const page = grok.shell.o as HTMLElement;
      expect(page?.querySelector?.('[data-u2-name="program-title"]')?.textContent, kin.name);
    } finally {
      view.close();
    }
  });
});
