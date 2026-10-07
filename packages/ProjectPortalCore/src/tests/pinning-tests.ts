import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {assembleTiles, checkGroupPin, describeTile, isProgramRecordName, NOTHING_PINNED, PinHolders, pinOf, pinsOf,
  programAddress, TileInput, tileTitle, unseenNews, whyNotPinnable} from '../model/pinning';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {myPins, pinByCode} from '../pinning';
import {pinnedFrame, pinnedTiles} from '../pinned-programs';
import {COL, portfolioTables} from '../portfolio-table';
import {ProgramRecordHandler} from '../program-record';
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
  test('the favorite is the program\'s record: no folder needed, only the record found', async () => {
    expect(whyNotPinnable({recordFound: true}), null);
    expect(whyNotPinnable({recordFound: false})!.startsWith('Its record is not found'), true);
    expect(isProgramRecordName('portal.program:657f|EXS-KIN1'), true);
    expect(isProgramRecordName('portal.sync_log:1'), false);
    expect(isProgramRecordName(null), false);
    expect(NOTHING_PINNED.includes('"Add to favorites"'), true, 'the widget names the button');
  });

  test('pinning for a group is the approvers\'', async () => {
    expect(checkGroupPin({isApprover: true}).length, 0);
    expect(checkGroupPin({isApprover: false})[0], 'Only the portfolio\'s approvers pin a program for a group');
  });

  test('pins are the program records and folders among the user\'s and their groups\' favorites', async () => {
    const pins = pinsOf([{id: 'p1', kind: 'program'}, {id: 'q', kind: null}],
      [{group: 'Chemists', entities: [{id: 'p1', kind: 'program'}, {id: 's2', kind: 'folder'}]},
        {group: 'Biology', entities: [{id: 's2', kind: 'folder'}]}]);
    expect([...pins.keys()].sort().join(), 'p1,s2', 'a favorite that is neither is no pin');
    expect(JSON.stringify(pins.get('p1')), '{"mine":true,"groups":["Chemists"]}');
    expect(JSON.stringify(pins.get('s2')), '{"mine":false,"groups":["Chemists","Biology"]}');
  });

  test('a program\'s pin: its record\'s and its old folder\'s together', async () => {
    const pins = new Map<string, PinHolders>([['p1', {mine: false, groups: ['Chemists']}], ['s1', {mine: true, groups: ['Chemists', 'Biology']}]]);
    expect(JSON.stringify(pinOf(pins, {programId: 'p1', spaceId: 's1'})), '{"mine":true,"groups":["Chemists","Biology"]}');
    expect(JSON.stringify(pinOf(pins, {programId: 'p1', spaceId: null})), '{"mine":false,"groups":["Chemists"]}');
    expect(pinOf(pins, {programId: 'p2', spaceId: 's2'}) === undefined, true, 'pinned by neither');
    const [t] = assembleTiles([input({programId: 'p1', spaceId: null})], mine('p1'));
    expect(t?.code, 'EXS-KIN1', 'a program without a folder has a tile from its record');
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
    const none = input({programId: 'd', code: 'N-0', spaceId: null});   // pinned by neither
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
  let record: DG.Entity;
  /** How the favorites stood (the record's and the folder's, the user's and the group's), to put them back. */
  let had: [DG.Entity, DG.Group | undefined, boolean][] = [];

  const favorite = async (group?: DG.Group, entity?: DG.Entity) =>
    (await grok.dapi.entities.getFavorites(group ?? DG.User.current().group)).some((e) => e.id === (entity ?? record).id);

  before(async () => {
    await loadAllPortfolios(_package);
    const pf = await portfolioByCode('EX-SMD');
    program = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "${CODE}"`}))[0];
    space = await grok.dapi.spaces.find(program.space_id!);
    record = await grok.dapi.entities.find(program.id);
    viewers = await grok.dapi.groups.find(pf.viewers!);
    had = [];
    for (const e of [record, space]) {
      for (const g of [undefined, viewers])
        had.push([e, g, await favorite(g, e)]);
    }
  });

  after(async () => {
    for (const [e, group, on] of had) {
      if (on !== await favorite(group, e))
        await (on ? DG.Favorites.add(e, group) : DG.Favorites.remove(e, group));
    }
    for (const r of await portalDb.programs.query({filter: `code = "${SCRATCH}"`})) {
      const e = await grok.dapi.entities.find(r.id).catch(() => null);
      if (e && await favorite(undefined, e))
        await DG.Favorites.remove(e);
    }
    for (const r of await portalDb.programs.query({filter: `code = "${SCRATCH}"`}))
      await portalDb.programs.delete(r.id);
  });

  test('add and remove: a favorite on the program\'s record, a tile, the table\'s pinned column', async () => {
    await pinByCode('EX-SMD', CODE, false);
    expect(await favorite(), false);
    expect(await pinByCode('EX-SMD', CODE, true), `${CODE}: added to favorites; it is in the Workspace and has a tile on Home`);
    expect(await favorite(), true, 'the platform favorite on the record is the pin');
    expect(await favorite(undefined, space), false, 'the folder is not favorited');
    const listed = (await grok.dapi.entities.getFavorites(DG.User.current().group)).find((e) => e.id === program.id);
    expect(`${listed?.entityType} ${listed?.friendlyName}`, `EntityRecord ${program.code}`, 'the Workspace lists it by the program\'s name');
    expect(await pinByCode('EX-SMD', CODE, true), `${CODE}: already in favorites`);
    const tile = (await pinnedTiles()).find((t) => t.code === CODE);
    expect(tile?.portfolioCode, 'EX-SMD');
    expect(tile?.pinnedByMe, true);
    const df = pinnedFrame(await pinnedTiles());
    expect(df.col('program')!.toList().includes(CODE), true);
    const table = (await portfolioTables('EX-SMD')).programs;
    const row = table.col(COL.program)!.toList().indexOf(CODE);
    expect(table.get(COL.pinned, row), true, 'pinned column');
    expect(await pinByCode('EX-SMD', CODE, false), `${CODE}: removed from favorites`);
    expect(await favorite(), false);
    expect((await pinnedTiles()).some((t) => t.code === CODE), false);
  });

  test('unpinning in the platform takes the tile away', async () => {
    await pinByCode('EX-SMD', CODE, true);
    expect((await pinnedTiles()).some((t) => t.code === CODE), true);
    await DG.Favorites.remove(record);   // what the Workspace's Unpin and the star do
    expect((await myPins()).get(program.id)?.mine ?? false, false);
    expect((await pinnedTiles()).some((t) => t.code === CODE), false);
  });

  test('an approver pins for a group: the group\'s favorite, a tile marked for it', async () => {
    const name = viewers.friendlyName || viewers.name;
    expect(await pinByCode('EX-SMD', CODE, true, name), `${CODE} for ${name}: added to favorites; it is in the Workspace and has a tile on Home`);
    expect(await favorite(viewers), true);
    expect((await pinnedTiles()).find((t) => t.code === CODE)?.pinnedFor.includes(name), true);
    await pinByCode('EX-SMD', CODE, false, name);
    expect(await favorite(viewers), false);
  });

  test('a program without a folder can be added: its record is the favorite', async () => {
    const pf = await portfolioByCode('EX-SMD');
    await portalDb.programs.insert({portfolio_id: pf.id, code: SCRATCH, name: SCRATCH, stage_model_id: pf.stage_model_id,
      source: 'manual'});
    expect(await pinByCode('EX-SMD', SCRATCH, true), `${SCRATCH}: added to favorites; it is in the Workspace and has a tile on Home`);
    expect((await pinnedTiles()).some((t) => t.code === SCRATCH), true);
    expect(await pinByCode('EX-SMD', SCRATCH, false), `${SCRATCH}: removed from favorites`);
    expect(await pinByCode('EX-SMD', SCRATCH, false), `${SCRATCH}: was not in favorites`);
  });

  test('a folder pinned before 0.3.0 still counts; removing takes it out too', async () => {
    await pinByCode('EX-SMD', CODE, false);
    await DG.Favorites.add(space);   // the old way
    expect((await pinnedTiles()).find((t) => t.code === CODE)?.pinnedByMe, true);
    expect(await pinByCode('EX-SMD', CODE, true), `${CODE}: already in favorites`);
    expect(await pinByCode('EX-SMD', CODE, false), `${CODE}: removed from favorites`);
    expect(`${await favorite()} ${await favorite(undefined, space)}`, 'false false');
  });

  test('a favorite program opens its page from the Workspace', async () => {
    await pinByCode('EX-SMD', CODE, true);
    const listed = (await grok.dapi.entities.getFavorites(DG.User.current().group)).find((e) => e.id === program.id)!;
    const handler = DG.ObjectHandler.forEntity(listed);
    // the package's own instance, registered by portalInit: the test bundle's class is another copy
    expect(handler?.type, new ProgramRecordHandler().type);
    // what the Workspace's Open does: the handler's preview, added as a view
    const view = await handler!.renderPreview(listed);
    try {
      expect(view.path.endsWith(`/Portal/EX-SMD/${CODE}`), true, view.path);
    } finally {
      view.close();
      await pinByCode('EX-SMD', CODE, false);
    }
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
