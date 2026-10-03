import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after as afterAll, before, category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {
  describeEntry, FeedRecords, filterFeed, mergeFeeds, newsUnit, pageFeed, parseKinds, programFeedEntries,
} from '../model/feed';
import {checkNewDataOf, follow, myFollows, unfollow} from '../new-data';
import {myFeedOf, programFeedOf} from '../feed';
import {_package} from '../package-test';

const records = (): FeedRecords => ({
  program: {id: 'p1', code: 'P1', portfolio: 'EX'},
  stages: [{id: 's1', code: 'H2L', name: 'Hit-to-Lead'}, {id: 's2', code: 'LO', name: 'Lead Optimization'}],
  decisions: [{id: 'd1', decision: 'advance', from: 's1', to: 's2', by: 'u1', on: '2026-09-01T10:00:00Z', source: 'approver',
    rationale: 'potent and soluble'}],
  evidence: [
    {id: 'e1', decisionId: 'd1', kind: 'model_run', name: 'Run A', by: 'u1', on: '2026-09-01T10:00:30Z'},
    {id: 'e2', decisionId: 'd1', kind: 'url', name: 'Minutes', by: 'u2', on: '2026-09-03T09:00:00Z'},
  ],
  runs: [{id: 'r1', title: 'Run A', model: 'FedBatchTiter', by: 'u3', on: '2026-08-30T08:00:00Z'}],
  selections: [{id: 'a1', run: 'Run A', criterion: 'Process description', by: 'u3', on: '2026-08-31T08:00:00Z', rationale: 'closest'}],
  documents: [{id: 'doc1', title: 'Report', version: '1.0', status: 'approved', linkedBy: 'u3', linkedOn: '2026-08-29T08:00:00Z',
    approvedBy: 'u1', approvedOn: '2026-09-02T08:00:00Z'}],
  snapshots: [{id: 'sn1', label: 'Before the gate', by: 'u1', on: '2026-09-01T09:59:00Z'}],
  roleChanges: [{id: 'c1', compound: 'C-2', from: 'hit', to: 'lead', by: null, on: '2026-08-20T00:00:00Z'}],
  syncs: [
    {id: 'y1', on: '2026-09-04T06:00:00Z', by: 'u0', newData: ['P1 series A: 2 new values'],
      change: {program: 'P1', fields: ['modality'], units: {created: ['series B'], updated: ['series A'], missing: []},
        compounds: {added: ['C-9'], roleChanges: [], missing: []}}},
  ],
  news: [
    {id: 'f1', line: 'P1 series A: 2 new values', on: '2026-09-04T06:01:00Z'},
    {id: 'f2', line: 'P1: 1 new row', on: '2026-09-05T06:00:00Z'},
  ],
});

category('Portal: feed rules', () => {
  test('entries come from the records: evidence with its decision folds in, later evidence is its own', async () => {
    const e = mergeFeeds(programFeedEntries(records()));
    const decision = e.find((x) => x.kind === 'decision')!;
    expect(decision.text, 'Advance: Hit-to-Lead → Lead Optimization, with 1 piece of evidence — potent and soluble');
    expect(e.filter((x) => x.kind === 'evidence').map((x) => x.text).join(),
      'Evidence added to advance Hit-to-Lead → Lead Optimization: Minutes (url)');
    expect(e.filter((x) => x.kind === 'document').map((x) => x.text).join(' | '), 'Document approved: Report v1.0 | Document linked: Report v1.0');
    expect(e.find((x) => x.kind === 'selection')!.text, 'Run selected for Process description: Run A — closest');
    expect(e.find((x) => x.kind === 'registry')!.text, 'Sync: registry changed modality');
    expect(e.filter((x) => x.kind === 'work_unit').map((x) => `${x.unit}=${x.text}`).sort().join(), 'series A=series A: updated,series B=series B: new');
  });

  test('newest first; the same news told twice within minutes is one entry', async () => {
    const e = mergeFeeds(programFeedEntries(records()));
    expect(e[0].text, 'New data: P1: 1 new row');
    expect(e.filter((x) => x.kind === 'new_data').length, 2);
    for (let i = 1; i < e.length; i++)
      expect(e[i - 1].at >= e[i].at, true, `${e[i - 1].at} before ${e[i].at}`);
    expect(mergeFeeds(e, e).length, e.length, 'merged twice counts once');
    expect(describeEntry(e[0], 'Ann', true), '2026-09-05 06:00 · P1 · New data: P1: 1 new row · Ann');
  });

  test('program roles and approved model cards are entries; a removed role says so', async () => {
    const e = programFeedEntries({...records(),
      programRoles: [{id: 'r1', group: 'Chem team', role: 'contributor', active: true, by: 'u1', on: '2026-09-06T08:00:00Z'},
        {id: 'r2', group: 'Ann', role: 'approver', active: false, by: 'u1', on: '2026-09-06T09:00:00Z'}],
      modelCards: [{id: 'c1', model: 'Solubility estimate', version: '1.0', by: 'u2', on: '2026-09-07T08:00:00Z'},
        {id: 'c2', model: 'Old model', version: '0.1', by: null, on: null}]});
    expect(e.filter((x) => x.kind === 'role').map((x) => x.text).sort().join(' | '),
      'Program role removed: Ann | Program role: Chem team is contributor');
    expect(e.filter((x) => x.kind === 'model_card').map((x) => `${x.text} · ${x.who}`).join(),
      'Model card approved: Solubility estimate 1.0 · u2', 'a card with no approval date is not an entry');
    expect(parseKinds('role, model_card').join(), 'role,model_card');
  });

  test('a series narrows to its work units, its news and its compounds; kinds filter; pages', async () => {
    const e = mergeFeeds(programFeedEntries(records()));
    const a = filterFeed(e, {series: 'A', seriesCompounds: ['C-2']});
    expect(a.map((x) => x.text).join(' | '), 'series A: updated | New data: P1 series A: 2 new values | C-2: hit → lead');
    expect(filterFeed(e, {kinds: ['run', 'snapshot']}).map((x) => x.kind).join(), 'snapshot,run');
    expect(filterFeed(e, {since: '2026-09-03T00:00:00Z'}).length, 7);
    const p = pageFeed(e, 5, 4);
    expect(`${p.page} ${p.pages} ${p.total} ${p.items.length}`, `${p.pages - 1} ${Math.ceil(e.length / 4)} ${e.length} ${e.length - 4 * (p.pages - 1)}`);
    expect(newsUnit('P1', 'P1 series A: 2 new values'), 'series A');
    expect(newsUnit('P1', 'P1: 1 new row'), null);
    expect(parseKinds(' run, new_data ').join(), 'run,new_data');
    let bad = '';
    try {
      parseKinds('runs');
    } catch (x: any) {
      bad = x.message;
    }
    expect(bad.startsWith('"runs" is not a feed kind'), true);
  });
});

category('Portal: activity feed', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;
  const lower = async (series: string | null, values: number) => {
    const f = (await myFollows(kin1.id)).find((x) => series ? x.work_unit_code === series : !x.work_unit_code)!;
    const s = f.checked as {rows: number, values: number};
    await portalDb.follows.update(f.id, {checked: {rows: s.rows, values: s.values - values}});
  };

  before(async () => {
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`}))[0];
    await unfollow(kin1);
    await unfollow(kin1, {kind: 'series', code: 'KIN-A'});
  });

  afterAll(async () => {
    await unfollow(kin1);
    await unfollow(kin1, {kind: 'series', code: 'KIN-A'});
  });

  test('new data a check found is in the program\'s feed and in the follower\'s', async () => {
    await follow(smd, kin1);
    await lower(null, 3);
    const r = await checkNewDataOf(smd, {notify: false});
    expect(r.news.filter((x) => x.user === DG.User.current().id).map((x) => x.line).join(), 'EXS-KIN1: 3 new values');
    const feed = await programFeedOf(smd, kin1);
    // found by what it says, not by being first: a newer record of EXS-KIN1 from elsewhere may lead
    const news = feed.find((x) => x.kind === 'new_data' && x.text === 'New data: EXS-KIN1: 3 new values');
    expect(!!news && Date.now() - Date.parse(news.at) < 10 * 60 * 1000, true, feed.slice(0, 3).map((x) => x.text).join('; '));
    const decisions = (await portalDb.stageTransitions.query({filter: `program_id = "${kin1.id}"`})).length;
    expect(feed.filter((x) => x.kind === 'decision').length, decisions, 'one entry per gate-log row');
    const mine = await myFeedOf();
    expect(mine.entries.some((x) => x.program === 'EXS-KIN1' && x.text === 'New data: EXS-KIN1: 3 new values'), true);
  });

  test('a series follow sees only its series; the series feed leaves the program\'s news out', async () => {
    await unfollow(kin1);
    await follow(smd, kin1, {kind: 'series', code: 'KIN-A'});
    await lower('KIN-A', 2);
    await checkNewDataOf(smd, {notify: false});
    const series = await programFeedOf(smd, kin1, {series: 'KIN-A'});
    expect(series.some((x) => x.text === 'New data: EXS-KIN1 series KIN-A: 2 new values'), true);
    expect(series.some((x) => x.text === 'New data: EXS-KIN1: 3 new values'), false);
    const mine = await myFeedOf();
    expect(mine.entries.filter((x) => x.program === 'EXS-KIN1').every((x) => x.unit === 'series KIN-A' ||
      (x.compounds ?? []).some((c) => c.startsWith('EXC-10'))), true, 'only the series');
  });

  test('a program role given on the program is in its feed, with who gave it', async () => {
    const me = DG.User.current();
    const stamp = {granted_by: me.id, granted_on: new Date().toISOString() as any};
    // one row per program and group: the user's own may be there already (inactive), and is put back
    const before = (await portalDb.programRoles.query({filter: `program_id = "${kin1.id}" and group_id = "${me.group.id}"`}))[0];
    const id = before ? before.id : (await portalDb.programRoles.insert({program_id: kin1.id, group_id: me.group.id, role: 'viewer',
      active: true, ...stamp}))[0].id;
    try {
      if (before)
        await portalDb.programRoles.update(id, {role: 'viewer', active: true, ...stamp});
      const feed = await programFeedOf(smd, kin1, {kinds: ['role']});
      const entry = feed.find((x) => x.key.startsWith(`program-role:${id}:`));
      expect(entry?.text.endsWith(' is viewer'), true, feed.map((x) => x.text).join('; '));
      expect(entry?.who, me.id);
    } finally {
      if (before) {
        await portalDb.programRoles.update(id, {role: before.role, active: before.active, granted_by: before.granted_by,
          granted_on: before.granted_on}).catch(() => undefined);
      } else
        await portalDb.programRoles.delete(id).catch(() => undefined);
    }
  });

  test('programFeed and myFeed return the feed as a table; kinds filter it', async () => {
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:programFeed',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', kinds: 'new_data'});
    expect(df.columns.names().join(), 'at,portfolio,program,kind,what,who,work unit,record');
    expect(df.rowCount >= 2, true);
    expect(df.col('kind')!.categories.join(), 'new_data');
    const mine: DG.DataFrame = await grok.functions.call('ProjectPortalCore:myFeed', {kinds: 'new_data', limit: 1});
    expect(mine.rowCount, 1);
    await expectExceptionAsync(() => grok.functions.call('ProjectPortalCore:programFeed',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', kinds: 'gossip'}));
  });
});
