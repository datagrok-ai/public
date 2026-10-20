import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after as afterAll, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programData, readPortfolioFiles} from '../portal';
import {validatePortfolio} from '../model/config';
import {after, compare, describeNews, followLabel, rowsOf, snapshot} from '../model/new-data';
import {checkNewDataOf, follow, myFollows, unfollow, whatsNew} from '../new-data';
import {recordedSync} from '../sync-runs';
import {_package} from '../package-test';

const rows = [
  {compound: 'C1', series: 'A', ic50: 0.1, sol: 10, run: '2026-09-01'},
  {compound: 'C2', series: 'A', ic50: 0.2, sol: null, run: '2026-09-15'},
  {compound: 'C3', series: 'B', ic50: null, sol: null, run: null},
];

category('Portal: new-data rules', () => {
  test('a snapshot counts rows, filled values and the watermark\'s largest value', async () => {
    expect(JSON.stringify(snapshot(rows)), '{"rows":3,"values":11}');
    expect(JSON.stringify(snapshot(rows, 'run')), '{"rows":3,"values":11,"mark":"2026-09-15"}');
    expect(JSON.stringify(snapshot(rowsOf(rows, 'B'))), '{"rows":1,"values":2}', 'one series');
    expect(after(10, 9) && !after(9, 10) && after('2026-10-01', '2026-09-30') && after(1, null) && !after(null, 1), true);
  });

  test('news is growth since the snapshot; the first look is a baseline; less is never news', async () => {
    const first = compare(null, rows);
    expect(first.baseline && first.newRows === 0 && first.newValues === 0, true);
    const more = [...rows, {compound: 'C4', series: 'B', ic50: 2, sol: null, run: '2026-10-01'}];
    const n = compare(first.now, more);
    expect(`${n.newRows} ${n.newValues}`, '1 4');
    expect(describeNews(followLabel('P1'), n), 'P1: 1 new row, 4 new values');
    const filled = compare(first.now, rows.map((r) => r.compound === 'C2' ? {...r, sol: 5} : r));
    expect(describeNews(followLabel('P1', {kind: 'series', code: 'A'}), filled), 'P1 series A: 1 new value',
      'a measurement in a blank cell');
    expect(describeNews('P1', compare(first.now, rows.slice(0, 2))), null, 'a row gone is not news');
    expect(describeNews('P1', first), null, 'a baseline is not news');
  });

  test('with a watermark, a row replaced by a newer one is new although the count is the same', async () => {
    const base = compare(null, rows, 'run').now;
    const replaced = [rows[0], rows[1], {compound: 'C5', series: 'B', ic50: 1, sol: 3, run: '2026-10-02'}];
    const n = compare(base, replaced, 'run');
    expect(`${n.newRows} ${n.now.mark}`, '1 2026-10-02');
    expect(compare(base, replaced).newRows, 0, 'counts alone see no new row');
  });

  test('a watermark must name a column', async () => {
    const smd = (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-SMD')!;
    const bad = {...smd, programData: {...smd.programData!, watermark: ' '}};
    expect(validatePortfolio(bad).some((e) => e.includes('watermark must name a column')), true);
    expect(validatePortfolio({...smd, programData: {...smd.programData!, watermark: 'run_date'}}).join('; '), '');
  });
});

category('Portal: new-data notification', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;
  const lower = async (series: string | null, by: {rows: number, values: number}, field: 'checked' | 'seen' = 'checked') => {
    const f = (await myFollows(kin1.id)).find((x) => series ? x.work_unit_code === series : !x.work_unit_code)!;
    const s = f[field] as {rows: number, values: number};
    await portalDb.followings.update(f.id, {[field]: {rows: s.rows - by.rows, values: s.values - by.values}});
  };

  const mine = <T extends {user: string}>(news: T[]) => news.filter((x) => x.user === DG.User.current().id);

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

  test('following takes the data now as the baseline; twice is once; a series is its own follow', async () => {
    const f = await follow(smd, kin1);
    expect(JSON.stringify(f.checked), '{"rows":11,"values":87}');
    expect((await follow(smd, kin1)).id, f.id);
    const s = await follow(smd, kin1, {kind: 'series', code: 'KIN-A'});
    expect((s.checked as any).rows, 8);
    expect((await myFollows(kin1.id)).length, 2);
  });

  test('nothing new: no notification; new rows and values: the followers are told, once', async () => {
    const quiet = await checkNewDataOf(smd);
    // other users' follows on the stand may add news and counts: only this user's are this test's
    expect(mine(quiet.news).length, 0, quiet.summary);
    expect(/^EX-SMD new data: .*\(\d+ follow\(s\)/.test(quiet.summary), true, quiet.summary);
    await lower(null, {rows: 1, values: 4});
    const r = await checkNewDataOf(smd);
    expect(mine(r.news).map((x) => x.line).join(), 'EXS-KIN1: 1 new row, 4 new values');
    expect(r.notified.includes(`${DG.User.current().login}: notified`), true, r.notified.join());
    expect((await myFollows(kin1.id)).find((f) => !f.work_unit_code)!.last_news, 'EXS-KIN1: 1 new row, 4 new values');
    expect(mine((await checkNewDataOf(smd)).news).length, 0, 'the snapshot moved on');
  });

  test('a series follower hears only about the series', async () => {
    await lower('KIN-A', {rows: 0, values: 2});
    const r = await checkNewDataOf(smd, {notify: false});
    expect(mine(r.news).map((x) => x.line).join(), 'EXS-KIN1 series KIN-A: 2 new values');
  });

  test('what\'s new since the last visit, then nothing', async () => {
    await lower(null, {rows: 2, values: 9}, 'seen');
    const rowsNow = await programData(smd, 'EXS-KIN1');
    const first = await whatsNew(smd, kin1, rowsNow);
    expect(first!.lines.join(), 'EXS-KIN1: 2 new rows, 9 new values');
    expect((await whatsNew(smd, kin1, rowsNow))!.lines.length, 0);
  });

  test('every sync checks the followed programs, and the function says what it found', async () => {
    await lower(null, {rows: 0, values: 1});
    const r = await recordedSync('EX-SMD', 'manual');
    expect(r.report.includes('EX-SMD new data: EXS-KIN1: 1 new value'), true, r.report);
    const text: string = await grok.functions.call('ProjectPortal:checkNewData', {code: 'EX-SMD'});
    expect(text.startsWith('EX-SMD new data: nothing new'), true, text);
    const biod: string = await grok.functions.call('ProjectPortal:checkNewData', {code: 'EX-BIOD'});
    // nobody follows EX-BIOD's programs unless another test or user left a follow there
    expect(biod.startsWith('EX-BIOD new data: '), true, biod);
  });

  test('unfollowing removes the follow; the functions follow and unfollow by code', async () => {
    expect(await grok.functions.call('ProjectPortal:unfollowProgram', {portfolio: 'EX-SMD', program: 'EXS-KIN1', series: 'KIN-A'}),
      'EXS-KIN1 series KIN-A: no longer followed');
    expect((await myFollows(kin1.id)).length, 1);
    expect(await grok.functions.call('ProjectPortal:followProgram', {portfolio: 'EX-SMD', program: 'EXS-KIN1'}),
      'EXS-KIN1: followed (11 rows now)');
    expect(await unfollow(kin1), 1);
    expect((await myFollows(kin1.id)).length, 0);
  });
});
