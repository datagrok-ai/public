import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {PortfolioRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
import {accessFrame, checkConnections, dataAccessChecks, dataAccessWarnings} from '../data-access';
import {
  accessFix, accessFrom, AccessCheck, accessNotice, accessStatus, connectionUses, describeAccessProblems,
  missingPermissions, QUERY, summarizeAccess, VIEW,
} from '../model/data-access';
import {withRestrictedUser} from './restricted-user';
import {_package} from '../package-test';

category('Portal: data access rules', () => {
  const check = (over: Partial<AccessCheck>): AccessCheck => ({portfolio: 'EX-SMD', connection: 'ExData',
    uses: ['dataset program'], needed: true, role: 'viewers', group: 'EX viewers', access: {view: true, query: true}, ...over});

  test('a connection is listed once with everything that uses it; only datasets need the groups', async () => {
    const uses = connectionUses(
      [{code: 'program', connection: 'ExData'}, {code: 'assays', connection: 'ExData'}, {code: 'sar', connection: 'ExSar'}],
      [{code: 'REG', connection: 'ExRegistry', target: 'program'}, {code: 'compounds', connection: 'ExData', target: 'compound'},
        {code: 'units', connection: 'ExRegistry', target: null}]);
    expect(uses.map((u) => `${u.connection}|${u.uses.join(',')}|${u.needed}`).join('\n'), [
      'ExData|dataset program,dataset assays,compound feed compounds|true',
      'ExSar|dataset sar|true',
      'ExRegistry|feed REG,feed units|false',
    ].join('\n'));
  });

  test('View and Query are read by name; Edit does not stand in for either', async () => {
    expect(JSON.stringify(accessFrom([VIEW, QUERY, 'DataQuery.Execute'])), '{"view":true,"query":true}');
    expect(JSON.stringify(accessFrom(['Edit', 'DataConnection.GetSchema'])), '{"view":false,"query":false}');
    expect(missingPermissions({view: true, query: false}).join(), 'Query');
    expect(missingPermissions({view: false, query: false}).join(), 'View,Query');
    expect(missingPermissions(null).length, 0);
  });

  test('status: ok, missing, unknown, and not needed for a feed-only connection', async () => {
    expect([
      check({}), check({access: {view: true, query: false}}), check({access: null, reason: 'not visible'}),
      check({needed: false, access: {view: false, query: false}}),
    ].map(accessStatus).join(), 'ok,missing,unknown,not needed');
  });

  test('a warning names the connection, the group, the missing permission and the fix', async () => {
    const lines = describeAccessProblems([
      check({}),
      check({role: 'approvers', group: 'EX approvers', access: {view: true, query: false}}),
      check({needed: false, access: {view: false, query: false}}),
      check({connection: 'ExSar', uses: ['dataset sar'], access: null, reason: 'you may not read the connection\'s sharing'}),
    ]);
    expect(lines.join('\n'), [
      'EX-SMD: connection ExData (dataset program): EX approvers (approvers) lacks Query; its members cannot open the data. ' +
        'Fix: share ExData with EX approvers: View + Query',
      'EX-SMD: connection ExSar (dataset sar): access of EX viewers (viewers) unknown: you may not read the connection\'s sharing',
    ].join('\n'));
    expect(accessFix('ExData', 'EX viewers'), 'share ExData with EX viewers: View + Query');
  });

  test('the board notice: none when every group can query, else who cannot and what was not checked', async () => {
    expect(accessNotice([check({}), check({needed: false, access: {view: false, query: false}})]), null);
    const n = accessNotice([check({access: {view: false, query: false}}), check({connection: 'ExSar', access: null})])!;
    expect(`${n.missing} ${n.unknown} ${n.text}`,
      '1 1 EX viewers cannot query ExData, so the program data does not open for them; access to ExSar could not be checked');
  });

  test('the text report lists each group per connection, then the fixes', async () => {
    expect(summarizeAccess([
      check({}), check({role: 'approvers', group: 'EX approvers', access: {view: false, query: false}}),
      check({connection: 'ExRegistry', uses: ['feed REG'], needed: false, access: {view: false, query: false}}),
    ]), [
      'EX-SMD: ExData (dataset program)',
      '  EX viewers (viewers): ok',
      '  EX approvers (approvers): missing: lacks View and Query',
      'EX-SMD: ExRegistry (feed REG) — feeds only, run as whoever syncs',
      '  EX viewers (viewers): not needed: lacks View and Query',
      'fixes:',
      '  share ExData with EX approvers: View + Query',
    ].join('\n'));
    expect(summarizeAccess([]), 'no data connections to check');
  });
});

category('Portal: data access', () => {
  let smd: PortfolioRow;
  let groups: {viewers?: string, contributors?: string, approvers?: string};

  before(async () => {
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    groups = {viewers: smd.viewers, contributors: smd.contributors, approvers: smd.approvers};
  });

  /** Runs [body] with a throwaway connection only this test knows, deleted after. */
  async function withScratchConnection<T>(body: (c: DG.DataConnection) => Promise<T>): Promise<T> {
    const name = `ExAccessAudit${Date.now()}`;
    const conn = await grok.dapi.connections.save(DG.DataConnection.create(name,
      {dataSource: 'Postgres', server: 'localhost', db: 'ex_access_audit'}));
    try {
      return await body(conn);
    } finally {
      await grok.dapi.connections.delete(conn).catch(() => undefined);
    }
  }

  test('the EX portfolios\' connections are queryable by all their groups', async () => {
    const checks = await dataAccessChecks('EX-SMD');
    expect(checks.length > 0, true, 'EX-SMD uses at least one connection');
    expect([...new Set(checks.map((c) => c.role))].sort().join(), 'approvers,contributors,viewers');
    expect(checks.filter((c) => accessStatus(c) !== 'ok').map((c) => `${c.connection} ${c.group}`).join('; '), '');
    const df = accessFrame(checks);
    expect(df.rowCount, checks.length);
    expect(df.col('status')!.toList().every((s) => s === 'ok'), true);
  });

  test('a dataset on a connection not shared with the groups: every group lacks View and Query, with the fix', async () => {
    await withScratchConnection(async (conn) => {
      const checks = await checkConnections('EX-SMD', groups,
        [{connection: conn.name, uses: ['dataset probe'], needed: true}]);
      expect(checks.map((c) => `${c.role}:${accessStatus(c)}:${missingPermissions(c.access).join('+')}`).join(),
        'viewers:missing:View+Query,contributors:missing:View+Query,approvers:missing:View+Query');
      const df = accessFrame(checks);
      expect(df.col('fix')!.get(0), `share ${conn.name} with ${checks[0].group}: View + Query`);
    });
  });

  test('the loader warns about such a dataset and grants nothing', async () => {
    await withScratchConnection(async (conn) => {
      const file = (await readPortfolioFiles(_package)).find((c) => c.code === 'EX-SMD')!;
      const probe = {...file, programDatasets: [...(file.programDatasets ?? []),
        {code: 'probe', name: 'Probe', connection: conn.name, sql: 'select @program as program'}]};
      const lines = await dataAccessWarnings(probe, groups);
      expect(lines.length, 3, lines.join('\n'));
      expect(lines.every((l) => l.startsWith(`warning: EX-SMD: connection ${conn.name} (dataset probe): `) &&
        l.includes(`Fix: share ${conn.name} with `)), true, lines.join('\n'));
      const after = await grok.dapi.permissions.get(conn);
      expect([...(after.view ?? []), ...(after.edit ?? [])].some((g) => Object.values(groups).includes(g.id)), false, 'no group was given the connection');
      expect((await dataAccessWarnings(file, groups)).join('\n'), '', 'the shipped EX-SMD file warns about nothing');
    });
  });

  test('a share with the viewers covers the contributors and approvers through the parent group', async () => {
    await withScratchConnection(async (conn) => {
      // only this test's own throwaway connection is shared. The js-api's View share
      // (permissions.grant, edit false) gives the group the whole Use set, Query included.
      await grok.dapi.permissions.grant(conn, await grok.dapi.groups.find(smd.viewers!), false);
      const checks = await checkConnections('EX-SMD', groups, [{connection: conn.name, uses: ['dataset probe'], needed: true}]);
      expect(checks.map((c) => `${c.role}:${accessStatus(c)}`).join(), 'viewers:ok,contributors:ok,approvers:ok');
      expect(describeAccessProblems(checks).length, 0);
    });
  });

  test('a user who cannot see the connection gets "unknown", not an error', async () => {
    await withScratchConnection(async (conn) => {
      const ran = await withRestrictedUser('ppaccess', async (user) => {
        const checks = await user.asUser(() => checkConnections('EX-SMD', groups,
          [{connection: conn.name, uses: ['dataset probe'], needed: true}]));
        expect(checks.map(accessStatus).join(), 'unknown,unknown,unknown');
        expect(checks[0].reason, 'the connection is not found or not shared with you');
        return true;
      });
      if (ran === null)
        console.log('skipped: no restricted user on this server');
    });
  });
});
