import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {ApprovalRow, portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {moveRecordsOutOfRoot, moveToSecuredTables, myPortfolioIds, placeApproval, placeSyncLog, recordFolderId, rootIdOf}
  from '../secured-rows';
import {
  compoundOfLink, copyValues, describeMove, describeRecordMoves, hasContent, inPortfolios, linkValuesOf, planLinks,
  planRecordMoves, portfoliosOfGroups, recordFolderFor, requestSpace, toCopy,
} from '../model/secured-rows';
import {folderLinks} from '../program-rows';
import {SYNC_HISTORY_FOLDER} from '../model/workspace-plan';
import {withStandUser} from './stand-users';
import {_package} from '../package-test';

category('Portal: who reads a portfolio\'s records (rules)', () => {
  test('a copy keeps the id and drops the server\'s columns and empty values', async () => {
    const v = copyValues({id: 'r1', version: 3, created_on: 'x', updated_on: 'y', author_id: 'u', title: 'T', rationale: null});
    expect(JSON.stringify(v), '{"id":"r1","title":"T"}');
    expect(toCopy([{id: 'a'}, {id: 'b'}], new Set(['a'])).map((r) => r.id).join(), 'b');
  });

  test('a request goes into its program\'s folder, else the portfolio\'s _Approvals', async () => {
    expect(requestSpace('folder', 'approvals'), 'folder');
    expect(requestSpace(null, 'approvals'), 'approvals', 'a model card, or a program without a folder');
    expect(requestSpace(undefined, null), null, 'no workspace: nowhere');
  });

  test('records out of the root: sync runs to _Sync history, requests to _Approvals, nothing else', async () => {
    expect(recordFolderFor('portal.sync_log:1'), '_Sync history');
    expect(recordFolderFor('portal.approval:2'), '_Approvals');
    expect(recordFolderFor('EXS-KIN1'), null, 'a program folder stays');
    expect(recordFolderFor('portal.program:x|EXS-KIN1'), null);
    const moves = planRecordMoves([{id: 'a', name: 'portal.sync_log:a'}, {id: 'b', name: 'portal.sync_log:b'},
      {id: 'c', name: 'portal.approval:c'}, {id: 'f', name: 'EXS-KIN1'}], new Map([['_Sync history', new Set(['b'])]]));
    expect(moves.map((m) => `${m.id}>${m.folder}${m.linked ? ' (linked)' : ''}`).join(', '),
      'a>_Sync history, b>_Sync history (linked), c>_Approvals');
    expect(describeRecordMoves('EX-SMD', moves, true)[0], 'EX-SMD: would move 3 record(s) out of the root: 2 into _Sync history, 1 into _Approvals');
    expect(describeRecordMoves('EX-SMD', [], false)[0], 'EX-SMD: no records in the root');
  });

  test('links: what to link, what is linked, what has nowhere to go', async () => {
    const plan = planLinks([
      {id: 'a', label: 'A', spaceId: 's1'}, {id: 'b', label: 'B', spaceId: 's1'},
      {id: 'c', label: 'C', spaceId: null}, {id: 'd', label: 'D', spaceId: 'gone'},
    ], new Map([['s1', new Set(['b'])]]));
    expect(`${plan.link.map((p) => p.id).join()} ${plan.linked} ${plan.unplaced.map((p) => p.id).join()}`, 'a 1 c,d');
  });

  test('the app keeps only the user\'s portfolios; an admin keeps all', async () => {
    const rows = [{portfolio_id: 'p1'}, {portfolio_id: 'p2'}];
    expect(inPortfolios(rows, new Set(['p2'])).map((r) => r.portfolio_id).join(), 'p2');
    expect(inPortfolios(rows, null).length, 2);
    const mine = portfoliosOfGroups([{id: 'p1', viewers: 'g1'}, {id: 'p2', contributors: 'g2'}, {id: 'p3', approvers: 'g3'}],
      new Set(['g2', 'g3']));
    expect([...mine].sort().join(), 'p2,p3');
  });

  test('a link carries its compound\'s values; an old link borrows the compound row\'s', async () => {
    expect(JSON.stringify(linkValuesOf({code: 'C1', name: 'n', smiles: 'CCO', aliases: null})),
      '{"compound_code":"C1","compound_name":"n","smiles":"CCO"}');
    const c = compoundOfLink({compound_id: 'id1', smiles: 'CCN'}, {code: 'C1', name: 'old', smiles: 'CCO'});
    expect(`${c.id} ${c.code} ${c.name} ${c.smiles}`, 'id1 C1 old CCN');
    expect(compoundOfLink({compound_id: 'id1'}).code, 'id1', 'nothing known: the id');
    expect(`${hasContent({code: 'C1', smiles: 'C'})} ${hasContent({code: 'C1', smiles: null, name: ''})}`, 'true false');
    expect(describeMove({portfolio: 'EX-SMD', copied: {approval: 1}, already: {approval: 2}, linked: {approval: 3},
      unplaced: ['approval "x"'], linkValues: 4, emptied: 5, errors: []}, true).join('\n'), [
      'EX-SMD approval: would copy 1, 2 copied before; would link 3',
      'EX-SMD compounds: would copy values onto 4 program link(s), would empty 5 compound row(s)',
      'EX-SMD approval "x": nowhere to link it (no workspace): its author and admins read it',
    ].join('\n'));
  });
});

const PROBE = 'read-scope probe';

/** A probe request in [portfolio] about [program], placed like a real one. */
async function probeRequest(portfolio: PortfolioRow, program: ProgramRow): Promise<ApprovalRow> {
  const [r] = await portalDb.approvals.insert({portfolio_id: portfolio.id, program_id: program.id, kind: 'decision',
    decision: 'hold', title: `${PROBE} ${program.code}`, rationale: `${PROBE}: private rationale`, source: 'person',
    proposed_by: DG.User.current().id, proposed_on: dayjs(), approvals_needed: 1, status: 'withdrawn'});
  const row = await portalDb.approvals.get(r.id);
  expect(await placeApproval(row, portfolio), null, `${program.code}'s probe request not linked into its folder`);
  return row;
}

/** Rows of [table] the session can read, as the domain API answers; a refusal reads as none. */
const readable = async (table: string, filter?: string): Promise<any[]> =>
  grok.dapi.domains.table(`portal.${table}`).query({...(filter ? {filter} : {}), limit: 10000}).catch(() => []);

category('Portal: what another portfolio cannot read', () => {
  let smd: PortfolioRow;
  let bioc: PortfolioRow;
  let biod: PortfolioRow;
  let others: Set<string>;
  let kin1: ProgramRow;
  let b101: ProgramRow;

  before(async () => {
    await loadAllPortfolios(_package);
    [smd, bioc, biod] = await Promise.all(['EX-SMD', 'EX-BIOC', 'EX-BIOD'].map(portfolioByCode));
    // every portfolio but EX-BIOD: what the outsider must not read, other customers' portfolios included
    others = new Set((await portalDb.portfolios.query({limit: 1000})).filter((p) => p.id !== biod.id).map((p) => p.id));
    for (const code of ['EX-SMD', 'EX-BIOD', 'EX-BIOC'])
      await moveToSecuredTables(code);
    const withFolder = async (pf: PortfolioRow) => (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}"`, limit: 1000}))
      .filter((p) => p.space_id).sort((a, b) => a.code.localeCompare(b.code))[0];
    kin1 = await withFolder(smd);
    b101 = await withFolder(biod);
    for (const r of await portalDb.approvals.query({filter: `title = "${PROBE} ${kin1?.code}" or title = "${PROBE} ${b101?.code}"`, limit: 100}))
      await portalDb.approvals.delete(r.id).catch(() => undefined);
  });

  test('nobody but admins reads the tables every user could read', async () => {
    const seen = await withStandUser('outsider', async (u) => {
      const out: string[] = [];
      for (const t of ['approval_request', 'sync_run', 'compound', 'follow', 'parked_folder'])
        out.push(`${t} ${(await u.asUser(() => readable(t))).length}`);
      return out.join(', ');
    });
    expect(seen, 'approval_request 0, sync_run 0, compound 0, follow 0, parked_folder 0');
  });

  test('requests: a portfolio\'s people read them, another portfolio\'s do not', async () => {
    expect(!!kin1 && !!b101, true, 'EX-SMD and EX-BIOD each need a program with a workspace folder');
    const mineRequest = await probeRequest(smd, kin1);
    const theirs = await probeRequest(biod, b101);
    try {
      const viewer = await withStandUser('viewer', (u) => u.asUser(() => portalDb.approvals.query({limit: 10000})));
      const outsider = await withStandUser('outsider', (u) => u.asUser(() => portalDb.approvals.query({limit: 10000})));
      expect(viewer !== null && outsider !== null, true, 'no plain user could be had on this stand');
      const has = (rows: ApprovalRow[] | null, id: string) => (rows ?? []).some((r) => r.id === id);
      expect(`${has(viewer, mineRequest.id)} ${has(viewer, theirs.id)}`, 'true false', 'the viewer: EX-SMD\'s, not EX-BIOD\'s');
      expect(`${has(outsider, mineRequest.id)} ${has(outsider, theirs.id)}`, 'false true', 'the outsider: EX-BIOD\'s only');
      expect(outsider!.filter((r) => others.has(r.portfolio_id)).map((r) => r.title).join('; '), '',
        'the outsider read another portfolio\'s request');
      // read through a link is View only
      const wrote = await withStandUser('viewer', (u) => u.asUser(() =>
        portalDb.approvals.update(mineRequest.id, {result: 'written by a viewer'} as any).then(() => 'allowed', () => 'refused')));
      expect(wrote, 'refused');
    } finally {
      await portalDb.approvals.delete(mineRequest.id).catch(() => undefined);
      await portalDb.approvals.delete(theirs.id).catch(() => undefined);
    }
  });

  test('sync history: the portfolio\'s people read it, nobody else', async () => {
    const [r] = await portalDb.syncLogs.insert({portfolio_id: smd.id, started_on: new Date().toISOString() as any,
      trigger: 'manual', run_by: DG.User.current().id, summary: PROBE});
    try {
      expect(await placeSyncLog(r.id, smd), null, 'the probe run was not linked into EX-SMD\'s _Sync history');
      const [rootId, historyId] = [await rootIdOf(smd), await recordFolderId(smd, SYNC_HISTORY_FOLDER, false)];
      expect(!!historyId && historyId !== rootId, true, 'EX-SMD has a _Sync history folder');
      const links = await folderLinks([rootId!, historyId!]);
      expect(`${links.get(historyId!)?.has(r.id)} ${links.get(rootId!)?.has(r.id)}`, 'true false', 'in _Sync history, not the root');
      const read = (kind: 'viewer' | 'outsider', pf: PortfolioRow) => withStandUser(kind, (u) => u.asUser(async () =>
        (await portalDb.syncLogs.query({filter: `portfolio_id = "${pf.id}"`, limit: 10000}))));
      const count = async (kind: 'viewer' | 'outsider', pf: PortfolioRow) => (await read(kind, pf))?.length;
      const viewerSmd = await read('viewer', smd);
      expect((viewerSmd ?? []).some((x) => x.id === r.id), true, 'the viewer reads EX-SMD\'s sync history through _Sync history');
      expect(await count('viewer', biod), 0, 'the viewer read EX-BIOD\'s sync history');
      expect(await count('outsider', smd), 0, 'the outsider read EX-SMD\'s sync history');
      expect(await count('outsider', bioc), 0, 'the outsider read EX-BIOC\'s sync history');
      const all = await withStandUser('outsider', (u) => u.asUser(() => portalDb.syncLogs.query({limit: 10000})));
      expect((all ?? []).filter((x) => others.has(x.portfolio_id)).length, 0, 'the outsider read another portfolio\'s run');
    } finally {
      await portalDb.syncLogs.delete(r.id).catch(() => undefined);
    }
  });

  test('the root keeps program folders only: no record is linked into it', async () => {
    for (const pf of [smd, bioc, biod])
      expect((await moveRecordsOutOfRoot(pf.code, true))[0], `${pf.code}: no records in the root`);
    // only the portfolio's approvers apply it
    const refused = await withStandUser('viewer', (u) => u.asUser(() =>
      moveRecordsOutOfRoot('EX-SMD', false).then(() => 'allowed', (e) => e?.message ?? String(e))));
    expect(refused, 'Only the approvers of EX-SMD move its records');
  });

  test('compounds: structures through the program, not the compound table', async () => {
    const withCompounds = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`, limit: 1}))[0];
    const links = (kind: 'viewer' | 'outsider') => withStandUser(kind, (u) => u.asUser(async () =>
      (await portalDb.programCompounds.query({filter: `program_id = "${withCompounds.id}"`, limit: 10000}))));
    const v = await links('viewer');
    expect((v ?? []).some((l) => l.smiles), true, 'the viewer sees EXS-KIN1\'s structures');
    expect((await links('outsider'))?.length, 0, 'the outsider read EXS-KIN1\'s compounds');
    const left = (await readable('compound')).filter(hasContent).map((c) => c.code);
    const smdCodes = new Set((v ?? []).map((l) => l.compound_code));
    expect(left.filter((c) => smdCodes.has(c)).join(), '', 'EX-SMD compounds still hold structures in the compound table');
  });

  test('follows: only their user reads them, admins included', async () => {
    const seen = await withStandUser('viewer', async (u) => {
      const [f] = await u.asUser(() => portalDb.followings.insert({user_id: u.id, program_id: kin1.id}));
      try {
        const outsider = await withStandUser('outsider', (o) => o.asUser(() => portalDb.followings.query({limit: 10000})));
        const admin = await portalDb.followings.query({limit: 10000});
        const own = await u.asUser(() => portalDb.followings.query({limit: 10000}));
        const has = (rows: any[] | null) => (rows ?? []).some((x) => x.id === f.id);
        return `${has(own)} ${has(outsider)} ${has(admin)}`;
      } finally {
        // only its author can remove it
        await u.asUser(() => portalDb.followings.delete(f.id)).catch(() => undefined);
      }
    });
    expect(seen, 'true false false', 'the viewer, the outsider, the admin');
  });

  test('the app\'s own check: the viewer\'s portfolios are EX-SMD and EX-BIOC', async () => {
    const mine = await withStandUser('viewer', (u) => u.asUser(() => myPortfolioIds()));
    expect(mine !== null, true, 'the viewer counted as an admin');
    expect(`${mine!.has(smd.id)} ${mine!.has(bioc.id)} ${mine!.has(biod.id)}`, 'true true false');
  });
});
