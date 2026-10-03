import * as grok from 'datagrok-api/grok';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode} from '../portal';
import {setProgramRole} from '../program-roles';
import {ALL_USERS, describeWriteGrants, OPEN_TABLES, planWriteGrants, writerGroup, WRITERS} from '../model/write-grants';
import {currentWriteGrants} from '../write-grants';
import {RestrictedUser, withRestrictedUser} from './restricted-user';
import {_package} from '../package-test';

category('Portal: write permission rules', () => {
  const smd = {code: 'EX-SMD', approvers: 'g-appr', contributors: 'g-contr'};

  test('approvers write roles and parked folders, contributors write approval requests', async () => {
    expect(`${writerGroup(smd, 'program_role')} ${writerGroup(smd, 'parked_folder')} ${writerGroup(smd, 'approval_request')}`,
      'g-appr g-appr g-contr');
    expect(writerGroup({code: 'X', approvers: 'g-appr'}, 'approval_request'), 'g-appr', 'no contributors group: its approvers');
    expect(writerGroup(smd, 'follow'), null, 'follow is not managed');
    expect(Object.keys(OPEN_TABLES).join(), 'follow');
  });

  test('the plan grants each writer group once and takes All users\' Edit and Delete away', async () => {
    const plan = planWriteGrants([smd, {code: 'EX-BIOD', approvers: 'g-appr', contributors: 'g-b'}, {code: 'NONE'}], {
      program_role: [[ALL_USERS, 'View'], [ALL_USERS, 'Edit']],
      approval_request: [[ALL_USERS, 'View'], [ALL_USERS, 'Edit'], ['g-contr', 'Edit']],
      parked_folder: [[ALL_USERS, 'View'], [ALL_USERS, 'Delete'], ['g-other', 'Edit']],
    });
    expect(describeWriteGrants(plan).join('\n'), [
      'program_role: Edit to g-appr (EX-SMD)',
      'parked_folder: Edit to g-appr (EX-SMD)',
      'approval_request: Edit to g-b (EX-BIOD)',
      'program_role: Edit taken from All users',
      'parked_folder: Delete taken from All users',
      'approval_request: Edit taken from All users',
    ].join('\n'));
    const done = planWriteGrants([smd], {program_role: [['g-appr', 'Edit']], approval_request: [['g-contr', 'Edit']],
      parked_folder: [['g-appr', 'Edit'], ['g-other', 'Edit']]});
    expect(done.grant.length + done.revoke.length, 0, 'in place: nothing to do; another group\'s grant is left alone');
    expect(Object.keys(WRITERS).join(), 'program_role,parked_folder,approval_request');
  });
});

/** Runs [f] and says whether the server refused it. */
async function outcome(f: () => Promise<any>): Promise<string> {
  try {
    await f();
    return 'allowed';
  } catch (e: any) {
    const m = String(e?.message ?? e);
    return /insufficient privileges|forbidden/i.test(m) ? 'refused' : `failed: ${m}`;
  }
}

category('Portal: write permissions', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;
  const PROBE = 'write-permission probe';

  before(async () => {
    // loading the portfolios sets the grants for EX-SMD, EX-BIOD and EX-BIOC
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`}))[0];
  });

  const request = (title = PROBE) => ({portfolio_id: smd.id, program_id: kin1.id, kind: 'decision', decision: 'hold', title,
    status: 'withdrawn'} as any);
  const parked = (name: string) => ({portfolio_id: smd.id, space_id: `probe-${name}`, name, kind: 'unclaimed',
    from_parent_id: 'probe', from_path: `/probe/${name}`, parked_in_id: 'probe', parked_path: `/probe/_Parked/${name}`,
    status: 'restored'} as any);

  /** Removes what a test wrote: the user's roles on EXS-KIN1, probe requests and probe parked folders. */
  async function cleanUp(user: RestrictedUser | null) {
    if (user) {
      for (const r of await portalDb.programRoles.query({filter: `program_id = "${kin1.id}" and group_id = "${user.group}"`}))
        await portalDb.programRoles.delete(r.id).catch(() => undefined);
    }
    for (const r of await portalDb.approvalRequests.query({filter: `portfolio_id = "${smd.id}" and title = "${PROBE}"`}))
      await portalDb.approvalRequests.delete(r.id).catch(() => undefined);
    for (const r of await portalDb.parkedFolders.query({filter: `portfolio_id = "${smd.id}" and from_parent_id = "probe"`}))
      await portalDb.parkedFolders.delete(r.id).catch(() => undefined);
  }

  /** Runs [body] with [user] in [groupId], taken out again after. */
  async function asMemberOf<T>(user: RestrictedUser, groupId: string, body: () => Promise<T>): Promise<T> {
    const group = await grok.dapi.groups.find(groupId);
    const personal = await grok.dapi.groups.find(user.group);
    await grok.dapi.groups.includeTo(personal, group);
    try {
      return await body();
    } finally {
      await grok.dapi.groups.excludeFrom(personal, group);
    }
  }

  test('the server grants: View to All users, Edit to EX-SMD\'s approvers or contributors, no Edit or Delete to All users', async () => {
    const grants = await currentWriteGrants();
    const has = (table: string, group: string, perm: string) => grants[table].some(([g, p]) => g === group && p === perm);
    for (const table of Object.keys(WRITERS)) {
      expect(`${table} ${has(table, ALL_USERS, 'View')} ${has(table, ALL_USERS, 'Edit')} ${has(table, ALL_USERS, 'Delete')}`,
        `${table} true false false`);
    }
    expect(`${has('program_role', smd.approvers!, 'Edit')} ${has('parked_folder', smd.approvers!, 'Edit')} ` +
      `${has('approval_request', smd.contributors!, 'Edit')}`, 'true true true');
  });

  test('a plain user can neither insert, update nor delete program roles, approval requests or parked folders', async () => {
    const ran = await withRestrictedUser('ppwrite', async (user) => {
      // rows the admin made, which the user then tries to change
      const [role] = await portalDb.programRoles.insert({program_id: kin1.id, group_id: user.group, role: 'viewer', active: true});
      const [req] = await portalDb.approvalRequests.insert(request());
      const [folder] = await portalDb.parkedFolders.insert(parked('plain-admin'));
      try {
        const tried = [
          await outcome(() => user.asUser(() => setProgramRole(kin1, user.group, 'approver'))),
          await outcome(() => user.asUser(() => portalDb.programRoles.update(role.id, {role: 'approver'}))),
          await outcome(() => user.asUser(() => portalDb.programRoles.delete(role.id))),
          await outcome(() => user.asUser(() => portalDb.approvalRequests.insert(request()))),
          await outcome(() => user.asUser(() => portalDb.approvalRequests.update(req.id, {status: 'approved'} as any))),
          await outcome(() => user.asUser(() => portalDb.approvalRequests.delete(req.id))),
          await outcome(() => user.asUser(() => portalDb.parkedFolders.insert(parked('plain')))),
          await outcome(() => user.asUser(() => portalDb.parkedFolders.update(folder.id, {status: 'parked'} as any))),
          await outcome(() => user.asUser(() => portalDb.parkedFolders.delete(folder.id))),
        ];
        // the app's own check answers first for the role; the console call under it is the server's
        expect(tried[0].includes('Only the portfolio\'s approvers manage program roles'), true, tried[0]);
        expect(tried.slice(1).join(), Array(8).fill('refused').join());
        const direct = await outcome(() => user.asUser(() => portalDb.programRoles.insert({program_id: kin1.id,
          group_id: user.group, role: 'approver', active: true})));
        expect(direct, 'refused', 'making oneself an approver from the console');
        expect(`${(await portalDb.programRoles.get(role.id)).role} ${(await portalDb.approvalRequests.get(req.id)).status} ` +
          `${(await portalDb.parkedFolders.get(folder.id)).status}`, 'viewer withdrawn restored', 'nothing changed');
        // still readable: the app lists them for everyone
        expect(await user.asUser(() => portalDb.programRoles.count({filter: `id = "${role.id}"`})), 1);
      } finally {
        await cleanUp(user);
      }
      return true;
    });
    if (ran === null)
      console.log('skipped: no self-signup on this stand');
  });

  test('an EX-SMD approver still manages roles, reviews requests and parks; a contributor proposes but does not manage roles', async () => {
    const ran = await withRestrictedUser('ppwrite', async (user) => {
      try {
        await asMemberOf(user, smd.approvers!, async () => {
          const changed = await user.asUser(() => setProgramRole(kin1, user.group, 'viewer'));
          expect(changed.line.endsWith(' is viewer'), true, changed.line);
          const role = (await portalDb.programRoles.query({filter: `program_id = "${kin1.id}" and group_id = "${user.group}"`}))[0];
          expect(await outcome(() => user.asUser(() => portalDb.programRoles.update(role.id, {active: false}))), 'allowed');
          const [req] = await user.asUser(() => portalDb.approvalRequests.insert({...request(), status: 'pending'}));
          expect(await outcome(() => user.asUser(() => portalDb.approvalRequests.update(req.id,
            {status: 'withdrawn', result: 'probe'} as any))), 'allowed', 'a review or a withdrawal writes the request');
          const [folder] = await user.asUser(() => portalDb.parkedFolders.insert(parked('approver')));
          expect(await outcome(() => user.asUser(() => portalDb.parkedFolders.update(folder.id, {status: 'restored'} as any))),
            'allowed', 'unparking');
          // Delete is no one's but the admins': the app never deletes these rows
          expect(await outcome(() => user.asUser(() => portalDb.parkedFolders.delete(folder.id))), 'refused');
        });
        await asMemberOf(user, smd.contributors!, async () => {
          expect(await outcome(() => user.asUser(() => portalDb.approvalRequests.insert(request()))), 'allowed', 'proposing');
          expect(await outcome(() => user.asUser(() => portalDb.programRoles.insert({program_id: kin1.id,
            group_id: user.group, role: 'approver', active: true}))), 'refused', 'a contributor making themselves approver');
        });
      } finally {
        await cleanUp(user);
      }
      return true;
    });
    if (ran === null)
      console.log('skipped: no self-signup on this stand');
  });
});
