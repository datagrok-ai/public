import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {isApprover, isContributor, loadAllPortfolios, portfolioByCode} from '../portal';
import {
  accessOf, checkRoleChange, describeRoles, grantChanges, grantsFor, programRoleOf, RoleGrant,
} from '../model/program-roles';
import {checkRunSelection} from '../model/run-selection';
import {namedRolesOf, programAccess, resolvePrincipal, rolesOf, ROW_GRANTS, setProgramRole} from '../program-roles';
import {buildPortfolioWorkspace} from '../workspace';
import {withRestrictedUser} from './restricted-user';
import {_package} from '../package-test';

const roles: RoleGrant[] = [{groupId: 'g-chem', role: 'viewer'}, {groupId: 'u-ann', role: 'contributor'},
  {groupId: 'g-chem', role: 'approver'}, {groupId: 'g-bio', role: 'viewer'}];

category('Portal: program role rules', () => {
  test('a user\'s program role is the strongest any of their groups holds', async () => {
    expect(programRoleOf(roles, ['g-chem']), 'approver');
    expect(programRoleOf(roles, ['u-ann', 'g-bio']), 'contributor');
    expect(programRoleOf(roles, ['g-bio']), 'viewer');
    expect(programRoleOf(roles, ['g-other']), null);
  });

  test('a program role adds to the portfolio\'s groups and never takes away', async () => {
    const none = {approver: false, contributor: false};
    expect(JSON.stringify(accessOf(none, roles, ['g-bio'])), '{"approver":false,"contributor":false,"programRole":"viewer"}');
    expect(JSON.stringify(accessOf(none, roles, ['u-ann'])), '{"approver":false,"contributor":true,"programRole":"contributor"}');
    expect(JSON.stringify(accessOf(none, roles, ['g-chem'])), '{"approver":true,"contributor":true,"programRole":"approver"}');
    expect(accessOf({approver: true, contributor: true}, roles, ['g-bio']).approver, true, 'a viewer role keeps a portfolio approver an approver');
    // what the access means where the app checks it
    const select = (a: {approver: boolean, contributor: boolean}) =>
      checkRunSelection({isApprover: a.approver, isContributor: a.contributor, programStatus: 'active'}).length === 0;
    expect(`${select(accessOf(none, roles, ['g-bio']))} ${select(accessOf(none, roles, ['u-ann']))}`, 'false true');
  });

  test('only the portfolio\'s approvers change program roles, to a known role or none', async () => {
    expect(checkRoleChange({isPortfolioApprover: true, role: 'viewer'}).length, 0);
    expect(checkRoleChange({isPortfolioApprover: true, role: null}).length, 0);
    expect(checkRoleChange({isPortfolioApprover: false, role: 'viewer'})[0], 'Only the portfolio\'s approvers manage program roles');
    expect(checkRoleChange({isPortfolioApprover: true, role: 'owner'}).length, 1);
  });

  test('grants: view to viewers, edit to contributors and approvers; only the groups roles manage change', async () => {
    expect(JSON.stringify(grantsFor(roles)), '[{"groupId":"g-chem","edit":true},{"groupId":"u-ann","edit":true},{"groupId":"g-bio","edit":false}]');
    const current = new Map([['g-chem', false], ['g-bio', false], ['g-gone', true], ['portfolio-viewers', false]]);
    const c = grantChanges(current, grantsFor(roles), ['g-chem', 'u-ann', 'g-bio', 'g-gone']);
    expect(JSON.stringify(c.grant), '[{"groupId":"g-chem","edit":true},{"groupId":"u-ann","edit":true}]', 'an upgrade and a new grant');
    expect(c.revoke.join(), 'g-gone', 'a removed role; the portfolio\'s group is not touched');
    expect(describeRoles([{name: 'Chemists', role: 'viewer'}, {name: 'ann', role: 'approver'}, {name: 'Bio', role: 'viewer'}]),
      'Approver: ann; Viewer: Bio, Chemists');
    expect(ROW_GRANTS, false, 'row grants wait on GROK-21041');
  });
});

category('Portal: program roles', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;

  before(async () => {
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`}))[0];
    await buildPortfolioWorkspace('EX-SMD', false, kin1.id);
    kin1 = await portalDb.programs.get(kin1.id);
  });

  test('a role on EXS-KIN1 gives a user access in the app and shares the program\'s folder; removing it takes both back', async () => {
    const ran = await withRestrictedUser('pprole', async (user) => {
      const folder = (await grok.dapi.getEntities([kin1.space_id!]))[0];
      const grantOf = async () => {
        const p = await grok.dapi.permissions.get(folder) as {view?: DG.Group[], edit?: DG.Group[]};
        return (p.edit ?? []).some((g) => g.id === user.group) ? 'edit' : (p.view ?? []).some((g) => g.id === user.group) ? 'view' : 'none';
      };
      const access = async () => {
        const a = accessOf({approver: false, contributor: false},
          (await rolesOf(kin1.id)).map((r) => ({groupId: r.group_id, role: r.role as any})), [user.group]);
        return `${a.programRole} ${a.contributor} ${a.approver}`;
      };
      const seesFolder = () => user.asUser(async () => !!(await grok.dapi.spaces.find(folder.id).catch(() => null)));
      try {
        expect(await seesFolder(), false, 'a new user does not see the folder');
        const viewer = await setProgramRole(kin1, user.group, 'viewer');
        expect(viewer.line.startsWith('EXS-KIN1: ') && viewer.line.endsWith(' is viewer'), true, viewer.line);
        expect(viewer.notes.some((n) => n.startsWith('workspace folder shared with ') && n.endsWith('(view)')), true, viewer.notes.join('; '));
        expect(`${await access()} ${await grantOf()}`, 'viewer false false view');
        expect(await seesFolder(), true, 'the folder is shared with the viewer');
        await setProgramRole(kin1, user.group, 'contributor');
        expect(`${await access()} ${await grantOf()}`, 'contributor true false edit');
        await setProgramRole(kin1, user.group, 'approver');
        expect(`${await access()} ${await grantOf()}`, 'approver true true edit');
        const row = (await namedRolesOf(kin1.id)).find((r) => r.group_id === user.group)!;
        expect(`${row.role} ${row.granted_by === DG.User.current().id}`, 'approver true', row.name);
        await setProgramRole(kin1, user.group, 'viewer');
        expect(await grantOf(), 'view', 'a downgrade takes edit away');
      } finally {
        await setProgramRole(kin1, user.group, null);
      }
      expect(`${await access()} ${await grantOf()}`, 'null false false none');
      expect(await seesFolder(), false, 'the folder is no longer shared');
      // the throwaway user's row, inactive now, would stay on EXS-KIN1 (and in its activity feed) for good
      for (const r of await portalDb.programRoles.query({filter: `program_id = "${kin1.id}" and group_id = "${user.group}"`}))
        await portalDb.programRoles.delete(r.id).catch(() => undefined);
      return true;
    });
    if (ran === null)
      console.log('skipped: no self-signup on this stand');
  });

  test('a user who is not a portfolio approver may not change roles', async () => {
    await withRestrictedUser('pprole', async (user) => {
      let refused = '';
      await user.asUser(() => setProgramRole(kin1, user.group, 'approver')).catch((e) => refused = e.message);
      expect(refused, 'Only the portfolio\'s approvers manage program roles');
      expect((await rolesOf(kin1.id)).some((r) => r.group_id === user.group), false);
    });
  });

  test('the app asks the program roles: an admin keeps what the portfolio gives', async () => {
    const a = await programAccess(smd, kin1.id);
    expect(`${a.approver} ${a.approver === await isApprover(smd)} ${a.contributor || !await isContributor(smd)}`, 'true true true');
  });

  test('the functions: a role by user login, the roles as a table, refused names', async () => {
    const admin = await resolvePrincipal('admin');
    const r = await grok.functions.call('ProjectPortalCore:setProgramRole',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', principal: 'admin', role: 'contributor'});
    try {
      expect(String(r).startsWith('EXS-KIN1: ') && String(r).includes(' is contributor'), true, String(r));
      const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:programRoles', {portfolio: 'EX-SMD', program: 'EXS-KIN1'});
      expect(df.col('role')!.toList().includes('contributor'), true);
    } finally {
      await setProgramRole(kin1, admin.id, null);
    }
    let refused = '';
    await resolvePrincipal('no such group or user 2026').catch((e) => refused = e.message);
    expect(refused, 'No group or user "no such group or user 2026"');
  });
});
