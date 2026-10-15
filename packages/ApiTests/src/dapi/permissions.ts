import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown} from './domain-lifecycle';

/* What a share gives the account it reaches, read as that account. A user cannot be deleted, so the
   grantee is one account kept on the stand and reused; it acts through a session minted from its
   developer key, which the JS API has no call for. */
const GRANTEE = 'apitestsgrantee';
const BASIC = ['View', 'Edit', 'Delete', 'Share'];

category('Dapi: permissions', () => {
  let grantee: _DG.User;
  let token: string;
  let allUsers: _DG.Group;
  let script: _DG.Script;
  let query: _DG.DataQuery;

  const asGrantee = async <T>(action: () => Promise<T>): Promise<T> => {
    grok.dapi.impersonationToken = token;
    try {
      return await action();
    }
    finally {
      grok.dapi.impersonationToken = null;
    }
  };

  const held = (e: _DG.Entity, privileges: string[] = BASIC): Promise<string[]> => asGrantee(async () => {
    const result: string[] = [];
    for (const p of privileges) {
      if (await grok.dapi.permissions.check(e, p))
        result.push(p);
    }
    return result;
  });

  before(async () => {
    grantee = await grok.dapi.users.filter(`login = "${GRANTEE}"`).include('group').first();
    if (!grantee) {
      const user = DG.User.create();
      user.login = GRANTEE;
      user.firstName = GRANTEE;
      user.status = DG.USER_STATUS.STATUS_ACTIVE;
      await grok.dapi.users.save(user);
      grantee = await grok.dapi.users.filter(`login = "${GRANTEE}"`).include('group').first();
    }
    const keyText = await (await fetch(`${grok.dapi.root}/users/${grantee.id}/dev_key`)).text();
    const key = keyText.startsWith('"') ? JSON.parse(keyText) : keyText;
    // credentials: 'omit' keeps the server's answer from replacing the runner's own session cookie
    token = (await (await fetch(`${grok.dapi.root}/users/login/dev`,
      {method: 'POST', credentials: 'omit', headers: {Authorization: `Dev ${key}`}})).json()).token;
    expect(!!token, true, `no session for ${GRANTEE}`);
    allUsers = await grok.dapi.groups.filter('friendlyName = "All users"').first();
    const name = `ApiTestsPerm${DG.Utils.randomString(6)}`;
    script = await grok.dapi.scripts.save(DG.Script.create(`//name: ${name}\n//language: javascript\nlet x = 1;`));
    const connection = await grok.dapi.connections.filter('name = "Datagrok"').first();
    query = await grok.dapi.queries.save(connection.query(name, 'select 1 as x'));
  });

  after(async () => {
    if (script)
      await grok.dapi.scripts.delete(script);
    if (query)
      await grok.dapi.queries.delete(query);
    expect(await grok.dapi.scripts.find(script.id) == null, true, 'the script is left on the server');
  });

  test('View and use: the grantee sees and uses, and may not change, delete or share', async () => {
    expect((await held(script)).join(), '', 'before the grant');
    try {
      await grok.dapi.permissions.grant(script, grantee.group, false);
      await grok.dapi.permissions.grant(query, grantee.group, false);
      expect((await held(script)).join(), 'View');
      expect((await held(query, ['View', 'Edit', 'DataQuery.Execute'])).join(), 'View,DataQuery.Execute');
      expect(await asGrantee(() => grok.dapi.scripts.find(script.id)) != null, true, 'the grantee cannot read the script');

      const edit = await asGrantee(async () => thrown(async () => {
        const s = await grok.dapi.scripts.find(script.id);
        s.description = 'edited by the grantee';
        await grok.dapi.scripts.save(s);
      }));
      expect(edit != null, true, 'the grantee saved an edit');
      expect(await asGrantee(() => thrown(() => grok.dapi.scripts.delete(script))) != null, true, 'the grantee deleted the script');
      const reshare = await asGrantee(() => thrown(() => grok.dapi.permissions.grant(script, allUsers, false)));
      expect(String(reshare?.message ?? reshare).includes('You don\'t have a permission to share this object'), true,
        `the grantee shared the script further: ${reshare}`);
    }
    finally {
      await grok.dapi.permissions.revoke(script, grantee.group);
      await grok.dapi.permissions.revoke(query, grantee.group);
    }
  });

  test('Full access: every privilege', async () => {
    try {
      await grok.dapi.permissions.grant(script, grantee.group, true);
      expect((await held(script)).join(), BASIC.join());
    }
    finally {
      await grok.dapi.permissions.revoke(script, grantee.group);
    }
  });

  test('a revoke takes the access back', async () => {
    await grok.dapi.permissions.grant(query, grantee.group, false);
    expect((await held(query, ['View', 'DataQuery.Execute'])).join(), 'View,DataQuery.Execute', 'after the grant');
    await grok.dapi.permissions.revoke(query, grantee.group);
    expect((await held(query, ['View', 'DataQuery.Execute'])).join(), '');
    expect(await asGrantee(() => grok.dapi.queries.find(query.id)) == null, true, 'the grantee still reads the query');
  });

  test('a grant to a group reaches its members, and leaves with them', async () => {
    const group = await grok.dapi.groups.createNew(`apitests-perm-${DG.Utils.randomString(6)}`);
    try {
      await grok.dapi.groups.addMember(group, grantee.group);
      await grok.dapi.permissions.grant(script, group, false);
      expect((await held(script)).join(), 'View', 'as a member');
      await grok.dapi.groups.removeMember(group, grantee.group);
      expect((await held(script)).join(), '', 'after leaving the group');
    }
    finally {
      await grok.dapi.permissions.revoke(script, group);
      await grok.dapi.groups.delete(group);
    }
  });

  test('a grant to All users reaches every account', async () => {
    try {
      await grok.dapi.permissions.grant(script, allUsers, false);
      expect((await held(script)).join(), 'View');
    }
    finally {
      await grok.dapi.permissions.revoke(script, allUsers);
    }
    expect((await held(script)).join(), '', 'after the revoke');
  });

  test('a space shared for viewing passes View and use to its content, and takes no new content from the grantee', async () => {
    const space = await grok.dapi.spaces.createRootSpace(`ApiTestsPerm${DG.Utils.randomString(6)}`);
    try {
      await grok.dapi.spaces.id(space.id).addEntity(query.id);
      await grok.dapi.permissions.grant(space, grantee.group, false);
      expect((await held(query, ['View', 'Edit', 'DataQuery.Execute'])).join(), 'View,DataQuery.Execute');
      expect(await asGrantee(() => thrown(() => grok.dapi.spaces.id(space.id).addEntity(script.id, true))) != null, true,
        'the grantee linked a script into the space');
      await grok.dapi.permissions.revoke(space, grantee.group);
      expect((await held(query)).join(), '', 'after the revoke');
    }
    finally {
      await grok.dapi.spaces.delete(space);
    }
  });
});
