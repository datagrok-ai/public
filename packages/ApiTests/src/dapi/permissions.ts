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

  test('the owner keeps every privilege after every grant is revoked', async () => {
    const own = await asGrantee(() => grok.dapi.scripts.save(DG.Script.create(
      `//name: ApiTestsPerm${DG.Utils.randomString(6)}\n//language: javascript\nlet x = 1;`)));
    try {
      await asGrantee(() => grok.dapi.permissions.grant(own, allUsers, true));
      await asGrantee(() => grok.dapi.permissions.revoke(own, allUsers));
      expect((await held(own)).join(), BASIC.join());
    }
    finally {
      await grok.dapi.scripts.delete(own);
    }
  });

  test('a space shared for viewing cannot be changed, deleted or shared further by the grantee', async () => {
    const space = await grok.dapi.spaces.createRootSpace(`ApiTestsPerm${DG.Utils.randomString(6)}`);
    try {
      await grok.dapi.permissions.grant(space, grantee.group, false);
      expect((await held(space)).join(), 'View');
      const edit = await asGrantee(() => thrown(async () => {
        const s = await grok.dapi.projects.find(space.id);
        s.description = 'edited by the grantee';
        await grok.dapi.projects.save(s);
      }));
      expect(edit != null, true, 'the grantee saved an edit');
      expect(await asGrantee(() => thrown(() => grok.dapi.spaces.delete(space))) != null, true, 'the grantee deleted the space');
      expect(await grok.dapi.projects.find(space.id) != null, true, 'the space is gone');
      const reshare = await asGrantee(() => thrown(() => grok.dapi.permissions.grant(space, allUsers, false)));
      expect(String(reshare?.message ?? reshare).includes('You don\'t have a permission to share this object'), true,
        `the grantee shared the space further: ${reshare}`);
    }
    finally {
      await grok.dapi.permissions.revoke(space, grantee.group);
      await grok.dapi.spaces.delete(space);
    }
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

  test('a database connection: View and use gives query/schema/files but no edit; Full access gives all', async () => {
    const name = `ApiTestsPerm${DG.Utils.randomString(6)}`;
    const connection = await grok.dapi.connections.save(DG.DataConnection.create(name,
      {dataSource: 'Postgres', server: 'localhost', port: 1, db: 'apitests', login: 'apitests', password: 'apitests'}));
    const use = ['View', 'Edit', 'Delete', 'Share',
      'DataConnection.Query', 'DataConnection.GetSchema', 'DataConnection.ListFiles'];
    try {
      await grok.dapi.permissions.grant(connection, grantee.group, false);
      expect((await held(connection, use)).join(),
        'View,DataConnection.Query,DataConnection.GetSchema,DataConnection.ListFiles');
      const edit = await asGrantee(() => thrown(async () => {
        const c = await grok.dapi.connections.find(connection.id);
        c.parameters['server'] = 'changed.example';
        await grok.dapi.connections.save(c);
      }));
      expect(edit != null, true, 'the grantee saved the connection');
      expect((await grok.dapi.connections.find(connection.id)).parameters['server'], 'localhost',
        'the grantee changed the connection parameters');
      await grok.dapi.permissions.revoke(connection, grantee.group);
      await grok.dapi.permissions.grant(connection, grantee.group, true);
      expect((await held(connection)).join(), BASIC.join());
    }
    finally {
      await grok.dapi.permissions.revoke(connection, grantee.group);
      await grok.dapi.connections.delete(connection);
    }
  });

  test('a database connection shared on its own shares none of its queries', async () => {
    const name = `ApiTestsPerm${DG.Utils.randomString(6)}`;
    const connection = await grok.dapi.connections.save(DG.DataConnection.create(name,
      {dataSource: 'Postgres', server: 'localhost', port: 1, db: 'apitests', login: 'apitests', password: 'apitests'}));
    const q = await grok.dapi.queries.save(connection.query(name, 'select 1 as x'));
    try {
      await grok.dapi.permissions.grant(connection, grantee.group, false);
      expect((await held(connection, ['View'])).join(), 'View', 'the connection');
      expect((await held(q, ['View', 'DataQuery.Execute'])).join(), '', 'its query');
    }
    finally {
      await grok.dapi.permissions.revoke(connection, grantee.group);
      await grok.dapi.queries.delete(q);
      await grok.dapi.connections.delete(connection);
    }
  });

  test('a query: the View-and-use grantee cannot edit its SQL or re-share it', async () => {
    const dg = await grok.dapi.connections.filter('name = "Datagrok"').first();
    const q = await grok.dapi.queries.save(dg.query(`ApiTestsPerm${DG.Utils.randomString(6)}`, 'select 1 as x'));
    try {
      await grok.dapi.permissions.grant(q, grantee.group, false);
      expect((await held(q, ['View', 'Edit', 'Delete', 'Share', 'DataQuery.Execute'])).join(), 'View,DataQuery.Execute');
      const edit = await asGrantee(() => thrown(async () => {
        const saved = await grok.dapi.queries.find(q.id);
        saved.query = 'select 2 as x';
        await grok.dapi.queries.save(saved);
      }));
      expect(edit != null, true, 'the grantee edited the query SQL');
      expect((await grok.dapi.queries.find(q.id)).query, 'select 1 as x', 'the SQL changed');
      const reshare = await asGrantee(() => thrown(() => grok.dapi.permissions.grant(q, allUsers, false)));
      expect(String(reshare?.message ?? reshare).includes('You don\'t have a permission to share this object'), true,
        `the grantee shared the query further: ${reshare}`);
    }
    finally {
      await grok.dapi.permissions.revoke(q, grantee.group);
      await grok.dapi.queries.delete(q);
    }
  });

  test('a table: View and use lets the grantee read the data but not change it; Full access gives all', async () => {
    const df = DG.DataFrame.fromCsv('x,y\n1,2\n3,4');
    df.name = `ApiTestsPerm${DG.Utils.randomString(6)}`;
    const id = await grok.dapi.tables.uploadDataFrame(df);
    const table = await grok.dapi.tables.find(id);
    try {
      await grok.dapi.permissions.grant(table, grantee.group, false);
      expect((await held(table, ['View', 'Edit', 'Delete', 'Share', 'Table.ReadData'])).join(), 'View,Table.ReadData');
      expect(await asGrantee(() => grok.dapi.tables.getTable(id).then((t) => t.rowCount)), 2,
        'the grantee cannot read the table data');
      await grok.dapi.permissions.revoke(table, grantee.group);
      await grok.dapi.permissions.grant(table, grantee.group, true);
      expect((await held(table)).join(), BASIC.join());
    }
    finally {
      await grok.dapi.permissions.revoke(table, grantee.group);
      await grok.dapi.tables.delete(await grok.dapi.tables.find(id));
    }
  });

  test('a layout: the View-and-use grantee sees it but cannot save or delete it; Full access gives all', async () => {
    const view = grok.shell.addTableView(grok.data.demo.demog(20));
    const draft = view.saveLayout();
    draft.name = `ApiTestsPerm${DG.Utils.randomString(6)}`;
    view.close();
    const layout = await grok.dapi.layouts.save(draft);
    try {
      await grok.dapi.permissions.grant(layout, grantee.group, false);
      expect((await held(layout)).join(), 'View');
      const edit = await asGrantee(() => thrown(async () => {
        const l = await grok.dapi.layouts.find(layout.id);
        l.name = `${l.name}X`;
        await grok.dapi.layouts.save(l);
      }));
      expect(edit != null, true, 'the grantee saved the layout');
      expect(await asGrantee(() => thrown(() => grok.dapi.layouts.delete(layout))) != null, true,
        'the grantee deleted the layout');
      await grok.dapi.permissions.revoke(layout, grantee.group);
      await grok.dapi.permissions.grant(layout, grantee.group, true);
      expect((await held(layout)).join(), BASIC.join());
    }
    finally {
      await grok.dapi.permissions.revoke(layout, grantee.group);
      await grok.dapi.layouts.delete(layout);
    }
  });

  test('a model: the View-and-use grantee sees it but cannot change, delete or share it', async () => {
    const name = `ApiTestsPerm${DG.Utils.randomString(6)}`;
    const created = await fetch(`${grok.dapi.root}/ml`, {method: 'POST', credentials: 'include',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({'#type': 'PredictiveModelInfo', name, friendlyName: name, input: [], output: []})});
    expect(created.status, 200, 'the model was not created');
    const model = await grok.dapi.models.filter(`name = "${name}"`).first();
    try {
      await grok.dapi.permissions.grant(model, grantee.group, false);
      expect((await held(model)).join(), 'View');
      const edit = await asGrantee(() => thrown(async () => {
        const m = await grok.dapi.models.find(model.id);
        m.friendlyName = 'edited by the grantee';
        await grok.dapi.models.save(m);
      }));
      expect(edit != null, true, 'the grantee edited the model');
      expect(await asGrantee(() => thrown(() => grok.dapi.models.delete(model))) != null, true,
        'the grantee deleted the model');
      const reshare = await asGrantee(() => thrown(() => grok.dapi.permissions.grant(model, allUsers, false)));
      expect(String(reshare?.message ?? reshare).includes('You don\'t have a permission to share this object'), true,
        `the grantee shared the model further: ${reshare}`);
      await grok.dapi.permissions.revoke(model, grantee.group);
      expect((await held(model, ['View'])).join(), '', 'after the revoke');
    }
    finally {
      await grok.dapi.permissions.revoke(model, grantee.group);
      await grok.dapi.models.delete(model);
    }
  });

  test('a share with notifications reaches the grantee; without them it does not', async () => {
    const viewId = (await (await fetch(`${grok.dapi.root}/privileges/types/Project/permissions`)).json())
      .find((p: any) => p.name === 'View').id;
    // No JS wrapper carries notify/message, so the dialog's endpoint is called directly.
    const makeScript = () => grok.dapi.scripts.save(DG.Script.create(
      `//name: ApiTestsPerm${DG.Utils.randomString(6)}\n//language: javascript\nlet x = 1;`));
    const share = async (entity: _DG.Entity, notify: boolean) => {
      const project = await (await fetch(`${grok.dapi.root}/projects/for_entity/${entity.id}`)).json();
      const res = await fetch(`${grok.dapi.root}/privileges/permissions/list?notify=${notify}`, {method: 'POST',
        credentials: 'include', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify([{entityId: project.id, userGroup: {id: grantee.group.id}, permission: {id: viewId}}])});
      expect(res.ok, true, `share failed: ${res.status}`);
    };
    const notified = () => asGrantee(() => grok.dapi.users.notifications.forCurrentUser().list({pageSize: 100}));
    const notifying = await makeScript();
    const silent = await makeScript();
    try {
      const before = (await notified()).map((n) => n.id);
      await share(notifying, true);
      const after = await notified();
      const fresh = after.filter((n) => !before.includes(n.id));
      expect(fresh.some((n) => (n.text ?? '').includes(notifying.id)), true,
        'the grantee got no notification naming the shared entity');

      const before2 = (await notified()).map((n) => n.id);
      await share(silent, false);
      expect((await notified()).filter((n) => !before2.includes(n.id)).length, 0,
        'a silent share still notified the grantee');
    }
    finally {
      await grok.dapi.permissions.revoke(notifying, grantee.group);
      await grok.dapi.permissions.revoke(silent, grantee.group);
      await grok.dapi.scripts.delete(notifying);
      await grok.dapi.scripts.delete(silent);
    }
  });
});
