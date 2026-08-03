import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Resolves to the thrown error (null when the action succeeded) — shared by both
// categories in this file for typed-error assertions.
async function thrown(action: () => Promise<any>): Promise<any> {
  try {
    await action();
    return null;
  } catch (e) {
    return e;
  }
}

// WO-4a (GROK-20601): schema lifecycle handle (grok.dapi.domains.schema(name):
// manifest/apply/audit/delete + createSchema) and table-scoped admin surface
// (grants/grant/revoke, shareColumn/restrictColumn/restoreColumnVisibility).
// The lifecycle round-trip needs the CreateDomainSchema privilege and skips
// cleanly without it; grants/column tests ride the package's own 'apitests'
// schema and grant only to 'All users' with finally-cleanup.
category('Dapi: domain lifecycle', () => {
  const allUsers = async () => (await grok.dapi.groups.filter('friendlyName = "All users"').first());

  test('schema lifecycle: create → apply (dryRun first) → data → audit → typed errors → delete', async () => {
    const name = `zzlc${`${Date.now()}`.slice(-8)}`;
    try {
      await grok.dapi.domains.createSchema(name, {friendlyName: 'Lifecycle probe'});
    } catch (e: any) {
      if (e instanceof DG.DomainError && (e.code === 'forbidden' || e.status === 403)) {
        console.log('skipped: no CreateDomainSchema privilege');
        return;
      }
      throw e;
    }
    const handle = grok.dapi.domains.schema(name);
    try {
      const tables = {things: {columns: {name: {type: 'string', required: true}, qty: {type: 'int'}}}};
      const plan = await handle.apply({tables}, {dryRun: true});
      expect(plan != null, true, 'dryRun must return the change plan');
      await handle.apply({tables});

      const manifest = await handle.manifest();
      expect(manifest['tables']?.['things'] != null, true, JSON.stringify(manifest));

      // Data flows through the standard table client immediately after apply.
      const things = grok.dapi.domains.table(`${name}.things`);
      const [ins] = await things.insert({name: 'w1', qty: 2});
      expect(ins.created, true);
      expect((await things.query({filter: 'name = "w1"'})).length, 1);

      // Whole-schema audit: ddl registry events + the row insert, newest first.
      const audit = await handle.audit({limit: 100});
      expect(audit.some((a) => a.op === 'ddl'), true, 'ddl events expected');
      expect(audit.some((a) => a.op === 'insert'), true, 'row events expected');

      // Stale ifVersion → typed version conflict.
      const stale = await thrown(() => handle.apply({tables, ifVersion: '999'}));
      expect(stale instanceof DG.DomainVersionConflictError, true,
        `expected DomainVersionConflictError, got ${stale?.constructor?.name}: ${stale?.message}`);

      // Destructive apply without confirmation → typed error carrying the plan.
      const destructive = await thrown(() => handle.apply({dropTables: ['things']}));
      expect(destructive instanceof DG.DomainError, true,
        `expected DomainError, got ${destructive?.constructor?.name}: ${destructive?.message}`);
      expect(destructive.code, 'destructive-confirmation-required');
      expect(destructive.body['plan'] != null, true, JSON.stringify(destructive?.body));
    } finally {
      await handle.delete();
    }
    const schemas = await grok.dapi.domains.schemas.list();
    expect(schemas.some((s) => s.name === name), false, 'deleted schema must not be listed');
  });

  test('table grants: list, grant, revoke round-trip', async () => {
    const items = grok.dapi.domains.table('apitests.item');
    const group = await allUsers();
    const before = (await items.grants()).length;
    await items.grant(group.id, 'View');
    try {
      const listed = await items.grants();
      expect(listed.some((g) => g.group.id === group.id && g.permission === 'View'), true,
        JSON.stringify(listed));
    } finally {
      await items.revoke(group.id, 'View');
    }
    const after = await items.grants();
    expect(after.some((g) => g.group.id === group.id && g.permission === 'View'), false,
      'revoked grant must disappear from a second read');
    expect(after.length, before);
  });

  test('column security: shareColumn, idempotent restrictColumn, restoreColumnVisibility', async () => {
    const items = grok.dapi.domains.table('apitests.item');
    const group = await allUsers();
    // Sharing to All users first means a mid-test failure never hides the column.
    const shared = await items.shareColumn('name', group.id);
    try {
      expect(shared.id != null, true, JSON.stringify(shared));
      expect(shared.name, 'apitests.item.name');
      // The grant half really landed. Per-column property schemas have NO
      // entities rows (registry-entities rule), so grok.dapi.getEntities /
      // permissions cannot read them and no typed JS accessor takes arbitrary
      // registry ids (the WO-4a security boundary) — read the domains-native
      // grants route directly, the ApiTests authenticated-fetch pattern
      // (cf. user-settings-storage.ts).
      const res = await fetch(`${grok.dapi.root}/domains/grants/${shared.id}`,
        {credentials: 'include'});
      const granted = await res.json();
      expect(Array.isArray(granted), true, JSON.stringify(granted));
      expect(granted.some((g: any) => g.group?.id === group.id && g.permission === 'View'), true,
        'shareColumn must grant View on the per-column schema');
      const again = await items.restrictColumn('name');
      expect(again.id, shared.id, 'restrict is idempotent — same per-column schema');
    } finally {
      await items.restoreColumnVisibility('name');
    }
    // Restored: the core schema serves the column again — a query NAMING it
    // succeeding is the assertion (a restricted column would 400, no-oracle).
    const rows = await items.query({limit: 1, columns: ['name']});
    expect(Array.isArray(rows), true);
  });

  test('schema grants: list, grant, revoke round-trip', async () => {
    const schema = grok.dapi.domains.schema('apitests');
    const group = await allUsers();
    const before = (await schema.grants()).length;
    await schema.grant(group.id, 'View');
    try {
      const listed = await schema.grants();
      expect(listed.some((g) => g.group.id === group.id && g.permission === 'View'), true,
        JSON.stringify(listed));
    } finally {
      await schema.revoke(group.id, 'View');
    }
    const after = await schema.grants();
    expect(after.some((g) => g.group.id === group.id && g.permission === 'View'), false,
      'revoked schema grant must disappear from a second read');
    expect(after.length, before);
  });
});

// Permission-driven capabilities (ui-js-api WO-2): DomainTableClient.capabilities()
// composes server-truth permission probes on the securing entity with the
// writable-column mirror of column security, cached until a grant change. The flip
// test signs up a throwaway restricted user (internal auth) and runs the probes
// under THAT session against a real grant round-trip — no mocks; the user is
// blocked in finally (users are not API-deletable). It skips cleanly where
// self-signup is unavailable (SSO-only or email-confirm setups).
category('Dapi: domain capabilities', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const hidden = () => grok.dapi.domains.table('apitests.hidden_item');

  test('admin sees full capabilities; cache survives reads, drops on invalidation', async () => {
    const caps = await items().capabilities();
    expect(caps.canView && caps.canInsert && caps.canEdit && caps.canDelete && caps.canShareTable,
      true, JSON.stringify(caps));
    expect(caps.securityMode, 'row');
    expect(caps.audit, true);
    expect(caps.hasBusinessKey, true);
    for (const c of ['sku', 'name', 'quantity'])
      expect(caps.writableColumns.includes(c), true, `writableColumns must include ${c}: ${JSON.stringify(caps.writableColumns)}`);
    grok.dapi.domains.invalidateUiCaches();
    expect((await items().capabilities()).canEdit, true, 'recompute after invalidation must succeed');
  });

  test('unknown table rejects with a typed validation error', async () => {
    const e = await thrown(() => grok.dapi.domains.table('apitests.nosuch').capabilities());
    expect(e instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${e?.constructor?.name}: ${e?.message}`);
  });

  test('capabilities flip on a real grant round-trip for a restricted user', async () => {
    // Same-origin sessions authenticate via the 'auth' COOKIE — DelegatingHttpClient
    // deliberately attaches no Authorization header when dapi.root starts with the
    // origin — so impersonating the restricted user means swapping the cookie
    // (grok.dapi.token alone would be a no-op here). Where the cookie is unreadable
    // (HttpOnly deployments) the restore path could only DELETE it and would take the
    // admin session down with it: skip before signing anyone up or swapping anything.
    const adminCookie = document.cookie.match(/(?:^|; )auth=([^;]*)/)?.[1] ?? null;
    if (adminCookie == null) {
      console.log('skipped: the auth cookie is not readable (HttpOnly) — impersonation is not restorable');
      return;
    }
    const adminToken = grok.dapi.token;
    const stamp = `${Date.now()}${Math.floor(Math.random() * 1e4)}`;
    const login = `wo2caps${stamp}`;
    // Self-signup issues the restricted session token the probes run under.
    const signup = await (await fetch(`${grok.dapi.root}/users/signup`, {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({login: login, email: `${login}@test.datagrok.ai`,
        password: btoa(`Wo2-${stamp}`), firstName: 'WO2', lastName: 'CapsProbe'}),
    })).json();
    if (!signup?.token) {
      console.log(`skipped: self-signup unavailable (${signup?.reason ?? 'no token'})`);
      return;
    }
    const setAuth = (token: string) => {
      document.cookie = `auth=${encodeURIComponent(token)}; path=/`;
      grok.dapi.token = token;
    };
    // Restore exactly what was captured — the cookie VALUE as it was written
    // (re-encoding a decoded value is not always the same string) and the token
    // the client held, which the cookie need not carry.
    const restoreAdmin = () => {
      document.cookie = `auth=${adminCookie}; path=/`;
      grok.dapi.token = adminToken;
    };
    const user = await grok.dapi.users.filter(`login = "${login}"`).include('group').first();
    expect(user != null, true, 'signed-up user not found');
    // Runs capabilities() under the restricted user's session: the permission probes
    // are per-call server requests, so they follow the swapped auth; the client-side
    // cache is dropped first because its user key tracks the SESSION user.
    const asUser = async () => {
      setAuth(signup.token);
      try {
        grok.dapi.domains.invalidateUiCaches();
        return await items().capabilities();
      } finally {
        restoreAdmin();
      }
    };
    try {
      const before = await asUser();
      expect(before.canInsert, false, `no grant yet, canInsert must deny: ${JSON.stringify(before)}`);
      expect(before.canEdit, false, 'no grant yet, canEdit must deny');
      // KNOWN LIMITATION, asserted so it cannot drift unnoticed: writableColumns is
      // admin-fast-pathed off the SESSION's Auth.adminMode, which a cookie swap does
      // not change — under an admin session it stays full even for the ungranted
      // user. So the flip below is proven by canInsert (canEdit only rides on it);
      // a genuinely restricted browser session is what makes writableColumns empty.
      expect(before.writableColumns.length > 0, true,
        `writableColumns is admin-fast-pathed here: ${JSON.stringify(before.writableColumns)}`);
      await items().grant(user!.group.id, 'Edit'); // drops the caches automatically
      const after = await asUser();
      expect(after.canInsert, true, `Edit grant must flip canInsert: ${JSON.stringify(after)}`);
      expect(after.canEdit, true, 'Edit grant must flip canEdit');
      await items().revoke(user!.group.id, 'Edit');
      expect((await asUser()).canEdit, false, 'revoke must flip canEdit back');
    } finally {
      restoreAdmin();
      try {
        await items().revoke(user!.group.id, 'Edit');
      } catch (_) { /* already revoked */ }
      // Users are not API-deletable — blocking is the platform cleanup for test
      // users (revokes their sessions too, including the signup session).
      const blocked = await fetch(`${grok.dapi.root}/users/block`, {
        method: 'POST', credentials: 'include', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({'#type': 'User', 'id': user!.id, 'login': login}),
      });
      expect(blocked.ok, true, `test user ${login} not blocked: ${blocked.status}`);
      grok.dapi.domains.invalidateUiCaches();
    }
  });

  test('a table grant does not reach the rows of a defaultRowVisibility:"none" table', async () => {
    // capabilities() ANDs canEdit/canDelete with "a grant on the securing TABLE
    // reaches its ROWS" — false exactly for this fixture. The client-side term
    // cannot be observed from here (Auth.adminMode short-circuits it, and a
    // cookie swap does not change the CLIENT's user), so what this pins is the
    // SERVER semantics the term mirrors, side by side with apitests.item whose
    // rows default to table visibility. STANDING LIMITATION: the client term
    // itself needs a genuinely restricted browser session (ui-js-api progress).
    const e = await thrown(() => hidden().capabilities());
    if (e instanceof DG.DomainValidationError) {
      console.log('skipped: apitests.hidden_item is not deployed (schema.json < 1.2.0)');
      return;
    }
    expect(e, null, `apitests.hidden_item capabilities failed: ${e?.message}`);
    const adminCookie = document.cookie.match(/(?:^|; )auth=([^;]*)/)?.[1] ?? null;
    if (adminCookie == null) {
      console.log('skipped: the auth cookie is not readable (HttpOnly) — impersonation is not restorable');
      return;
    }
    const adminToken = grok.dapi.token;
    const stamp = `${Date.now()}${Math.floor(Math.random() * 1e4)}`;
    const login = `wo4fvis${stamp}`;
    const signup = await (await fetch(`${grok.dapi.root}/users/signup`, {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({login: login, email: `${login}@test.datagrok.ai`,
        password: btoa(`Wo4f-${stamp}`), firstName: 'WO4F', lastName: 'VisProbe'}),
    })).json();
    if (!signup?.token) {
      console.log(`skipped: self-signup unavailable (${signup?.reason ?? 'no token'})`);
      return;
    }
    const restoreAdmin = () => {
      document.cookie = `auth=${adminCookie}; path=/`;
      grok.dapi.token = adminToken;
    };
    const asUser = async <T>(action: () => Promise<T>): Promise<T> => {
      document.cookie = `auth=${encodeURIComponent(signup.token)}; path=/`;
      grok.dapi.token = signup.token;
      try {
        return await action();
      } finally {
        restoreAdmin();
      }
    };
    // Every server call is labelled: this test touches two tables, two sessions
    // and four grants, and a bare 'forbidden' would not say which.
    const step = async <T>(label: string, action: () => Promise<T>): Promise<T> => {
      try {
        return await action();
      } catch (x: any) {
        throw new Error(`${label}: ${x?.message ?? x} (${JSON.stringify(x?.body ?? {})})`);
      }
    };
    const user = await grok.dapi.users.filter(`login = "${login}"`).include('group').first();
    expect(user != null, true, 'signed-up user not found');
    const group = user!.group.id;
    const sku = `SKU-VIS-${stamp}`;
    let visibleId: string | undefined;
    let hiddenId: string | undefined;
    try {
      visibleId = (await step('insert into apitests.item', () => items().insert({sku, name: 'Visibility probe'})))[0].id;
      hiddenId = (await step('insert into apitests.hidden_item',
        () => hidden().insert({sku, name: 'Visibility probe'})))[0].id;
      await step('grant View on apitests.item', () => items().grant(group, 'View'));
      await step('grant View on apitests.hidden_item', () => hidden().grant(group, 'View'));
      // Same grant, same row-mode: only the table whose rows default to table
      // visibility lets it through.
      const seenVisible = await step('read apitests.item as the restricted user',
        () => asUser(() => items().query({filter: `sku = "${sku}"`})));
      const seenHidden = await step('read apitests.hidden_item as the restricted user',
        () => asUser(() => hidden().query({filter: `sku = "${sku}"`})));
      expect(seenVisible.length, 1, 'a table View grant must reach rows of a table-visibility table');
      expect(seenHidden.length, 0,
        `a table View grant must NOT reach rows of a defaultRowVisibility:"none" table: ${JSON.stringify(seenHidden)}`);
    } finally {
      restoreAdmin();
      // Best-effort, one report per failure: a throw here would mask the result
      // AND skip the rest of the cleanup.
      const cleanup: [string, () => Promise<any>][] = [
        ['revoke View on apitests.item', () => items().revoke(group, 'View')],
        ['revoke View on apitests.hidden_item', () => hidden().revoke(group, 'View')],
        ['delete the apitests.item probe row', async () => visibleId && await items().delete(visibleId)],
        ['delete the apitests.hidden_item probe row', async () => hiddenId && await hidden().delete(hiddenId)],
      ];
      for (const [label, action] of cleanup)
        try {
          await action();
        } catch (x) {
          console.error(`cleanup — ${label} failed: ${x}`);
        }
      const blocked = await fetch(`${grok.dapi.root}/users/block`, {
        method: 'POST', credentials: 'include', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({'#type': 'User', 'id': user!.id, 'login': login}),
      });
      expect(blocked.ok, true, `test user ${login} not blocked: ${blocked.status}`);
      grok.dapi.domains.invalidateUiCaches();
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
