/* What only the Connections features need: the tree rows and schema boxes the core does not name
   yet, a connection's chats, a catalog's Database meta comment, a table dropped in a shared test
   database, the hidden providers of the tree, and the balloon a connection test answers with.
   Connections on the server, secrets, the local file and the reload are the library's. */
import type {Page} from '@playwright/test';
import {Given, Then, When, kind} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* The Schemas / Catalogs row under a connection carries the connection's own tree name; only its
   wrapper is distinct (`div-<Provider>-<ConnName>-Schemas`). Switch to
   `Databases---<Provider>---<Conn>---Schemas tree node` once the core names the row. */
kind('tree group', {
  selector: '.d4-tree-view-group[name^="div-"]',
  match: ['dart'],
  dartNames: ['div-{q}'],
  description: 'a Databases tree group the core names only by its wrapper: "Postgres-Chembl-Schemas" tree group',
});

/* The schema view draws one box per table with no name; switch to the core's name once it lands. */
kind('schema table', {
  selector: '.d4-sketch-item',
  match: ['label'],
  labelSelector: '.d4-sketch-item-header',
  description: 'a table box of the schema view (Schemas > <schema> > Browse), by its header',
});

type Conn = {id: string; name: string; friendlyName: string};

function connections(page: Page): Promise<Conn[]> {
  return page.evaluate(async () => (await grok.dapi.connections.list({pageSize: 5000}))
    .map((c: any) => ({id: c.id, name: c.name, friendlyName: c.friendlyName})));
}

/** The session's own endpoints for what the JS API does not wrap (a connection's chats). */
async function api(page: Page) {
  const {root, token} = await page.evaluate(() => ({root: new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''),
    token: String(grok.dapi.token)}));
  return {
    get: async <T>(path: string): Promise<T> => (await page.request.get(`${root}${path}`, {headers: {Authorization: token}})).json(),
    remove: async (path: string): Promise<void> => {
      const done = await page.request.delete(`${root}${path}`, {headers: {Authorization: token}});
      if (!done.ok())
        throw new Error(`DELETE ${path}: HTTP ${done.status()}`);
    },
  };
}

/* A chat about an entity is found by the entity (`/chats?entityId=`); a group's chat lives in a
   hidden group made for it (`/chats/with_groups`). Either goes before the entity does, or the chat
   outlives what it was about (libraries/bdd CLAUDE.md, groups). */
async function chatIdsOf(page: Page, id: string): Promise<string[]> {
  const a = await api(page);
  const about = await a.get<{id: string}[]>(`/chats?entityId=${id}`);
  const grouped = await a.get<{id: string}[]>(`/chats/with_groups?ids=${id}`);
  return [...new Set([...(Array.isArray(about) ? about : []), ...(Array.isArray(grouped) ? grouped : [])].map((c) => c.id))];
}

async function deleteChatsOf(page: Page, id: string): Promise<void> {
  const a = await api(page);
  for (const chat of await chatIdsOf(page, id))
    await a.remove(`/chats/${chat}`);
}

const chatsOf = async (page: Page, id: string): Promise<number> => (await chatIdsOf(page, id)).length;

export const connectionHasChat = Then('the {string} connection should have a chat on the server', async (page: Page, name: string) => {
  const c = (await connections(page)).find((x) => x.friendlyName === name);
  expect(c, `the "${name}" connection`).toBeTruthy();
  await expect.poll(() => chatsOf(page, c!.id), {message: `chats on "${name}"`, timeout: pollMs(15000)}).toBeGreaterThan(0);
}, {tier: 'api'});

export const connectionHasNoChat = Then('the {string} connection should have no chat on the server', async (page: Page, name: string) => {
  const c = (await connections(page)).find((x) => x.friendlyName === name);
  expect(c, `the "${name}" connection`).toBeTruthy();
  await expect.poll(() => chatsOf(page, c!.id), {message: `chats on "${name}"`, timeout: pollMs(15000)}).toBe(0);
}, {tier: 'api'});

export const deleteChatOfConnection = When('user deletes the chat of the {string} connection', async (page: Page, name: string) => {
  const c = (await connections(page)).find((x) => x.friendlyName === name);
  expect(c, `the "${name}" connection`).toBeTruthy();
  await deleteChatsOf(page, c!.id);
}, {tier: 'api', description: 'the chat goes first — a chat must never outlive the entity it is about'});

/* A catalog's Database meta comment, read and written through the connection's database info —
   what the pane's SAVE writes. The connection is named by its name or its friendly name. */
async function catalogComment(page: Page, catalog: string, connection: string, set?: string): Promise<string> {
  return page.evaluate(async ([cat, conn, value]) => {
    const c = (await grok.dapi.connections.list({pageSize: 5000})).find((x: any) => x.name === conn || x.nqName === conn);
    if (!c)
      throw new Error(`no connection named ${conn}`);
    const info = (await grok.dapi.connections.getDatabaseInfo(c, cat)).find((d: any) => d.name === cat);
    if (!info)
      throw new Error(`no catalog ${cat} in ${conn}`);
    if (value !== null)
      await info.setComment(value);
    return info.comment ?? '';
  }, [catalog, connection, set ?? null] as [string, string, string | null]);
}

export const catalogHasNoComment = Given('the {string} catalog of the {string} connection has no comment', async (page: Page, catalog: string, connection: string) => {
  const clear = async () => {
    await catalogComment(page, catalog, connection, '');
    await expect.poll(() => catalogComment(page, catalog, connection), {message: `the comment on ${connection}/${catalog}`}).toBe('');
  };
  atFeatureEnd(page, clear);
  await clear();
}, {tier: 'api', description: 'clears the Database meta comment now and again at feature end — a shared connection keeps no trace'});

export const catalogCommentIs = Then('the {string} catalog of the {string} connection should have the comment {string}', async (page: Page, catalog: string, connection: string, text: string) => {
  await expect.poll(() => catalogComment(page, catalog, connection), {message: `the comment on ${connection}/${catalog}`, timeout: pollMs(15000)}).toBe(text);
}, {tier: 'api', description: 'what the server holds for the catalog, not what the pane still shows'});

/* A connection test answers on a balloon, and a test that cannot log in answers only when the
   connector gives up — long after its task-bar entry has gone, so the balloon checks' 5 s is too
   short. The connector's own socket timeout is 180 s (grok_connect, `socketTimeout|180`), and a
   stand whose network drops the packets rather than refusing them takes all of it: ~40 s on dev,
   the full timeout on a local stand. The test's own answer is awaited up to 200 s. */
export const connectionTestEnded = Then('the connection test should have ended on an {word} balloon containing {string}',
  async (page: Page, type: string, text: string) => {
    let shown: string[] = [];
    await expect.poll(async () => {
      shown = shown.concat((await viewers.takeBalloons(page)).map((b: {type: string; message: string}) => `${b.type}: ${b.message}`));
      return shown.some((s) => s.startsWith(`${type}: `) && s.includes(text));
    }, {timeout: pollMs(200000), message: `an ${type} balloon containing "${text}"; balloons so far: ${shown.join(' | ') || 'none'}`}).toBe(true);
  }, {description: 'the balloon a connection test answers with (info "connected successfully" or error "failed to connect: …"), awaited up to 200 s, past the connector 180 s socket timeout'});

/* A table a scenario creates in a shared test database (the Dbtests connections) is dropped before
   and after the feature, so a run that failed half-way leaves nothing behind. The connection is named
   by its name. */
export const noTableOnConnection = Given('no table {string} is in the database of the {string} connection', async (page: Page, table: string, connection: string) => {
  if (!/^[a-z_][a-z0-9_]*$/.test(table))
    throw new Error(`"${table}" is not a plain table name`);
  const drop = async () => {
    const error = await page.evaluate(async ([t, conn]) => {
      const c = (await grok.dapi.connections.list({pageSize: 5000})).find((x: any) => x.name === conn || x.nqName === conn);
      if (!c)
        return `no connection named ${conn}`;
      try {
        await grok.data.db.query(c.nqName, `drop table if exists ${t}`);
        return '';
      }
      catch (e: any) {
        return String(e?.message ?? e);
      }
    }, [table, connection] as [string, string]);
    // a DDL statement returns no result set, which some connectors report as an error
    if (error && !/result|no results|resultset/i.test(error))
      throw new Error(`dropping ${table} on ${connection}: ${error}`);
  };
  atFeatureEnd(page, drop);
  await drop();
}, {tier: 'api', description: 'drop table if exists, now and at feature end, through the connection itself'});

/* The Databases tree hides the rarer providers behind "Show more" until its icon is clicked, and a
   page keeps them shown after that: an idempotent state for the scenarios after the one that claims
   the reveal itself. */
export const hiddenProvidersShown = Given('the hidden providers of the Databases tree are shown', async (page: Page) => {
  const more = page.locator('[name="tree-Databases---Show-more"]').filter({visible: true});
  if (await more.count() > 0)
    await more.first().locator('[name="icon-ellipsis-h"]').click();
  await expect(page.locator('[name="tree-Databases---Show-more"]').filter({visible: true}), 'the "Show more" row, gone once the providers show')
    .toHaveCount(0, {timeout: pollMs(15000)});
}, {tier: 'ui', description: 'clicks "Show more" when it is there; the row goes once the hidden providers show'});

/* What only the Connections features need: the tree rows and schema boxes the core does not name
   yet, a connection's chats, a catalog's Database meta comment, a table dropped in a shared test
   database, the hidden providers of the tree, and the balloon a connection test answers with.
   Connections on the server, secrets, the local file and the reload are the library's. */
