/* What only the Connections features need: the schema boxes the library's kinds do not reach, a
   connection's chats and the hidden providers of the tree. Connections on the server are the
   library's. */
import type {Page} from '@playwright/test';
import {Given, Then, When, kind} from '@datagrok-libraries/bdd';
import {chatIdsOf, deleteChatsOf, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* The schema view names each table box's host (`div-table-<table>`, db_views.dart), which is neither
   a u2 element nor a widget, so the library's generic `element` kind does not reach it. */
kind('schema table', {
  selector: '.d4-host',
  match: ['dart'],
  dartNames: ['div-table-{q}'],
  description: 'a table box of the schema view (Schemas > <schema> > Browse), by its table name',
});

type Conn = {id: string; name: string; friendlyName: string};

function connections(page: Page): Promise<Conn[]> {
  return page.evaluate(async () => (await grok.dapi.connections.list({pageSize: 5000}))
    .map((c: any) => ({id: c.id, name: c.name, friendlyName: c.friendlyName})));
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
