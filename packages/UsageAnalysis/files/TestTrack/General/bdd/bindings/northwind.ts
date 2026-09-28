/* What only the dev-only NorthwindTest features need: the md's "second user who can access the
   NorthwindTest connection". The connection (Dbtests:PostgresTest) is not shared with the second
   account on dev, and a project's Share dialog does not share the connection the project reads, so
   the feature grants the second account "View and use" on the connection for its run and takes the
   grant back at the end, reading the server back both times. */
import type {Page} from '@playwright/test';
import {Given} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect} from '@datagrok-libraries/bdd/runtime';
import {sharingLogin} from '@datagrok-libraries/bdd/bindings/platform/steps';

declare const grok: any;

/* The connection's id and the id of the sharing user's personal group. A package connection is not
   found by `name = …` on dev; `shortName` finds it. */
function ids(page: Page, nqName: string): Promise<{connection: string; group: string}> {
  return page.evaluate(async ([nq, login]) => {
    const name = nq.split(':').pop();
    const listed = await grok.dapi.connections.filter(`shortName = "${name}"`).list();
    const connection = listed.find((e: any) => e.nqName.toLowerCase() === nq.toLowerCase());
    if (connection == null)
      throw new Error(`no connection ${nq} on the server; with that short name: ${listed.map((e: any) => e.nqName).join(', ') || 'none'}`);
    const user = await grok.dapi.users.filter(`login = "${login}"`).first();
    if (user == null)
      throw new Error(`no user with login ${login} on the server`);
    return {connection: connection.id as string, group: user.group.id as string};
  }, [nqName, sharingLogin()] as const);
}

/* The privileges the group holds on the entity, as the Share dialog lists them (the route walks the
   entity's project links; `grok.dapi.permissions.get` answers nothing for a connection on dev). Read
   in the page: the library's `serverRequests` is a Playwright request, whose failure log prints the
   session's Authorization header. */
function privilegesOf(page: Page, entity: string, group: string): Promise<string[]> {
  return page.evaluate(async ([e, g]) => {
    const response = await fetch(`${grok.dapi.root}/privileges/permissions?entityId=${e}&all=true`,
      {headers: {Authorization: String(grok.dapi.token)}});
    if (!response.ok)
      throw new Error(`the privileges of ${e} were not read: HTTP ${response.status}`);
    const rows: any[] = await response.json();
    return rows.filter((r) => r?.userGroup?.id === g).map((r) => String(r?.permission?.name ?? '?')).sort();
  }, [entity, group] as const);
}

function changeGrant(page: Page, connection: string, group: string, grant: boolean): Promise<void> {
  return page.evaluate(async ([c, g, on]) => {
    const entity = await grok.dapi.connections.find(c);
    const userGroup = await grok.dapi.groups.find(g);
    if (on)
      await grok.dapi.permissions.grant(entity, userGroup, false);
    else
      // the (group, entity) order is the one every server version's JS API reads
      await grok.dapi.permissions.revoke(userGroup, entity);
  }, [connection, group, grant] as const);
}

export const sharingUserUsesConnection = Given('the sharing user may use the {string} connection until the feature ends',
  async (page: Page, nqName: string) => {
    const {connection, group} = await ids(page, nqName);
    atFeatureEnd(page, async () => {
      if ((await privilegesOf(page, connection, group)).length > 0)
        await changeGrant(page, connection, group, false);
      await expect.poll(() => privilegesOf(page, connection, group),
        {message: `the privileges the sharing user still holds on ${nqName}`}).toEqual([]);
    });
    if (!(await privilegesOf(page, connection, group)).includes('DataConnection.Query'))
      await changeGrant(page, connection, group, true);
    await expect.poll(() => privilegesOf(page, connection, group),
      {message: `the privileges the sharing user holds on ${nqName}`}).toContain('DataConnection.Query');
  }, {tier: 'api', description: 'grants the sharing user "View and use" on the connection (by nqName) and reads the Query privilege back; at feature end revokes it and reads back that none is left'});
