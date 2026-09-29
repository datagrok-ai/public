/* Endpoints the JS API does not wrap (chats, global permissions), called with the page's session —
   shared by the platform bindings and by a package's own bindings, so a suite never hand-rolls a
   fetch against the server. */
import {Page} from '@playwright/test';
import {atFeatureEnd} from './harness.js';

declare const grok: any;

export interface ServerApi {
  get<T>(path: string): Promise<T>;
  post(path: string, data: unknown): Promise<void>;
  remove(path: string): Promise<void>;
}

export async function serverRequests(page: Page): Promise<ServerApi> {
  const {root, token} = await page.evaluate(() => ({root: new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''),
    token: String(grok.dapi.token)}));
  const headers = {Authorization: token};
  const checked = async (method: string, path: string, response: Promise<{ok(): boolean; text(): Promise<string>}>) => {
    const done = await response;
    const body = await done.text();
    if (!done.ok() || body.includes('ApiError'))
      throw new Error(`${method} ${path}: ${body.slice(0, 200)}`);
  };
  return {
    async get<T>(path: string): Promise<T> {
      const got = await page.request.get(`${root}${path}`, {headers});
      if (!got.ok())
        throw new Error(`GET ${path} failed: HTTP ${got.status()}`);
      return got.json();
    },
    post: (path: string, data: unknown) => checked('POST', path, page.request.post(`${root}${path}`, {headers, data})),
    remove: (path: string) => checked('DELETE', path, page.request.delete(`${root}${path}`, {headers})),
  };
}

/** The chats about an entity: its own (`/chats?entityId=`) and, for a group, the one kept in the
 * hidden group made for it (`/chats/with_groups`). */
export async function chatIdsOf(page: Page, id: string): Promise<string[]> {
  const api = await serverRequests(page);
  const listed = await Promise.all([api.get<{id: string}[]>(`/chats/with_groups?ids=${id}`),
    api.get<{id: string}[]>(`/chats?entityId=${id}`)]);
  return [...new Set(listed.flat().filter(Boolean).map((c) => c.id))];
}

/* Deleting the entity first leaves a chat that throws in every profile's chat listing (forum.dart),
   so the chat goes first, whichever way it is held. */
export async function deleteChatsOf(page: Page, id: string): Promise<void> {
  const api = await serverRequests(page);
  for (const chat of await chatIdsOf(page, id))
    await api.remove(`/chats/${chat}`);
}

/** The server's clock, from the Date header of an API answer (whole seconds): what the `createdOn` of an
 * entity the server stamps counts in. The JS API has no getter for it. */
export async function serverNow(page: Page): Promise<number> {
  const root = await page.evaluate(() => new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''));
  const response = await page.request.get(`${root}/info/server`);
  if (!response.ok())
    throw new Error(`server clock request failed: HTTP ${response.status()}`);
  const time = Date.parse(response.headers()['date'] ?? '');
  if (!Number.isFinite(time))
    throw new Error('server clock response has no valid Date header');
  return time;
}

/* A fixture name ends in its run's {run} or {time}. A run that was killed never reached its
   feature-end cleanup, so the fixtures of the same family that are older than any live feature go too. */
export const RUN_SUFFIX = /-(\d{13,}|[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})$/;
export const STALE_AFTER_MS = 60 * 60 * 1000;

export const fixtureFamilies = (names: string[]): string[] =>
  names.filter((name) => RUN_SUFFIX.test(name)).map((name) => name.replace(RUN_SUFFIX, ''));

export function isStaleFixture(entity: {name: string; friendlyName: string; createdOn: number}, families: string[],
  now = Date.now()): boolean {
  if (entity.createdOn === 0 || now - entity.createdOn < STALE_AFTER_MS)
    return false;
  return [entity.friendlyName, entity.name]
    .some((name) => RUN_SUFFIX.test(name) && families.includes(name.replace(RUN_SUFFIX, '')));
}

/** A view that saves a layout leaves the layout and the project that wraps it on the server; the
 *  entity the layout belongs to does not take them with it. The signed-in account's layouts whose name
 *  matches `pattern` and that were made since this call go at feature end, each with the wrapper of its
 *  name, and the listing is read again to see them gone. Those matching `stale` and made over an hour
 *  before it are a killed run's: they go now and at the end. */
export async function deleteLayoutsAtEnd(page: Page, pattern: string, stale?: string): Promise<void> {
  // the client stamps a layout's createdOn (layout.dart), so the browser's clock is the one to compare
  const now: number = await page.evaluate(() => Date.now());
  const sweep = (recent: boolean): Promise<string[]> => page.evaluate(async ([source, old, since, before]) => {
    const me = String((await grok.dapi.users.current()).id);
    // the grok name drops what the friendly name keeps ("BDD-Q-layout-1" is "BDDQLayout1")
    const names = (x: any): string[] => [String(x.friendlyName ?? ''), String(x.name ?? '')].filter(Boolean);
    const created = (x: any): number => x.createdOn ? new Date(x.createdOn.toString()).getTime() : 0;
    const named = (x: any, re: string): boolean => names(x).some((n) => new RegExp(re, 'i').test(n));
    const ours = (x: any): boolean => String(x.author?.id ?? '') === me &&
      (since !== null && created(x) >= since && named(x, source) ||
        old !== null && created(x) > 0 && created(x) < before && named(x, old));
    const listed = async (): Promise<any[]> => {
      const all: any[] = [];
      for (let pageNumber = 1; ; pageNumber++) {
        const layouts = await grok.dapi.layouts.list({pageSize: 1000, pageNumber});
        all.push(...layouts);
        if (layouts.length < 1000)
          return all.filter(ours);
      }
    };
    for (const layout of await listed()) {
      // the wrapper is found by name, so it is taken only when it is ours and of the same run too:
      // "Df" is a name another account may hold on a shared stand
      for (const n of names(layout)) {
        const filter = `name = ${JSON.stringify(n)} or friendlyName = ${JSON.stringify(n)}`;
        const project = (await grok.dapi.projects.filter(filter).list())
          .find((p: any) => ours(p) && names(p).some((pn) => names(layout).includes(pn)));
        if (project)
          await grok.dapi.projects.delete(project);
      }
      await grok.dapi.layouts.delete(layout);
    }
    return (await listed()).map((l) => names(l)[0]);
  }, [pattern, stale ?? null, recent ? now : null, now - STALE_AFTER_MS] as const);
  if (stale)
    await sweep(false);
  atFeatureEnd(page, async () => {
    const left = await sweep(true);
    if (left.length)
      throw new Error(`layouts still on the server after the delete: ${left.join(', ')}`);
  });
}
