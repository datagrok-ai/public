/* Endpoints the JS API does not wrap (chats, global permissions), called with the page's session —
   shared by the platform bindings and by a package's own bindings, so a suite never hand-rolls a
   fetch against the server. */
import {Page} from '@playwright/test';
import {atFeatureEnd} from './harness.js';
import {expect, pollMs} from './patience.js';

declare const grok: any;

export interface ServerApi {
  get<T>(path: string): Promise<T>;
  post(path: string, data: unknown): Promise<void>;
  remove(path: string): Promise<void>;
}

/** With `session`, the requests go as that session's account instead of the page's. */
export async function serverRequests(page: Page, session?: string): Promise<ServerApi> {
  const {root, token} = await page.evaluate((s) => ({root: new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''),
    token: s ?? String(grok.dapi.token)}), session ?? null);
  const headers = {Authorization: token};
  // a request that fails in transport throws with Playwright's call log, which prints the headers:
  // the session token must not reach the run's log, report or trace
  const scrubbed = async <T>(method: string, path: string, request: Promise<T>): Promise<T> => {
    try {
      return await request;
    }
    catch (error) {
      const message = String((error as Error)?.message ?? error).split(token).join('<token>')
        .replace(/(authorization["']?\s*[:=]\s*["']?)[^\s"',}]+/gi, '$1<token>');
      throw new Error(`${method} ${path}: ${message.split('\n')[0]}`);
    }
  };
  const checked = async (method: string, path: string, response: Promise<{ok(): boolean; text(): Promise<string>}>) => {
    const done = await scrubbed(method, path, response);
    const body = await done.text();
    if (!done.ok() || body.includes('ApiError'))
      throw new Error(`${method} ${path}: ${body.slice(0, 200)}`);
  };
  return {
    async get<T>(path: string): Promise<T> {
      const got = await scrubbed('GET', path, page.request.get(`${root}${path}`, {headers}));
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

/** The script the script view's Save made in a feature: its id for the cleanup, its name and grok name for
 * the claims that look for it elsewhere (a gallery card links to the grok name). */
export type SavedScript = {id: string; name: string; nqName: string};
const savedScripts = new WeakMap<Page, SavedScript>();

/** The page outlives the feature: the script is forgotten when the feature that saved it ends. */
export function rememberSavedScript(page: Page, script: SavedScript | null): void {
  if (script)
    savedScripts.set(page, script);
  else
    savedScripts.delete(page);
}

export function savedScriptOf(page: Page): SavedScript {
  const script = savedScripts.get(page);
  if (!script)
    throw new Error('no script has been saved in this feature yet: "user saves the script" comes first');
  return script;
}

/** The picture the Save dialog or the Layouts pane stores for an entity (`<pictureId>.png`): the server
 * keeps it when it deletes the entity, so the feature that made the entity deletes it, read back gone.
 * The thumbnails the server cuts from it (`<pictureId>_<width>.png`) have no delete of their own. A copy
 * saved with "Save a copy" shares its original's picture, which an earlier sweep may have taken, so only
 * a picture still stored is deleted (a storage backend may refuse the delete of a missing file). */
export async function deletePictures(page: Page, ids: string[]): Promise<void> {
  if (ids.length === 0)
    return;
  const api = await serverRequests(page);
  const root = await page.evaluate(() => new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''));
  // a picture the server does not hold is answered 200 with a JSON error, not a 404: stored is an image answer
  const stored = async (id: string) => {
    const got = await page.request.get(`${root}/entities/picture/${id}`);
    return got.ok() && (got.headers()['content-type'] ?? '').startsWith('image/');
  };
  for (const id of ids)
    if (await stored(id))
      await api.remove(`/entities/picture/${id}`);
  for (const id of ids)
    await expect.poll(() => stored(id), {message: `the picture ${id} still on the server`, timeout: pollMs(15000)}).toBe(false);
}

/** The picture id an entity's pictureUrl names ("…/entities/picture/<id>"); null for the default picture. */
export const pictureIdOf = (pictureUrl: unknown): string | null =>
  /\/entities\/picture\/([^/?]+)/.exec(String(pictureUrl ?? ''))?.[1] ?? null;

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

/** What a capability gate makes of the service health a stand reports (`grok.dapi.admin.getServiceInfos()`):
 * '' when the named service is enabled and Running, else why not. A stand that reports no health at all (a
 * dev stack whose datlas runs without `checkHealth`) says nothing about the service: a lenient gate lets the
 * test go on, a strict one (a tutorial, which refuses to start on such a stand) counts it as absent. */
export function serviceGap(services: {name: string; enabled: boolean; status: string}[], name: string, strict = false): string {
  if (services.length === 0)
    return strict ? 'the stand reports no service health' : '';
  const service = services.find((s) => s.name === name);
  if (service == null)
    return 'absent';
  return service.enabled && service.status === 'Running' ? '' : `${service.enabled ? '' : 'disabled, '}${service.status}`;
}

/** The service health the page's stand reports, for `serviceGap`. */
export const reportedServices = (page: Page): Promise<{name: string; enabled: boolean; status: string}[]> =>
  page.evaluate(async () => (await grok.dapi.admin.getServiceInfos())
    .map((s: any) => ({name: String(s.name), enabled: !!s.enabled, status: String(s.status)})));

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
 *  before it are a killed run's: they go now and at the end. With no `pattern`, only those go, now. */
export async function deleteLayoutsAtEnd(page: Page, pattern: string | null, stale?: string): Promise<void> {
  // the server stamps createdOn when it saves an entity (dinq repository_query.dart `save`), whatever the
  // client set, so its clock is the one to compare; the Date header comes from the proxy in front of it,
  // which may run on a clock of its own, hence the margin
  const now = await serverNow(page) - 5000;
  const sweep = async (recent: boolean): Promise<string[]> => {
    const {left, pictures} = await sweepPage(recent);
    await deletePictures(page, [...new Set(pictures.map(pictureIdOf).filter((p): p is string => p !== null))]);
    return left;
  };
  const sweepPage = (recent: boolean): Promise<{left: string[]; pictures: string[]}> => page.evaluate(async ([source, old, since, before]) => {
    const pictureOf = (window as any).grok_PictureMixin_Get_PictureUrl;
    if (typeof pictureOf !== 'function')
      throw new Error('grok_PictureMixin_Get_PictureUrl is gone from the client: a layout\'s picture cannot be found to delete');
    const me = String((await grok.dapi.users.current()).id);
    // the grok name drops what the friendly name keeps ("BDD-Q-layout-1" is "BDDQLayout1")
    const names = (x: any): string[] => [String(x.friendlyName ?? ''), String(x.name ?? '')].filter(Boolean);
    // valueOf keeps the milliseconds a toString() of the dayjs value drops
    const created = (x: any): number => x.createdOn ? Number(x.createdOn.valueOf()) : 0;
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
    const pictures: string[] = [];
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
      // the Layouts pane saves a picture of the view with the layout
      pictures.push(String(pictureOf(layout.dart) ?? ''));
      await grok.dapi.layouts.delete(layout);
    }
    return {left: (await listed()).map((l) => names(l)[0]), pictures};
  }, [pattern ?? '(?!)', stale ?? null, recent ? now : null, now - STALE_AFTER_MS] as const);
  if (stale)
    await sweep(false);
  if (pattern === null)
    return;
  atFeatureEnd(page, async () => {
    const left = await sweep(true);
    if (left.length)
      throw new Error(`layouts still on the server after the delete: ${left.join(', ')}`);
  });
}
