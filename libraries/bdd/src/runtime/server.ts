/* Endpoints the JS API does not wrap (chats, global permissions), called with the page's session —
   shared by the platform bindings and by a package's own bindings, so a suite never hand-rolls a
   fetch against the server. */
import {Page} from '@playwright/test';

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
