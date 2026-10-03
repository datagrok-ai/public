/** Display names of platform users, by id, for every view that says who did something. One lookup
 * per distinct id, all at once, kept for the session: the same people appear on every program. A
 * user that cannot be read (deleted, not visible to the current user) keeps the id, and is asked
 * again next time. */
import * as grok from 'datagrok-api/grok';

const cache = new Map<string, Promise<string>>();

function nameOf(id: string): Promise<string> {
  let name = cache.get(id);
  if (!name) {
    // a failed read is not kept: the next view tries again
    name = grok.dapi.users.find(id).then((u) => u?.friendlyName || u?.login || id).catch(() => {
      cache.delete(id);
      return id;
    });
    cache.set(id, name);
  }
  return name;
}

/** Display names of [ids] (nulls skipped), by id. Never fails: an unreadable user maps to its id. */
export async function userNames(ids: (string | null | undefined)[]): Promise<Map<string, string>> {
  const distinct = [...new Set(ids.filter((x): x is string => !!x))];
  const names = await Promise.all(distinct.map(nameOf));
  return new Map(distinct.map((id, i) => [id, names[i]]));
}

/** Forgets the names read so far: after a rename, or in tests. */
export function clearUserNames(): void {
  cache.clear();
}

/** How many ids have a name read or being read; for tests. */
export function cachedUserNames(): number {
  return cache.size;
}
