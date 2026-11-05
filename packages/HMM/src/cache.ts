/** A small IndexedDB store for downloaded model data (InterPro Pfam families),
 * so families are fetched once per browser. Every failure (private mode,
 * blocked storage, quota) degrades to "not cached" instead of an error. */

const DB_NAME = 'datagrok-hmm';
const STORE = 'pfam';
const VERSION = 1;

let opening: Promise<IDBDatabase | null> | null = null;

function open(): Promise<IDBDatabase | null> {
  return opening ??= new Promise<IDBDatabase | null>((resolve) => {
    try {
      const request = indexedDB.open(DB_NAME, VERSION);
      request.onupgradeneeded = () => request.result.createObjectStore(STORE);
      request.onsuccess = () => resolve(request.result);
      request.onerror = () => resolve(null);
      request.onblocked = () => resolve(null);
    } catch {
      resolve(null);
    }
  });
}

function run<T>(mode: IDBTransactionMode, action: (store: IDBObjectStore) => IDBRequest): Promise<T | null> {
  return open().then((db) => new Promise<T | null>((resolve) => {
    if (!db) return resolve(null);
    try {
      const request = action(db.transaction(STORE, mode).objectStore(STORE));
      request.onsuccess = () => resolve((request.result ?? null) as T | null);
      request.onerror = () => resolve(null);
    } catch {
      resolve(null);
    }
  }));
}

export function cacheGet<T>(key: string): Promise<T | null> {
  return run<T>('readonly', (store) => store.get(key));
}

export async function cachePut(key: string, value: unknown): Promise<void> {
  await run('readwrite', (store) => store.put(value, key));
}

export async function cacheDelete(key: string): Promise<void> {
  await run('readwrite', (store) => store.delete(key));
}
