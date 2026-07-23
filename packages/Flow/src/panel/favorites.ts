/** Toolbox favorites — node types the user starred, persisted in localStorage.
 *  A favorite is `{type, label}` so the Favorites tab can render entries even
 *  when the underlying function is temporarily absent (package not loaded);
 *  `type` is the registered node type name (e.g. `DG Functions/…/Join Tables`,
 *  `Inputs/Table Input`, `Viewers/Scatter Plot`) — the same string every
 *  toolbox gesture (double-click, drag) already uses to create nodes. */

export interface FavoriteEntry {
  type: string;
  label: string;
}

const LS_KEY = 'funcflow.favorites.v1';
const listeners = new Set<() => void>();
let cache: FavoriteEntry[] | null = null;

function load(): FavoriteEntry[] {
  if (cache) return cache;
  try {
    const raw = localStorage.getItem(LS_KEY);
    const parsed = raw ? JSON.parse(raw) as unknown : [];
    cache = Array.isArray(parsed) ?
      parsed.filter((e): e is FavoriteEntry =>
        typeof (e as FavoriteEntry)?.type === 'string' && typeof (e as FavoriteEntry)?.label === 'string') :
      [];
  } catch {
    cache = []; // corrupt/blocked storage — start empty
  }
  return cache;
}

function save(items: FavoriteEntry[]): void {
  cache = items;
  try {
    localStorage.setItem(LS_KEY, JSON.stringify(items));
  } catch {/* storage blocked/full — favorites stay in-memory for the session */}
  for (const cb of [...listeners]) cb();
}

export function getFavorites(): FavoriteEntry[] {
  return [...load()];
}

export function isFavorite(type: string): boolean {
  return load().some((e) => e.type === type);
}

/** Stars/unstars a node type; returns the new state (true = now a favorite). */
export function toggleFavorite(entry: FavoriteEntry): boolean {
  const items = load();
  const idx = items.findIndex((e) => e.type === entry.type);
  if (idx >= 0) {
    save([...items.slice(0, idx), ...items.slice(idx + 1)]);
    return false;
  }
  save([...items, {type: entry.type, label: entry.label}]);
  return true;
}

/** Subscribes to favorites changes; returns the unsubscribe function. */
export function onFavoritesChanged(cb: () => void): () => void {
  listeners.add(cb);
  return () => listeners.delete(cb);
}

/** Test helper — wipes the store (and storage) and notifies subscribers. */
export function clearFavorites(): void {
  save([]);
}
