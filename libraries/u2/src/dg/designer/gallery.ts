/* The sample gallery's save medium (F12) — localStorage, explicitly the P3.5 stopgap: P5's
   dg-form asset (WO-18) imports these entries once and deletes this path. Corrupt JSON, a full
   quota and a missing storage all answer empty — the gallery is a convenience, never a store
   of record. */

export const GALLERY_KEY = 'u2.designer.gallery';

function storage(): Storage | null {
  try {
    return globalThis.localStorage ?? null;
  } catch {
    return null;
  }
}

/** Saved specs by name. */
export function loadGallery(): Record<string, object> {
  try {
    const raw = storage()?.getItem(GALLERY_KEY);
    if (!raw)
      return {};
    const parsed = JSON.parse(raw) as unknown;
    return parsed !== null && typeof parsed === 'object' && !Array.isArray(parsed) ?
      parsed as Record<string, object> : {};
  } catch {
    return {};
  }
}

/** An existing `name` is overwritten — the save dialog says so. False when there is no storage or
 * it refused the write (quota, privacy mode), so the caller never reports a save that never was. */
export function saveToGallery(name: string, spec: object): boolean {
  const store = storage();
  if (!store)
    return false;
  const entries = loadGallery();
  entries[name] = spec;
  try {
    store.setItem(GALLERY_KEY, JSON.stringify(entries));
    return true;
  } catch {
    return false;
  }
}

export function listGallery(): string[] {
  return Object.keys(loadGallery());
}
