/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)
import {UUID, UUID_RE} from './registry';

const TYPED_RE = new RegExp(`^(.+):(${UUID})$`, 'i');

/**
 * Depth-first copy with every mapped UUID replaced, in values and in keys, bare or as
 * `<Type>:<uuid>`. UUIDs with no mapping are left in place and collected into `orphans`.
 */
export function rewrite<T>(json: T, idmap: Record<string, string>, orphans?: Set<string>): T {
  return walk(json, idmap, orphans);
}

function walk(v: any, idmap: Record<string, string>, orphans?: Set<string>): any {
  if (typeof v === 'string') return mapString(v, idmap, orphans);
  if (Array.isArray(v)) return v.map((x) => walk(x, idmap, orphans));
  if (v === null || typeof v !== 'object') return v;
  const out: any = {};
  for (const k of Object.keys(v))
    out[mapString(k, idmap, orphans)] = walk(v[k], idmap, orphans);
  return out;
}

function mapString(s: string, idmap: Record<string, string>, orphans?: Set<string>): string {
  if (UUID_RE.test(s)) return mapId(s, idmap, orphans);
  const typed = TYPED_RE.exec(s);
  return typed ? `${typed[1]}:${mapId(typed[2], idmap, orphans)}` : s;
}

function mapId(id: string, idmap: Record<string, string>, orphans?: Set<string>): string {
  if (idmap[id]) return idmap[id];
  orphans?.add(id);
  return id;
}

const isRef = (o: any): boolean => Object.keys(o).every((k) => k === 'id' || k === '#type');

/**
 * Ids of nested rows that own a primary key (`relations[]`, `entityTags[]`, query `params[]`,
 * view/layout child rows) — reusing them on the same server re-points the source's rows at
 * the copy, so a re-id or an adoption has to mint fresh ones. `{id}` / `{'#type', id}`
 * references are not rows and keep pointing at whatever the idmap says.
 */
export function nestedIds(json: any): string[] {
  const ids: string[] = [];
  const collect = (v: any, top: boolean) => {
    if (Array.isArray(v)) {
      for (const x of v) collect(x, false);
      return;
    }
    if (v === null || typeof v !== 'object') return;
    if (!top && typeof v.id === 'string' && UUID_RE.test(v.id) && !isRef(v))
      ids.push(v.id);
    for (const k of Object.keys(v))
      collect(v[k], false);
  };
  collect(json, true);
  return ids;
}
