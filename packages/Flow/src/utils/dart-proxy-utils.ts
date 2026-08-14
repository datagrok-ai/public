/** Utilities for safely accessing Dart proxy objects (tags, options) from JS side */
import * as DG from 'datagrok-api/dg';
import {propertyNameToFriendly} from './naming';

export function safeGet(obj: any, key: string): any {
  try {
    if (!obj) return undefined;
    // Direct key access works for MapProxy.
    const val = obj[key];
    if (val !== undefined) return val;
    const entries = Object.entries(obj);
    for (const [k, v] of entries)
      if (k === key) return v;

    return undefined;
  } catch {
    return undefined;
  }
}

export function getRole(func: DG.Func): string | null {
  try {
    const opts = func.options;
    if (!opts) return null;
    const role = safeGet(opts, 'role');
    return role ? String(role) : null;
  } catch {
    return null;
  }
}

export function getTags(func: DG.Func): string[] {
  try {
    const tags: unknown = func.tags;
    if (!tags) return [];
    if (Array.isArray(tags)) return tags.map(String);
    if (typeof tags === 'string') return tags.split(',').map((t: string) => t.trim()).filter(Boolean);
    return Object.values(tags).map(String);
  } catch {
    return [];
  }
}

export function getPackageName(func: DG.Func): string {
  try {
    if (func.package && func.package.name)
      return (func.package as DG.Package).name! as string;
    return '';
  } catch {
    return '';
  }
}

export function getFuncQualifiedName(func: DG.Func): string {
  const pkg = getPackageName(func);
  const name = func.name || '';
  return pkg ? `${pkg}:${name}` : name;
}

/** Reads, in order: `Property.isOptional`, `nullable`, and the options map's `optional` flag. */
export function isInputOptional(prop: DG.Property): boolean {
  try {
    if (prop.isOptional) return true;
  } catch {/* older platform without the getter — fall through */}
  if (prop.nullable)
    return true;
  try {
    const opt = safeGet((prop as unknown as {options?: unknown}).options, 'optional');
    return opt === true || opt === 'true';
  } catch {
    return false;
  }
}

/** Parameter description from the Dart-proxy options map, falling back to `caption`. */
export function getParamDescription(prop: DG.Property): string {
  try {
    const opts = (prop as unknown as {options?: unknown}).options;
    const desc = safeGet(opts, 'description') ?? safeGet(opts, 'caption') ??
      (prop as unknown as {description?: unknown}).description;
    return desc ? String(desc).trim() : '';
  } catch {
    return '';
  }
}

/** Strip one pair of wrapping quotes from a default that arrives double-encoded from the annotation. */
export function unquoteDefault(s: string): string {
  const t = s.trim();
  if (t.length >= 2 && ((t.startsWith('\'') && t.endsWith('\'')) || (t.startsWith('"') && t.endsWith('"'))))
    return t.slice(1, -1);
  return t;
}

/** Declared default (`defaultValue ?? initialValue`), read defensively; strings unquoted. */
export function getParamDefault(prop: DG.Property): unknown {
  let v: unknown;
  try {v = prop.defaultValue;} catch {/* proxy read failed */}
  if (v === undefined || v === null)
    try {v = (prop as unknown as {initialValue?: unknown}).initialValue;} catch {/* proxy read failed */}
  if (v === undefined || v === null) return undefined;
  return typeof v === 'string' ? unquoteDefault(v) : v;
}

/** Display label: the declared caption, else a `friendlyName` differing from the raw name
 *  (an equal value is just the Dart fallback), else the humanized name. Identity stays `prop.name`. */
export function getParamDisplayName(prop: DG.Property): string {
  try {
    const cap = safeGet((prop as unknown as {options?: unknown}).options, 'caption');
    if (typeof cap === 'string' && cap.trim() !== '') return cap.trim();
    const fn = (prop as unknown as {friendlyName?: unknown}).friendlyName;
    if (typeof fn === 'string' && fn.trim() !== '' && fn.trim() !== prop.name) return fn.trim();
  } catch {/* fall through to name */}
  return propertyNameToFriendly(prop.name);
}

/** Node-header display name: the friendlyName's last `|` segment; humanized when it equals the
 *  raw name — a `//friendlyName:` annotation does not reliably survive package publishing. */
export function getFuncDisplayName(func: DG.Func): string {
  const raw = func.friendlyName || func.name || '';
  const parts = raw.split('|');
  const last = parts[parts.length - 1].trim();
  const declared = raw.trim() !== (func.name ?? '').trim() || parts.length > 1;
  return declared ? last : propertyNameToFriendly(last);
}
