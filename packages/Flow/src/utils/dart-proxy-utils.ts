/** Utilities for safely accessing Dart proxy objects (tags, options) from JS side */
import * as DG from 'datagrok-api/dg';

export function safeGetEntries(obj: any): [string, any][] {
  try {
    if (!obj) return [];
    if (typeof obj !== 'object') return [];
    return Object.entries(obj);
  } catch {
    return [];
  }
}

export function safeGet(obj: any, key: string): any {
  try {
    if (!obj) return undefined;
    // Try direct key access first (works for MapProxy)
    const val = obj[key];
    if (val !== undefined) return val;
    // Fallback to iteration
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
    const tags = (func as any).tags ?? (!!func.dart ? (window as any).grok_Script_Get_Tags?.(func.dart) : null);
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

/** Whether a function input parameter is optional (declared `{optional: true}`).
 *  Read defensively from the Dart-proxy `options` map. */
export function isInputOptional(prop: DG.Property): boolean {
  try {
    const opt = safeGet((prop as unknown as {options?: unknown}).options, 'optional');
    return opt === true || opt === 'true';
  } catch {
    return false;
  }
}

/** The human description of a function parameter, read defensively from the
 *  Dart-proxy `options` map (`description` set via `@grok.decorators.param`),
 *  falling back to a `caption`. Empty string when none is declared. */
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

/** Returns the display name for a function node header.
 * Prefers friendlyName over name, then splits by '|' and takes the last segment. */
export function getFuncDisplayName(func: DG.Func): string {
  const raw = func.friendlyName || func.name || '';
  const parts = raw.split('|');
  return parts[parts.length - 1].trim();
}
