/// Argument parsing and table formatting shared by `grok s observe alerts|errors|logger|capture|timeline`.

const UNIT_MS: Record<string, number> = {m: 60000, h: 3600000, d: 86400000, w: 604800000};
const BLOCKS = '▁▂▃▄▅▆▇█';

const given = (value: any) => hasValue(value) ? `got '${value}'` : 'got no value';

/** `30m`, `2h`, `7d`, `1w` → milliseconds; `m` is minutes. A leading `-` is accepted, as `--since -2w` in migrate. */
export function parseDuration(value: any, flag: string): number {
  const m = /^-?(\d+)\s*([mhdw])$/.exec(String(value ?? '').trim());
  if (!m || Number(m[1]) <= 0)
    throw new Error(`${flag} expects a duration such as 30m, 2h, 7d or 1w (m is minutes), ${given(value)}`);
  return Number(m[1]) * UNIT_MS[m[2]];
}

/** A `--since` value as the server takes it: validated, without the optional leading `-`. */
export function sinceArg(value: any): string {
  parseDuration(value, '--since');
  return String(value).trim().replace(/^-/, '');
}

/** ISO (UTC when it has no offset), a date (UTC midnight), relative `-7d`, or `HH:MM` UTC today. */
export function parseTime(value: any, flag: string, now: Date = new Date()): Date {
  const s = String(value ?? '').trim();
  const rel = /^-(\d+)([mhdw])$/.exec(s);
  if (rel)
    return new Date(now.getTime() - Number(rel[1]) * UNIT_MS[rel[2]]);
  const hm = /^(\d{1,2}):(\d{2})$/.exec(s);
  if (hm && Number(hm[1]) < 24 && Number(hm[2]) < 60) {
    const d = new Date(now.getTime());
    d.setUTCHours(Number(hm[1]), Number(hm[2]), 0, 0);
    return d;
  }
  if (/^\d{4}-\d{2}-\d{2}/.test(s)) {
    const d = new Date(/^\d{4}-\d{2}-\d{2}$/.test(s) ? `${s}T00:00Z` : /(Z|[+-]\d{2}:?\d{2})$/i.test(s) ? s : `${s}Z`);
    if (!isNaN(d.getTime()))
      return d;
  }
  throw new Error(`${flag} expects an ISO time (2026-10-04T06:00), a relative time (-7d) or HH:MM, ${given(value)}`);
}

function toDate(value: any): Date | null {
  if (value === null || value === undefined || value === '') return null;
  const d = new Date(value);
  return isNaN(d.getTime()) ? null : d;
}

/** `HH:MMZ` (`HH:MM:SS.mmmZ` with [seconds]) when [value] is today, else prefixed with `MM-DD`, in UTC. */
export function fmtTime(value: any, now: Date = new Date(), seconds: boolean = false): string {
  const d = toDate(value);
  if (!d) return '';
  const t = `${d.toISOString().slice(11, seconds ? 23 : 16)}Z`;
  return d.toISOString().slice(0, 10) === now.toISOString().slice(0, 10) ? t : `${fmtDate(d)} ${t}`;
}

export function fmtDateTime(value: any): string {
  const d = toDate(value);
  return d ? `${d.toISOString().slice(0, 16).replace('T', ' ')}Z` : '';
}

export function fmtDate(value: any): string {
  const d = toDate(value);
  return d ? d.toISOString().slice(5, 10) : '';
}

/** A span as `2 d`, `5 h` or `30 min`. */
export function fmtSpan(ms: number): string {
  if (ms >= UNIT_MS.d) return `${Math.round(ms / UNIT_MS.d)} d`;
  if (ms >= UNIT_MS.h) return `${Math.round(ms / UNIT_MS.h)} h`;
  return `${Math.max(0, Math.round(ms / UNIT_MS.m))} min`;
}

/** `190` → `3 h 10 min`. */
export function fmtMinutes(minutes: any): string {
  if (minutes === null || minutes === undefined || isNaN(Number(minutes))) return '—';
  const total = Math.round(Number(minutes));
  const h = Math.floor(total / 60);
  return h ? `${h} h ${total % 60} min` : `${total} min`;
}

export function sparkline(values: any): string {
  const nums: number[] = Array.isArray(values) ? values.map((v) => Number(v) || 0) : [];
  const max = Math.max(0, ...nums);
  return nums.map((v) => BLOCKS[max ? Math.round((v / max) * (BLOCKS.length - 1)) : 0]).join('');
}

/** `<id>.<n>` → `…X7K2QM.3`: the last six characters of the id (its leading ones are time), suffix kept. */
export function shortRequestId(id: any): string {
  if (!id) return '';
  const s = String(id);
  const dot = s.indexOf('.');
  const action = dot < 0 ? s : s.slice(0, dot);
  return (action.length > 6 ? `…${action.slice(-6)}` : action) + (dot < 0 ? '' : s.slice(dot));
}

/** A stack signature as the examples print it: its first six hex characters. */
export function shortSig(sig: any): string {
  return sig ? String(sig).replace(/-/g, '').slice(0, 6) : '';
}

export function truncate(value: any, max: number): string {
  const s = value === null || value === undefined ? '' : String(value);
  return s.length > max ? `${s.slice(0, max - 1)}…` : s;
}

/** A flag that may be repeated or comma-separated, as a flat list. */
export function listArg(value: any): string[] {
  if (value === undefined || value === null || value === true || value === false) return [];
  return (Array.isArray(value) ? value : [value]).flatMap((v) => String(v).split(',')).map((s) => s.trim()).filter(Boolean);
}

/** `queries` and `files` are the names the worked examples use for the platform's debug flags. */
export function normalizeFlag(name: string): string {
  const s = name.toLowerCase();
  return ({queries: 'query', files: 'storage'} as Record<string, string>)[s] ?? s;
}

/** A string option minimist left empty (`--signature` with no value) counts as not given. */
export function hasValue(value: any): boolean {
  return value !== undefined && value !== null && value !== true && value !== false && value !== '';
}

export function optString(value: any): string | undefined {
  return hasValue(value) ? String(value) : undefined;
}

export function rows<T>(list: T[], output: string, fn: (x: T) => Record<string, any>): any[] {
  return output === 'json' || output === 'quiet' ? list : list.map(fn);
}

/**
 * `a,b` replaces the list; `+a,-b` adds and removes against [current]. Mixing both forms is
 * refused: `a,+b` has no single meaning.
 */
export function applyListSpec(current: string[] | null | undefined, spec: any, normalize: (s: string) => string, flag: string): string[] {
  const items = listArg(spec);
  if (!items.length)
    throw new Error(`${flag} expects a comma list (a,b) or signed items (+a,-b)`);
  const signed = items.filter((i) => /^[+-]/.test(i));
  if (signed.length && signed.length !== items.length)
    throw new Error(`${flag}: use either a plain list (a,b) or signed items (+a,-b), not both`);
  if (!signed.length)
    return [...new Set(items.map(normalize))];
  const result = [...(current ?? [])];
  for (const item of items) {
    const name = normalize(item.slice(1));
    const at = result.indexOf(name);
    if (item[0] === '+' && at < 0) result.push(name);
    if (item[0] === '-' && at >= 0) result.splice(at, 1);
  }
  return result;
}

export function valueText(v: any): string {
  if (v === null || v === undefined) return '(none)';
  if (Array.isArray(v)) return v.length ? v.join(', ') : '(none)';
  if (typeof v === 'object') return JSON.stringify(v);
  return String(v);
}

export function printBlock(lines: [string, string][]): void {
  const width = Math.max(0, ...lines.map(([k]) => k.length)) + 2;
  for (const [k, v] of lines)
    console.log(`${k.padEnd(width)}${v}`.trimEnd());
}
