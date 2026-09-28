/// Argument parsing and table formatting shared by `grok s alerts|errors|logger|capture|timeline`.

/** Names the worked examples use for the platform's debug flags; the server validates the rest. */
export const DEBUG_FLAG_ALIASES: Record<string, string> = {queries: 'query', files: 'storage'};

const UNIT_MS: Record<string, number> = {m: 60000, h: 3600000, d: 86400000, w: 604800000};
const DAYS: Record<string, string> = {SUN: '0', MON: '1', TUE: '2', WED: '3', THU: '4', FRI: '5', SAT: '6', DAILY: '*', WEEKDAYS: '1-5'};
const BLOCKS = '▁▂▃▄▅▆▇█';

const pad2 = (n: number) => String(n).padStart(2, '0');

/** `30m`, `2h`, `7d`, `1w` → milliseconds; `m` is minutes. A leading `-` is accepted, as `--since -2w` in migrate. */
export function parseDuration(value: any, flag: string): number {
  const m = /^-?(\d+)\s*([mhdw])$/.exec(String(value ?? '').trim());
  if (!m || Number(m[1]) <= 0)
    throw new Error(`${flag} expects a duration such as 30m, 2h, 7d or 1w (m is minutes), got '${value}'`);
  return Number(m[1]) * UNIT_MS[m[2]];
}

/** A `--since` value as the server takes it: validated, without the optional leading `-`. */
export function sinceArg(value: any, flag: string = '--since'): string {
  parseDuration(value, flag);
  return String(value).trim().replace(/^-/, '');
}

/** ISO (local when it has no offset), a date (local midnight), relative `-7d`, or `HH:MM` today. */
export function parseTime(value: any, flag: string, now: Date = new Date()): Date {
  const s = String(value ?? '').trim();
  const rel = /^-(\d+)([mhdw])$/.exec(s);
  if (rel)
    return new Date(now.getTime() - Number(rel[1]) * UNIT_MS[rel[2]]);
  const hm = /^(\d{1,2}):(\d{2})$/.exec(s);
  if (hm && Number(hm[1]) < 24 && Number(hm[2]) < 60) {
    const d = new Date(now.getTime());
    d.setHours(Number(hm[1]), Number(hm[2]), 0, 0);
    return d;
  }
  if (/^\d{4}-\d{2}-\d{2}/.test(s)) {
    const d = new Date(/^\d{4}-\d{2}-\d{2}$/.test(s) ? `${s}T00:00` : s);
    if (!isNaN(d.getTime()))
      return d;
  }
  throw new Error(`${flag} expects an ISO time (2026-10-04T06:00), a relative time (-7d) or HH:MM, got '${value}'`);
}

function toDate(value: any): Date | null {
  if (value === null || value === undefined || value === '') return null;
  const d = new Date(value);
  return isNaN(d.getTime()) ? null : d;
}

const sameDay = (a: Date, b: Date) => a.getFullYear() === b.getFullYear() && a.getMonth() === b.getMonth() && a.getDate() === b.getDate();

/** `HH:MM` when [value] is today, else `MM-DD HH:MM`, local time. */
export function fmtTime(value: any, now: Date = new Date()): string {
  const d = toDate(value);
  if (!d) return '';
  const hm = `${pad2(d.getHours())}:${pad2(d.getMinutes())}`;
  return sameDay(d, now) ? hm : `${pad2(d.getMonth() + 1)}-${pad2(d.getDate())} ${hm}`;
}

/** `HH:MM:SS.mmm`, prefixed with `MM-DD` when not today: timeline rows of one click share a minute. */
export function fmtClock(value: any, now: Date = new Date()): string {
  const d = toDate(value);
  if (!d) return '';
  const t = `${pad2(d.getHours())}:${pad2(d.getMinutes())}:${pad2(d.getSeconds())}.${String(d.getMilliseconds()).padStart(3, '0')}`;
  return sameDay(d, now) ? t : `${pad2(d.getMonth() + 1)}-${pad2(d.getDate())} ${t}`;
}

export function fmtDateTime(value: any): string {
  const d = toDate(value);
  return d ? `${d.getFullYear()}-${pad2(d.getMonth() + 1)}-${pad2(d.getDate())} ${pad2(d.getHours())}:${pad2(d.getMinutes())}` : '';
}

export function fmtDate(value: any): string {
  const d = toDate(value);
  return d ? `${pad2(d.getMonth() + 1)}-${pad2(d.getDate())}` : '';
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

/** `<id>.<n>` → `01J9…7K`: the action part, first four and last two characters. */
export function shortRequestId(id: any): string {
  if (!id) return '';
  const s = String(id);
  const action = s.includes('.') ? s.slice(0, s.indexOf('.')) : s;
  return action.length > 8 ? `${action.slice(0, 4)}…${action.slice(-2)}` : action;
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

export function normalizeLevel(name: string): string {
  return name.toLowerCase();
}

export function normalizeFlag(name: string): string {
  return DEBUG_FLAG_ALIASES[name.toLowerCase()] ?? name.toLowerCase();
}

/** A string option minimist left empty (`--signature` with no value) counts as not given. */
export function hasValue(value: any): boolean {
  return value !== undefined && value !== null && value !== true && value !== false && value !== '';
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

/** `MON 07:00`, `DAILY 07:00`, `WEEKDAYS 07:00` → cron; a five-field cron string passes through. */
export function cronFromSchedule(schedule: any): string {
  const s = String(schedule ?? '').trim();
  const m = /^([A-Za-z]+)\s+(\d{1,2}):(\d{2})$/.exec(s);
  if (m) {
    const day = DAYS[m[1].toUpperCase()];
    const hour = Number(m[2]);
    const minute = Number(m[3]);
    if (day === undefined || hour > 23 || minute > 59)
      throw new Error(`--schedule '${s}': expected MON..SUN, DAILY or WEEKDAYS and HH:MM`);
    return `${minute} ${hour} * * ${day}`;
  }
  if (s.split(/\s+/).length === 5)
    return s;
  throw new Error(`--schedule '${s}': expected "MON 07:00", "DAILY 07:00", "WEEKDAYS 07:00" or a five-field cron`);
}

export function slug(name: string): string {
  return String(name).toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-+|-+$/g, '');
}

export function valueText(v: any): string {
  if (v === null || v === undefined) return '(none)';
  if (Array.isArray(v)) return v.length ? v.join(', ') : '(none)';
  if (typeof v === 'object') return JSON.stringify(v);
  return String(v);
}

/** A key/value block: labels padded to one column. */
export function printBlock(lines: [string, string][]): void {
  const width = Math.max(0, ...lines.map(([k]) => k.length)) + 2;
  for (const [k, v] of lines)
    console.log(`${k.padEnd(width)}${v}`.trimEnd());
}
