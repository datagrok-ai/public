/** New-data notification (F-11, FR-83): what landed in a program's dataset since a follower last
 * looked. Interim until the platform can push events (PLAT-03/06): each sync reads the program
 * dataset once per followed program and compares it with a snapshot kept on each follow. Pure: the
 * caller reads the rows and keeps the snapshots.
 *
 * A snapshot is cheap and conservative: the number of rows, the number of filled cells (so a
 * measurement landing in a blank cell of a wide SAR table counts), and, when the portfolio names a
 * watermark column, its largest value. Only growth is news; a row or value that disappears, or a
 * value overwritten by another, is not reported. */
import {columnValue, DataRow, isBlank} from './conditions';

export interface DataSnapshot {
  rows: number;
  /** Filled cells over every column. */
  values: number;
  /** The largest watermark value, when the portfolio declares a watermark column. */
  mark?: number | string | null;
}

export interface NewData {
  now: DataSnapshot;
  /** No snapshot before: this one is the baseline, nothing is news yet. */
  baseline: boolean;
  newRows: number;
  newValues: number;
}

const num = (v: unknown): number | null => typeof v === 'number' && !Number.isNaN(v) ? v :
  typeof v === 'string' && v.trim() !== '' && !Number.isNaN(Number(v)) ? Number(v) : null;

/** Whether watermark value [a] comes after [b]: numbers by value, anything else (ISO dates) as text. */
export function after(a: unknown, b: unknown): boolean {
  if (isBlank(a))
    return false;
  if (isBlank(b))
    return true;
  const x = num(a);
  const y = num(b);
  if (x !== null && y !== null)
    return x > y;
  const t = (v: unknown) => v instanceof Date ? v.toISOString() : typeof (v as any)?.toISOString === 'function' ?
    (v as any).toISOString() : String(v);
  return t(a) > t(b);
}

/** The rows a follow is about: every row, or one series' (by the `series` column). */
export const rowsOf = (rows: DataRow[], series?: string | null): DataRow[] =>
  series ? rows.filter((r) => String(columnValue(r, 'series') ?? '') === series) : rows;

export function snapshot(rows: DataRow[], watermark?: string | null): DataSnapshot {
  let values = 0;
  let mark: unknown = null;
  for (const r of rows) {
    for (const v of Object.values(r)) {
      if (!isBlank(v))
        values++;
    }
    if (watermark) {
      const m = columnValue(r, watermark);
      if (after(m, mark))
        mark = m;
    }
  }
  const out: DataSnapshot = {rows: rows.length, values};
  if (watermark)
    out.mark = isBlank(mark) ? null : typeof mark === 'number' ? mark : String((mark as any)?.toISOString?.() ?? mark);
  return out;
}

/** What is new in [rows] since [before]. With a watermark and a mark before, new rows are those
 * past it, so a row replaced by a newer one counts although the count is the same. */
export function compare(before: DataSnapshot | null | undefined, rows: DataRow[], watermark?: string | null): NewData {
  const now = snapshot(rows, watermark);
  if (!before || typeof before.rows !== 'number')
    return {now, baseline: true, newRows: 0, newValues: 0};
  const newRows = watermark && !isBlank(before.mark) ?
    rows.filter((r) => after(columnValue(r, watermark), before.mark)).length :
    Math.max(0, now.rows - before.rows);
  return {now, baseline: false, newRows, newValues: Math.max(0, now.values - (before.values ?? 0))};
}

export const isNews = (n: NewData): boolean => !n.baseline && (n.newRows > 0 || n.newValues > 0);

const plural = (n: number, one: string) => `${n} ${n === 1 ? one : `${one}s`}`;

/** One line for a notification or the program page; null when nothing is new. */
export function describeNews(label: string, n: NewData): string | null {
  if (!isNews(n))
    return null;
  const parts = [...(n.newRows ? [plural(n.newRows, 'new row')] : []), ...(n.newValues ? [plural(n.newValues, 'new value')] : [])];
  return `${label}: ${parts.join(', ')}`;
}

/** What a follow is called in a message: the program, or the program and its series. */
export const followLabel = (program: string, unit?: {kind?: string | null, code?: string | null}): string =>
  unit?.code ? `${program} ${String(unit.kind ?? 'series').replace('_', ' ')} ${unit.code}` : program;
