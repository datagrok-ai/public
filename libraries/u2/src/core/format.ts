import {backends} from '../sources/backends.js';
import {Dates} from './dates.js';

export interface FormatOptions {
  /** A platform format — a date pattern ('yyyy-MM-dd', 'MMM d, yyyy') or a number format ('#,##0.00') — as a
   * column or a property carries it. */
  format?: string;
  /** A Date on midnight UTC is a date, not a moment, and reads in UTC (`Dates.isDateOnly`): a domain `datetime`
   * that carries a date. Without it every Date reads in the local zone, as the platform's grid shows it. */
  utcDates?: boolean;
}

interface Parts {
  year: number;
  month: number;
  day: number;
  weekday: number;
  hours: number;
  minutes: number;
  seconds: number;
  ms: number;
}

const MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
  'November', 'December'];
const DAYS = ['Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday'];

const pad = (n: number, width = 2) => String(n).padStart(width, '0');
const hour12 = (p: Parts) => p.hours % 12 || 12;

/** The platform's date tokens (`date_time_parser.dart` `dateTimeFormatters`): a run of one letter is a token, any
 * other run prints as it is. */
const TOKENS: Record<string, (p: Parts) => string> = {
  d: (p) => String(p.day),
  dd: (p) => pad(p.day),
  ddd: (p) => DAYS[p.weekday].slice(0, 3),
  dddd: (p) => DAYS[p.weekday],
  h: (p) => String(hour12(p)),
  hh: (p) => pad(hour12(p)),
  H: (p) => String(p.hours),
  HH: (p) => pad(p.hours),
  m: (p) => String(p.minutes),
  mm: (p) => pad(p.minutes),
  M: (p) => String(p.month),
  MM: (p) => pad(p.month),
  MMM: (p) => MONTHS[p.month - 1].slice(0, 3),
  MMMM: (p) => MONTHS[p.month - 1],
  s: (p) => String(p.seconds),
  ss: (p) => pad(p.seconds),
  f: (p) => String(Math.floor(p.ms / 100)),
  ff: (p) => pad(Math.floor(p.ms / 10)),
  fff: (p) => pad(p.ms, 3),
  tt: (p) => p.hours >= 12 ? 'PM' : 'AM',
  yy: (p) => pad(p.year % 100),
  yyyy: (p) => String(p.year),
};

/** The one formatter: the text u2 shows for a value — a table cell, a group header, a filter's value, a stat — so a
 * value reads one way wherever it shows. The platform's own formatting (the user's date format, number formats) comes
 * through `backends.format`, which u2/dg fills; headless, the platform's defaults stand in. */
export class Format {
  /** The platform's date format while the user has not chosen one (`dtu.defaultDateFormat`). */
  static readonly DATE = 'M/d/yyyy';

  /** A value as u2 shows it: empty for null, undefined and NaN; a boolean `true` or `false`; a number as it is, or
   * in `format`; a Date (or anything with `toDate()`, a dayjs) in the platform's date format — the user's, else
   * {@link DATE} — with the time of day where it carries one, as the platform's grid shows a datetime column; a
   * list item by item; a record of no known shape empty, never "[object Object]". */
  static value(v: unknown, options: FormatOptions = {}): string {
    if (v === null || v === undefined || (typeof v === 'number' && isNaN(v)))
      return '';
    if (typeof v === 'string')
      return v;
    if (typeof v === 'number')
      return options.format ? backends.format?.number(v, options.format) ?? String(v) : String(v);
    if (v instanceof Date)
      return Format.date(v, options);
    if (Array.isArray(v))
      return v.map((x) => Format.value(x, options)).filter((t) => t !== '').join(', ');
    if (typeof v === 'object') {
      if (typeof (v as {toDate?: unknown}).toDate === 'function')
        return Format.date((v as {toDate(): Date}).toDate(), options);
      const proto = Object.getPrototypeOf(v);
      if (proto === Object.prototype || proto === null)
        return '';
    }
    return String(v);
  }

  /** `date` in `format`, a pattern of the platform's tokens (d dd ddd dddd, M MM MMM MMMM, yy yyyy, h hh H HH, m mm,
   * s ss, f ff fff, tt); without one, the platform's date format and the time of day the value carries. An invalid
   * date is empty. */
  static date(date: Date, options: FormatOptions = {}): string {
    if (isNaN(date.getTime()))
      return '';
    const p = Format._parts(date, options.utcDates === true && Dates.isDateOnly(date));
    const pattern = options.format || Format._pattern(p);
    return pattern.replace(/(.)\1*/g, (run) => TOKENS[run]?.(p) ?? run);
  }

  /** A count or a quantity as people read it, its thousands grouped by the locale: 1,234,567. */
  static count(n: number): string {
    return n.toLocaleString();
  }

  /** A date spelled out to read beside an item (`timestamp`): the locale's short month, the day and the year —
   * Oct 2, 2026 — and with `time` the hour and the minute; `utc` reads it in UTC. A cell's text is {@link value}. */
  static spelled(date: Date, options: {time?: boolean, utc?: boolean} = {}): string {
    return date.toLocaleString(undefined, {month: 'short', day: 'numeric', year: 'numeric',
      ...(options.time ? {hour: '2-digit', minute: '2-digit'} : {}), ...(options.utc ? {timeZone: 'UTC'} : {})});
  }

  // the platform's `initColumnFormat`, decided per value rather than per column
  private static _pattern(p: Parts): string {
    const day = backends.format?.dateFormat() || Format.DATE;
    return p.ms !== 0 ? `${day} HH:mm:ss.fff` : p.seconds !== 0 ? `${day} HH:mm:ss` :
      p.hours !== 0 || p.minutes !== 0 ? `${day} HH:mm` : day;
  }

  private static _parts(d: Date, utc: boolean): Parts {
    return utc ?
      {year: d.getUTCFullYear(), month: d.getUTCMonth() + 1, day: d.getUTCDate(), weekday: d.getUTCDay(),
        hours: d.getUTCHours(), minutes: d.getUTCMinutes(), seconds: d.getUTCSeconds(), ms: d.getUTCMilliseconds()} :
      {year: d.getFullYear(), month: d.getMonth() + 1, day: d.getDate(), weekday: d.getDay(),
        hours: d.getHours(), minutes: d.getMinutes(), seconds: d.getSeconds(), ms: d.getMilliseconds()};
  }
}
