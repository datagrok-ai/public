/* One answer to "is this value a DATE rather than a moment?" and one set of UTC parts to render and
   edit it by (Astra B6). The schema has no `date` type, so a domain `datetime` that carries a date
   is stamped at UTC midnight; reading it in the local zone moves it a day back west of Greenwich.
   Three surfaces used to decide this three different ways — `timestamp()` on whole seconds,
   `DomainGrid` on microsecond divisibility, `DateField` on local parts. */

export interface UtcParts {
  year: number;
  /** 1-based, as an ISO day is spelled — `fromUtcParts` is the inverse. */
  month: number;
  day: number;
  hours: number;
  minutes: number;
}

export class Dates {
  /** Microseconds in a day — the unit a DG datetime column's raw data is in. */
  static readonly MICROSECONDS_PER_DAY = 86400000000;

  /** Whether the instant falls exactly on midnight UTC. */
  static isDateOnly(value: Date | number | string): boolean {
    const time = value instanceof Date ? value.getTime() : new Date(value).getTime();
    return Number.isFinite(time) && time % (Dates.MICROSECONDS_PER_DAY / 1000) === 0;
  }

  /** {@link isDateOnly} for a raw DG datetime cell (microseconds since the epoch); a null
   * sentinel and a non-finite value answer true — they say nothing either way. The sentinel is
   * `Filters.FLOAT_NULL`, which no count of microseconds can be: it is not an integer. */
  static isDateOnlyMicros(value: number): boolean {
    return !Number.isFinite(value) || !Number.isInteger(value) ||
      value % Dates.MICROSECONDS_PER_DAY === 0;
  }

  /** The UTC calendar parts of an instant — what a UTC date is FORMATTED and EDITED by. */
  static utcParts(date: Date): UtcParts {
    return {year: date.getUTCFullYear(), month: date.getUTCMonth() + 1, day: date.getUTCDate(),
      hours: date.getUTCHours(), minutes: date.getUTCMinutes()};
  }

  /** An instant from UTC calendar parts — the inverse of {@link utcParts}. */
  static fromUtcParts(parts: {year: number, month: number, day: number, hours?: number, minutes?: number}): Date {
    return new Date(Date.UTC(parts.year, parts.month - 1, parts.day, parts.hours ?? 0, parts.minutes ?? 0));
  }
}
