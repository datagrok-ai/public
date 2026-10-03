/** Times as read back from the domain store. A datetime column comes back as a dayjs object, whose
 * string form starts with the weekday (`Fri, 02 Oct 2026 …`), so comparing those strings sorts by
 * weekday. Compare milliseconds instead. */
import dayjs from 'dayjs';

/** Milliseconds since the epoch of a dayjs object, a Date or an ISO string; 0 when absent or unreadable. */
export function timeOf(d: unknown): number {
  if (d === null || d === undefined || d === '')
    return 0;
  const t = dayjs(d as any).valueOf();
  return Number.isNaN(t) ? 0 : t;
}

/** A comparator for newest first by [when]. */
export const newestFirstBy = <T>(when: (x: T) => unknown) => (a: T, b: T): number => timeOf(when(b)) - timeOf(when(a));

/** A comparator for oldest first by [when]. */
export const oldestFirstBy = <T>(when: (x: T) => unknown) => (a: T, b: T): number => timeOf(when(a)) - timeOf(when(b));
