/* One value just measured, checked: Monitor's rule, shared by Monitor (monitor.ts), the day's entry (day-entry.ts)
   and the study table's corrections (study-table.ts). A module of its own, importing nothing of the app: study-table.ts
   reading it from monitor.ts made an import cycle (monitor.ts imports study-table.ts). */
import {minus} from './format';

/** A volume this many times (or this fraction of) the animal's last trusted one, between two sessions: a slipped
 *  decimal or a missed key, not growth (a tumor at its fastest doubles in about a day, eightfold in three). */
export const TYPO_FOLD = 7;

/** How `v` stands against `before` when it is a likely typo ("10 times", "1/10 of"), else null. */
export function likelyTypo(before: number, v: number): string | null {
  if (!(before > 0 && v > 0))
    return null;
  const r = v / before;
  return r >= TYPO_FOLD ? `${Math.round(r)} times` : r <= 1 / TYPO_FOLD ? `1/${Math.round(1 / r)} of` : null;
}

/** How one value just measured reads, the one rule for Monitor's checks and the day's entry: impossible (a negative
 *  volume, a body weight not above 0), a likely typo against the animal's last value (`before`), outside the range
 *  that asks to measure again (`lower`–`upper`, volumes only), else ok; `why` in words. Null: no value. */
export type ReadingCheck = {kind: 'ok' | 'again' | 'typo' | 'impossible', why: string};
export function checkReading(v: number | null | undefined, o: {weight?: boolean, before?: number | null,
  lower?: number | null, upper?: number | null}): ReadingCheck | null {
  if (v === null || v === undefined || !Number.isFinite(v))
    return null;
  const unit = o.weight ? 'g' : 'mm³';
  if (o.weight ? v <= 0 : v < 0)
    return {kind: 'impossible', why: o.weight ? 'A body weight is above 0 g' : 'A volume cannot be negative'};
  if (v > (o.weight ? MAX_WEIGHT : MAX_VOLUME)) {
    return {kind: 'impossible', why: `${mm3Text(v)} ${unit} is not a possible ${o.weight ? 'body weight' :
      'tumor volume'}`};
  }
  // a tumor gone after a value: asked, never refused (a complete response reads 0; QA 0.1.3 round 4 chaos N8)
  if (!o.weight && v === 0 && (o.before ?? 0) > 0)
    return {kind: 'again', why: '0 mm³: tumor gone? Check it; a second reading that agrees stands'};
  const t = o.before ? likelyTypo(o.before, v) : null;
  if (t)
    return {kind: 'typo', why: `Likely a typo: ${t} the last value (${mm3Text(o.before!)} ${unit}): check it`};
  if (!o.weight && o.lower != null && o.upper != null && (v < o.lower || v > o.upper)) {
    return {kind: 'again', why: `Outside the expected ${mm3Text(o.lower)}–${mm3Text(o.upper)} mm³: measure again; ` +
      'a second reading that agrees stands'};
  }
  return {kind: 'ok', why: ''};
}

/** Above these a value is no measurement at all (impossible, not a typo to check): a tumor of 20 cm³, far past any
 *  humane endpoint (QA 0.1.3 round 4 technician R4-4: 50,000 mm³ passed as a "likely typo"); a body weight beyond any
 *  rat's (a heavier one than usual is still only asked, as a typo). ponytail: fixed bounds, a lab setting if a
 *  species ever needs others. */
export const MAX_VOLUME = 20_000;
export const MAX_WEIGHT = 2_000;

/** What a typed value that is no number reads as (NaN: "abc", "1e9"). */
export const NOT_A_NUMBER: ReadingCheck = {kind: 'impossible', why: 'Not a number: write it as a plain number'};

/** A typed value, read as the day's entry and a correction both read it (null: empty; NaN: not a number), checked:
 *  the one check for both (QA 0.1.3 round 3: technician N3, N4; chaos 8). */
export const checkTyped = (v: number | null, o: Parameters<typeof checkReading>[1]): ReadingCheck | null =>
  v !== null && Number.isNaN(v) ? NOT_A_NUMBER : checkReading(v, o);

/** A volume in words, thousands separated, its sign as a minus (10,790; −50, not -50). */
export const mm3Text = (v: number) => minus(Math.round(v).toLocaleString('en-US'));
