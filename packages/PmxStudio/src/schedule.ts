/* The one reader of dosing schedules as protocols write them, and the days they dose: the arm's name (arm-name.ts
   finds the schedule in it with SCHEDULE), the study calendar, the modelling dataset, the doses a day the exposure
   multiple reads and the dose totals all come from here. QD, BID, TID, QOD, Q3D, QW, BIW, Q2W, Q12H, "every 3 days",
   "twice weekly"; a count (×14, 5 doses), a first day (from day 3, d3–d16, days 0–4), 5 on / 2 off; continuous
   exposure (pump, drinking water, chow, diet) every day. The study's rules (no dosing at weekends) move the days,
   never the schedule as written. */

import {sig} from './format';

// the parts of a schedule as a name writes them (lower case); SCHEDULE finds one in a name
// (the unit right after the number: "Q1 w/ PD-1" is not Q1W)
const RATE = 'q\\d+(?:d|h|w)|qd|bid|tid|qid|qod|qw|biw|tiw|' +
  '(?:once|twice|three times|\\d+ times) (?:a |per )?(?:day|daily|week|weekly)|daily|weekly|' +
  'every (?:other|\\d+) (?:day|days|week|weeks)|\\d+\\s*(?:days?\\s*)?on\\s*[/,]?\\s*\\d+\\s*(?:days?\\s*)?off';
const COUNT = '[x×]\\s*\\d+(?:\\s*doses)?|\\d+\\s*doses';
const FIRST = '(?:from|starting(?:\\s+on)?|start(?:ing)?)\\s+day\\s*\\d+|d\\d+\\s*[-–]\\s*d\\d+|' +
  'days?\\s*\\d+\\s*(?:[-–]|to)\\s*\\d+';
/** A schedule in an arm's name: a rate with its count and first day, a count alone ("x3") or days alone ("d0–d4").
 *  A repeat count is a word of its own, not the X1 of "Cpd-X1". */
export const SCHEDULE = new RegExp(`\\b(?:(?:${RATE})(?:\\s*(?:${COUNT}))?(?:\\s*,?\\s*(?:${FIRST}))?|` +
  `(?<![a-z-])(?:${COUNT})(?:\\s*,?\\s*(?:${FIRST}))?|(?<![a-z-])(?:${FIRST}))(?![a-z])`, 'i');

/** Exposure that does not come in doses: every day, never moved. Read from the schedule or the route only, never
 *  from a group's name ("High-fat diet + Cpd-A Q3D" is dosed Q3D, a diet control is not dosed). */
const CONTINUOUS = /\b(?:(?:mini)?pump|drinking water|in (?:the )?(?:water|diet|chow|food|feed))\b/i;

/** A schedule as read: the days of a cycle it doses (QD [0] of 1, BIW [0, 3] of 7, 5 on / 2 off [0…4] of 7), the
 *  doses on a dosed day (BID 2), how many dosed days in all (null: until its last day or the study's), its first day
 *  and its last (null: none written). */
export interface Schedule {
  kind: 'daily' | 'interval' | 'continuous';
  pattern: number[];
  cycle: number;
  perDay: number;
  count: number | null;
  first: number;
  end: number | null;
  /** written as days on and off: no dosing at weekends doses its weekdays */
  onOff?: boolean;
}

/** A group's dose changed on a study day, in Plan or in Monitor: skipped, or added on a day its schedule does not
 *  dose; why, who and when. The schedule as written stays as it is. */
export interface DoseException {day: number, group: string, change: 'skip' | 'add', note: string, who: string,
  when: string}

/** What the study adds to every schedule: day 0's date (YYYY-MM-DD; null: none yet), whether no one doses at
 *  weekends, and its dated exceptions (applied per group by excepted). */
export interface DosingRules {start: string | null; weekendsOff?: boolean; exceptions?: DoseException[]}

/** A group's dosed days (`days`, as dosedDays gives them) with the study's exceptions for it: a skipped day out, an
 *  added one in (from day 0 up to `last`; null: no end). Every reader of dosing days goes through it. */
export function excepted(days: number[], group: string, rules: DosingRules | null, last: number | null): number[] {
  const own = (rules?.exceptions ?? []).filter((e) => e.group === group);
  if (!own.length)
    return days;
  const set = new Set(days);
  for (const e of own) {
    if (e.change === 'skip')
      set.delete(e.day);
    else if (e.day >= 0 && (last === null || e.day <= last))
      set.add(e.day);
  }
  return [...set].sort((a, b) => a - b);
}

const spread = (n: number) => n === 2 ? [0, 3] : n === 3 ? [0, 2, 4] :
  Array.from({length: n}, (_, i) => Math.floor(i * 7 / n));
const word = (w: string) => w === 'once' ? 1 : w === 'twice' ? 2 : w === 'three times' ? 3 : parseInt(w);

/** The rate of a schedule's rate word: its pattern, cycle and doses a dosed day; null: not one. */
function rateOf(r: string): Pick<Schedule, 'pattern' | 'cycle' | 'perDay'> | null {
  const s = r.replace(/\s+/g, ' ').trim();
  const fixed: {[s: string]: [number[], number, number]} = {qd: [[0], 1, 1], daily: [[0], 1, 1], bid: [[0], 1, 2],
    tid: [[0], 1, 3], qid: [[0], 1, 4], qod: [[0], 2, 1], qw: [[0], 7, 1], weekly: [[0], 7, 1], biw: [[0, 3], 7, 1],
    tiw: [[0, 2, 4], 7, 1]};
  if (s in fixed)
    return {pattern: fixed[s][0], cycle: fixed[s][1], perDay: fixed[s][2]};
  let m = /^q ?(\d+) ?(h|d|w)$/.exec(s);
  if (m) {
    const k = Number(m[1]);
    if (!(k > 0))
      return null;
    if (m[2] !== 'h')
      return {pattern: [0], cycle: k * (m[2] === 'w' ? 7 : 1), perDay: 1};
    return 24 % k === 0 ? {pattern: [0], cycle: 1, perDay: 24 / k} :
      k % 24 === 0 ? {pattern: [0], cycle: k / 24, perDay: 1} : null;
  }
  m = /^(once|twice|three times|(\d+) times) (?:a |per )?(day|daily|week|weekly)$/.exec(s);
  if (m) {
    const n = word(m[1]);
    if (!(n > 0))
      return null;
    // ponytail: 7 a week or more only in whole doses a day (14 a week is BID); 10 a week is not read
    return m[3].startsWith('d') ? {pattern: [0], cycle: 1, perDay: n} : n < 7 ?
      {pattern: spread(n), cycle: 7, perDay: 1} : n % 7 === 0 ? {pattern: [0], cycle: 1, perDay: n / 7} : null;
  }
  m = /^every (other|\d+) (day|days|week|weeks)$/.exec(s);
  if (m) {
    const k = m[1] === 'other' ? 2 : Number(m[1]);
    return k > 0 ? {pattern: [0], cycle: k * (m[2].startsWith('week') ? 7 : 1), perDay: 1} : null;
  }
  return null;
}

/** A schedule as protocols write it, read (null: no rate in it; "×3" or "5 doses" alone give none). `route`: the
 *  group's route as its details give it, which may say continuous exposure (pump, drinking water) when the schedule
 *  does not. */
export function readSchedule(text: string | null | undefined, route = ''): Schedule | null {
  let s = ` ${(text ?? '').toLowerCase().replace(/[–—]/g, '-')} `;
  const cut = (re: RegExp): RegExpExecArray | null => {
    const m = re.exec(s);
    if (m)
      s = s.slice(0, m.index) + ' ' + s.slice(m.index + m[0].length);
    return m;
  };
  const range = cut(/\bd(\d+)\s*-\s*d(\d+)/) ?? cut(/\bdays?\s*(\d+)\s*(?:-|to)\s*(\d+)/);
  const from = cut(/(?:from|starting(?:\s+on)?|start(?:ing)?)\s+day\s*(\d+)/);
  const first = range ? Number(range[1]) : from ? Number(from[1]) : 0;
  const end = range && Number(range[2]) >= first ? Number(range[2]) : null;
  if (CONTINUOUS.test(s) || (!s.trim() && CONTINUOUS.test(route)))
    return {kind: 'continuous', pattern: [0], cycle: 1, perDay: 1, count: null, first, end};
  const onOff = cut(/(\d+)\s*(?:days?\s*)?on\s*[/,]?\s*(\d+)\s*(?:days?\s*)?off/);
  // (QDx14: the count may follow the rate with no space)
  const rateWord = onOff ? null : cut(new RegExp(`\\b(?:${RATE})(?=$|[^a-z]|x\\s*\\d)`));
  const count = cut(/[x×*]\s*(\d+)(?:\s*doses)?|(\d+)\s*doses/);
  const n = count ? Number(count[1] ?? count[2]) : null;
  let rate: Pick<Schedule, 'pattern' | 'cycle' | 'perDay'> | null;
  if (onOff) {
    const [on, off] = [Number(onOff[1]), Number(onOff[2])];
    rate = on > 0 ? {pattern: Array.from({length: on}, (_, i) => i), cycle: on + off, perDay: 1} : null;
  } else
    rate = rateWord ? rateOf(rateWord[0]) : range ? {pattern: [0], cycle: 1, perDay: 1} : null;
  if (!rate)
    return null;
  return {kind: rate.cycle === 1 || onOff ? 'daily' : 'interval', ...rate,
    // 5 on / 2 off ×3: three cycles
    count: n === null || !(n > 0) ? null : onOff ? n * rate.pattern.length : n, first, end,
    ...onOff ? {onOff: true} : {}};
}

const gcd = (a: number, b: number): number => b ? gcd(b, a % b) : a;
/** The days a schedule spans when nothing ends it: whole weeks and whole cycles (QD 7, Q3D 21, Q2W 14). */
const spanOf = (s: Schedule) => s.cycle * 7 / gcd(s.cycle, 7);

/** The weekday of day 0 (Monday 0) when the study doses no weekends and day 0 has a valid date; else null. */
function weekdayZero(rules: DosingRules | null): number | null {
  const iso = rules?.weekendsOff ? validDay(rules.start) : null;
  return iso ? (new Date(`${iso}T00:00:00Z`).getUTCDay() + 6) % 7 : null;
}

/** `iso` when it is a date that exists (YYYY-MM-DD; 2026-02-31 is not), else null. */
export function validDay(iso: string | null | undefined): string | null {
  if (!iso || !/^\d{4}-\d{2}-\d{2}$/.test(iso))
    return null;
  const t = Date.parse(`${iso}T00:00:00Z`);
  return !Number.isNaN(t) && new Date(t).toISOString().slice(0, 10) === iso ? iso : null;
}

/** The days a schedule doses, from its first: up to `last` (the study's last day; null: its own span when nothing
 *  else ends it), under the study's rules. With no dosing at weekends and day 0's date known, daily schedules skip
 *  Saturdays and Sundays (a count, QD×14, is of doses: they go on past the weekends); others move a weekend dose to
 *  the nearest weekday (Saturday to Friday, Sunday to Monday, the other way when the arm is dosed that day already or
 *  the move would leave its days or the study), keeping the count; continuous exposure stays. */
export function dosedDays(s: Schedule | null, last: number | null = null, rules: DosingRules | null = null): number[] {
  if (!s)
    return [];
  const weekday0 = s.kind === 'continuous' ? null : weekdayZero(rules);
  const weekend = (d: number) => weekday0 !== null && ((weekday0 + d) % 7 + 7) % 7 >= 5;
  // 5 on / 2 off at weekends off: weekdays, as many as its cycles dose
  const sched: Schedule = s.onOff && weekday0 !== null ? {...s, pattern: [0], cycle: 1} : s;
  const skip = weekday0 !== null && sched.kind === 'daily';
  const stop = sched.end ?? (sched.count === null ? last ?? sched.first + spanOf(sched) - 1 : Infinity);
  const nominal: number[] = [];
  for (let k = 0; nominal.length < (sched.count ?? Infinity) && k < 5000; k++) {
    const p = sched.pattern.length;
    const d = sched.first + Math.floor(k / p) * sched.cycle + sched.pattern[k % p];
    if (d > stop)
      break;
    if (!(skip && weekend(d)))
      nominal.push(d);
  }
  let days = nominal;
  if (weekday0 !== null && sched.kind === 'interval') {
    // a moved dose stays within the schedule's days and the study's
    const ceiling = Math.min(sched.end ?? Infinity, last ?? Infinity);
    const taken = new Set(nominal.filter((d) => !weekend(d)));
    days = [];
    for (const d of nominal) {
      if (!weekend(d)) {
        days.push(d);
        continue;
      }
      // Saturday: Friday, else Monday, Thursday, Tuesday…; Sunday: Monday, else Friday, Tuesday, Thursday…
      const sat = ((weekday0 + d) % 7 + 7) % 7 === 5;
      const tries = sat ? [-1, 2, -2, 3, -3, 4] : [1, -2, 2, -3, 3, -4];
      const to = tries.map((t) => d + t)
        .find((x) => x >= sched.first && x <= ceiling && !weekend(x) && !taken.has(x));
      if (to !== undefined) {
        taken.add(to);
        days.push(to);
      }
    }
    days.sort((a, b) => a - b);
  }
  return last === null ? days : days.filter((d) => d <= last);
}

/** Doses a day of a schedule on average (QD 1, BID 2, QOD ½, Q3D ⅓, QW 1/7, BIW 2/7, Q12H 2), as dosedDays places
 *  them: with no dosing at weekends, a daily schedule doses 5 days in 7 whatever day 0's weekday (a count goes on past
 *  the weekends); a moved dose is still a dose. null: no rate ("×3", "5 doses"), or continuous exposure (no doses to
 *  count). */
export function dosesPerDay(schedule: string | Schedule | null, rules: DosingRules | null = null): number | null {
  const s = typeof schedule === 'string' ? readSchedule(schedule) : schedule;
  if (!s || s.kind === 'continuous')
    return null;
  return s.kind === 'daily' && weekdayZero(rules) !== null ? s.perDay * 5 / 7 : s.perDay * s.pattern.length / s.cycle;
}

/** The doses of a schedule up to day `last` (dosed days × doses a dosed day). */
export const doseCount = (s: Schedule | null, last: number | null, rules: DosingRules | null = null) =>
  s ? dosedDays(s, last, rules).length * s.perDay : 0;

/** What a group needs of its compound: losses in preparation, and a mouse when no weight is known (kg). */
export const OVERAGE = 0.2;
export const DEFAULT_WEIGHT_KG = 0.025;

/** A group's doses and what they add up to: the dose in all per kg and the compound to prepare (mg, with the overage)
 *  for `animals` of `weightKg` (null: a 25 g mouse); only for doses in mg/kg (others: the count alone). */
export interface DoseTotal {doses: number, perKg: number | null, mg: number | null, assumedWeight: boolean}
export function doseTotal(doses: number, dose: {value: number, unit: string} | null, animals: number | null,
  weightKg: number | null): DoseTotal {
  const mgkg = dose && dose.value > 0 && /^mg\s*\/\s*kg$/i.test(dose.unit.trim());
  const perKg = mgkg ? dose!.value * doses : null;
  const known = weightKg !== null && weightKg > 0;
  return {doses, perKg, mg: perKg !== null && animals !== null && animals > 0 ?
    perKg * animals * (known ? weightKg : DEFAULT_WEIGHT_KG) * (1 + OVERAGE) : null, assumedWeight: !known};
}

const amount = (x: number) => x >= 100 ? sig(Math.round(x), 15) : sig(x, 2);

/** "14 doses · 280 mg/kg in total · about 84 mg of compound" (`given`: the doses so far, "6 of 14 doses · 120 of
 *  280 mg/kg"), and its tooltip. */
export function doseLine(t: DoseTotal, animals: number | null, given?: DoseTotal,
  withAmount = true): {text: string, tip: string} {
  const n = (x: number) => `${amount(x)} dose${x === 1 ? '' : 's'}`;
  const parts = [given ? `${amount(given.doses)} of ${n(t.doses)} so far` : n(t.doses)];
  if (t.perKg !== null) {
    parts.push(given && given.perKg !== null ? `${amount(given.perKg)} of ${amount(t.perKg)} mg/kg` :
      `${amount(t.perKg)} mg/kg in total`);
  }
  const compound = t.mg === null ? '' :
    `about ${amount(t.mg)} mg of compound for the group${t.assumedWeight ? ', assuming 25 g mice' : ''}`;
  // Monitor: the compound is ordering detail, so it stays in the tooltip
  if (compound && withAmount)
    parts.push(compound);
  const tip = t.mg === null ? 'The doses the schedule gives over the study.' :
    `The doses the schedule gives over the study. ${withAmount ? '' : `You would prepare ${compound}. `}` +
    `Compound: dose × doses × ${animals} animals × ` +
    `${t.assumedWeight ? '25 g' : 'their mean body weight'}, plus ${OVERAGE * 100}% for losses in preparation. ` +
    'The amount is for the compound as dosed; if a salt is weighed, scale it by the salt to free base ratio.';
  return {text: parts.join(' · '), tip};
}

/** A row of the study for groupWeights: its group, animal, day, body weight and whether it is excluded. */
export interface WeightRow {group: string, animal: string, day: number, weight: number | null, excluded: boolean}

/** Each group's animals still in the study (measured, not excluded, on the group's last day: one gone before then
 *  is dead or removed) and their mean last body weight in kg (null: none weighed). `unit`: the weight column's unit;
 *  none: a study mean of 10 or more is grams, below 10 kilograms. */
export function groupWeights(rows: WeightRow[], unit: string | null):
  Map<string, {animals: number, kg: number | null}> {
  const kept = rows.filter((r) => !r.excluded);
  const lastOf = new Map<string, number>();
  for (const r of kept)
    lastOf.set(r.group, Math.max(lastOf.get(r.group) ?? -Infinity, r.day));
  // per group, each animal's last day and its last weight
  const by = new Map<string, Map<string, {day: number, w: number | null, wDay: number}>>();
  for (const r of kept) {
    const animals = by.get(r.group) ?? new Map();
    by.set(r.group, animals);
    const a = animals.get(r.animal) ?? {day: -Infinity, w: null, wDay: -Infinity};
    a.day = Math.max(a.day, r.day);
    if (r.weight !== null && Number.isFinite(r.weight) && r.weight > 0 && r.day >= a.wDay)
      [a.w, a.wDay] = [r.weight, r.day];
    animals.set(r.animal, a);
  }
  const alive = [...by].map(([g, as]) => [g, [...as.values()].filter((a) => a.day === lastOf.get(g))] as const);
  const all = alive.flatMap(([, as]) => as.map((a) => a.w)).filter((w): w is number => w !== null);
  const u = (unit ?? '').trim().toLowerCase();
  const toKg = u === 'kg' ? 1 : u === 'g' ? 1e-3 : u === 'mg' ? 1e-6 :
    all.length && all.reduce((s, w) => s + w, 0) / all.length >= 10 ? 1e-3 : 1;
  return new Map(alive.map(([g, as]) => {
    const ws = as.map((a) => a.w).filter((w): w is number => w !== null);
    return [g, {animals: as.length, kg: ws.length ? ws.reduce((s, w) => s + w, 0) / ws.length * toKg : null}];
  }));
}
