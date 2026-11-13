import {SUMMARY, SUMMARY_LIKE, dayOf, dayWord, letter, numberOf, readWide} from './wide-table';

/** The answers given to the reader's questions, by question id. */
export type Answers = {[id: string]: string};
/** A question the reader cannot answer from the table: one choice among `options` ([value, label]). With
 *  `preset`, a choice already made that the user may change; without, nothing is read until it is answered. */
export interface Question {id: string, text: string, options: [string, string][], preset?: string}

export type VolumeFormula = 'LW2/2' | 'pi/6*LW2' | 'pi/6*LWH' | 'LWH/2';
/** Where a study's volumes come from (study-file.ts keeps it in the study's provenance): as the table gave them,
 *  computed from length and width with `formula`, or given and found to follow `formula`. */
export interface VolumeSource {formula: VolumeFormula | null, source: 'given' | 'computed' | 'deduced'}

interface Formula {label: string, height: boolean, v: (l: number, w: number, h: number) => number}
export const FORMULAS: {[f in VolumeFormula]: Formula} = {
  'LW2/2': {label: 'length × width² / 2', height: false, v: (l, w) => l * w * w / 2},
  'pi/6*LW2': {label: 'π/6 × length × width²', height: false, v: (l, w) => Math.PI / 6 * l * w * w},
  'pi/6*LWH': {label: 'π/6 × length × width × height', height: true, v: (l, w, h) => Math.PI / 6 * l * w * h},
  'LWH/2': {label: 'length × width × height / 2', height: true, v: (l, w, h) => l * w * h / 2},
};

export interface StudyRow {animal: string, group: string, day: number, volume: number | null,
  weight: number | null}
export interface TableRead {
  shape: 'long' | 'by animal' | 'by day',
  rows: StudyRow[],
  /** Animals read per group, in the table's order. */
  groups: [string, number][],
  /** What was read and how, in plain words: the volumes' source, the units, day 0. */
  notes: string[],
  /** What was left out on purpose, each with its reason (never a group's name). */
  ignored: string[],
  /** The table's other columns, kept as they are (a value per row): numbers, or text. */
  extra: [string, (number | string | null)[]][],
  volume: VolumeSource | null,
  questions: Question[],
  /** Why nothing was read; null when read (or when a question waits). */
  error: string | null,
}

// ---- roles ---------------------------------------------------------------------------------------------------------

export type Role = 'animal' | 'group' | 'day' | 'date' | 'start' | 'volume' | 'length' | 'width' | 'height' | 'weight';
/** The names each role is known by, lowercased, a unit and the separators taken out ("Tumor Volume (mm3)" is
 *  tumor_volume): closed lists, as in NCA's column-autodetect; a role is filled only by exactly one column. */
const SYNONYMS: {readonly [r in Role]: readonly string[]} = Object.freeze({
  animal: ['animal', 'animal_id', 'animalid', 'animal_no', 'animal_number', 'mouse', 'mouse_id', 'mouseid', 'mouse_no',
    'rat', 'rat_id', 'subject', 'subject_id', 'subjid', 'id', 'tag', 'ear_tag', 'tumor_id', 'tumour_id'],
  group: ['group', 'groups', 'group_name', 'arm', 'treatment', 'treatment_group', 'cohort', 'tx', 'dose_group'],
  day: ['day', 'days', 'study_day', 'studyday', 'day_of_study', 'time', 'timepoint', 'time_point', 'dpi',
    'days_post_implant', 'days_post_inoculation', 'week', 'weeks', 'wk', 'wks'],
  date: ['date', 'measurement_date', 'date_of_measurement', 'meas_date', 'visit_date', 'assessment_date'],
  start: ['randomization_date', 'randomisation_date', 'rand_date', 'randomization', 'randomisation',
    'date_of_randomization', 'date_of_randomisation', 'treatment_start', 'treatment_start_date', 'start_date',
    'dosing_start', 'dosing_start_date', 'first_dose_date'],
  volume: ['volume', 'vol', 'tv', 'tvol', 'tumor_volume', 'tumour_volume', 'tumor_vol', 'tumour_vol', 'tumorvolume'],
  length: ['length', 'l', 'tumor_length', 'tumour_length', 'long_axis', 'major_axis', 'longest_diameter'],
  width: ['width', 'w', 'tumor_width', 'tumour_width', 'short_axis', 'minor_axis', 'shortest_diameter'],
  height: ['height', 'h', 'depth', 'thickness', 'tumor_height', 'tumour_height'],
  weight: ['weight', 'wt', 'bw', 'body_weight', 'bodyweight', 'body_wt', 'animal_weight', 'mouse_weight'],
});
const ROLE_WORD: {[r in Role]: string} = {animal: 'animal', group: 'group', day: 'day', date: 'measurement date',
  start: 'randomization or treatment start date', volume: 'tumor volume', length: 'tumor length',
  width: 'tumor width', height: 'tumor height', weight: 'body weight'};

type Unit = 'mm' | 'cm' | 'mm3' | 'cm3' | 'ml' | 'ul' | 'g' | 'kg' | 'mg';
/** Each unit in the study's own (mm, mm³, g), and the roles it may be read for. */
const UNITS: {[u in Unit]: {factor: number, of: 'length' | 'volume' | 'weight', name: string}} = {
  mm: {factor: 1, of: 'length', name: 'mm'}, cm: {factor: 10, of: 'length', name: 'cm'},
  mm3: {factor: 1, of: 'volume', name: 'mm³'}, ul: {factor: 1, of: 'volume', name: 'µL'},
  cm3: {factor: 1000, of: 'volume', name: 'cm³'}, ml: {factor: 1000, of: 'volume', name: 'mL'},
  g: {factor: 1, of: 'weight', name: 'g'}, kg: {factor: 1000, of: 'weight', name: 'kg'},
  mg: {factor: 0.001, of: 'weight', name: 'mg'},
};
const KIND: {[r in Role]?: 'length' | 'volume' | 'weight'} = {length: 'length', width: 'length', height: 'length',
  volume: 'volume', weight: 'weight'};

const unitOf = (s: string): Unit | null => {
  const u = s.toLowerCase().replace(/[µμ]/g, 'u').replace(/³|\^3/g, '3');
  return u in UNITS ? u as Unit : null;
};

/** A column's name as its role's key and its unit: the unit in brackets ("Volume (mm3)", "[cm]") or after a space or
 *  an underscore ("weight_g"); a bracket that is not a unit stays in the key, so "Body weight (% change)" is no
 *  body weight. */
export function nameOf(header: string): {key: string, unit: Unit | null} {
  let name = header.trim().toLowerCase();
  let unit: Unit | null = null;
  const m = /\s*[([]\s*([^)\]]+?)\s*[)\]]$/.exec(name) ?? /[\s_]+([a-zµμ]{1,2}(?:3|³|\^3)?)$/.exec(name);
  if (m && unitOf(m[1])) {
    unit = unitOf(m[1]);
    name = name.slice(0, m.index);
  }
  return {key: name.replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, ''), unit};
}

/** The role a column's name gives it, with its unit (a unit of another kind than the role's: none). */
export function roleOf(header: string): {role: Role, unit: Unit | null} | null {
  const {key, unit} = nameOf(header);
  const role = (Object.keys(SYNONYMS) as Role[]).find((r) => SYNONYMS[r].includes(key));
  if (!role || unit && UNITS[unit].of !== KIND[role])
    return null;
  return {role, unit};
}

// ---- dates ---------------------------------------------------------------------------------------------------------

/** A date cell: ISO (2026-03-02, a time after it or not), or two numbers and a four-digit year with / . or - between
 *  them, whose order (day or month first) the column says. */
const ISO = /^(\d{4})-(\d{1,2})-(\d{1,2})(?:[T ][\d:.]+(?:Z|[+-]\d{2}:?\d{2})?)?$/;
const DMY = /^(\d{1,2})[/.-](\d{1,2})[/.-](\d{4})(?:[T ][\d:.]+)?$/;
const utcDay = (y: number, m: number, d: number): number | null => {
  const t = Date.UTC(y, m - 1, d);
  const x = new Date(t);
  return x.getUTCFullYear() === y && x.getUTCMonth() === m - 1 && x.getUTCDate() === d ? t / 864e5 : null;
};
export const isoOf = (day: number) => new Date(day * 864e5).toISOString().slice(0, 10);
const isDate = (s: string) => ISO.test(s) || DMY.test(s);

/** The day (days since 1970) of a date cell, `order` saying how a d/m/y or m/d/y is written; null when it is none. */
export function dateOf(s: string, order: 'dmy' | 'mdy'): number | null {
  const iso = ISO.exec(s);
  if (iso)
    return utcDay(+iso[1], +iso[2], +iso[3]);
  const m = DMY.exec(s);
  return m ? order === 'dmy' ? utcDay(+m[3], +m[2], +m[1]) : utcDay(+m[3], +m[1], +m[2]) : null;
}

/** How a column's dates are written, from all of them: day first when one starts above 12, month first when one has
 *  its second number above 12, unknown when none does (then the user is asked), mixed when both. */
export function dateOrder(values: string[]): 'dmy' | 'mdy' | 'unknown' | 'mixed' {
  let dmy = false;
  let mdy = false;
  for (const s of values) {
    const m = DMY.exec(s);
    if (m && +m[1] > 12)
      dmy = true;
    if (m && +m[2] > 12)
      mdy = true;
  }
  return dmy && mdy ? 'mixed' : dmy ? 'dmy' : mdy ? 'mdy' : values.some((s) => DMY.test(s)) ? 'unknown' : 'dmy';
}

// ---- volumes -------------------------------------------------------------------------------------------------------

/** The formula a lab's volumes follow from its lengths and widths (and heights): the one whose volumes come
 *  within 1% of the given ones at the median and within 2% in 90% of the rows with all of them, closest first; null
 *  when none does, or with fewer than five rows to tell by. */
export function deduceFormula(sets: {l: number, w: number, h: number | null, v: number}[]): VolumeFormula | null {
  const usable = sets.filter((x) => x.v > 0);
  // five rows at least: fewer cannot tell a formula from a coincidence
  if (usable.length < 5)
    return null;
  // with a height column the formulas with a height come first; one is chosen only if it beats the others
  const tall = usable.some((x) => x.h !== null);
  const all = Object.keys(FORMULAS) as VolumeFormula[];
  const order = tall ? [...all.filter((f) => FORMULAS[f].height), ...all.filter((f) => !FORMULAS[f].height)] : all;
  let best: VolumeFormula | null = null;
  let bestError = Infinity;
  for (const f of order) {
    const rows = FORMULAS[f].height ? usable.filter((x) => x.h !== null) : usable;
    if (rows.length === 0 || FORMULAS[f].height && rows.length < usable.length)
      continue;
    const errors = rows.map((x) => Math.abs(FORMULAS[f].v(x.l, x.w, x.h ?? 0) - x.v) / x.v).sort((a, b) => a - b);
    const median = errors[Math.floor((errors.length - 1) / 2)];
    if (median < 0.01 && errors.filter((e) => e <= 0.02).length >= 0.9 * errors.length && median < bestError) {
      best = f;
      bestError = median;
    }
  }
  return best;
}

const median = (xs: number[]) => {
  const s = [...xs].sort((a, b) => a - b);
  return s.length === 0 ? NaN : (s[Math.floor((s.length - 1) / 2)] + s[Math.ceil((s.length - 1) / 2)]) / 2;
};

// ---- the reader ----------------------------------------------------------------------------------------------------

const plural = (n: number, word: string) => `${n} ${word}${n === 1 ? '' : 's'}`;
/** At most six cells named, the rest counted. */
const fewOf = (xs: string[]) =>
  xs.length > 6 ? `${xs.slice(0, 6).join(', ')}, … and ${xs.length - 6} more` : xs.join(', ');
const listed = (xs: string[]) => `${xs.slice(0, 6).join(', ')}${xs.length > 6 ? ', …' : ''}`;
const andList = (xs: string[]) => xs.length > 1 ? `${xs.slice(0, -1).join(', ')} and ${xs[xs.length - 1]}` : xs[0];
const empty = (shape: TableRead['shape']): TableRead => ({shape, rows: [], groups: [], notes: [], ignored: [],
  extra: [], volume: null, questions: [], error: null});
const refused = (shape: TableRead['shape'], error: string): TableRead => ({...empty(shape), error});
const has = (x: number | null | undefined): x is number => x !== null && x !== undefined;

function counted(rows: StudyRow[]): [string, number][] {
  const seen = new Map<string, Set<string>>();
  rows.forEach((x) => seen.set(x.group, (seen.get(x.group) ?? new Set()).add(x.animal)));
  return [...seen].map(([g, a]) => [g, a.size]);
}

/** The questions on a unit assumed, asked when the volumes come out as no study's. */
const UNIT_QUESTIONS: {[kind in 'volume' | 'length' | 'weight' | 'weightBig']: Question} = {
  volume: {id: 'unit:volume', text: 'The volumes are very small or very large for mm³. In which unit are they?',
    options: [['mm3', 'mm³ (or µL), as they are'], ['cm3', 'cm³ (or mL): multiplied by 1000']]},
  length: {id: 'unit:length', text: 'The volumes come out very small or very large from lengths in mm. In which ' +
    'unit are the lengths and widths?', options: [['mm', 'mm, as they are'], ['cm', 'cm: multiplied by 10']]},
  weight: {id: 'unit:weight', text: 'The body weights are very small for grams. In which unit are they?',
    options: [['g', 'g, as they are'], ['kg', 'kg: multiplied by 1000']]},
  weightBig: {id: 'unit:weight', text: 'The body weights are very large for grams. In which unit are they?',
    options: [['g', 'g, as they are'], ['mg', 'mg: divided by 1000']]},
};
const STANDARD: {[kind in 'length' | 'volume' | 'weight']: Unit} = {length: 'mm', volume: 'mm3', weight: 'g'};
const KIND_WORD = {length: 'Lengths', volume: 'Volumes', weight: 'Body weights'};

/** The units of the measurements read (`units`: each kind's units as named or answered, null when assumed): what
 *  was converted and what assumed, in words; or the question on the unit the volumes come from (`from`), when it was
 *  assumed and they are not a study's (a median below 1 or above 50 000 mm³), or the weights'
 *  (below 1 or above 1000 g). */
function withUnits(read: TableRead, units: [kind: 'length' | 'volume' | 'weight', unit: Unit | null][],
  from: 'length' | 'volume'): TableRead {
  const vs = read.rows.flatMap((r) => has(r.volume) && r.volume > 0 ? [r.volume] : []);
  const m = median(vs);
  const ws = read.rows.flatMap((r) => has(r.weight) && r.weight > 0 ? [r.weight] : []);
  const asked: Question[] = [];
  if (units.some(([k, u]) => k === from && u === null) && (m < 1 || m > 50000))
    asked.push(UNIT_QUESTIONS[from]);
  if (units.some(([k, u]) => k === 'weight' && u === null) && (median(ws) < 1 || median(ws) > 1000))
    asked.push(median(ws) < 1 ? UNIT_QUESTIONS.weight : UNIT_QUESTIONS.weightBig);
  if (asked.length)
    return {...read, rows: [], notes: [], ignored: [], extra: [], questions: [...read.questions, ...asked]};
  const converted = units.filter(([, u]) => u && UNITS[u].factor !== 1)
    .filter(([k], i, all) => all.findIndex(([x]) => x === k) === i)
    .map(([k, u]) => `${KIND_WORD[k]} converted from ${UNITS[u!].name} to ${UNITS[STANDARD[k]].name}.`);
  const assumed = [...new Set(units.filter(([, u]) => u === null).map(([k]) => UNITS[STANDARD[k]].name))];
  const day0 = read.notes.filter((n) => n.startsWith('Day 0'));
  return {...read, groups: counted(read.rows), notes: [...read.notes.filter((n) => !n.startsWith('Day 0')),
    ...converted, ...assumed.length ? [`No unit in the column names: assumed ${andList(assumed)}.`] : [], ...day0]};
}

/** A table of measurements read, in steps: the columns' roles by their names, the shape (one row per measurement,
 *  or wide: wide-table.ts), the volumes (given, computed from the caliper, or the lab's own formula found from both),
 *  the units, the days (numbers, or dates from a day 0), and the check that every measurement cell is placed at one
 *  animal and day or left out with its reason, else nothing is read. `answers`: to the questions it asked before (it
 *  asks only what it cannot tell); `decimalComma` where semicolons separate the cells. Pure, so tested in node.
 *  The seam for spreadsheets: it takes one sheet's cells as text, so an .xlsx sheet (PowerPack's xlsx handler, a
 *  soft dependency) plugs in as its cells, numbers written with a dot and dates as ISO text. */
export function readStudy(cells: string[][], answers: Answers = {}, decimalComma = false): TableRead {
  const wide = readWide(cells, decimalComma);
  if (!wide)
    return readLong(cells, answers, decimalComma);
  if (wide.error)
    return refused(wide.shape, wide.error);
  // a wide table holds volumes only; their unit, when a day's heading gives it (Day 7 (mm3)), else as answered
  const named = cells.some((row) => row.some((s) => /\((mm3|mm³|mm\^3)\)/i.test(s))) ? 'mm3' : null;
  const unit: Unit | null = named ?? unitOf(answers['unit:volume'] ?? '');
  const f = UNITS[unit ?? 'mm3'].factor;
  return withUnits({...empty(wide.shape), ignored: wide.ignored, volume: {formula: null, source: 'given'},
    notes: ['Volumes as given.', ...wide.notes], rows: wide.rows.map((x) => ({animal: x.animal, group: x.group,
      day: x.day, volume: x.volume * f, weight: null}))}, [['volume', unit]], 'volume');
}

/** One row per measurement: a header row naming the columns, then a row per animal and day. */
function readLong(cells: string[][], answers: Answers, decimalComma: boolean): TableRead {
  const filled = cells.map((_, r) => r).filter((r) => cells[r].some((s) => s !== ''));
  if (filled.length < 2)
    return refused('long', 'Not read: the table has no rows of measurements under its header.');
  const [h, ...rest] = filled;
  const head = cells[h];
  const at = (r: number, c: number) => cells[r]?.[c] ?? '';
  const name = (c: number) => head[c] ? `column ${letter(c)} ("${head[c]}")` : `column ${letter(c)}`;

  // a column for each role by its name, exactly one; between two the user is asked, never a guess
  const named = head.map(roleOf);
  const col: {[r in Role]?: number} = {};
  const questions: Question[] = [];
  const ask = (role: Role, among: number[]): number | undefined => {
    const id = `role:${role}`;
    if (among.includes(Number(answers[id])) && answers[id] !== '')
      return Number(answers[id]);
    if (among.length) {
      questions.push({id, text: `Which column is the ${ROLE_WORD[role]}?`,
        options: among.map((c) => [String(c), head[c] || `Column ${letter(c)}`])});
    }
    return undefined;
  };
  for (const role of Object.keys(SYNONYMS) as Role[]) {
    const among = named.flatMap((x, c) => x?.role === role ? [c] : []);
    col[role] = among.length === 1 ? among[0] : among.length > 1 ? ask(role, among) : undefined;
  }
  const asked = (re: RegExp) => questions.some((q) => re.test(q.id));
  if ((Object.keys(col) as Role[]).filter((r) => r !== 'start' && col[r] !== undefined).length + questions.length < 2) {
    return refused('long', 'Not read: the columns are not named as a study\'s. Name them Animal, Group, Day (or ' +
      'Date) and Volume (or Length and Width), one row per animal and day; or open a wide table with a column per ' +
      'day or per animal.');
  }
  // a role the study needs that no column is named for is asked among the columns no role took
  const free = () => head.map((_, c) => c).filter((c) => head[c] !== '' && !Object.values(col).includes(c));
  for (const role of ['animal', 'group'] as Role[]) {
    if (col[role] === undefined && !asked(new RegExp(role)))
      col[role] = ask(role, free());
  }
  if (col.day === undefined && col.date === undefined && !asked(/day|date/))
    col.day = ask('day', free());
  const caliper = col.length !== undefined && col.width !== undefined;
  if (col.volume === undefined && !caliper && !asked(/volume|length|width/))
    col.volume = ask('volume', free());
  if (questions.length)
    return {...empty('long'), questions};
  const none = (['animal', 'group'] as Role[]).find((r) => col[r] === undefined) ??
    (col.day === undefined && col.date === undefined ? 'day' : col.volume === undefined && !caliper ? 'volume' : null);
  if (none) {
    return refused('long', `Not read: no column is the ${ROLE_WORD[none]}${none === 'volume' ? ', nor its length and ' +
      'width' : ''}.`);
  }

  // the days as numbers, or dates counted from day 0
  const values = (c: number | undefined) => c === undefined ? [] : rest.map((r) => at(r, c)).filter((s) => s !== '');
  // a Day column with nothing in it, beside a filled date column, is no day: the dates are used
  const emptyDay = col.day !== undefined && col.date !== undefined && values(col.day).length === 0;
  let dayCol = emptyDay ? undefined : col.day;
  // a table with both counts its days as given, its dates kept as they are
  let dateCol = dayCol === undefined ? col.date : undefined;
  if (dayCol !== undefined && !values(dayCol).every((s) => dayOf(s) !== null)) {
    if (!values(dayCol).every(isDate)) {
      const bad = rest.find((r) => at(r, dayCol!) !== '' && dayOf(at(r, dayCol!)) === null)!;
      return refused('long', `Not read: ${name(dayCol)} holds neither days (0, 7, D7 or Day 7) nor dates ` +
        `throughout (row ${bad + 1}).`);
    }
    [dateCol, dayCol] = [dayCol, undefined];
  }
  const notes: string[] = [];
  let weeks = 1;
  if (dayCol !== undefined) {
    const big = values(dayCol).map((s) => dayOf(s)!).find((d) => d > 1000);
    if (big !== undefined)
      return refused('long', `Not read: ${name(dayCol)} holds ${dayWord(big)}.`);
    const key = nameOf(head[dayCol]).key;
    if (/^(weeks?|wks?)$/.test(key)) {
      weeks = 7;
      notes.push(`Weeks in ${name(dayCol)} read as days (× 7).`);
    } else if (!['day', 'days', 'study_day', 'studyday', 'd'].includes(key))
      notes.push(`Days read from ${name(dayCol)}.`);
  }
  let order: 'dmy' | 'mdy' = 'dmy';
  // a start-date column with no date in it is no start date: it stays among the other columns
  const startCol = dayCol === undefined && col.start !== undefined && values(col.start).some(isDate) ? col.start :
    undefined;
  if (dateCol !== undefined) {
    const dated = [dateCol, ...startCol === undefined ? [] : [startCol]];
    const bad = dated.find((c) => values(c).some((s) => !isDate(s)));
    if (bad !== undefined)
      return refused('long', `Not read: ${name(bad)} holds something that is not a date (2026-03-02 or 02/03/2026).`);
    const all = dated.flatMap(values);
    const found = dateOrder(all);
    if (found === 'mixed')
      return refused('long', 'Not read: the dates are written two ways, day first and month first.');
    if (found === 'unknown' && answers['dates'] !== 'dmy' && answers['dates'] !== 'mdy') {
      const ex = all.find((s) => DMY.test(s))!;
      const [a, b, y] = DMY.exec(ex)!.slice(1, 4).map(Number);
      const month = (k: number) => new Date(Date.UTC(2000, k - 1, 1)).toLocaleString('en', {month: 'long',
        timeZone: 'UTC'});
      return {...empty('long'), questions: [{id: 'dates', text: 'How are the dates written?', options: [
        ['dmy', `Day first: ${ex} is ${a} ${month(b)} ${y}`],
        ['mdy', `Month first: ${ex} is ${month(a)} ${b}, ${y}`]]}]};
    }
    order = found === 'unknown' ? answers['dates'] as 'dmy' | 'mdy' : found;
  }

  // units as named, else as answered, else assumed (mm, mm³, g)
  const kindOf = (role: Role) => KIND[role]!;
  const unitFor = (role: Role): Unit | null => named[col[role]!]?.unit ?? unitOf(answers[`unit:${kindOf(role)}`] ?? '');
  const factor = (role: Role) => UNITS[unitFor(role) ?? STANDARD[kindOf(role)]].factor;

  // each row one animal on one day, or left out with its reason
  const measured = (['volume', 'length', 'width', 'height', 'weight'] as Role[]).filter((r) => col[r] !== undefined);
  const ignored: string[] = [];
  const notNumbers: string[] = [];
  const placed = new Set<string>();
  const skipped = new Set<number>();
  interface Raw {r: number, animal: string, group: string, when: number, start: number | null,
    m: {[r in Role]?: number | null}}
  const raws: Raw[] = [];
  for (const r of rest) {
    const id = at(r, col.animal!);
    const leave = (why: string) => {
      // a number anywhere in the row, in an extra column too, is told
      if (cells[r].some((s) => /\d/.test(s))) {
        ignored.push(`row ${r + 1}, ${why}`);
        skipped.add(r);
      }
    };
    const day = dayCol !== undefined ? dayOf(at(r, dayCol)) : null;
    const when = dayCol !== undefined ? day === null ? null : day * weeks : dateOf(at(r, dateCol!), order);
    if (SUMMARY.test(id))
      leave(`a summary row (${id.replace(/[\s=]+$/, '')})`);
    else if (SUMMARY_LIKE.test(id))
      leave(`a summary row (${/^[a-z]+/i.exec(id)![0]})`);
    else if (id === '')
      leave('with no animal');
    else if (at(r, col.group!) === '')
      leave('with no group');
    else if (when === null)
      leave('with no day');
    else {
      const m: Raw['m'] = {};
      for (const role of measured) {
        const c = col[role]!;
        const x = numberOf(at(r, c), decimalComma);
        if (x !== null && Number.isNaN(x))
          notNumbers.push(`${letter(c)}${r + 1}`);
        else if (x !== null)
          placed.add(`${r},${c}`);
        m[role] = x !== null && Number.isNaN(x) ? null : x === null ? null : x * factor(role);
      }
      if (measured.some((role) => has(m[role]))) {
        const start = startCol !== undefined && at(r, startCol) !== '' ? dateOf(at(r, startCol), order) : null;
        raws.push({r, animal: id, group: at(r, col.group!), when, start, m});
      }
    }
  }
  if (notNumbers.length) {
    ignored.push(`${plural(notNumbers.length, 'cell')} that ${notNumbers.length === 1 ? 'is' : 'are'} not a plain ` +
      `number (${fewOf(notNumbers)})`);
  }
  // numbers in a column with no heading (inside the header's width or past it) are no one's measurement
  const width = Math.max(head.length, ...rest.map((r) => cells[r].length));
  Array.from({length: width}, (_, c) => c).filter((c) => head[c] === undefined || head[c] === '')
    .filter((c) => rest.some((r) => /\d/.test(at(r, c))))
    .forEach((c) => ignored.push(`column ${letter(c)}, with no heading`));

  // the check, cell by cell, that every number in a measurement column is placed, or in a row or cell left out
  const lost: string[] = [];
  for (const r of rest) {
    for (const role of measured) {
      const c = col[role]!;
      const cell = `${letter(c)}${r + 1}`;
      if (/\d/.test(at(r, c)) && !placed.has(`${r},${c}`) && !skipped.has(r) && !notNumbers.includes(cell))
        lost.push(cell);
    }
  }
  if (lost.length)
    return refused('long', `Not read: ${plural(lost.length, 'number')} could not be placed (${listed(lost)}).`);
  if (raws.length === 0)
    return refused('long', 'Not read: no measurements under the header.');

  // day 0 is each animal's own randomization date when the table has them, else the study's first measurement
  const first = Math.min(...raws.map((x) => x.when));
  const starts = new Map<string, Raw>();
  if (dateCol !== undefined && startCol !== undefined) {
    for (const x of raws.filter((y) => y.start !== null)) {
      const other = starts.get(x.animal);
      if (other && other.start !== x.start) {
        return refused('long', `Not read: rows ${other.r + 1} and ${x.r + 1} give one animal two different dates ` +
          `in ${name(startCol)}.`);
      }
      starts.set(x.animal, x);
    }
    const without = raws.find((x) => !starts.has(x.animal));
    if (without) {
      return refused('long', `Not read: ${name(startCol)} has no date for the animal in row ${without.r + 1}. ` +
        'Give each animal its date, or leave the column out to count the days from the first measurement.');
    }
    const one = new Set([...starts.values()].map((x) => x.start!));
    notes.push(one.size === 1 ? `Day 0 is the date in ${name(startCol)}, ${isoOf([...one][0])}.` :
      `Day 0 is each animal's date in ${name(startCol)}.`);
  } else if (dateCol !== undefined)
    notes.push(`Day 0 is the first measurement, ${isoOf(first)}.`);
  const dayOfRaw = (x: Raw) => dateCol === undefined ? x.when : x.when - (startCol !== undefined ?
    starts.get(x.animal)!.start! : first);
  if (dateCol !== undefined && startCol !== undefined && raws.some((x) => dayOfRaw(x) < 0))
    notes.push('Measurements before day 0 are at negative days.');

  // one row per animal and day; an id found in two groups is two animals
  const groupsOf = new Map<string, Set<string>>();
  raws.forEach((x) => groupsOf.set(x.animal, (groupsOf.get(x.animal) ?? new Set()).add(x.group)));
  const shared = new Set([...groupsOf].filter(([, g]) => g.size > 1).map(([id]) => id));
  // only the shared ids are told apart, by the order their group first appears (never by its name: a blinded study)
  const groupNo = new Map<string, number>();
  raws.forEach((x) => groupNo.has(x.group) || groupNo.set(x.group, groupNo.size + 1));
  const animal = (x: Raw) => shared.has(x.animal) ? `${x.animal} (g${groupNo.get(x.group)})` : x.animal;
  if (shared.size)
    notes.push(`Animal ids repeat across groups: read as different animals (${plural(shared.size, 'id')}).`);
  const seen = new Map<string, number>();
  for (const x of raws) {
    const key = `${animal(x)}\u0000${dayOfRaw(x)}`;
    if (seen.has(key)) {
      return refused('long', `Not read: rows ${seen.get(key)! + 1} and ${x.r + 1} are one animal on one day ` +
        `(${x.animal}, day ${dayOfRaw(x)}); keep one row per animal and day.`);
    }
    seen.set(key, x.r);
  }

  // volumes as given, computed from the caliper, or given and found to follow a formula
  const height = col.height !== undefined;
  let volume: VolumeSource;
  if (col.volume !== undefined && caliper) {
    const formula = deduceFormula(raws.flatMap((x) => has(x.m.length) && has(x.m.width) && has(x.m.volume) ?
      [{l: x.m.length, w: x.m.width, h: x.m.height ?? null, v: x.m.volume}] : []));
    volume = {formula, source: formula ? 'deduced' : 'given'};
    const tested = raws.filter((x) => has(x.m.length) && has(x.m.width) && has(x.m.volume) && x.m.volume > 0).length;
    notes.unshift(formula ? `Volumes match ${FORMULAS[formula].label}.` : tested < 5 ? 'Volumes as given.' :
      'The volumes given do not follow a known formula from length and width: they are used as given.');
  } else if (col.volume !== undefined) {
    volume = {formula: null, source: 'given'};
    notes.unshift('Volumes as given.');
  } else {
    const offered = (Object.keys(FORMULAS) as VolumeFormula[]).filter((f) => height || !FORMULAS[f].height);
    const formula = offered.find((f) => f === answers['formula']) ?? 'LW2/2';
    volume = {formula, source: 'computed'};
    notes.unshift(`Volumes computed as ${FORMULAS[formula].label}.`);
    questions.push({id: 'formula', text: 'Volume from length and width', preset: formula,
      options: offered.map((f) => [f, FORMULAS[f].label])});
  }
  // a volume missing where the lab's formula is known is computed with it
  const volumeOf = (x: Raw): number | null => {
    if (volume.source !== 'computed' && has(x.m.volume))
      return x.m.volume;
    const f = volume.formula ? FORMULAS[volume.formula] : null;
    return f && has(x.m.length) && has(x.m.width) && (!f.height || has(x.m.height)) ?
      f.v(x.m.length, x.m.width, x.m.height ?? 0) : null;
  };
  const rows: StudyRow[] = raws.map((x) => ({animal: animal(x), group: x.group, day: dayOfRaw(x), volume: volumeOf(x),
    weight: x.m.weight ?? null}));
  if (volume.source === 'deduced') {
    const filled = raws.filter((x) => !has(x.m.volume) && has(volumeOf(x))).length;
    if (filled)
      notes[0] = `${notes[0].slice(0, -1)}; ${plural(filled, 'missing volume')} computed with it.`;
  }
  // length and width but no volume, and no formula to compute it: left out, said
  if (col.volume !== undefined && caliper && !volume.formula) {
    const no = raws.filter((x) => !has(x.m.volume) && has(x.m.length) && has(x.m.width)).map((x) => x.r + 1);
    const why = 'length and width but no volume, and no formula matches the volumes given';
    if (no.length > 3)
      ignored.push(`${no.length} rows (${listed(no.map(String))}), ${why}`);
    else
      no.forEach((r) => ignored.push(`row ${r}, ${why}`));
  }

  // the other columns, kept as they are (the caliper's readings, a date, a cage, a comment)
  const taken = [col.animal, col.group, dayCol, col.volume, col.weight, ...emptyDay ? [col.day] : []];
  const canonical = ['animal', 'group', 'day', 'volume', 'body_weight'];
  const extraRaw: TableRead['extra'] = head.flatMap((title, c) => {
    if (taken.includes(c) || title === '')
      return [];
    const texts = raws.map((x) => at(x.r, c));
    const nums = texts.map((s) => numberOf(s, decimalComma));
    const numeric = texts.some((s) => s !== '') && nums.every((x) => x === null || !Number.isNaN(x));
    return [[canonical.includes(title.toLowerCase()) ? `${title} (as given)` : title,
      numeric ? nums : texts.map((s) => s === '' ? null : s)]];
  });
  // two columns that read alike (Datagrok's names ignore case): the later one numbered, once, here
  const used = new Set(['animal', 'group', 'day', 'volume',
    ...rows.some((x) => x.weight !== null) ? ['body_weight'] : []]);
  const extra = extraRaw.map(([title, xs]): TableRead['extra'][number] => {
    let named = title;
    for (let k = 2; used.has(named.toLowerCase()); k++)
      named = `${title} (${k})`;
    used.add(named.toLowerCase());
    return [named, xs];
  });

  // the units the volumes come from are the lengths' when computed, the volumes' when given
  const roles = measured.filter((r) => r !== 'volume' || volume.source !== 'computed')
    .filter((r) => KIND[r] !== 'length' || volume.source !== 'given');
  return withUnits({...empty('long'), rows, notes, ignored, extra, volume, questions},
    roles.map((r) => [kindOf(r), unitFor(r)]), volume.source === 'computed' ? 'length' : 'volume');
}

/** The plain-words account of a read, as Open a study shows it: animals per group and the days, how the volumes,
 *  units and day 0 were read, the columns kept, then what was left out and why. Blinded, the counts alone, in the
 *  table's order. */
export function studySummary(read: TableRead, blinded = false): string[] {
  if (read.error)
    return [read.error];
  if (read.rows.length === 0)
    return [];
  const counts = read.groups.map(([, n]) => String(n));
  const groups = blinded ? `${andList(counts)} animals` :
    read.groups.map(([g, n]) => `${g} ${plural(n, 'animal')}`).join(', ');
  const days = read.rows.map((r) => r.day);
  return [`Read ${plural(read.groups.length, 'group')}: ${groups}, days ${Math.min(...days)}–${Math.max(...days)}.`,
    ...read.notes, ...read.extra.length ? [`Also kept: ${read.extra.map(([n]) => n).join(', ')}.`] : [],
    ...read.ignored.map((l) => `Not read: ${l}.`)];
}

/** Whether the read waits on an answer: a question with no choice made yet. */
export const waiting = (read: TableRead): boolean => read.questions.some((q) => q.preset === undefined);

/** A volume source read from a file someone may have edited: as it was saved, or undefined. */
export const volumeSourceOf = (x: any): VolumeSource | undefined =>
  ['given', 'computed', 'deduced'].includes(x?.source) && (x.formula === null || x.formula in FORMULAS) ?
    {formula: x.formula, source: x.source} : undefined;

// ---- an Excel sheet as cells ---------------------------------------------------------------------------------------

/** What a sheet's frame needs to be read as cells: DG.DataFrame fits it. */
export interface SheetFrame {name: string, rowCount: number,
  columns: {length: number, byIndex(i: number): {name: string, type?: string, getString(i: number): string,
    get?(i: number): any, isNone?(i: number): boolean}}}

// an Excel date as the sheet reader hands it over (Date.toString: Mon Nov 13 2026 00:00:00 GMT+0100 (...)), as written
const JS_DATE = /^[A-Z][a-z]{2} [A-Z][a-z]{2} \d{2} \d{4} \d{2}:\d{2}:\d{2} GMT[+-]\d{4}(?: \(.*\))?$/;
// Excel dates are UTC midnight: printing them in local time would move the day west of UTC
const dayText = (s: string): string => JS_DATE.test(s) ? new Date(s).toISOString().slice(0, 10) : s;

/** A sheet's cells as text: the column names, then a row per line; a date cell as yyyy-mm-dd, a blank as ''.
 *  A date column (PowerPack's reader makes one of the cells it recognises as dates) is read from its stored UTC
 *  instant; every other cell as written (getString: PowerPack sets each number column's format to its decimals). */
export function cellsOfFrame(df: SheetFrame): string[][] {
  const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
  const rows = Array.from({length: df.rowCount}, (_, r) => cols.map((c) => c.type === 'datetime' && c.get ?
    (c.isNone?.(r) ? '' : c.get(r).toISOString().slice(0, 10)) : dayText((c.getString(r) ?? '').trim())));
  return [cols.map((c) => dayText(c.name.trim())), ...rows];
}

/** Why a sheet cannot be read, or null: PowerPack sizes the columns by the sheet's first row, so a sheet that does
 *  not start with its column names loses cells. */
export function sheetProblem(df: SheetFrame): string | null {
  const names = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
  // PowerPack's own names for what the first row left blank: colN at position N, and "Name (1)" beside "Name"
  const invented = (n: string, i: number) => n === `col${i + 1}` || (/ \(\d+\)$/.test(n) &&
    names.includes(n.replace(/ \(\d+\)$/, '')));
  if (names.length < 2 || names.some(invented)) {
    return `Sheet «${df.name}»: the first row must hold the column names (Animal, Day, …). Delete the rows above ` +
      'them, start the table in row 1, or save the sheet as CSV.';
  }
  return null;
}

/** The sheet's animal column read as decimal numbers, which cannot tell 1.1 from 1.10: a refusal, or null. */
export function animalProblem(df: SheetFrame): string | null {
  for (let i = 0; i < df.columns.length; i++) {
    const c = df.columns.byIndex(i);
    if (c.type === 'double' && roleOf(c.name.trim())?.role === 'animal') {
      return `The animal IDs in «${df.name}» were read as decimal numbers, so «1.1» and «1.10» can't be told apart. ` +
        'Format the animal column as text in Excel, or save the sheet as CSV.';
    }
  }
  return null;
}

/** Cells as CSV text, as a file's would be; a cell with a separator, quote or line break is quoted. */
export const csvOf = (cells: string[][]): string =>
  cells.map((r) => r.map((c) => /[",;\t\n\r]/.test(c) ? `"${c.replace(/"/g, '""')}"` : c).join(',')).join('\n');

/** The sheets (with data rows) that read as a study. */
export function studySheets(sheets: {name: string, cells: string[][]}[]): string[] {
  return sheets.filter((s) => {
    try {
      const r = readStudy(s.cells);
      return r.rows.length > 0 || r.questions.length > 0;
    } catch {
      return false;
    }
  }).map((s) => s.name);
}
