/* The study as a modelling dataset, in NONMEM's layout: one row per event, tumour volumes in compartment 1, the doses
   from the study's schedule and NCA's concentrations of the study's own animals in each drug's compartment (2 for the
   first drug, 3 for the second...: a combination's drugs dosed each into its own). A copy of the data for
   a modeller: nothing here changes the analysis or the study. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import type {StudyInputs} from './explore';
import {download} from './report';
import {EXCLUDED, FLAGGED, ticked} from './study-file';
import {dosedDays, excepted, readSchedule} from './schedule';
import type {Schedule} from './schedule';
import {StudyMeta, armParts} from './study-meta';
import {drugKey} from './arm-name';
import {okNamed, showCentred} from './ui';
import {plural} from './wide-table';

/** The concentrations NCA read: its table and mapping (subject, time, concentration, blq, lloq and their units), and
 *  the column that names the group, if any. */
export interface PkSource {table: DG.DataFrame, mapping: any, groupCol: string | null}

/** What the dialog lets the user set: each group's first dosing day (default 0; a first day its schedule writes wins),
 *  a dose in mg/kg its details and name do not give (by line label: the group, or "group, drug" for a combination's
 *  drug; a written one wins), and the day of the dose the concentrations' hours count from (default the first
 *  dosing day). */
export interface DatasetOptions {firstDose?: Map<string, number>, dose?: Map<string, number>, pkDoseDay?: number}

/** How a group (a combination: each of its drugs) goes into the file. */
export type LineStatus = 'control' | 'dosed' | 'no dose' | 'no schedule' | 'unread';
export interface DatasetLine {group: string, label: string, dose: number | null, schedule: string,
  status: LineStatus,
  /** a dose per animal (100 µg) with no mg/kg: it over the group's mean body weight at its first dose, and whence
   *  (what the dialog offers in its box) */
  proposed?: {dose: number, from: string}}

// a dose per animal in mg, by its unit (arm-name.ts writes µ as u)
const MG_PER: {[unit: string]: number} = {g: 1000, mg: 1, ug: 1e-3, ng: 1e-6};

export const MODEL_COLUMNS = ['ID', 'TIME', 'DV', 'MDV', 'EVID', 'CMT', 'AMT', 'GROUP', 'DOSE', 'COMBO', 'BW', 'FLAG'];
/** A missing covariate: NONMEM reads "." as 0, so a missing DOSE or BW is written -99. */
export const MISSING = -99;
type Cell = number | '.';
export type ModelRow = {[c: string]: Cell};

export interface ModellingDataset {
  rows: ModelRow[];
  /** each animal's number, as written in the study, and its group's number and name */
  ids: {id: number, animal: string, group: number, groupName: string}[];
  /** what the file holds beyond its volumes: each group's doses, the concentrations read */
  held: string[];
  /** what it lacks, and why */
  notes: string[];
  /** whether any dosing record was written */
  dosed: boolean;
  /** the groups with dosing records */
  dosedGroups: string[];
  /** each drug's compartment (2, 3...), in the order the groups first give them */
  compartments: {drug: string, cmt: number}[];
  /** each group (a combination: each of its drugs), in order, and how it goes in */
  lines: DatasetLine[];
}

const HOURS_PER: {[u: string]: number} = {h: 1, hr: 1, hrs: 1, hour: 1, hours: 1, min: 1 / 60, d: 24, day: 24,
  days: 24};
const keyOf = (group: string, animal: string) => `${group}\u0000${animal}`;

/** The schedule of drug `i` of a group (a combination's drugs each have their own) as its details write it, else
 *  its name, and as read. */
export function scheduleOfGroup(g: string, meta: StudyMeta | null, i = 0): {text: string, read: Schedule | null} {
  const parts = armParts(g, meta);
  const p = parts[i];
  // the details' row of that drug as typed (a route the name reader does not know, "drinking water", kept)
  const rows = meta?.groups.filter((x) => x.group === g) ?? [];
  const row = parts.length > 1 ? rows.find((x) => x.compound.trim() && drugKey(x.compound) === p?.key) ?? rows[i] :
    rows[0];
  const text = row?.schedule?.trim() || (p?.schedule?.trim() ?? '');
  return {text, read: readSchedule(text, row?.route || (p?.route ?? ''))};
}

/** The study as NONMEM rows: ID (the animals, by group and name, numbered in the order they first appear), TIME
 *  (days from randomization), DV, MDV, EVID, CMT, AMT, and the covariates GROUP (the groups numbered in order), DOSE
 *  (the group's dose in mg/kg as its name gives it; 0 for the control, -99 when unknown or a combination: its doses
 *  are in AMT), COMBO (1 for a combination) and BW (-99 when unknown), with FLAG 1 for a flagged measurement.
 *  Excluded measurements are left out. Each drug has its own compartment (2, 3... in the order the groups first give
 *  them; a combination's drugs each dosed into theirs). Doses come from each drug's schedule in the study's details
 *  (else the group's name) as schedule.ts reads it
 *  (QD, QD×14, Q3D×4, BIW, from day 3...; BID's second dose half a day later), under the study's rules (no dosing
 *  at weekends), from the group's first dosing day, up to each animal's own last measurement day (excluded ones
 *  included); AMT is the dose in mg/kg. At one time the observation comes before the dose (measured before dosing). */
export function modellingDataset(s: StudyInputs, meta: StudyMeta | null, pk: PkSource | null = null,
  opts: DatasetOptions = {}): ModellingDataset {
  const df = s.study;
  const excluded = df.col(EXCLUDED);
  const flagged = df.col(FLAGGED);
  const numAt = (c: string | null, i: number): number | null => {
    const col = c ? df.col(c) : null;
    if (!col || col.isNone(i))
      return null;
    const v = Number(col.get(i));
    return Number.isFinite(v) ? v : null;
  };
  const ids = new Map<string, number>();
  const groups = new Map<string, number>();
  const animals = new Map<string, {animal: string, group: string, last: number}>();
  // each animal's body weight by day, for its doses
  const weights = new Map<string, [number, number][]>();
  const rows: ModelRow[] = [];
  // a drug's dose in mg/kg as its details (else its name) give it, else as typed in the dialog
  const mgPerKg = (p: ReturnType<typeof armParts>[number], label: string): number =>
    p.dose && p.dose.unit.toLowerCase() === 'mg/kg' && p.dose.value > 0 ? p.dose.value :
      (opts.dose?.get(label) ?? 0) > 0 ? opts.dose!.get(label)! : MISSING;
  // the group's dose: 0 for a group with no drug, unknown for a combination
  const doseOf = (g: string): number => {
    const parts = armParts(g, meta);
    return parts.length === 0 ? 0 : parts.length === 1 ? mgPerKg(parts[0], g) : MISSING;
  };
  const covariates = (key: string) => {
    const g = animals.get(key)!.group;
    return {GROUP: groups.get(g)!, DOSE: doseOf(g), COMBO: armParts(g, meta).length > 1 ? 1 : 0};
  };
  for (let i = 0; i < df.rowCount; i++) {
    const animal = String(df.get(s.animal, i) ?? '').trim();
    const group = String(df.get(s.group, i) ?? '').trim();
    const day = numAt(s.day, i);
    if (!animal || day === null)
      continue;
    const key = keyOf(group, animal);
    if (!ids.has(key)) {
      ids.set(key, ids.size + 1);
      animals.set(key, {animal, group, last: day});
    }
    if (!groups.has(group))
      groups.set(group, groups.size + 1);
    // its last day, excluded rows counted: its doses run to it
    animals.get(key)!.last = Math.max(animals.get(key)!.last, day);
    if (ticked(excluded, i))
      continue;
    const bw = numAt(s.bodyWeight, i);
    if (bw !== null)
      weights.set(key, [...weights.get(key) ?? [], [day, bw]]);
    const v = numAt(s.volume, i);
    if (v === null)
      continue;
    rows.push({ID: ids.get(key)!, TIME: day, DV: v, MDV: 0, EVID: 0, CMT: 1, AMT: '.', ...covariates(key),
      BW: bw ?? MISSING, FLAG: ticked(flagged, i) ? 1 : 0});
  }
  weights.forEach((w) => w.sort((a, b) => a[0] - b[0]));
  // the animal's body weight on a day: that day's, else the last before it, else its first
  const weightOn = (key: string, day: number): number => {
    const w = weights.get(key) ?? [];
    let out = w.length ? w[0][1] : MISSING;
    for (const [d, bw] of w) {
      if (d <= day)
        out = bw;
    }
    return out;
  };
  // a dose per animal (100 µg) as mg/kg: over the group's mean body weight (g) on its first dosing day; none
  // without body weights
  const proposedDose = (p: ReturnType<typeof armParts>[number], g: string, day: number):
    {proposed?: DatasetLine['proposed']} => {
    const perMg = p.dose && p.dose.value > 0 ? MG_PER[p.dose.unit.toLowerCase().replace(/[µμ]/g, 'u')] : undefined;
    const bws = [...animals].filter(([, a]) => a.group === g).map(([k]) => weightOn(k, day)).filter((w) => w > 0);
    if (!perMg || !bws.length)
      return {};
    const mean = bws.reduce((a, b) => a + b, 0) / bws.length;
    return {proposed: {dose: Math.round(100 * p.dose!.value * perMg / (mean / 1000)) / 100,
      from: `from ${p.dose!.value} ${p.dose!.unit.replace(/^u/, 'µ')} at ${mean.toFixed(1)} g mean weight`}};
  };
  // the doses, group by group
  const notes: string[] = [];
  const held: string[] = [];
  const undosed: string[] = [];
  const dosedGroups: string[] = [];
  const lines: DatasetLine[] = [];
  const lastDay = Math.max(0, ...[...animals.values()].map((a) => a.last));
  // each drug's doses, group by group, then a compartment for each drug dosed (2, 3... in the order the groups give
  // them; a drug never dosed, such as radiation in Gy, has none)
  const plans: {g: string, label: string, key: string, amt: number, days: number[], schedule: string}[] = [];
  for (const g of groups.keys()) {
    const parts = armParts(g, meta);
    if (parts.length === 0)
      lines.push({group: g, label: g, dose: 0, schedule: '', status: 'control'});
    parts.forEach((p, i) => {
      const label = parts.length > 1 ? `${g}, ${p.drug}` : g;
      const amt = mgPerKg(p, label);
      const {text: schedule, read} = scheduleOfGroup(g, meta, i);
      // one first dosing day: the schedule's when it writes one ("from day 3"), else the dialog's (default 0)
      // ponytail: "from day 0" written reads as none written; the same day unless the dialog moves it
      const start = read && read.first > 0 ? read.first : opts.firstDose?.get(g) ?? 0;
      // the schedule on that day; continuous exposure has no doses to write
      const moved = read && read.kind !== 'continuous' ? {...read, first: start,
        end: read.end === null ? null : read.end - read.first + start} : null;
      // a day's doses spread over it (BID: day and day + 0.5), as doseCount counts them; TIME to 4 decimals; the
      // study's exceptions for the group (a skipped day out, an added one in: every drug of the arm), and with no
      // schedule read no doses, whatever exceptions it kept
      const days = (moved ? excepted(dosedDays(moved, lastDay, meta), g, meta, lastDay) : []).flatMap((d) =>
        Array.from({length: moved!.perDay}, (_, k) => Math.round((d + k / moved!.perDay) * 1e4) / 1e4));
      lines.push({group: g, label, dose: amt === MISSING ? null : amt, schedule, status: amt === MISSING ?
        'no dose' : days.length ? 'dosed' : schedule ? 'unread' : 'no schedule',
      ...amt === MISSING ? proposedDose(p, g, start) : {}});
      if (amt === MISSING || days.length === 0) {
        undosed.push(`${label} (${amt === MISSING ? 'no dose in mg/kg in its details or name' : schedule ?
          `its schedule, ${schedule}, gives no dosing days (write it as QD, Q3D×4 or QW×3)` :
          'no schedule in the study details or its name'})`);
        return;
      }
      plans.push({g, label, key: p.key, amt, days, schedule});
    });
  }
  const cmts = new Map<string, {drug: string, cmt: number}>();
  for (const pl of plans) {
    if (!cmts.has(pl.key))
      cmts.set(pl.key, {drug: armParts(pl.g, meta).find((p) => p.key === pl.key)!.drug, cmt: cmts.size + 2});
  }
  const several = cmts.size > 1;
  for (const {g, label, key: drug, amt, days, schedule} of plans) {
    const cmt = cmts.get(drug)!.cmt;
    if (!dosedGroups.includes(g))
      dosedGroups.push(g);
    held.push(`${label}: ${amt} mg/kg on day${days.length === 1 ? '' : 's'} ${days.join(', ')} (${schedule})` +
      `${several ? `, compartment ${cmt}` : ''}, each animal up to its last measurement.`);
    for (const [key, a] of animals) {
      if (a.group !== g)
        continue;
      // the dose's day against the animal's last (BID's second dose on that day kept)
      for (const day of days.filter((d) => Math.floor(d) <= a.last)) {
        rows.push({ID: ids.get(key)!, TIME: day, DV: '.', MDV: 1, EVID: 1, CMT: cmt, AMT: amt, ...covariates(key),
          BW: weightOn(key, day), FLAG: 0});
      }
    }
  }
  if (undosed.length)
    notes.push(`No dosing records for ${undosed.join('; ')}: their observations only.`);
  // NCA's concentrations of the study's own animals (by group and animal; the group NCA names, else the animal's
  // only group in the study)
  if (pk) {
    const m = pk.mapping ?? {};
    const t = pk.table;
    const per = HOURS_PER[String(m.units?.time ?? 'h').toLowerCase()];
    const [subj, time, conc, blq, lloq] = [m.subject, m.time, m.concentration, m.blq, m.lloq]
      .map((c) => c ? t.col(c) : null);
    const doseDay = opts.pkDoseDay ?? 0;
    // a concentration's compartment: its drug's, as the table names it (a column NCA groups by), else the drug of
    // the animal's group when it has one; the only one when the study gives one drug
    const drugCols = (m.groupBy ?? []).filter((c: string) => c !== pk.groupCol).map((c: string) => t.col(c))
      .filter((c: DG.Column | null): c is DG.Column => !!c);
    const cmtOf = (i: number, key: string): number | null => {
      if (!several)
        return 2;
      for (const c of drugCols) {
        const hit = cmts.get(drugKey(String(c.get(i) ?? '')));
        if (hit)
          return hit.cmt;
      }
      const own = armParts(animals.get(key)!.group, meta);
      return own.length === 1 ? cmts.get(own[0].key)?.cmt ?? null : null;
    };
    let unplaced = 0;
    const byName = new Map<string, string[]>();
    animals.forEach((a, key) => byName.set(a.animal, [...byName.get(a.animal) ?? [], key]));
    if (!subj || !time || !conc)
      notes.push(`The concentrations in ${t.name} could not be read: no PK rows.`);
    else if (per === undefined)
      notes.push(`The concentrations' time unit (${m.units?.time}) is not read: no PK rows.`);
    else {
      let n = 0;
      let low = 0;
      const seen = new Set<string>();
      for (let i = 0; i < t.rowCount; i++) {
        const animal = String(subj.get(i) ?? '').trim();
        const keys = pk.groupCol ? [keyOf(String(t.get(pk.groupCol, i) ?? '').trim(), animal)] :
          byName.get(animal) ?? [];
        // without a group column, an animal name two groups share cannot be placed
        const key = keys.length === 1 && ids.has(keys[0]) ? keys[0] : null;
        if (key === null || time.isNone(i) || conc.isNone(i))
          continue;
        const h = Number(time.get(i));
        const c = Number(conc.get(i));
        if (!Number.isFinite(h) || !Number.isFinite(c))
          continue;
        const limit = lloq && !lloq.isNone(i) ? Number(lloq.get(i)) : NaN;
        if (c <= 0 || ticked(blq, i) || c < limit) {
          low++;
          continue;
        }
        const cmt = cmtOf(i, key);
        if (cmt === null) {
          unplaced++;
          continue;
        }
        const day = doseDay + h * per / 24;
        rows.push({ID: ids.get(key)!, TIME: day, DV: c, MDV: 0, EVID: 0, CMT: cmt, AMT: '.', ...covariates(key),
          BW: weightOn(key, day), FLAG: 0});
        seen.add(key);
        n++;
      }
      (n ? held : notes).push(n ? `${n} concentrations of ${plural(seen.size, 'animal')} from ${t.name} in ` +
        `${several ? 'their drug\'s compartment' : 'compartment 2'} (${m.units?.concentration ?? 'their unit'}), ` +
        `their hours counted from the dose of day ${doseDay}.` :
        `No animal of ${t.name} is one of the study's (PK animals apart, as satellites): no PK rows.`);
      if (unplaced) {
        notes.push(`${plural(unplaced, 'concentration')} of animals given several drugs left out: ` +
          `${t.name} does not say which drug each is.`);
      }
      if (low) {
        notes.push(`${plural(low, 'concentration')} at or below 0 or below the limit of quantification ` +
          'left out.');
      }
    }
  }
  // by animal and time; at one time the observations first, the tumour before the concentration
  rows.sort((a, b) => (a.ID as number) - (b.ID as number) || (a.TIME as number) - (b.TIME as number) ||
    (a.EVID as number) - (b.EVID as number) || (a.CMT as number) - (b.CMT as number));
  return {rows, dosed: dosedGroups.length > 0, dosedGroups, held, notes, compartments: [...cmts.values()], lines,
    ids: [...animals].map(([key, a]) =>
      ({id: ids.get(key)!, animal: a.animal, group: groups.get(a.group)!, groupName: a.group}))};
}

const csvCell = (x: string) => /[",\n\r]/.test(x) ? `"${x.replace(/"/g, '""')}"` : x;
// 7 significant digits: what a 32-bit Datagrok column holds (0.3 comes back as 0.30000001192…)
const numText = (x: Cell) => x === '.' ? '.' : String(Number(x.toPrecision(7)));

/** The column dictionary, as comment lines at the top of the file: each starts with @, so IGNORE=@ skips them with the
 *  header line, and the dialog need not carry it. */
export const CSV_DICTIONARY = ['@ Modelling dataset, one row per event (NONMEM layout). Read with IGNORE=@.',
  '@ TIME: days from randomization. CMT 1: tumor volume (mm3); CMT 2, 3...: each drug dose and concentration.',
  '@ EVID 0: observation (DV); EVID 1: dose (AMT, mg/kg). MDV 1: DV not read. "." only in DV and AMT where' +
    ' MDV and EVID say so.',
  '@ -99: missing dose or body weight (NONMEM reads "." as 0). COMBO 1: a combination group (DOSE -99). FLAG 1: a ' +
    'flagged measurement, kept.',
  '@ ID is the animal number: the second file (IDs) gives the animal and group as written in the study.'];

/** The dataset as CSV: the dictionary lines, then a header line (both skipped by IGNORE=@ in NONMEM); "." only in DV
 *  and AMT where MDV and EVID say they are not read, -99 for a missing covariate. */
export function modellingCsv(d: ModellingDataset): string {
  return [...CSV_DICTIONARY, MODEL_COLUMNS.join(','),
    ...d.rows.map((r) => MODEL_COLUMNS.map((c) => numText(r[c])).join(','))].join('\n') + '\n';
}

/** The animals' numbers as CSV: ID, the animal as written in the study, GROUP and the group's name. */
export function idMapCsv(d: ModellingDataset): string {
  return ['ID,ANIMAL,GROUP,GROUP_NAME', ...d.ids.map((x) => [String(x.id), csvCell(x.animal), String(x.group),
    csvCell(x.groupName)].join(','))].join('\n') + '\n';
}

/** What the dialog says when no concentrations go into the dataset, and how to add them; `recorded`: the table NCA
 *  read when the study kept its exposure (not open now). */
export const noConcentrations = (recorded: string | null): string => 'No concentrations are included: ' +
  (recorded ? `${recorded}, the table the study's exposure was computed from, is not open. Open it, or pick a ` +
    'concentration table in Exposure (NCA), ' : 'pick a concentration table in Exposure (NCA), ') +
  'then export again.';

/** A group's status from its lines (a combination: one per drug): the first thing missing, in the order the dialog
 *  asks for it (a dose, then a schedule that reads), else dosed. */
export function groupStatus(lines: DatasetLine[]): LineStatus {
  const order: LineStatus[] = ['control', 'no dose', 'unread', 'no schedule', 'dosed'];
  return order.find((st) => lines.some((l) => l.status === st)) ?? 'dosed';
}

/** A group's status as the dialog shows it, and its detail on hover (a combination's, drug by drug). */
function statusOf(lines: DatasetLine[], d: ModellingDataset): [string, string] {
  const tip = (l: DatasetLine) => {
    switch (l.status) {
    case 'control': return 'No drug: its volumes only';
    case 'dosed': return d.held.find((h) => h.startsWith(`${l.label}:`)) ?? 'Dosed from its schedule';
    case 'no dose': return 'No dose in mg/kg in its details or name: type it to add its doses';
    case 'no schedule': return 'No schedule in the study details or its name (QD, Q3D×4, QW×3...): its volumes only';
    default: return `"${l.schedule}" gives no dosing days: write it as QD, Q3D×4 or QW×3 in the study details`;
    }
  };
  const text = {'control': 'control', 'dosed': '✓', 'no dose': '⚠ no dose in mg/kg',
    'no schedule': '⚠ no schedule', 'unread': '⚠ schedule not read'}[groupStatus(lines)];
  return [text, lines.map((l) => lines.length > 1 ? `${l.label.slice(l.group.length + 2)}: ${tip(l)}` : tip(l))
    .join('\n')];
}

/** What the dataset holds, group by group (its dose, first dosing day and whether it is dosed; a missing dose typed
 *  there), and the concentrations, then its download: the dataset and the animals' numbers, two CSV files.
 *  `recorded`: the concentration table the study kept its exposure from (said when it is not used). */
export function modellingDatasetDialog(s: StudyInputs, meta: StudyMeta | null, pk: PkSource | null,
  recorded: string | null = null): void {
  const name = s.study.name || 'study';
  const lines = modellingDataset(s, meta, null).lines;
  const narrow = (x: DG.InputBase<number | null>) => {
    x.input.style.width = '64px';
    return x;
  };
  // a dose typed where the group's details and name give none; a dose per animal offered over the group's weight
  const doseIn = new Map(lines.filter((l) => l.status === 'no dose').map((l) => [l.label,
    narrow(ui.input.float('', {value: l.proposed?.dose, nullable: true, min: 0,
      tooltipText: 'mg/kg per dose'}))]));
  // each group with a drug: its first dosing day as its schedules write it (shown; set in the study details), else
  // typed (day 0 by default)
  const written = new Map<string, number>();
  const starts = new Map<string, DG.InputBase<number | null>>();
  for (const g of new Set(lines.filter((l) => l.status === 'dosed' || l.status === 'no dose').map((l) => l.group))) {
    const firsts = armParts(g, meta).map((_, i) => scheduleOfGroup(g, meta, i).read?.first ?? 0);
    if (firsts.every((f) => f > 0))
      written.set(g, Math.min(...firsts));
    else {
      starts.set(g, narrow(ui.input.int('', {value: 0, min: 0,
        tooltipText: 'The study day of this group\'s first dose (0: randomization)'})));
    }
  }
  const pkDay = pk ? ui.input.int('Hours counted from the dose of day', {value: 0, min: 0,
    tooltipText: 'The study day of the dose the concentration table\'s times count from'}) : null;
  let pkDayTouched = false;
  let setting = false;
  const optsOf = (): DatasetOptions => ({
    firstDose: new Map([...starts].map(([g, x]) => [g, x.value ?? 0])),
    dose: new Map([...doseIn].map(([l, x]) => [l, x.value ?? 0])),
    pkDoseDay: pkDay?.value ?? 0});
  // one row per group, a combination's drugs' doses side by side; each group's status, its detail on hover (bound
  // once, read as it changes)
  const groupNames = [...new Set(lines.map((l) => l.group))];
  const tips = new Map(groupNames.map((g) => [g, '']));
  const status = new Map(groupNames.map((g) => [g, ui.tooltip.bind(ui.span([], {style: {whiteSpace: 'nowrap'}}),
    () => tips.get(g)!)]));
  const doseCell = (l: DatasetLine) => {
    const x = doseIn.get(l.label)?.input;
    const e = x ?? ui.span([l.dose === null || l.status === 'control' ? '—' : String(l.dose)]);
    // a combination's drug by name; a dose offered, whence ("from 100 µg at 20.1 g mean weight")
    const tip = [l.label !== l.group ? l.label.slice(l.group.length + 2) : '', x && l.proposed ?
      `mg/kg ${l.proposed.from}; change it if needed` : ''].filter((t) => t).join(': ');
    if (tip)
      ui.tooltip.bind(e, tip);
    return e;
  };
  const table = ui.table(groupNames, (g) => {
    const own = lines.filter((l) => l.group === g);
    // a combination: one line per drug, its name first ("1 + 10" did not say which dose is whose)
    const doses = own.length < 2 ? doseCell(own[0]) : ui.divV(own.map((l) => ui.divH([
      ui.span([l.label.slice(l.group.length + 2)]), doseCell(l)], {style: {alignItems: 'center', gap: '4px'}})));
    const day = written.get(g);
    const fixed = ui.span([day === undefined ? '—' : String(day)]);
    if (day !== undefined)
      ui.tooltip.bind(fixed, 'From its schedule; change it in the study details');
    return [g, doses, starts.get(g)?.input ?? fixed, status.get(g)!];
  }, ['Group', 'Dose (mg/kg)', 'First dose (day)', '']);
  const footer = ui.divText('');
  let concTip = '';
  const conc = ui.tooltip.bind(ui.divText(''), () => concTip);
  const pkNotes = ui.divV([]);
  let d = modellingDataset(s, meta, pk, optsOf());
  const refresh = () => {
    // the concentrations' dose day follows the first dosing day until it is set
    if (pkDay && !pkDayTouched) {
      setting = true;
      const days = [...written.values(), ...[...starts.values()].map((x) => x.value ?? 0)];
      pkDay.value = days.length ? Math.min(...days) : 0;
      setting = false;
    }
    d = modellingDataset(s, meta, pk, optsOf());
    for (const g of groupNames) {
      const [text, tip] = statusOf(d.lines.filter((l) => l.group === g), d);
      const e = status.get(g)!;
      e.textContent = text;
      e.style.color = text.startsWith('⚠') ? 'var(--orange-2, #b26a00)' : text === '✓' ?
        'var(--green-2, #2e7d32)' : 'var(--grey-4, #8a959f)';
      tips.set(g, tip);
    }
    footer.textContent = `${d.ids.length} animals · ${d.rows.filter((r) => r.CMT === 1).length} volumes`;
    const nConc = d.rows.filter((r) => r.EVID === 0 && r.CMT !== 1).length;
    conc.textContent = `Concentrations: ${pk ? `${nConc} from ${pk.table.name}` : 'none'}`;
    concTip = pk ? d.held.filter((h) => !lines.some((l) => h.startsWith(`${l.label}:`))).join('\n') ||
      'None of its animals is one of the study\'s' : noConcentrations(recorded);
    // what else the file lacks (the concentrations'), beside them
    pkNotes.replaceChildren(...d.notes.filter((n) => !n.startsWith('No dosing records')).map((n) => ui.divText(`⚠ ${n}`,
      {style: {color: 'var(--orange-2, #b26a00)'}})));
  };
  [...starts.values(), ...doseIn.values()].forEach((x) => x.onChanged.subscribe(refresh));
  pkDay?.onChanged.subscribe(() => {
    if (setting)
      return;
    pkDayTouched = true;
    refresh();
  });
  refresh();
  ui.tooltip.bind(footer, 'Excluded measurements left out; flagged ones kept with FLAG = 1. Doses in mg/kg, one ' +
    'compartment per drug (2, 3...). The columns are explained at the top of the file.');
  footer.style.color = 'var(--grey-5, #68727f)';
  const dlg = ui.dialog('Modelling dataset (NONMEM)')
    .add(ui.divV([table, conc, ...pkDay ? [ui.inputs([pkDay])] : [], pkNotes, footer], {style: {gap: '8px'}}))
    .onOK(() => {
      // one toast for the two files
      const [data, ids] = [`${name} - modelling dataset.csv`, `${name} - modelling dataset IDs.csv`];
      download(data, modellingCsv(d), 'text/csv', []);
      download(ids, idMapCsv(d), 'text/csv', [data, ids]);
    });
  okNamed(dlg, 'Export');
  showCentred(dlg);
}
