/* The rules the study table flags measurements to review by, and where they are kept. Each rule is a check a person
   would make by eye, with its threshold: a value that cannot be, a jump a caliper typo makes, a tumor far from its
   group on its day, a body weight that falls fast. A finding only asks for a look: nothing leaves the analysis unless
   a person excludes it (study-table.ts). The rules in force are the user's own, else the lab's (set by its group in
   the package settings), else these defaults. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from './package';
import {isFiniteNumber as ok, median} from './engine';

export interface OutlierRules {
  /** a negative volume or body weight, an animal measured twice on one day */
  impossible: boolean;
  /** a volume that changed since the animal's previous one this many times more, or less, than its group's did over
   *  the same days (no group measured on both days: than no change) */
  jump: {on: boolean, fold: number};
  /** a volume this many robust deviations (MAD of log volume) from its group's median on its day, in groups of at
   *  least `minAnimals` */
  group: {on: boolean, mads: number, minAnimals: number};
  /** a body weight this percent below the animal's previous weighing */
  weightDrop: {on: boolean, percent: number};
}

export const DEFAULT_RULES: OutlierRules = {
  impossible: true,
  jump: {on: true, fold: 3},
  group: {on: true, mads: 3.5, minAnimals: 5},
  weightDrop: {on: true, percent: 10},
};

/** The study's columns as the rules read them (an empty cell as null). */
export interface StudyColumns {
  animal: (string | null)[]; group: (string | null)[]; day: (number | null)[]; volume: (number | null)[];
  weight: (number | null)[] | null;
}

/** One finding: the row, the value it is about and why, in words. */
export interface Finding {row: number, on: 'volume' | 'weight' | 'row', reason: string}

const mm3 = (v: number) => `${Math.round(v)} mm³`;
/** The robust standard deviation of normal data from its median absolute deviation. */
const MAD_SD = 1.4826;
/** A group's spread is not taken below a caliper's own error on a large tumor (sd of log volume about 0.1): a group
 *  tighter than that is tight by its randomization window or by chance, and would flag every animal a little off. */
const SD_FLOOR = 0.1;
/** The animals of a group measured on the same two days that it takes to say what the group did between them. */
const MIN_SAME_DAYS = 3;

/** What the rules find in a study, row by row. */
export function findOutliers(s: StudyColumns, r: OutlierRules): Finding[] {
  const out: Finding[] = [];
  const rows = s.animal.map((_, i) => i);
  // each animal's rows by day, the rows with no animal or day left to the impossible rule
  const byAnimal = new Map<string, number[]>();
  for (const i of rows) {
    const a = s.animal[i];
    if (a !== null && a !== '' && ok(s.day[i]))
      byAnimal.set(a, [...(byAnimal.get(a) ?? []), i]);
  }
  byAnimal.forEach((is) => is.sort((p, q) => s.day[p]! - s.day[q]!));

  if (r.impossible) {
    for (const i of rows) {
      if (ok(s.volume[i]) && s.volume[i]! < 0)
        out.push({row: i, on: 'volume', reason: 'A negative volume'});
      if (s.weight && ok(s.weight[i]) && s.weight[i]! <= 0)
        out.push({row: i, on: 'weight', reason: 'A body weight of zero or less'});
    }
    byAnimal.forEach((is, a) => is.forEach((i, k) => {
      if (k > 0 && s.day[is[k - 1]] === s.day[i])
        out.push({row: i, on: 'row', reason: `${a} measured twice on day ${s.day[i]}`});
    }));
  }

  // each value against the animal's previous one of the same kind
  const previous = (is: number[], k: number, v: (number | null)[]) => {
    for (let j = k - 1; j >= 0; j--) {
      if (ok(v[is[j]]) && v[is[j]]! > 0)
        return is[j];
    }
    return null;
  };
  if (r.jump.on && r.jump.fold > 1) {
    // each animal's change since its previous volume, on the log scale, and the days it spans
    const steps: {row: number, from: number, change: number, days: string}[] = [];
    byAnimal.forEach((is) => is.forEach((i, k) => {
      const v = s.volume[i];
      const p = previous(is, k, s.volume);
      if (ok(v) && v > 0 && p !== null) {
        steps.push({row: i, from: p, change: Math.log(v / s.volume[p]!),
          days: `${s.group[i]}\u0000${s.day[p]}\u0000${s.day[i]}`});
      }
    }));
    // against what its group did over the same days: a tumor that triples in a week is growth where its group's do,
    // a typo where they do not. Too few of the group measured on both days: against no change
    const together = new Map<string, number[]>();
    for (const t of steps)
      together.set(t.days, [...(together.get(t.days) ?? []), t.change]);
    for (const t of steps) {
      const same = together.get(t.days)!;
      const usual = same.length >= MIN_SAME_DAYS ? median(same) : 0;
      if (Math.abs(t.change - usual) < Math.log(r.jump.fold))
        continue;
      const ratio = Math.exp(t.change);
      const times = (x: number) => x >= 1 ? `${x.toFixed(1)} times` : `1/${(1 / x).toFixed(1)} of`;
      out.push({row: t.row, on: 'volume', reason: `${times(ratio)} its volume on day ${s.day[t.from]} ` +
        `(${mm3(s.volume[t.from]!)} to ${mm3(s.volume[t.row]!)})` +
        (same.length >= MIN_SAME_DAYS ? `; its group's changed ${Math.exp(usual).toFixed(1)}-fold` : '')});
    }
  }
  if (r.weightDrop.on && s.weight) {
    const w = s.weight;
    byAnimal.forEach((is) => is.forEach((i, k) => {
      const p = previous(is, k, w);
      if (!ok(w[i]) || w[i]! <= 0 || p === null)
        return;
      const drop = 100 * (1 - w[i]! / w[p]!);
      if (drop >= r.weightDrop.percent)
        out.push({row: i, on: 'weight', reason: `${drop.toFixed(0)}% lighter than on day ${s.day[p]}`});
    }));
  }

  // each volume against its group on its day, on the log scale tumors grow on
  if (r.group.on) {
    const cells = new Map<string, number[]>();
    for (const i of rows) {
      const [g, d, v] = [s.group[i], s.day[i], s.volume[i]];
      if (g !== null && ok(d) && ok(v) && v > 0)
        cells.set(`${g}\u0000${d}`, [...(cells.get(`${g}\u0000${d}`) ?? []), i]);
    }
    cells.forEach((is) => {
      if (is.length < r.group.minAnimals)
        return;
      const logs = is.map((i) => Math.log(s.volume[i]!));
      const med = median(logs);
      const sd = Math.max(MAD_SD * median(logs.map((x) => Math.abs(x - med))), SD_FLOOR);
      is.forEach((i, k) => {
        if (Math.abs(logs[k] - med) / sd > r.group.mads) {
          out.push({row: i, on: 'volume', reason: `Far from its group on day ${s.day[i]}: ${mm3(s.volume[i]!)}, ` +
            `the group's median ${mm3(Math.exp(med))}`});
        }
      });
    });
  }
  return out;
}

// --- where the rules are kept: the user's own, else the lab's, else the defaults

const SETTINGS = 'PmxStudio';
const USER_KEY = 'outlierRules';
/** The package settings: the lab's rules (JSON) and the group that may change them. */
const LAB_RULES = 'outlierRules';
const LAB_GROUP = 'labGroup';

/** Rules read from storage, each field from the defaults when missing (a rule added later keeps its default). */
function parse(text: unknown): OutlierRules | null {
  if (typeof text !== 'string' || text.trim() === '')
    return null;
  try {
    const x = JSON.parse(text);
    return {impossible: x.impossible ?? DEFAULT_RULES.impossible, jump: {...DEFAULT_RULES.jump, ...x.jump},
      group: {...DEFAULT_RULES.group, ...x.group}, weightDrop: {...DEFAULT_RULES.weightDrop, ...x.weightDrop}};
  } catch {
    return null;
  }
}

export type RulesSource = 'user' | 'lab' | 'default';

export const labRules = (): OutlierRules | null => parse(_package.settings?.[LAB_RULES]);
export const userRules = (): OutlierRules | null => parse(grok.userSettings.getValue(SETTINGS, USER_KEY));

/** The rules in force and whose they are. */
export function rulesInForce(): {rules: OutlierRules, from: RulesSource} {
  const user = userRules();
  if (user)
    return {rules: user, from: 'user'};
  const lab = labRules();
  return lab ? {rules: lab, from: 'lab'} : {rules: DEFAULT_RULES, from: 'default'};
}

/** The user's own rules; null: back to the lab's. */
export function saveUserRules(r: OutlierRules | null): void {
  if (r)
    grok.userSettings.add(SETTINGS, USER_KEY, JSON.stringify(r));
  else
    grok.userSettings.delete(SETTINGS, USER_KEY);
  void grok.userSettings.flush();
}

/** The lab's rules, for every user (who has none of their own). */
export async function saveLabRules(r: OutlierRules): Promise<void> {
  // every setting kept, as panels.ts does: the lab's group and its layouts live there too
  // @ts-ignore (setSettings is declared with a Map; the platform takes the plain object, as ApiTests does)
  await _package.setSettings({..._package.settings, [LAB_RULES]: JSON.stringify(r)}, DG.Group.allUsers);
}

/** Whether this user may change the lab's rules: an Administrator, or a member of the lab's group (by name, set in
 *  the package settings). */
export async function canEditLabRules(): Promise<boolean> {
  const lab = String(_package.settings?.[LAB_GROUP] ?? '').trim().toLowerCase();
  const admins = DG.Group.defaultGroupsIds['Administrators'];
  const groups = await grok.dapi.groups.currentUserGroups();
  return groups.some((g) => g.id === admins || (lab !== '' &&
    [g.name, g.friendlyName].some((n) => String(n ?? '').toLowerCase() === lab)));
}
