/* The user's earlier studies Plan may scale an arm's effect from (study-effect.ts): those open in the workspace, as
   Analyze has them, and those saved in their Home. Apart from study-effect.ts, which the arms' cards read: these read
   the library and the study files, which lead back to the cards (plan-file.ts). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {parseArmName} from './arm-name';
import {LibraryEntry} from './library';
import {entryText, myEntries} from './store';
import {patientsDoseOf} from './exposure';
import {Anchor, PreviousStudy, anchorsIn, nearestAnchor, undosedIn} from './study-effect';
import {exposureOfTag, groupsShown, lineOf, parseStudy, roleColumns, SavedStudy} from './study-file';
import {metaOf, parseMeta} from './study-meta';
import type {DosingRules} from './schedule';

/** A study open in the workspace (as Analyze has it), as Plan reads it; null: no study table, or one whose groups
 *  are hidden (blinded). */
export function previousOfTable(df: DG.DataFrame): PreviousStudy | null {
  try {
    if (!groupsShown(df))
      return null;
    const {animal, group, day, volume} = roleColumns(df);
    const exposure = exposureOfTag(df);
    return {name: df.name, line: lineOf(df), meta: metaOf(df), ...exposure ? {exposure} : {},
      rows: Array.from({length: df.rowCount}, (_, i) =>
        ({animal: String(animal.get(i)), group: group.isNone(i) ? null : String(group.get(i)),
          day: day.isNone(i) ? null : day.get(i), volume: volume.isNone(i) ? null : volume.get(i)}))};
  } catch {
    return null; // no study table
  }
}

export const previousOfSaved = (s: SavedStudy): PreviousStudy => ({name: s.name,
  line: s.plan?.predicted.line ?? s.line ?? null, meta: parseMeta(s.meta) ?? parseMeta(s.plan?.study), rows: s.rows,
  ...s.exposure ? {exposure: s.exposure} : {}});

/** The study whose plan is open (its ID, its name): never an earlier study of its own arms. */
export interface OwnStudy {id: string, name: string}

/** Whether `s` is the study `own` plans: the same ID, or the same name (a followed study's table takes its plan's
 *  name, and has no ID until saved or followed). */
export const isOwnStudy = (s: Pick<PreviousStudy, 'name' | 'meta'>, own: OwnStudy | null | undefined): boolean =>
  !!own && (!!own.id && s.meta?.id === own.id || !!own.name && s.name === own.name);

// the library listed once a minute at most, and each study file read once a version (its date)
let listed: {at: number, entries: Promise<LibraryEntry[]>} | null = null;
const read = new Map<string, Promise<PreviousStudy | null>>();

/** The earlier studies that may hold `armName`'s drug: those open in the workspace, then the user's saved ones whose
 *  listed compounds name it (a blinded study names none). Never throws: a library that cannot be read gives none. */
export function previousStudies(armName: string, own: OwnStudy | null = null): Promise<PreviousStudy[]> {
  // the card's lookups on one name change (anchor, undosed, patients' dose) share one read
  if (!last || last.name !== armName || Date.now() - last.at > 1000)
    last = {name: armName, at: Date.now(), studies: readPrevious(armName)};
  // a study's own plan never offers that study: its arms would be estimated from themselves
  return last.studies.then((all) => all.filter((s) => !isOwnStudy(s, own)));
}
let last: {name: string, at: number, studies: Promise<PreviousStudy[]>} | null = null;

async function readPrevious(armName: string): Promise<PreviousStudy[]> {
  const parts = parseArmName(armName);
  if (parts.length !== 1)
    return [];
  const key = parts[0].key;
  const open = grok.shell.tables.map(previousOfTable).filter((s): s is PreviousStudy => s !== null);
  if (!listed || Date.now() - listed.at > 60000)
    listed = {at: Date.now(), entries: myEntries().catch(() => [])};
  // (a blinded study's row names no compound until revealed, so it is never an anchor before then)
  const entries = (await listed.entries).filter((e) => e.kind === 'study' && (e.file || e.row) &&
    parseArmName(e.compounds.replace(/,/g, ' + ')).some((p) => p.key === key));
  const saved = await Promise.all(entries.map((e) => {
    const path = `${e.file?.fullPath ?? e.row!.id} ${e.updated}`;
    if (!read.has(path))
      read.set(path, entryText(e).then((t) => previousOfSaved(parseStudy(t))).catch(() => null));
    return read.get(path)!;
  }));
  const names = new Set(open.map((s) => s.name));
  return [...open, ...saved.filter((s): s is PreviousStudy => s !== null && !names.has(s.name))];
}

/** The first earlier study with a group of `armName`'s drug alone that says no dose (study-effect.ts undosedIn). */
export async function findUndosed(armName: string, own: OwnStudy | null = null):
  Promise<{study: string, group: string} | null> {
  for (const s of await previousStudies(armName, own)) {
    const group = undosedIn(s, armName);
    if (group)
      return {study: s.name, group};
  }
  return null;
}

/** The patients' exposure as a mouse dose for `armName`, from the earlier study that kept the drug's exposure with
 *  the most doses around the patients' value (exposure.ts patientsDoseOf); on a tie the most recent (open tables
 *  first, then the library newest first); or the hint to add the schedule that would place it. null: none. Doses a day
 *  under `rules`, the plan's own (no dosing at weekends), not the earlier study's. */
export const findPatientsDose = async (armName: string, rules: DosingRules | null = null,
  own: OwnStudy | null = null): Promise<{text: string, tip: string} | null> =>
  patientsDoseOf(armName, await previousStudies(armName, own), rules);

/** The group that can anchor `armName` in the first earlier study that has one, a study of the plan's tumor model
 *  `line` first (an effect seen in another model is not this one's, and Plan scales only within one): of its groups,
 *  the one with the dose nearest the arm's (study-effect.ts nearestAnchor), not the first in the file. */
export async function findAnchor(armName: string, own: OwnStudy | null = null, line: string | null = null):
  Promise<Anchor | null> {
  const all = await previousStudies(armName, own);
  for (const s of [...all.filter((x) => x.line === line), ...all.filter((x) => x.line !== line)]) {
    const a = nearestAnchor(anchorsIn(s, armName), armName);
    if (a)
      return a;
  }
  return null;
}
