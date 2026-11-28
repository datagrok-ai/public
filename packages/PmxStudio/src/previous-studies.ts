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
import {FINISHED_TAG, exposureOfTag, groupsShown, lineOf, parseStudy, roleColumns, SavedStudy} from './study-file';
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
    // (not marked finished: still running, its effect interim)
    return {name: df.name, line: lineOf(df), meta: metaOf(df), ...exposure ? {exposure} : {},
      ...df.getTag(FINISHED_TAG) === 'true' ? {} : {interim: true as const},
      rows: Array.from({length: df.rowCount}, (_, i) =>
        ({animal: String(animal.get(i)), group: group.isNone(i) ? null : String(group.get(i)),
          day: day.isNone(i) ? null : day.get(i), volume: volume.isNone(i) ? null : volume.get(i)}))};
  } catch {
    return null; // no study table
  }
}

export const previousOfSaved = (s: SavedStudy): PreviousStudy => ({name: s.name,
  // (an excluded measurement is no longer the study's: as in the analysis, every reader leaves it out)
  line: s.plan?.predicted.line ?? s.line ?? null, meta: parseMeta(s.meta) ?? parseMeta(s.plan?.study),
  rows: s.rows.filter((r) => !r.excluded),
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

/** Where a plan's earlier studies come from in place of the workspace and the library (the demo plan: the demo's own
 *  finished study, never the user's library). */
export type StudySource = () => Promise<PreviousStudy[]>;

/** The earlier studies that may hold `armName`'s drug: those open in the workspace, then the user's saved ones whose
 *  listed compounds name it (a blinded study names none); with `from`, only those it gives. Never throws: a library
 *  that cannot be read gives none. */
export function previousStudies(armName: string, own: OwnStudy | null = null, from?: StudySource):
  Promise<PreviousStudy[]> {
  // the card's lookups on one name change (anchor, undosed, patients' dose) share one read
  if (!last || last.name !== armName || last.from !== from || Date.now() - last.at > 1000)
    last = {name: armName, from, at: Date.now(), studies: from ? from().catch(() => []) : readPrevious(armName)};
  // a study's own plan never offers that study: its arms would be estimated from themselves
  return last.studies.then((all) => all.filter((s) => !isOwnStudy(s, own)));
}
let last: {name: string, from?: StudySource, at: number, studies: Promise<PreviousStudy[]>} | null = null;

async function readPrevious(armName: string): Promise<PreviousStudy[]> {
  const parts = parseArmName(armName);
  if (parts.length !== 1)
    return [];
  const key = parts[0].key;
  const open = grok.shell.tables.map(previousOfTable).filter((s): s is PreviousStudy => s !== null);
  if (!listed || Date.now() - listed.at > 60000)
    listed = {at: Date.now(), entries: myEntries().catch(() => [])};
  // (a blinded study's row names no compound until revealed, so it is never an anchor before then). Only finished
  // studies (a running one's effect is interim), and only the user's own or one shared with them to read, never a
  // copy of someone else's (QA 0.1.3 final designer #9: the demo plan offered another user's interim QA copy); the
  // user's own first, each newest first (the library's order)
  const me = [grok.shell.user?.login, grok.shell.user?.friendlyName].filter((n): n is string => !!n)
    .map((n) => n.toLowerCase());
  const mine = (e: LibraryEntry) => !e.row || me.includes(e.owner.toLowerCase());
  const entries = (await listed.entries).filter((e) => e.kind === 'study' && (e.file || e.row) && e.stage === 3 &&
    !e.copyOf && (mine(e) || e.readOnly) && parseArmName(e.compounds.replace(/,/g, ' + ')).some((p) => p.key === key))
    .sort((a, b) => Number(mine(b)) - Number(mine(a)));
  const saved = await Promise.all(entries.map((e) => {
    const path = `${e.file?.fullPath ?? e.row!.id} ${e.updated}`;
    if (!read.has(path))
      read.set(path, entryText(e).then((t) => previousOfSaved(parseStudy(t))).catch(() => null));
    return read.get(path)!;
  }));
  return withSaved(open, saved.filter((s): s is PreviousStudy => s !== null));
}

/** The studies open in the workspace, then the saved ones not open: a saved study is left out for an open one of its
 *  name only when that one is finished too (an open copy still running, interim, never offered, hid the finished
 *  study of the plan's own model, and another model's was offered: QA 0.1.3 round 4 designer D6). */
export function withSaved(open: PreviousStudy[], saved: PreviousStudy[]): PreviousStudy[] {
  const names = new Set(open.filter((s) => !s.interim).map((s) => s.name));
  return [...open, ...saved.filter((s) => !names.has(s.name))];
}

/** The first earlier finished study with a group of `armName`'s drug alone that says no dose (study-effect.ts
 *  undosedIn). */
export async function findUndosed(armName: string, own: OwnStudy | null = null, from?: StudySource):
  Promise<{study: string, group: string} | null> {
  for (const s of (await previousStudies(armName, own, from)).filter((x) => !x.interim)) {
    const group = undosedIn(s, armName);
    if (group)
      return {study: s.name, group};
  }
  return null;
}

/** The patients' exposure as a mouse dose for `armName`, from the earlier finished study that kept the drug's exposure
 *  with the most doses around the patients' value (exposure.ts patientsDoseOf); on a tie the most recent (open tables
 *  first, then the library newest first); or the hint to add the schedule that would place it. null: none. Doses a day
 *  under `rules`, the plan's own (no dosing at weekends), not the earlier study's. */
export const findPatientsDose = async (armName: string, rules: DosingRules | null = null,
  own: OwnStudy | null = null, from?: StudySource): Promise<{text: string, tip: string} | null> =>
  patientsDoseOf(armName, (await previousStudies(armName, own, from)).filter((x) => !x.interim), rules);

/** The earlier studies offered, in order: finished ones only (a running one's effect is interim: never offered, QA
 *  round 3 designer W3), a study of the plan's tumor model `line` first (an effect seen in another model is not this
 *  one's, and Plan scales only within one). */
export const anchorOrder = (all: PreviousStudy[], line: string | null): PreviousStudy[] =>
  all.filter((s) => !s.interim).sort((a, b) => Number(a.line !== line) - Number(b.line !== line));

/** The group that can anchor `armName` in the first earlier study that has one (anchorOrder): of its groups, the one
 *  with the dose nearest the arm's (study-effect.ts nearestAnchor), not the first in the file. */
export async function findAnchor(armName: string, own: OwnStudy | null = null, line: string | null = null,
  from?: StudySource): Promise<Anchor | null> {
  for (const s of anchorOrder(await previousStudies(armName, own, from), line)) {
    const a = nearestAnchor(anchorsIn(s, armName), armName);
    if (a)
      return a;
  }
  return null;
}
