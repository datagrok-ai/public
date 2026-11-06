/* The user's studies, as the Overview lists them: each saved plan and followed study in their Home, read for its
   details (ID, line, what the groups receive, owner, tags) and how far it has gone. A study followed from a plan
   is listed once, at its furthest stage. Datagrok projects are not listed: a study shared as a project keeps its
   tags and marks in its table, and opens from Browse. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {myPlans, parsePlan} from './plan-file';
import {parseStudy, studiesFolder, STUDY_EXT} from './study-file';
import {compoundsOf, parseMeta, treatmentOf} from './study-meta';
import {futureDays} from './monitor';

export interface LibraryEntry {
  kind: 'plan' | 'study';
  /** 1 planned, 2 followed (Monitor) */
  stage: 1 | 2;
  file: DG.FileInfo | null;
  text: string;
  id: string;
  name: string;
  line: string;
  compounds: string;
  owner: string;
  tags: string[];
  /** ms since the epoch */
  updated: number;
  /** what each group receives, with its dose, route and schedule (else the compounds, else the arms) */
  treatment: string;
  /** day 0, YYYY-MM-DD (null: not fixed yet) */
  start: string | null;
  /** a followed study's last day measured and its next (from its usual gap); null for a plan */
  day: number | null;
  next: number | null;
  /** the study's last day as planned (null: no plan) */
  end: number | null;
}

/** A plan or study file's entry (null: not one, or damaged). */
export function entryOf(text: string, updated: number, file: DG.FileInfo | null = null): LibraryEntry | null {
  try {
    const x = JSON.parse(text);
    const study = x?.format === 'pmx-study' ? parseStudy(text) : null;
    const plan = study ? study.plan : parsePlan(text);
    const meta = parseMeta(study ? x.meta : plan?.study) ?? parseMeta(plan?.study);
    const days = (study?.rows ?? []).map((r) => r.day).filter((d) => Number.isFinite(d));
    const compounds = compoundsOf(meta) || (plan?.predicted.arms ?? []).map((a) => a.name).join(', ');
    return {kind: study ? 'study' : 'plan', stage: study ? 2 : 1, file, text, id: meta?.id ?? '', name: x.name ?? '',
      line: plan?.predicted.line ?? '', compounds, owner: meta?.owner || plan?.author || '', tags: meta?.tags ?? [],
      updated, treatment: treatmentOf(meta) || compounds, start: meta?.start ?? null,
      day: days.length ? days.reduce((a, b) => Math.max(a, b)) : null, next: days.length ? futureDays(days)[0] : null,
      end: plan?.predicted.studyEnd ?? null};
  } catch {
    return null;
  }
}

/** One entry per study, the furthest stage kept (a plan and the study followed from it share an ID), newest first. */
export function libraryOf(entries: (LibraryEntry | null)[]): LibraryEntry[] {
  const kept = new Map<string, LibraryEntry>();
  entries.forEach((e, i) => {
    if (!e)
      return;
    const key = e.id || `#${i}`;
    const was = kept.get(key);
    if (!was || e.stage > was.stage || (e.stage === was.stage && e.updated > was.updated))
      kept.set(key, e);
  });
  return [...kept.values()].sort((a, b) => b.updated - a.updated);
}

/** Whether an entry has every word of `query` in its ID, name, line, treatment, owner or tags. */
export function matches(e: LibraryEntry, query: string): boolean {
  const text = [e.id, e.name, e.line, e.treatment, e.owner, ...e.tags].join(' ').toLowerCase();
  return query.toLowerCase().split(/\s+/).filter((w) => w).every((w) => text.includes(w));
}

/** The user's plans and studies, read. ponytail: every file is read to list it; an index file once there are
 *  hundreds. */
export async function myLibrary(): Promise<LibraryEntry[]> {
  const studies = grok.dapi.files.list(studiesFolder(), false, STUDY_EXT).catch(() => [] as DG.FileInfo[]);
  const files = [...await myPlans().catch(() => [] as DG.FileInfo[]), ...await studies];
  return libraryOf(await Promise.all(files.map((f) => grok.dapi.files.readAsText(f)
    .then((text) => entryOf(text, f.updatedOn?.valueOf() ?? 0, f), () => null))));
}
