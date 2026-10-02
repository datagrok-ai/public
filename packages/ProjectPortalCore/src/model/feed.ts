/** Activity feed (F-21, FR-81): what happened to a program, assembled from the portal's own records
 * rather than a separate event store. Each record kind already says who did what and when: the gate
 * log, evidence, saved and selected runs, documents, snapshots, compound role changes, the sync
 * history (registry and work-unit changes, new data) and the follows' last news. Pure: the caller
 * reads the records; this turns them into entries, merges, filters and pages them. */
import type {ProgramChange} from './sync-runs';

export type FeedKind = 'decision' | 'evidence' | 'run' | 'selection' | 'document' | 'snapshot' |
  'compound' | 'work_unit' | 'registry' | 'new_data';

export const FEED_KINDS: FeedKind[] = ['decision', 'evidence', 'run', 'selection', 'document', 'snapshot', 'compound',
  'work_unit', 'registry', 'new_data'];

export const FEED_KIND_LABELS: Record<FeedKind, string> = {
  decision: 'Gate decisions', evidence: 'Evidence', run: 'Runs saved', selection: 'Runs selected',
  document: 'Documents', snapshot: 'Snapshots', compound: 'Compound roles',
  work_unit: 'Work units', registry: 'Registry', new_data: 'New data',
};

/** The record an entry opens, by the kinds provenance starts from. */
export interface FeedRef {
  kind: 'decision' | 'evidence' | 'model_run' | 'snapshot' | 'document' | 'program';
  id: string;
}

export interface FeedEntry {
  /** Stable across reads: the record and what happened to it. */
  key: string;
  /** ISO time. */
  at: string;
  kind: FeedKind;
  programId: string;
  program: string;
  portfolio?: string | null;
  /** The user's id; null for the sync or the data. */
  who: string | null;
  text: string;
  /** The work unit the entry is about, as `series KIN-A`. */
  unit?: string | null;
  /** The compounds it is about. */
  compounds?: string[];
  ref?: FeedRef | null;
}

type When = string | null | undefined;

/** The records one program's feed is read from; times are ISO strings. */
export interface FeedRecords {
  program: {id: string, code: string, portfolio?: string | null};
  stages: {id: string, code: string, name: string}[];
  decisions: {id: string, decision: string, from: string | null, to: string | null, by: string | null, on: When,
    source?: string | null, rationale?: string | null}[];
  evidence: {id: string, decisionId: string | null, kind: string, name: string, by: string | null, on: When}[];
  runs: {id: string, title: string, model: string, by: string | null, on: When}[];
  selections: {id: string, run: string, criterion: string, by: string | null, on: When, rationale?: string | null}[];
  documents: {id: string, title: string, version?: string | null, status: string, linkedBy: string | null, linkedOn: When,
    approvedBy?: string | null, approvedOn?: When}[];
  snapshots: {id: string, label: string, by: string | null, on: When}[];
  roleChanges: {id: string, compound: string, from: string | null, to: string, by: string | null, on: When}[];
  /** Sync runs of the portfolio: this program's change in each, and the new-data lines about it. */
  syncs: {id: string, on: When, by: string | null, change: ProgramChange | null, newData: string[]}[];
  /** The follows' last news about the program, with when it was told. */
  news: {id: string, line: string, on: When}[];
}

/** Evidence attached this close to its decision belongs to the decision; later, it is an entry. */
export const EVIDENCE_WINDOW_MS = 2 * 60 * 1000;
/** The same news line told this close together (to several followers, or by the sync and a check)
 * is one entry. */
export const NEWS_WINDOW_MS = 15 * 60 * 1000;

const ms = (s: When): number => s ? Date.parse(s) : NaN;
const plural = (n: number, one: string) => `${n} ${n === 1 ? one : `${one}s`}`;

/** The work unit a new-data line is about: `EXS-KIN1 series KIN-A: 2 new values` → `series KIN-A`. */
export function newsUnit(program: string, line: string): string | null {
  const head = line.split(':')[0];
  return head.startsWith(`${program} `) ? head.slice(program.length + 1) : null;
}

/** One program's entries, unsorted. */
export function programFeedEntries(r: FeedRecords): FeedEntry[] {
  const out: FeedEntry[] = [];
  const base = {programId: r.program.id, program: r.program.code, portfolio: r.program.portfolio ?? null};
  const stage = (id: string | null) => id ? r.stages.find((s) => s.id === id)?.name ?? '?' : '—';
  const push = (e: Omit<FeedEntry, 'programId' | 'program' | 'portfolio' | 'at'> & {at: When}) => {
    if (e.at && !Number.isNaN(ms(e.at)))
      out.push({...base, ...e, at: new Date(ms(e.at)).toISOString()});
  };

  const decisionLabel = new Map<string, string>();
  for (const d of r.decisions) {
    const how = d.source === 'data' ? ' (by data)' : d.source === 'registry' || d.decision === 'sync' ? ' (registry)' : '';
    decisionLabel.set(d.id, `${d.decision} ${stage(d.from)} → ${stage(d.to)}`);
    const withIt = r.evidence.filter((e) => e.decisionId === d.id && Math.abs(ms(e.on) - ms(d.on)) <= EVIDENCE_WINDOW_MS);
    push({key: `decision:${d.id}`, at: d.on, kind: 'decision', who: d.by,
      text: `${d.decision[0].toUpperCase()}${d.decision.slice(1)}: ${stage(d.from)} → ${stage(d.to)}${how}` +
        (withIt.length ? `, with ${plural(withIt.length, 'piece')} of evidence` : '') + (d.rationale ? ` — ${d.rationale}` : ''),
      ref: {kind: 'decision', id: d.id}});
  }
  for (const e of r.evidence) {
    const d = r.decisions.find((x) => x.id === e.decisionId);
    if (d && Math.abs(ms(e.on) - ms(d.on)) <= EVIDENCE_WINDOW_MS)
      continue;
    push({key: `evidence:${e.id}`, at: e.on, kind: 'evidence', who: e.by,
      text: `Evidence added${d ? ` to ${decisionLabel.get(d.id)}` : ''}: ${e.name} (${e.kind.replace('_', ' ')})`,
      ref: {kind: 'evidence', id: e.id}});
  }
  for (const x of r.runs)
    push({key: `run:${x.id}`, at: x.on, kind: 'run', who: x.by, text: `Run saved: ${x.title} (${x.model})`, ref: {kind: 'model_run', id: x.id}});
  for (const s of r.selections) {
    push({key: `selection:${s.id}`, at: s.on, kind: 'selection', who: s.by,
      text: `Run selected for ${s.criterion}: ${s.run}${s.rationale ? ` — ${s.rationale}` : ''}`});
  }
  for (const d of r.documents) {
    const title = `${d.title}${d.version ? ` v${d.version}` : ''}`;
    push({key: `document:${d.id}`, at: d.linkedOn, kind: 'document', who: d.linkedBy, text: `Document linked: ${title}`,
      ref: {kind: 'document', id: d.id}});
    if (d.approvedOn) {
      push({key: `document-approved:${d.id}`, at: d.approvedOn, kind: 'document', who: d.approvedBy ?? null,
        text: `Document approved: ${title}`, ref: {kind: 'document', id: d.id}});
    }
  }
  for (const s of r.snapshots)
    push({key: `snapshot:${s.id}`, at: s.on, kind: 'snapshot', who: s.by, text: `Snapshot taken: ${s.label}`, ref: {kind: 'snapshot', id: s.id}});
  for (const c of r.roleChanges) {
    push({key: `role:${c.id}`, at: c.on, kind: 'compound', who: c.by, compounds: [c.compound],
      text: `${c.compound}: ${c.from ?? 'new'} → ${c.to}`});
  }

  for (const s of r.syncs) {
    const c = s.change;
    if (c) {
      // stage moves and role changes are in the gate log and the role history already
      const reg: string[] = [];
      if (c.created)
        reg.push('new in the registry');
      if (c.fields?.length)
        reg.push(`registry changed ${c.fields.join(', ')}`);
      if (c.missing)
        reg.push('no longer in the registry feed');
      if (reg.length)
        push({key: `sync:${s.id}:registry`, at: s.on, kind: 'registry', who: s.by, text: `Sync: ${reg.join('; ')}`});
      const u = c.units;
      for (const [what, list] of [['new', u?.created], ['updated', u?.updated], ['missing from its feed', u?.missing]] as const) {
        for (const unit of list ?? [])
          push({key: `sync:${s.id}:unit:${what}:${unit}`, at: s.on, kind: 'work_unit', who: s.by, unit, text: `${unit}: ${what}`});
      }
      const k = c.compounds;
      if (k?.added.length) {
        push({key: `sync:${s.id}:compounds-added`, at: s.on, kind: 'compound', who: s.by, compounds: k.added,
          text: `${plural(k.added.length, 'compound')} added: ${k.added.join(', ')}`});
      }
      if (k?.missing.length) {
        push({key: `sync:${s.id}:compounds-missing`, at: s.on, kind: 'compound', who: s.by, compounds: k.missing,
          text: `${plural(k.missing.length, 'compound')} no longer in the feed: ${k.missing.join(', ')}`});
      }
    }
    s.newData.forEach((line, i) => push({key: `sync:${s.id}:news:${i}`, at: s.on, kind: 'new_data', who: null,
      unit: newsUnit(r.program.code, line), text: `New data: ${line}`}));
  }
  for (const n of r.news) {
    push({key: `news:${n.id}:${n.on}`, at: n.on, kind: 'new_data', who: null, unit: newsUnit(r.program.code, n.line),
      text: `New data: ${n.line}`});
  }
  return dedupeNews(out);
}

/** The same new-data line within NEWS_WINDOW_MS is one entry, the earliest. */
function dedupeNews(entries: FeedEntry[]): FeedEntry[] {
  const kept: FeedEntry[] = [];
  const news = entries.filter((e) => e.kind === 'new_data').sort((a, b) => ms(a.at) - ms(b.at));
  for (const n of news) {
    if (!kept.some((k) => k.programId === n.programId && k.text === n.text && Math.abs(ms(k.at) - ms(n.at)) <= NEWS_WINDOW_MS))
      kept.push(n);
  }
  return [...entries.filter((e) => e.kind !== 'new_data'), ...kept];
}

const KIND_ORDER = new Map(FEED_KINDS.map((k, i) => [k, i]));

/** Several feeds as one, newest first; an entry read twice (a program in two lists) counts once. */
export function mergeFeeds(...lists: FeedEntry[][]): FeedEntry[] {
  const by = new Map<string, FeedEntry>();
  for (const e of lists.flat())
    by.set(`${e.programId}|${e.key}`, e);
  return [...by.values()].sort((a, b) => ms(b.at) - ms(a.at) || KIND_ORDER.get(a.kind)! - KIND_ORDER.get(b.kind)! ||
    a.program.localeCompare(b.program) || a.text.localeCompare(b.text));
}

export interface FeedFilter {
  /** Only these kinds; empty or absent = all. */
  kinds?: FeedKind[] | null;
  /** Only what is about this series: entries tagged with it, and those about its compounds. */
  series?: string | null;
  /** The series' compounds (from the program dataset's `series` column). */
  seriesCompounds?: string[] | null;
  /** Only entries after this ISO time. */
  since?: string | null;
}

export function filterFeed(entries: FeedEntry[], f: FeedFilter = {}): FeedEntry[] {
  const kinds = f.kinds?.length ? new Set(f.kinds) : null;
  const unit = f.series ? `series ${f.series}` : null;
  const compounds = new Set(f.seriesCompounds ?? []);
  const since = f.since ? ms(f.since) : NaN;
  return entries.filter((e) => (!kinds || kinds.has(e.kind)) &&
    (!unit || e.unit === unit || (e.compounds ?? []).some((c) => compounds.has(c))) &&
    (Number.isNaN(since) || ms(e.at) > since));
}

/** Kinds as the filter offers them: a comma-separated list, each a FEED_KINDS code. */
export function parseKinds(text: string | null | undefined): FeedKind[] {
  const out: FeedKind[] = [];
  for (const k of String(text ?? '').split(',').map((s) => s.trim()).filter((s) => s.length > 0)) {
    if (!(FEED_KINDS as string[]).includes(k))
      throw new Error(`"${k}" is not a feed kind; one of ${FEED_KINDS.join(', ')}`);
    out.push(k as FeedKind);
  }
  return out;
}

export interface FeedPage {
  items: FeedEntry[];
  /** 0-based. */
  page: number;
  pages: number;
  total: number;
}

export function pageFeed(entries: FeedEntry[], page: number, size: number): FeedPage {
  const pages = Math.max(1, Math.ceil(entries.length / size));
  const p = Math.min(Math.max(0, page), pages - 1);
  return {items: entries.slice(p * size, (p + 1) * size), page: p, pages, total: entries.length};
}

/** Entries per kind, for the filter's labels. */
export function countByKind(entries: FeedEntry[]): Map<FeedKind, number> {
  const out = new Map<FeedKind, number>();
  for (const e of entries)
    out.set(e.kind, (out.get(e.kind) ?? 0) + 1);
  return out;
}

/** One line: date and time (UTC), program, what happened, who. */
export function describeEntry(e: FeedEntry, who?: string | null, withProgram = false): string {
  return `${e.at.slice(0, 16).replace('T', ' ')} · ${withProgram ? `${e.program} · ` : ''}${e.text}${who ? ` · ${who}` : ''}`;
}
