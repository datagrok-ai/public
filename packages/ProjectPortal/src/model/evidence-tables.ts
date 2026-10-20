/** Evidence tables generated into documents (F-40, FR-72). For a document criterion, the evidence
 * behind it is gathered into one table: the run selected for the criterion with its inputs and
 * outputs and their units, the newest program snapshot with its frozen values, the verdicts of the
 * gate's data conditions, and the decisions that cite the criterion, its run or its documents. The
 * table is written as a CSV and a Markdown file and linked as a new draft version of a document;
 * the document records which evidence it was generated from, each source with a content hash.
 *
 * "Evidence changed since written" (F-41, FR-73): later the same sources are read again and compared
 * with the record: another run selected, a newer snapshot, a data verdict that changed, a decision
 * added or gone. Any difference makes the document stale. Pure: the caller reads the records and
 * writes the files. */
import {formatRunValue} from './run-selection';

/** Where a row of the table comes from. */
export type EvidenceSourceKind = 'run' | 'snapshot' | 'data_condition' | 'decision';

export const SOURCE_LABELS: Record<EvidenceSourceKind, string> = {
  run: 'Run', snapshot: 'Snapshot', data_condition: 'Data condition', decision: 'Decision',
};

/** A parameter of a run: its name, caption and units, and the value the saved call holds. */
export interface RunParam {
  name: string;
  caption?: string | null;
  units?: string | null;
  value: unknown;
}

/** Everything one evidence table is assembled from, read by the caller. */
export interface EvidenceInputs {
  program: {code: string, name?: string | null};
  criterion: {id: string, code: string, name: string, stage: string};
  /** The run selected for the criterion, with the saved call's values; null when none is. */
  run?: {runId: string, title: string, model: string, selectedOn?: string | null, selectedBy?: string | null,
    rationale?: string | null, inputs: RunParam[], outputs: RunParam[]} | null;
  /** The program's newest snapshot; null when it has none. */
  snapshot?: {id: string, label: string, takenOn: string, sha256?: string | null, rows?: number | null,
    /** The criterion's state as the snapshot froze it. */
    criterionState?: string | null,
    /** Frozen rows of the program dataset, each as column → value; the caller may cut them short. */
    values?: Record<string, unknown>[]} | null;
  /** The data conditions of the criterion's gate, as readiness evaluates them now. */
  conditions: {id: string, code: string, name: string, state: string, summary?: string | null}[];
  /** The decisions that cite the criterion, its run or its documents. */
  decisions: {id: string, decision: string, on: string, from?: string | null, to?: string | null, rationale?: string | null}[];
}

/** One source as the document records it: enough to say later what changed. */
export interface RecordedSource {
  /** The slot it fills: `run`, `snapshot`, `condition:<code>`, `decision:<id>`. */
  key: string;
  kind: EvidenceSourceKind;
  /** The soft reference: a run_id, a snapshot id, a criterion id, a gate-log entry id. */
  ref: string;
  name: string;
  /** For a data condition, its verdict; for a run, when it was selected. */
  state?: string | null;
  hash: string;
}

/** A 53-bit string hash (cyrb53) as 14 hex digits. Not cryptographic: it tells content apart. */
export function contentHash(text: string, seed = 0): string {
  let h1 = 0xdeadbeef ^ seed;
  let h2 = 0x41c6ce57 ^ seed;
  for (let i = 0; i < text.length; i++) {
    const ch = text.charCodeAt(i);
    h1 = Math.imul(h1 ^ ch, 2654435761);
    h2 = Math.imul(h2 ^ ch, 1597334677);
  }
  h1 = Math.imul(h1 ^ (h1 >>> 16), 2246822507) ^ Math.imul(h2 ^ (h2 >>> 13), 3266489909);
  h2 = Math.imul(h2 ^ (h2 >>> 16), 2246822507) ^ Math.imul(h1 ^ (h1 >>> 13), 3266489909);
  return (4294967296 * (2097151 & h2) + (h1 >>> 0)).toString(16).padStart(14, '0');
}

const source = (key: string, kind: EvidenceSourceKind, ref: string, name: string, basis: unknown[],
  state?: string | null): RecordedSource => ({key, kind, ref, name, ...(state ? {state} : {}),
  hash: contentHash([kind, ref, ...basis.map((b) => b === null || b === undefined ? '' : String(b))].join('|'))});

/** The sources of an evidence table, each with its hash. What a value of the run or the snapshot
 * says is not hashed: a saved call and a snapshot are written once, so their id (and the snapshot's
 * SHA-256) stands for their content. */
export function evidenceSources(i: Pick<EvidenceInputs, 'run' | 'snapshot' | 'conditions' | 'decisions'>): RecordedSource[] {
  const out: RecordedSource[] = [];
  if (i.run)
    out.push(source('run', 'run', i.run.runId, i.run.title, [i.run.selectedOn], i.run.selectedOn ?? null));
  if (i.snapshot)
    out.push(source('snapshot', 'snapshot', i.snapshot.id, i.snapshot.label, [i.snapshot.takenOn, i.snapshot.sha256]));
  for (const c of i.conditions)
    out.push(source(`condition:${c.code}`, 'data_condition', c.id, c.name, [c.state, c.summary], c.state));
  for (const d of i.decisions)
    out.push(source(`decision:${d.id}`, 'decision', d.id, `${d.decision} ${d.on.slice(0, 10)}`, [d.decision, d.on]));
  return out;
}

/** One hash over all sources, the same whatever their order. */
export const evidenceHash = (sources: RecordedSource[]): string =>
  contentHash(sources.map((s) => `${s.key}=${s.hash}`).sort().join(';'));

/** One row of the table. */
export interface EvidenceRow {
  /** `run`, `run input`, `run output`, `snapshot`, `snapshot value`, `data condition`, `decision`. */
  section: string;
  /** The record the row comes from: the run's title, the snapshot's label, the decision. */
  source: string;
  /** Its soft reference. */
  ref: string;
  item: string;
  value: string;
  units: string;
  state: string;
  /** ISO date or time, when the row has one. */
  on: string;
}

export const EVIDENCE_COLUMNS: (keyof EvidenceRow)[] = ['section', 'source', 'ref', 'item', 'value', 'units', 'state', 'on'];

/** How many frozen rows of the snapshot the table shows; the rest are counted. */
export const MAX_SNAPSHOT_ROWS = 25;

const row = (r: Partial<EvidenceRow> & Pick<EvidenceRow, 'section' | 'source' | 'item'>): EvidenceRow =>
  ({ref: '', value: '', units: '', state: '', on: '', ...r});

/** The table: the run (its inputs, then outputs, with units), the snapshot (the criterion as frozen,
 * then its rows), the data conditions, the decisions. A part with nothing in it says so in one row. */
export function evidenceTableRows(i: EvidenceInputs): EvidenceRow[] {
  const out: EvidenceRow[] = [];
  const r = i.run;
  if (r) {
    out.push(row({section: 'run', source: r.title, ref: r.runId, item: 'selected', value: r.model.split(':').pop()!,
      state: r.selectedBy ?? '', on: r.selectedOn ?? '', units: ''}));
    if (r.rationale)
      out.push(row({section: 'run', source: r.title, ref: r.runId, item: 'rationale', value: r.rationale}));
    for (const [section, params] of [['run input', r.inputs], ['run output', r.outputs]] as const) {
      for (const p of params)
        out.push(row({section, source: r.title, ref: r.runId, item: p.caption || p.name, value: formatRunValue(p.value), units: p.units ?? ''}));
    }
  } else
    out.push(row({section: 'run', source: '', item: 'no run selected for the criterion'}));
  const s = i.snapshot;
  if (s) {
    out.push(row({section: 'snapshot', source: s.label, ref: s.id, item: i.criterion.code, state: s.criterionState ?? '',
      value: s.rows === null || s.rows === undefined ? 'no table' : `${s.rows} rows`, on: s.takenOn}));
    const values = s.values ?? [];
    for (const v of values.slice(0, MAX_SNAPSHOT_ROWS)) {
      const [first, ...rest] = Object.entries(v);
      out.push(row({section: 'snapshot value', source: s.label, ref: s.id, item: first ? formatRunValue(first[1]) : '',
        value: rest.map(([k, x]) => `${k}: ${formatRunValue(x)}`).join('; ')}));
    }
    const total = s.rows ?? values.length;
    if (total > MAX_SNAPSHOT_ROWS && values.length)
      out.push(row({section: 'snapshot value', source: s.label, ref: s.id, item: `${total - MAX_SNAPSHOT_ROWS} more rows in the snapshot`}));
  } else
    out.push(row({section: 'snapshot', source: '', item: 'no snapshot of the program'}));
  for (const c of i.conditions) {
    out.push(row({section: 'data condition', source: c.name, ref: c.id, item: c.code, value: c.summary ?? '',
      state: c.state.replace(/_/g, ' ')}));
  }
  if (!i.conditions.length)
    out.push(row({section: 'data condition', source: '', item: `the ${i.criterion.stage} gate has no data condition`}));
  for (const d of i.decisions) {
    out.push(row({section: 'decision', source: `${d.decision} ${d.on.slice(0, 10)}`, ref: d.id, item: d.decision,
      value: [d.from && d.to ? `${d.from} → ${d.to}` : '', d.rationale ?? ''].filter((x) => x).join(' — '), on: d.on}));
  }
  if (!i.decisions.length)
    out.push(row({section: 'decision', source: '', item: 'no decision cites it yet'}));
  return out;
}

const csvCell = (v: string) => /[",\n\r]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v;

export function evidenceCsv(rows: EvidenceRow[]): string {
  return [EVIDENCE_COLUMNS.join(','), ...rows.map((r) => EVIDENCE_COLUMNS.map((c) => csvCell(r[c])).join(','))].join('\n') + '\n';
}

const mdCell = (v: string) => v.replace(/\|/g, '\\|').replace(/\r?\n/g, ' ');

/** The document: a heading, what it was generated from and when, then one Markdown table per section. */
export function evidenceMarkdown(i: Pick<EvidenceInputs, 'program' | 'criterion'>, rows: EvidenceRow[],
  meta: {version: string, generatedOn: string, generatedBy?: string | null, hash: string}): string {
  const lines = [`# ${evidenceTableTitle(i.criterion.name)}`, '',
    `Program ${i.program.code}${i.program.name && i.program.name !== i.program.code ? ` (${i.program.name})` : ''}` +
    ` · criterion ${i.criterion.code} of the ${i.criterion.stage} gate · version ${meta.version}`, '',
    `Generated ${meta.generatedOn.slice(0, 16).replace('T', ' ')} UTC${meta.generatedBy ? ` by ${meta.generatedBy}` : ''}` +
    ` from evidence ${meta.hash}. Generated by the Project Portal; regenerate it rather than edit it.`];
  const sections = [...new Set(rows.map((r) => r.section))];
  for (const s of sections) {
    const mine = rows.filter((r) => r.section === s);
    lines.push('', `## ${s[0].toUpperCase()}${s.slice(1)}`, '', '| Source | Item | Value | Units | State | On |', '|---|---|---|---|---|---|');
    for (const r of mine)
      lines.push(`| ${[r.source, r.item, r.value, r.units, r.state, r.on.slice(0, 16).replace('T', ' ')].map(mdCell).join(' | ')} |`);
  }
  return lines.join('\n') + '\n';
}

export const evidenceTableTitle = (criterionName: string) => `Evidence table, ${criterionName}`;

/** The files of version [version]: `evidence-process_description-EXB-201-v2-20261002T103105Z.md` and `.csv`. */
export function evidenceFileName(criterionCode: string, programCode: string, version: string, at: string, ext: 'md' | 'csv'): string {
  const safe = (s: string) => s.replace(/[^A-Za-z0-9_.-]+/g, '_');
  const stamp = at.slice(0, 19).replace(/[-:]/g, '');
  return `evidence-${safe(criterionCode)}-${safe(programCode)}-v${safe(version)}-${stamp}Z.${ext}`;
}

/** What a generated document records about how it was made; stored on `document.generated`. */
export interface GeneratedRecord {
  criterion: string;
  criterionCode: string;
  generatedOn: string;
  generatedBy?: string | null;
  hash: string;
  sources: RecordedSource[];
  /** The CSV beside the Markdown document, with its SHA-256. */
  csv?: string | null;
  csvSha256?: string | null;
  rows: number;
}

/** The version a new generation takes: one more than the highest of the criterion's generated
 * documents, counted as whole numbers; `1` for the first. */
export function nextGeneratedVersion(versions: (string | null | undefined)[]): string {
  const n = versions.map((v) => parseInt(String(v ?? ''), 10)).filter((x) => !Number.isNaN(x));
  return String(n.length ? Math.max(...n) + 1 : 1);
}

/** Which older generated documents of the same criterion a new generation supersedes now: drafts and
 * rejected ones. An approved one stays until the new version is approved (directly, or through the
 * approval workflow where the portfolio asks for it), which supersedes it then; a submitted one waits
 * for its review. */
export function supersededOnRegenerate(docs: {id: string, status: string, generated: boolean, criterionId?: string | null}[],
  criterionId: string): string[] {
  return docs.filter((d) => d.generated && d.criterionId === criterionId && (d.status === 'draft' || d.status === 'rejected'))
    .map((d) => d.id);
}

/** One difference between the evidence a document was generated from and the evidence now. */
export interface EvidenceChange {
  key: string;
  kind: EvidenceSourceKind;
  change: 'changed' | 'added' | 'removed';
  text: string;
}

/** What changed since the document was written; empty when nothing did. */
export function compareEvidence(recorded: RecordedSource[], current: RecordedSource[]): EvidenceChange[] {
  const out: EvidenceChange[] = [];
  const now = new Map(current.map((s) => [s.key, s]));
  const then = new Map(recorded.map((s) => [s.key, s]));
  for (const r of recorded) {
    const c = now.get(r.key);
    if (!c) {
      out.push({key: r.key, kind: r.kind, change: 'removed', text: r.kind === 'run' ? `the run selection was removed (was ${r.name})` :
        r.kind === 'snapshot' ? `the snapshot ${r.name} is gone` : r.kind === 'decision' ? `decision ${r.name} no longer cites it` :
          `data condition ${r.key.slice('condition:'.length)} no longer applies`});
    } else if (c.hash !== r.hash) {
      out.push({key: r.key, kind: r.kind, change: 'changed', text: r.kind === 'run' ?
        (c.ref === r.ref ? `run ${c.name} selected again` : `another run selected: ${c.name} (was ${r.name})`) :
        r.kind === 'snapshot' ? (c.ref === r.ref ? `snapshot ${c.name} changed` : `a newer snapshot: ${c.name} (was ${r.name})`) :
          r.kind === 'data_condition' ? `${r.key.slice('condition:'.length)} is now ${(c.state ?? '').replace(/_/g, ' ')}` +
            ` (was ${(r.state ?? '').replace(/_/g, ' ')})` : `decision ${c.name} changed`});
    }
  }
  for (const c of current.filter((x) => !then.has(x.key))) {
    out.push({key: c.key, kind: c.kind, change: 'added', text: c.kind === 'run' ? `a run selected: ${c.name}` :
      c.kind === 'snapshot' ? `a snapshot taken: ${c.name}` : c.kind === 'decision' ? `decision added: ${c.name}` :
        `data condition ${c.key.slice('condition:'.length)} added: ${(c.state ?? '').replace(/_/g, ' ')}`});
  }
  return out;
}

/** `Evidence changed since written: another run selected: Feed 2 (was Feed 4); a snapshot taken: …` */
export const describeChanges = (changes: EvidenceChange[]): string =>
  changes.length ? `Evidence changed since written: ${changes.map((c) => c.text).join('; ')}` : 'Evidence unchanged since written';

/** The last check of a generated document, stored on `document.evidence_check`. */
export interface EvidenceCheck {
  state: 'current' | 'stale';
  /** The evidence hash now. */
  hash: string;
  changes: string[];
  checkedOn: string;
  /** The evidence hash the author and followers were last told about; they are told once per hash. */
  notified?: string | null;
  notifiedOn?: string | null;
}

/** Whether a check should tell people: the document is stale and this state of the evidence was not told yet. */
export const shouldNotify = (check: Pick<EvidenceCheck, 'state' | 'hash'>, previous: Pick<EvidenceCheck, 'notified'> | null | undefined) =>
  check.state === 'stale' && previous?.notified !== check.hash;

/** The notification text for one user: one line per stale document. */
export function staleMessage(portfolio: {code: string, name: string}, lines: string[]): string {
  return [`Evidence changed in ${portfolio.name} (${portfolio.code})`, ...lines].join('\n');
}
