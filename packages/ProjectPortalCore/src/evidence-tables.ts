/** Evidence tables generated into documents against the platform (F-40): reads the evidence behind
 * a document criterion (the run selected for it and its saved call, the newest snapshot and its
 * frozen table, the gate's data conditions as readiness evaluates them, the decisions citing it),
 * writes the table as a CSV and a Markdown file into the program's workspace folder, and links the
 * Markdown file as a new draft version of a document over the criterion, recording which evidence
 * it was generated from. Then (F-41) compares a generated document's record with the evidence now,
 * marks it stale and tells its author and the program's followers. The rules live in
 * ./model/evidence-tables. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {DocumentRow, portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {ContextCriterion, ProgramContext, programContext, programData, readinessFor} from './portal';
import {programAccess} from './program-roles';
import {selectedRunsOf} from './evidence';
import {frozenCriteria, sha256, snapshotsOf, snapshotTable} from './snapshots';
import {documentsOf, linkDocumentTo, uploadDocument, vocabularyOf} from './documents';
import {checkDocumentChange} from './model/documents';
import {DataRow} from './model/conditions';
import {
  evidenceCsv, evidenceFileName, evidenceHash, EvidenceInputs, evidenceMarkdown, EvidenceRow, evidenceSources,
  evidenceTableRows, evidenceTableTitle, EVIDENCE_COLUMNS, GeneratedRecord, MAX_SNAPSHOT_ROWS, nextGeneratedVersion,
  RunParam, supersededOnRegenerate, compareEvidence, describeChanges, EvidenceChange, EvidenceCheck, RecordedSource,
  shouldNotify, staleMessage,
} from './model/evidence-tables';
import {notifyUser} from './notify';
import {filterValue, ROW_LIMIT} from './rows';

/** The kind a generated table takes: `evidence_table` when the portfolio declares it, else `other`. */
export const EVIDENCE_KIND = 'evidence_table';

/** What a document records when it is a generated evidence table; null for a document a person wrote. */
export const generatedOf = (d: Pick<DocumentRow, 'generated'>): GeneratedRecord | null =>
  (d.generated as GeneratedRecord | undefined)?.hash ? d.generated as GeneratedRecord : null;

/** The criterion by id or code; it must apply to the program and be of kind document. */
export function documentCriterion(ctx: Pick<ProgramContext, 'criteria' | 'program'>, ref: string): ContextCriterion {
  const c = ctx.criteria.find((x) => x.id === ref || x.code === ref);
  if (!c)
    throw new Error(`no criterion ${ref} applies to ${ctx.program.code}`);
  if (c.kind !== 'document')
    throw new Error(`criterion ${c.code} is a ${c.kind.replace('_', ' ')}; evidence tables are generated for document criteria`);
  return c;
}

/** A saved call's parameters with their values; a table parameter shows as `saved table`. */
async function runParams(runId: string): Promise<{inputs: RunParam[], outputs: RunParam[]}> {
  const call = await grok.dapi.functions.calls.include('inputs, outputs').find(runId).catch(() => null) as DG.FuncCall | null;
  if (!call)
    return {inputs: [], outputs: []};
  const read = (which: 'inputs' | 'outputs'): RunParam[] => call.func[which].map((p) => {
    const v = call[which][p.name];
    return {name: p.name, caption: p.caption || null, units: p.units || p.options?.['units'] || null,
      value: p.propertyType === DG.TYPE.DATA_FRAME && !(v instanceof DG.DataFrame) ? (v ? 'saved table' : null) : v};
  });
  return {inputs: read('inputs'), outputs: read('outputs')};
}

/** The evidence of one document criterion now. With [values], the run's parameters and the
 * snapshot's frozen rows are read too (generating); without, only what the sources' hashes need. */
export async function evidenceInputsOf(ctx: ProgramContext, c: ContextCriterion,
  opts: {values?: boolean, rows?: DataRow[] | null, docs?: DocumentRow[]} = {}): Promise<EvidenceInputs> {
  const {program} = ctx;
  const stage = ctx.stages.find((s) => s.id === c.stageId);
  const selected = (await selectedRunsOf(program.id).catch(() => [])).find((r) => r.criterionId === c.id) ?? null;
  const run = selected ? {runId: selected.runId, title: selected.title, model: selected.model, selectedOn: selected.selectedOn ?? null,
    selectedBy: selected.selectedBy ?? null, rationale: selected.rationale ?? null,
    ...(opts.values ? await runParams(selected.runId) : {inputs: [], outputs: []})} : null;

  const newest = (await snapshotsOf(program.id).catch(() => []))[0] ?? null;
  let values: Record<string, unknown>[] | undefined;
  if (newest && opts.values) {
    const t = await snapshotTable(newest).catch(() => null);
    if (t?.df) {
      values = [];
      for (let i = 0; i < Math.min(t.df.rowCount, MAX_SNAPSHOT_ROWS); i++)
        values.push(Object.fromEntries(t.df.columns.toList().map((col) => [col.name, col.get(i)])));
    }
  }
  const snapshot = newest ? {id: newest.id, label: newest.label,
    takenOn: newest.taken_on ? dayjs(newest.taken_on as any).toISOString() : '', sha256: newest.table_sha256 ?? null,
    rows: newest.table_rows ?? null, criterionState: frozenCriteria(newest).find((x) => x.code === c.code)?.state ?? null,
    ...(values ? {values} : {})} : null;

  // the gate's data conditions, as readiness evaluates them now
  const gate = ctx.criteria.filter((x) => x.stageId === c.stageId && x.kind === 'data_condition');
  let conditions: EvidenceInputs['conditions'] = [];
  if (gate.length) {
    let readiness = ctx.readiness?.stageId === c.stageId ? ctx.readiness : null;
    if (!readiness) {
      let rows = opts.rows ?? null;
      let error: string | undefined;
      if (rows === null && gate.some((x) => x.condition)) {
        try {
          rows = await programData(ctx.portfolio, program.code);
        } catch (e: any) {
          error = e?.message ?? String(e);
        }
      }
      readiness = readinessFor(ctx, c.stageId, rows, error);
    }
    conditions = readiness.criteria.filter((x) => x.kind === 'data_condition')
      .map((x) => ({id: x.id, code: x.code, name: x.name, state: x.state, summary: x.result?.summary ?? x.gap ?? null}));
  }

  // decisions citing the criterion, its run, or a document over it that a person wrote
  const docs = opts.docs ?? await documentsOf(program.id).catch(() => [] as DocumentRow[]);
  const docIds = new Set(docs.filter((d) => d.gate_criterion_id === c.id && !generatedOf(d)).map((d) => d.id));
  const evidence = await portalDb.evidences.query({filter: `program_id = ${filterValue(program.id)}`, limit: ROW_LIMIT}).catch(() => []);
  const citing = new Set(evidence.filter((e) => (e.kind === 'criterion' && e.target === c.id) ||
    (e.kind === 'model_run' && run && e.target === run.runId) ||
    (e.kind === 'file' && docIds.has(String((e.context as any)?.document ?? ''))))
    .map((e) => e.stage_transition_id).filter((x): x is string => !!x));
  const transitions = citing.size ? await portalDb.stageTransitions.query({filter: `program_id = ${filterValue(program.id)}`, limit: ROW_LIMIT}) : [];
  const code = (id?: string | null) => id ? ctx.stages.find((s) => s.id === id)?.code ?? null : null;
  const decisions = transitions.filter((t) => citing.has(t.id) && t.decision !== 'sync')
    .map((t) => ({id: t.id, decision: t.decision, on: dayjs((t.decided_on ?? t.created_on) as any).toISOString(),
      from: code(t.from_stage_id), to: code(t.to_stage_id), rationale: t.rationale ?? null}))
    .sort((a, b) => a.on.localeCompare(b.on));

  return {program: {code: program.code, name: program.name}, criterion: {id: c.id, code: c.code, name: c.name,
    stage: stage?.name ?? ''}, run, snapshot, conditions, decisions};
}

export function evidenceFrame(rows: EvidenceRow[], name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns(EVIDENCE_COLUMNS.map((c) => DG.Column.fromStrings(c, rows.map((r) => r[c]))));
  df.name = name;
  return df;
}

export interface GeneratedTable {
  document: DocumentRow;
  df: DG.DataFrame;
  version: string;
  /** Older generated drafts of the criterion this one superseded. */
  superseded: number;
  notes: string[];
}

/** Generates the evidence table of [criterion] (id or code) for [program]: the CSV and the Markdown
 * file in the program's workspace folder, the Markdown linked as a new draft version over the
 * criterion, older generated drafts superseded. Contributors and approvers; needs the folder. */
export async function generateEvidenceTableFor(program: ProgramRow, criterion: string): Promise<GeneratedTable> {
  const ctx = await programContext(program.id);
  const c = documentCriterion(ctx, criterion);
  const access = await programAccess(ctx.portfolio, program.id);
  const refused = checkDocumentChange({isApprover: access.approver, isContributor: access.contributor, status: 'draft',
    programStatus: program.status});
  if (refused.length)
    throw new Error(refused.join('; '));
  const docs = await documentsOf(program.id);
  const mine = docs.filter((d) => generatedOf(d)?.criterion === c.id);
  const version = nextGeneratedVersion(mine.map((d) => d.doc_version));
  const inputs = await evidenceInputsOf(ctx, c, {values: true, docs});
  const sources = evidenceSources(inputs);
  const hash = evidenceHash(sources);
  const rows = evidenceTableRows(inputs);
  const at = new Date().toISOString();
  const me = DG.User.current();
  const title = evidenceTableTitle(c.name);
  const csv = evidenceCsv(rows);
  const md = evidenceMarkdown(inputs, rows, {version, generatedOn: at, generatedBy: me.friendlyName, hash});

  const csvFile = await uploadDocument(program, evidenceFileName(c.code, program.code, version, at, 'csv'), csv);
  let mdFile: {path: string, folder: string} | null = null;
  const removeFiles = async () => {
    for (const f of [csvFile, mdFile].filter((x): x is {path: string, folder: string} => !!x))
      await grok.dapi.spaces.id(f.folder).files.delete(f.path.split('/').pop()!).catch(() => undefined);
  };
  let document: DocumentRow;
  const notes: string[] = [];
  try {
    mdFile = await uploadDocument(program, evidenceFileName(c.code, program.code, version, at, 'md'), md);
    const generated: GeneratedRecord = {criterion: c.id, criterionCode: c.code, generatedOn: at, generatedBy: me.id, hash, sources,
      csv: csvFile.path, csvSha256: await sha256(csv), rows: rows.length};
    const kind = vocabularyOf(ctx).some((k) => k.code === EVIDENCE_KIND) ? EVIDENCE_KIND : 'other';
    const linked = await linkDocumentTo(program, {path: mdFile.path, folder: mdFile.folder, criterionId: c.id, kind, title,
      version, status: 'draft', generated: generated as any});
    document = linked.document;
    notes.push(...linked.notes);
  } catch (e) {
    // a table nobody points at would only confuse the folder's readers
    await removeFiles();
    throw e;
  }
  let superseded = 0;
  const older = supersededOnRegenerate(mine.map((d) => ({id: d.id, status: d.status, generated: true,
    criterionId: generatedOf(d)!.criterion})), c.id);
  for (const id of older) {
    await portalDb.documents.update(id, {status: 'superseded'});
    superseded++;
  }
  const still = mine.filter((d) => d.status === 'approved' || d.status === 'submitted');
  if (still.length)
    notes.push(`v${still.map((d) => d.doc_version).join(', v')} stays ${still[0].status} until this version is approved`);
  const df = evidenceFrame(rows, `${program.code} · ${title} v${version}`);
  df.setTag('portal.document', document.id);
  return {document, df, version, superseded, notes};
}

/** The table of a generated document, read back from its CSV. */
export async function evidenceTableOf(d: DocumentRow): Promise<DG.DataFrame> {
  const g = generatedOf(d);
  if (!g?.csv)
    throw new Error(`"${d.title}" is not a generated evidence table`);
  const text = await grok.dapi.files.readAsText(g.csv).catch(() => null);
  if (text === null)
    throw new Error(`the table of "${d.title}" is gone: ${g.csv}`);
  const df = DG.DataFrame.fromCsv(text);
  df.name = `${d.title}${d.doc_version ? ` v${d.doc_version}` : ''}`;
  df.setTag('portal.document', d.id);
  return df;
}

export function summarizeGenerated(program: ProgramRow, r: GeneratedTable): string {
  return [`${program.code}: ${r.document.title} v${r.version} generated as a draft, ${r.df.rowCount} rows` +
    (r.superseded ? `; ${r.superseded} older draft(s) superseded` : ''), ...r.notes.map((n) => `  ${n}`)].join('\n');
}

// ---------------------------------------------------------------------------------------------
// Evidence changed since written (F-41)

/** The last check stored on a generated document, if any. */
export const checkOf = (d: Pick<DocumentRow, 'evidence_check'>): EvidenceCheck | null =>
  (d.evidence_check as EvidenceCheck | undefined)?.state ? d.evidence_check as EvidenceCheck : null;

/** What changed since each generated document of the program was written, by document id, with the
 * evidence hash now. Superseded documents are not checked. Each criterion's evidence is read once. */
export async function evidenceChangesOf(ctx: ProgramContext, docs?: DocumentRow[]):
  Promise<Map<string, {changes: EvidenceChange[], hash: string}>> {
  const rows = docs ?? await documentsOf(ctx.program.id);
  const out = new Map<string, {changes: EvidenceChange[], hash: string}>();
  const now = new Map<string, RecordedSource[] | null>();
  for (const d of rows) {
    const g = generatedOf(d);
    if (!g || d.status === 'superseded')
      continue;
    if (!now.has(g.criterion)) {
      const c = ctx.criteria.find((x) => x.id === g.criterion);
      now.set(g.criterion, c ? evidenceSources(await evidenceInputsOf(ctx, c, {docs: rows})) : null);
    }
    const current = now.get(g.criterion);
    // a criterion that no longer applies: nothing to compare with
    if (current)
      out.set(d.id, {changes: compareEvidence(g.sources, current), hash: evidenceHash(current)});
  }
  return out;
}

/** The check of one document as it is stored; null when the stored one already says the same. */
export function nextCheck(d: DocumentRow, r: {changes: EvidenceChange[], hash: string}): EvidenceCheck | null {
  const before = checkOf(d);
  const state = r.changes.length ? 'stale' : 'current';
  if (before && before.state === state && before.hash === r.hash)
    return null;
  return {state, hash: r.hash, changes: r.changes.map((c) => c.text), checkedOn: new Date().toISOString(),
    notified: before?.notified ?? null, notifiedOn: before?.notifiedOn ?? null};
}

/** Stores what the program page found, without telling anyone: the feed then shows a stale document
 * before the next sync. A user who may not write documents leaves it to the sync. */
export async function storeChecks(docs: DocumentRow[], found: Map<string, {changes: EvidenceChange[], hash: string}>): Promise<void> {
  for (const d of docs) {
    const r = found.get(d.id);
    const check = r ? nextCheck(d, r) : null;
    if (check)
      await portalDb.documents.update(d.id, {evidence_check: check as any}).catch(() => undefined);
  }
}

export interface EvidenceChangeCheck {
  /** Generated documents checked. */
  documents: number;
  stale: {program: string, document: string, line: string}[];
  notified: string[];
  errors: string[];
  summary: string;
}

/** Checks every generated document of a portfolio's programs against the evidence now, stores the
 * result on each, and tells each stale document's author and the program's followers once per state
 * of the evidence, in the platform's notifications. */
export async function checkEvidenceChangesOf(portfolio: PortfolioRow, opts: {notify?: boolean} = {}): Promise<EvidenceChangeCheck> {
  const programs = await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)}`, limit: ROW_LIMIT});
  const out: EvidenceChangeCheck = {documents: 0, stale: [], notified: [], errors: [], summary: ''};
  const ids = programs.map((p) => `"${p.id}"`).join(', ');
  const all = ids ? await portalDb.documents.query({filter: `program_id in (${ids})`, limit: ROW_LIMIT}) : [];
  const docs = all.filter((d) => generatedOf(d) && d.status !== 'superseded');
  const tell = new Map<string, string[]>();
  for (const program of programs) {
    const mine = docs.filter((d) => d.program_id === program.id);
    if (!mine.length)
      continue;
    try {
      const ctx = await programContext(program.id);
      // every document of the program: the decisions citing a person-written one count as evidence too
      const found = await evidenceChangesOf(ctx, all.filter((d) => d.program_id === program.id));
      const followers = opts.notify === false ? [] : (await portalDb.followings.query({filter: `program_id = ${filterValue(program.id)}`, limit: ROW_LIMIT}))
        .filter((f) => f.active !== false).map((f) => f.user_id);
      for (const d of mine) {
        const r = found.get(d.id);
        if (!r)
          continue;
        out.documents++;
        let check = nextCheck(d, r);
        const state = check ?? checkOf(d)!;
        if (state.state === 'stale') {
          const line = `${program.code}: ${d.title}${d.doc_version ? ` v${d.doc_version}` : ''}: ${describeChanges(r.changes)}`;
          out.stale.push({program: program.code, document: d.id, line});
          if (opts.notify !== false && shouldNotify(state, checkOf(d))) {
            for (const u of new Set([d.written_by, ...followers].filter((x): x is string => !!x)))
              tell.set(u, [...tell.get(u) ?? [], line]);
            check = {...state, notified: state.hash, notifiedOn: new Date().toISOString()};
          }
        }
        if (check)
          await portalDb.documents.update(d.id, {evidence_check: check as any});
      }
    } catch (e: any) {
      out.errors.push(`${program.code} (${e?.message ?? e})`);
    }
  }
  for (const [u, lines] of tell) {
    const r = await notifyUser(u, staleMessage(portfolio, lines));
    const who = await grok.dapi.users.find(u).then((x) => x.login).catch(() => u);
    out.notified.push(`${who}: ${r.detail}`);
  }
  out.summary = `${portfolio.code} evidence: ` + (out.documents === 0 ? 'no generated documents' :
    out.stale.length ? out.stale.map((s) => s.line).join('; ') : `${out.documents} generated document(s), none stale`) +
    (out.notified.length ? `; told: ${out.notified.join(', ')}` : '') +
    (out.errors.length ? `; not checked: ${out.errors.join(', ')}` : '');
  return out;
}
