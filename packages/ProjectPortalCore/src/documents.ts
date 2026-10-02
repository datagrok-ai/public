/** Documents over files against the platform (F-27): linking a file as a document of a program,
 * uploading one into the program's workspace folder, approving it (which meets its criterion
 * unless waived, and supersedes the older approved version), the document cards of the program
 * page, and the example documents of EX-BIOC. The rules live in ./model/documents. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {DocumentRow, portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {ProgramContext, programContext, portfolioByCode} from './portal';
import {programAccess} from './program-roles';
import {existingSpaces} from './workspace';
// evidence.ts imports from here too; the cycle is safe because neither uses the other at load
import {selectedRunsOf} from './evidence';
import {CriterionState, nextStage} from './model/gates';
import {
  checkDocumentChange, criterionAfter, defaultKindFor, describeDocument, DocumentCard, documentCards, DocumentInput,
  DocumentKindConfig, DocumentRef, DocumentStatus, documentVocabulary, supersededBy, validateDocument,
} from './model/documents';

/** Where the package's own files are read from by path. */
export const PACKAGE_FILES = 'System:AppData/ProjectPortalCore/';

export const toRef = (d: DocumentRow): DocumentRef => ({id: d.id, kind: d.kind, title: d.title, version: d.doc_version ?? null,
  status: d.status as DocumentStatus, criterionId: d.gate_criterion_id ?? null, path: d.path, author: d.written_by ?? null,
  on: (d.approved_on ?? d.linked_on ?? d.created_on) ? dayjs((d.approved_on ?? d.linked_on ?? d.created_on) as any).toISOString() : null});

/** The kinds the portfolio's documents take: its `documentKinds`, else one per document criterion. */
export function vocabularyOf(ctx: Pick<ProgramContext, 'portfolio' | 'criteria'>): DocumentKindConfig[] {
  const declared = (ctx.portfolio.document_kinds as {kinds?: DocumentKindConfig[]} | undefined)?.kinds;
  const seen = new Set<string>();
  const criteria = ctx.criteria.filter((c) => c.kind === 'document' && !seen.has(c.code) && !!seen.add(c.code));
  return documentVocabulary(declared, criteria);
}

export async function documentsOf(programId: string): Promise<DocumentRow[]> {
  return (await portalDb.documents.query({filter: `program_id = "${programId}"`}))
    .sort((a, b) => String(b.created_on).localeCompare(String(a.created_on)));
}

export interface LinkedDocument {
  document: DocumentRow;
  /** Linked new, or an existing document of the same file changed. */
  created: boolean;
  /** What the criterion became, when the document changed it. */
  criterion: {code: string, state: CriterionState} | null;
  superseded: number;
  notes: string[];
}

/** Links a file to [program] as a document, or changes the document already over that file.
 * Contributors and approvers link drafts; approvers approve. Approving meets the criterion unless
 * it was waived, and supersedes the older approved document of the same kind and criterion. */
export async function linkDocumentTo(program: ProgramRow, input: Partial<DocumentInput> & {path: string,
  folder?: string | null, note?: string | null, criterion?: string | null}): Promise<LinkedDocument> {
  const ctx = await programContext(program.id);
  const access = await programAccess(ctx.portfolio, program.id);
  const vocabulary = vocabularyOf(ctx);
  const path = input.path.trim();
  const existing = (await documentsOf(program.id)).find((d) => d.path === path) ?? null;
  // the criterion by id or code; one that does not apply to the program is refused
  const named = input.criterionId ?? input.criterion ?? (existing ? existing.gate_criterion_id ?? null : null);
  const criterion = named ? ctx.criteria.find((c) => c.id === named || c.code === named) ?? null : undefined;
  const status: DocumentStatus = input.status ?? (existing?.status as DocumentStatus | undefined) ?? 'draft';
  const doc: DocumentInput = {path, criterionId: criterion?.id ?? null, status,
    kind: input.kind || existing?.kind || defaultKindFor(criterion?.code, vocabulary),
    title: input.title?.trim() || existing?.title || path.split(/[/\\:]/).filter((s) => s).pop() || path,
    version: input.version?.trim() || existing?.doc_version || undefined};
  const refused = [...checkDocumentChange({isApprover: access.approver, isContributor: access.contributor, status,
    programStatus: program.status}), ...validateDocument(doc, vocabulary, named ? criterion ?? null : undefined)];
  if (refused.length)
    throw new Error(refused.join('; '));
  const notes: string[] = [];
  const exists = await grok.dapi.files.exists(path).catch(() => undefined);
  if (exists === false)
    throw new Error(`no file at ${path}`);
  if (exists === undefined)
    notes.push(`${path} could not be checked`);

  const me = DG.User.current().id;
  const now = dayjs();
  const approving = status === 'approved' && existing?.status !== 'approved';
  const values = {kind: doc.kind, title: doc.title, doc_version: doc.version ?? null, status,
    gate_criterion_id: doc.criterionId ?? null, folder: input.folder ?? existing?.folder ?? null,
    note: input.note ?? existing?.note ?? null,
    ...(approving ? {approved_by: me, approved_on: now} : {})};
  let row: DocumentRow;
  if (existing) {
    await portalDb.documents.update(existing.id, values as any);
    row = await portalDb.documents.get(existing.id);
  } else {
    const [r] = await portalDb.documents.insert({program_id: program.id, path, written_by: me, linked_by: me,
      linked_on: now as any, ...values} as any);
    row = await portalDb.documents.get(r.id);
  }
  // a newer approved version supersedes the older
  let superseded = 0;
  if (status === 'approved') {
    for (const id of supersededBy((await documentsOf(program.id)).map(toRef), toRef(row))) {
      await portalDb.documents.update(id, {status: 'superseded'});
      superseded++;
    }
  }
  let changed: LinkedDocument['criterion'] = null;
  if (criterion) {
    const next = criterionAfter(ctx.recorded[criterion.id], status);
    if (next) {
      await portalDb.criterionStatuses.upsert({program_id: program.id, gate_criterion_id: criterion.id, status: next,
        evidence: `document "${describeDocument(toRef(row))}" ${status === 'approved' ? 'approved' : 'linked'} ${now.format('YYYY-MM-DD')}`});
      changed = {code: criterion.code, state: next};
    } else if (ctx.recorded[criterion.id] === 'waived' && status === 'approved')
      notes.push(`${criterion.code} stays waived`);
  }
  return {document: row, created: !existing, criterion: changed, superseded, notes};
}

/** Changes a document's status (approve, mark superseded), with the same checks as linking. */
export async function setDocumentStatus(program: ProgramRow, documentId: string, status: DocumentStatus): Promise<LinkedDocument> {
  const d = await portalDb.documents.get(documentId);
  if (d.program_id !== program.id)
    throw new Error(`${program.code} has no document ${documentId}`);
  return linkDocumentTo(program, {path: d.path, status});
}

/** Whether the program's workspace folder exists. */
async function liveFolder(program: ProgramRow): Promise<string | null> {
  return program.space_id && (await existingSpaces([program.space_id]).catch(() => new Set<string>())).has(program.space_id) ?
    program.space_id : null;
}

/** The files at the top of the program's workspace folder: name and full path. */
export async function folderFiles(program: ProgramRow): Promise<{name: string, path: string}[]> {
  const folder = await liveFolder(program);
  if (!folder)
    return [];
  const files = await grok.dapi.spaces.id(folder).children.ofTypes('FileInfo', false).list({pageSize: 500}).catch(() => []);
  return (files as any[]).filter((f) => f instanceof DG.FileInfo && !f.isDirectory)
    .map((f: DG.FileInfo) => ({name: f.name, path: f.fullPath})).sort((a, b) => a.name.localeCompare(b.name));
}

/** Writes a file into the program's workspace folder and returns its full path. An existing file
 * of the same name is not overwritten: a new version is a new file. */
export async function uploadDocument(program: ProgramRow, name: string, bytes: Uint8Array | string):
  Promise<{path: string, folder: string}> {
  const folder = await liveFolder(program);
  if (!folder)
    throw new Error(`${program.code} has no workspace folder: build it first, or link a file by its path`);
  const file = name.trim().replace(/[\\/:]+/g, '_');
  if (!file)
    throw new Error('the file needs a name');
  const space = grok.dapi.spaces.id(folder);
  if (await space.files.exists(file).catch(() => false))
    throw new Error(`${file} is already in ${program.code}'s folder: name the new version differently, or link that file`);
  if (typeof bytes === 'string')
    await space.files.writeString(file, bytes);
  else
    await space.files.write(file, bytes);
  const found = (await folderFiles(program)).find((f) => f.name === file);
  if (!found)
    throw new Error(`${file} was written but cannot be found in ${program.code}'s folder`);
  return {path: found.path, folder};
}

/** The decisions that cite each document, by document id: `advance 2026-10-02`. */
async function citationsOf(programId: string): Promise<Record<string, string[]>> {
  const evidence = (await portalDb.evidences.query({filter: `program_id = "${programId}" and kind = "file"`}).catch(() => []))
    .filter((e) => (e.context as any)?.document);
  if (!evidence.length)
    return {};
  const transitions = await portalDb.stageTransitions.query({filter: `program_id = "${programId}"`});
  const out: Record<string, string[]> = {};
  for (const e of evidence) {
    const t = transitions.find((x) => x.id === e.stage_transition_id);
    if (t) {
      const id = String((e.context as any).document);
      out[id] = [...out[id] ?? [], `${t.decision} ${(t.decided_on ?? t.created_on)?.format?.('YYYY-MM-DD') ?? ''}`];
    }
  }
  return out;
}

/** The document cards of the program page: the document criteria of the current and the next gate. */
export async function documentCardsOf(ctx: ProgramContext, docs?: DocumentRow[]): Promise<DocumentCard[]> {
  const stages = ctx.stages.map((s) => ({id: s.id, code: s.code, name: s.name, ordinal: s.ordinal, isTerminal: s.is_terminal}));
  const current = stages.find((s) => s.id === ctx.program.stage_id) ?? null;
  const next = nextStage(stages, ctx.program.stage_id ?? null);
  const rows = docs ?? await documentsOf(ctx.program.id);
  const runs = Object.fromEntries((await selectedRunsOf(ctx.program.id).catch(() => [])).map((r) => [r.criterionId, r.title]));
  return documentCards(ctx.criteria, {current, next}, rows.map(toRef), ctx.recorded, runs,
    await citationsOf(ctx.program.id).catch(() => ({})));
}

export function documentsFrame(rows: DocumentRow[], criteria: {id: string, code: string}[], name: string): DG.DataFrame {
  const code = (id?: string) => criteria.find((c) => c.id === id)?.code ?? '';
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('id', rows.map((r) => r.id)),
    DG.Column.fromStrings('criterion', rows.map((r) => code(r.gate_criterion_id))),
    DG.Column.fromStrings('kind', rows.map((r) => r.kind)),
    DG.Column.fromStrings('title', rows.map((r) => r.title)),
    DG.Column.fromStrings('version', rows.map((r) => r.doc_version ?? '')),
    DG.Column.fromStrings('status', rows.map((r) => r.status)),
    DG.Column.fromStrings('path', rows.map((r) => r.path)),
    DG.Column.fromStrings('linked_on', rows.map((r) => r.linked_on ? dayjs(r.linked_on as any).toISOString() : '')),
    DG.Column.fromStrings('approved_on', rows.map((r) => r.approved_on ? dayjs(r.approved_on as any).toISOString() : '')),
  ]);
  df.name = name;
  return df;
}

/** One example document, as `files/documents/examples.json` lists it. */
interface ExampleDocument {
  portfolio: string;
  program: string;
  criterion?: string;
  kind?: string;
  title: string;
  version?: string;
  status?: DocumentStatus;
  /** The package file, under files/. */
  file: string;
}

/** Links the example documents the package ships to the EX-* programs they name. Only portfolios
 * whose code starts with EX- are touched. Idempotent: a document already linked is left as it is. */
export async function loadExampleDocuments(pkg: DG.Package): Promise<string> {
  const list: ExampleDocument[] = JSON.parse(await pkg.files.readAsText('documents/examples.json'));
  const lines: string[] = [];
  const portfolios = new Map<string, PortfolioRow>();
  for (const e of list) {
    if (!e.portfolio.startsWith('EX-')) {
      lines.push(`${e.portfolio} ${e.program}: skipped, examples go to EX-* portfolios only`);
      continue;
    }
    try {
      if (!portfolios.has(e.portfolio))
        portfolios.set(e.portfolio, await portfolioByCode(e.portfolio));
      const pf = portfolios.get(e.portfolio)!;
      const program = (await portalDb.programs.query({filter: `portfolio_id = "${pf.id}" and code = "${e.program}"`}))[0];
      if (!program) {
        lines.push(`${e.portfolio} ${e.program}: no such program (sync the portfolio first)`);
        continue;
      }
      const path = `${PACKAGE_FILES}${e.file}`;
      if ((await documentsOf(program.id)).some((d) => d.path === path)) {
        lines.push(`${e.program}: ${e.title} already linked`);
        continue;
      }
      const r = await linkDocumentTo(program, {path, criterion: e.criterion ?? null, kind: e.kind, title: e.title,
        version: e.version, status: e.status ?? 'draft'});
      lines.push(`${e.program}: ${describeDocument(toRef(r.document))} linked${r.criterion ? `; ${r.criterion.code} ${r.criterion.state}` : ''}`);
    } catch (err: any) {
      lines.push(`${e.program}: ${e.title} not linked: ${err?.message ?? err}`);
    }
  }
  return lines.join('\n');
}
