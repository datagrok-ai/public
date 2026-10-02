/** Documents over files (F-27, FR-70). A document is a record over one file, in the program's
 * workspace folder or at any file path the platform reads: its kind (from the portfolio's
 * vocabulary), title, version, status and author, and the gate criterion of kind `document` it
 * answers. Approving a document meets its criterion unless the criterion was waived; a newer
 * approved version of the same kind for the same criterion supersedes the older one.
 *
 * The document cards on the program page: for each document criterion of the current and the next
 * gate, what exists, what is missing, its state, and the evidence linked to it (the run selected
 * for the criterion, the decisions that cite a document). Pure: the caller reads the records. */
import type {CriterionState} from './gates';
import type {EvidenceItem} from './evidence';

export type DocumentStatus = 'draft' | 'submitted' | 'approved' | 'rejected' | 'superseded';
/** The statuses a document is linked or set to directly. `submitted` and `rejected` come only from
 * the approval workflow (F-30), in portfolios that ask for approval of documents. */
export const DOCUMENT_STATUSES: DocumentStatus[] = ['draft', 'approved', 'superseded'];
export const REVIEW_STATUSES: DocumentStatus[] = ['submitted', 'rejected'];

/** One kind of document a portfolio knows: `process_description`, `ind_quality`. */
export interface DocumentKindConfig {
  code: string;
  name: string;
  description?: string;
}

/** The kind every vocabulary ends with, for a document that fits none of the others. */
export const OTHER_KIND: DocumentKindConfig = {code: 'other', name: 'Other'};

const CODE = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;

/** Problems with a portfolio's `documentKinds`, as sentences. */
export function validateDocumentKinds(kinds: DocumentKindConfig[] | undefined): string[] {
  if (kinds === undefined)
    return [];
  if (!Array.isArray(kinds))
    return ['documentKinds must be a list'];
  const errors: string[] = [];
  const seen = new Set<string>();
  for (const k of kinds) {
    if (!k?.code || !CODE.test(k.code))
      errors.push(`documentKinds: code "${k?.code ?? ''}" may hold only letters, digits, '_', '.', '-'`);
    else if (seen.has(k.code))
      errors.push(`documentKinds: "${k.code}" is declared twice`);
    else if (!k.name?.trim())
      errors.push(`documentKinds: ${k.code} has no name`);
    seen.add(k?.code);
  }
  return errors;
}

/** The kinds a portfolio's documents take: its `documentKinds` when declared, else one kind per
 * document criterion of its stage model (code and name of the criterion); always with `other`. */
export function documentVocabulary(declared: DocumentKindConfig[] | null | undefined,
  documentCriteria: {code: string, name: string}[]): DocumentKindConfig[] {
  const base = declared?.length ? declared : documentCriteria.map((c) => ({code: c.code, name: c.name}));
  const out: DocumentKindConfig[] = [];
  for (const k of [...base, OTHER_KIND]) {
    if (!out.some((x) => x.code === k.code))
      out.push(k);
  }
  return out;
}

/** The kind a document for [criterionCode] gets when none is given: the kind of that code, else `other`. */
export const defaultKindFor = (criterionCode: string | null | undefined, vocabulary: DocumentKindConfig[]): string =>
  vocabulary.find((k) => k.code === criterionCode)?.code ?? OTHER_KIND.code;

/** A document as it is linked or changed. */
export interface DocumentInput {
  /** The file's full path as the platform reads it: `System:AppData/...`, `Space:Program:Files/x.pdf`. */
  path: string;
  kind: string;
  title: string;
  version?: string;
  status?: DocumentStatus;
  /** The gate criterion it answers: the criterion's id, which must be of kind document. */
  criterionId?: string | null;
}

/** Problems with one document, as sentences. [criterion] is the criterion it names, when it names one. */
export function validateDocument(d: Partial<DocumentInput>, vocabulary: DocumentKindConfig[],
  criterion?: {code: string, kind: string} | null): string[] {
  const errors: string[] = [];
  if (!d.path?.trim())
    errors.push('a document needs a file');
  if (!d.title?.trim())
    errors.push('a document needs a title');
  if (!d.kind || !vocabulary.some((k) => k.code === d.kind))
    errors.push(`document kind "${d.kind ?? ''}" is not one of ${vocabulary.map((k) => k.code).join(', ')}`);
  if (d.status && !DOCUMENT_STATUSES.includes(d.status))
    errors.push(`document status "${d.status}" is not one of ${DOCUMENT_STATUSES.join(', ')}`);
  if (d.criterionId && criterion === null)
    errors.push('the criterion does not apply to this program');
  else if (criterion && criterion.kind !== 'document')
    errors.push(`criterion ${criterion.code} is a ${criterion.kind.replace('_', ' ')}; documents answer document criteria`);
  return errors;
}

/** Who may do what: contributors and approvers link drafts; approving (which meets a criterion)
 * and marking superseded are for approvers. A program that ended takes no documents. */
export function checkDocumentChange(c: {isApprover: boolean, isContributor: boolean, status: DocumentStatus,
  programStatus?: string | null}): string[] {
  const refused: string[] = [];
  if (c.programStatus && ['terminated', 'completed', 'transitioned'].includes(c.programStatus))
    refused.push(`the program is ${c.programStatus}`);
  if (!c.isApprover && !c.isContributor)
    refused.push('Only the portfolio\'s or the program\'s contributors and approvers link documents');
  else if (c.status !== 'draft' && !c.isApprover)
    refused.push(`Only approvers mark a document ${c.status}`);
  return refused;
}

/** What a document's status does to its criterion: approval meets it, a draft puts a criterion
 * not started in progress. A waiver stands; a met criterion stays met. Null: leave it. */
export function criterionAfter(current: CriterionState | null | undefined, status: DocumentStatus): CriterionState | null {
  if (current === 'waived')
    return null;
  if (status === 'approved')
    return current === 'met' ? null : 'met';
  if ((status === 'draft' || status === 'submitted') && (!current || current === 'not_started'))
    return 'in_progress';
  return null;
}

/** A document as the rules see it. */
export interface DocumentRef {
  id: string;
  kind: string;
  title: string;
  version?: string | null;
  status: DocumentStatus;
  criterionId?: string | null;
  path: string;
  author?: string | null;
  /** When it was linked or last changed, ISO; ordering only. */
  on?: string | null;
}

/** Versions compare as numbers where they are numbers: 1.10 after 1.9. */
export const compareVersions = (a?: string | null, b?: string | null): number =>
  String(a ?? '').localeCompare(String(b ?? ''), undefined, {numeric: true});

/** The approved documents [approved] supersedes: the others of its kind for the same criterion. */
export function supersededBy(docs: DocumentRef[], approved: DocumentRef): string[] {
  return docs.filter((d) => d.id !== approved.id && d.status === 'approved' && d.kind === approved.kind &&
    (d.criterionId ?? null) === (approved.criterionId ?? null)).map((d) => d.id);
}

/** Approved first, then submitted, drafts, rejected and superseded; newest version first within each. */
export function orderDocuments<T extends Pick<DocumentRef, 'status' | 'version' | 'title' | 'on'>>(docs: T[]): T[] {
  const rank = (d: T) => ['approved', 'submitted', 'draft', 'rejected', 'superseded'].indexOf(d.status);
  return [...docs].sort((a, b) => rank(a) - rank(b) || compareVersions(b.version, a.version) ||
    String(b.on ?? '').localeCompare(String(a.on ?? '')) || a.title.localeCompare(b.title));
}

/** `Process description · v0.3 · draft` */
export const describeDocument = (d: Pick<DocumentRef, 'title' | 'version' | 'status'>): string =>
  [d.title, d.version ? `v${d.version}` : '', d.status].filter((x) => x).join(' · ');

/** A document as evidence: kind `file` (the evidence kinds are kept), its path as the target, and
 * the record it is in the context, so the history and provenance can find it again. */
export function documentEvidence(d: DocumentRef, criterionCode?: string | null): EvidenceItem {
  return {kind: 'file', target: d.path, name: d.title, role: 'supporting',
    context: {document: d.id, documentKind: d.kind, version: d.version ?? '', status: d.status, path: d.path,
      ...(criterionCode ? {criterion: criterionCode} : {})}};
}

/** One document criterion, as its card shows it. */
export interface DocumentCard {
  criterionId: string;
  code: string;
  name: string;
  required: boolean;
  /** The stage whose gate it belongs to, and whether that is the current or the next gate. */
  stage: string;
  gate: 'current' | 'next';
  state: CriterionState;
  documents: DocumentRef[];
  /** What is missing, in one sentence; empty when an approved document answers it or it is waived. */
  missing: string;
  /** The run selected for the criterion, when it takes one: its title. */
  run: string | null;
  /** Decisions that cite one of its documents as evidence: one line each. */
  citedBy: string[];
}

export interface CardCriterion {
  id: string;
  stageId: string;
  code: string;
  name: string;
  kind: string;
  required: boolean;
}

/** The cards: every document criterion of the current gate (the program's stage) and the next,
 * in that order and in declaration order. */
export function documentCards(criteria: CardCriterion[], gates: {current: {id: string, name: string} | null,
  next: {id: string, name: string} | null}, docs: DocumentRef[], states: Record<string, CriterionState>,
  runs: Record<string, string> = {}, citations: Record<string, string[]> = {}): DocumentCard[] {
  const out: DocumentCard[] = [];
  for (const [gate, stage] of [['current', gates.current], ['next', gates.next]] as const) {
    if (!stage)
      continue;
    for (const c of criteria.filter((x) => x.stageId === stage.id && x.kind === 'document')) {
      const mine = orderDocuments(docs.filter((d) => d.criterionId === c.id));
      const state = states[c.id] ?? 'not_started';
      const approved = mine.some((d) => d.status === 'approved');
      const missing = state === 'waived' || approved ? '' :
        mine.some((d) => d.status === 'submitted') ? 'submitted: waiting for approval' :
          mine.some((d) => d.status === 'draft') ? 'only drafts: approving one meets the criterion' :
            mine.some((d) => d.status === 'rejected') ? 'rejected: link a new version, or submit it again' :
          mine.length ? 'only superseded documents' : 'no document yet';
      out.push({criterionId: c.id, code: c.code, name: c.name, required: c.required, stage: stage.name, gate, state,
        documents: mine, missing, run: runs[c.id] ?? null,
        citedBy: [...new Set(mine.flatMap((d) => citations[d.id] ?? []))]});
    }
  }
  return out;
}
