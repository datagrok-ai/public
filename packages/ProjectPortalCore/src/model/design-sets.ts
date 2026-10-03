/** Design sets and selection sets (F-23, FR-101/102), discovery. A design set is what a chemist
 * proposes for a series: structures enumerated, generated or drawn, each scored and ranked, and the
 * ones chosen for synthesis marked. A selection set is a hit list: compounds picked from a
 * campaign, with the rationale. Both are records of decisions about compounds, never folders.
 * Pure: the caller reads the program's series, campaigns and compounds. */

export type DesignSource = 'enumerated' | 'generative' | 'manual';
export const DESIGN_SOURCES: DesignSource[] = ['enumerated', 'generative', 'manual'];

export interface ProposalInput {
  label: string;
  structure: string;
  score?: number | null;
  note?: string | null;
}

export interface DesignSetInput {
  name: string;
  series: string;
  source: DesignSource;
  scoreName?: string | null;
  higherIsBetter?: boolean;
  description?: string | null;
  proposals: ProposalInput[];
}

/** Problems with a design set before it is written, as sentences. [series] are the program's own. */
export function validateDesignSet(d: Partial<DesignSetInput>, series: string[]): string[] {
  const errors: string[] = [];
  if (!d.name?.trim())
    errors.push('a design set needs a name');
  if (!d.series?.trim())
    errors.push('a design set extends a series');
  else if (!series.includes(d.series.trim()))
    errors.push(`the program has no series ${d.series.trim()}`);
  if (!d.source || !DESIGN_SOURCES.includes(d.source))
    errors.push(`source is one of ${DESIGN_SOURCES.join(', ')}`);
  const proposals = d.proposals ?? [];
  if (!proposals.length)
    errors.push('a design set proposes at least one structure');
  const seen = new Set<string>();
  proposals.forEach((p, i) => {
    if (!p.label?.trim())
      errors.push(`proposal ${i + 1} needs a label`);
    else if (seen.has(p.label.trim()))
      errors.push(`label ${p.label.trim()} is used twice`);
    else
      seen.add(p.label.trim());
    if (!p.structure?.trim())
      errors.push(`proposal ${p.label || i + 1} needs a structure`);
    if (p.score !== undefined && p.score !== null && !Number.isFinite(p.score))
      errors.push(`proposal ${p.label || i + 1}: the score is not a number`);
  });
  return errors;
}

/** Ranks by score, best first (descending unless [higherIsBetter] is false); ties share a rank
 * (1, 2, 2, 4); unscored proposals come last, unranked. Returns the rank per label. */
export function rankProposals(proposals: {label: string, score?: number | null}[], higherIsBetter = true): Map<string, number | null> {
  const scored = proposals.filter((p) => p.score !== null && p.score !== undefined && Number.isFinite(p.score));
  const sorted = [...scored].sort((a, b) => higherIsBetter ? b.score! - a.score! : a.score! - b.score!);
  const out = new Map<string, number | null>();
  sorted.forEach((p, i) => out.set(p.label, i > 0 && sorted[i - 1].score === p.score ? out.get(sorted[i - 1].label)! : i + 1));
  for (const p of proposals)
    if (!out.has(p.label))
      out.set(p.label, null);
  return out;
}

const STRUCTURE_COLUMNS = ['smiles', 'structure', 'molecule', 'mol', 'canonical_smiles'];
const LABEL_COLUMNS = ['label', 'id', 'name', 'code'];
const SCORE_COLUMNS = ['score'];

/** Which columns of a table hold the structure, the label and the score: a Molecule column or one
 * named smiles/structure/molecule; label/id/name/code; score, or [scoreColumn] when given. */
export function proposalColumns(columns: {name: string, semType?: string | null, numeric?: boolean}[], scoreColumn?: string | null):
  {structure: string | null, label: string | null, score: string | null} {
  const by = (names: string[]) => columns.find((c) => names.includes(c.name.toLowerCase()))?.name ?? null;
  return {
    structure: columns.find((c) => c.semType === 'Molecule')?.name ?? by(STRUCTURE_COLUMNS),
    label: by(LABEL_COLUMNS),
    score: scoreColumn ? columns.find((c) => c.name === scoreColumn)?.name ?? null : by(SCORE_COLUMNS),
  };
}

/** Proposals from table rows; a row without a label is named `<prefix>-<n>` (01, 02, …). */
export function proposalsFromRows(rows: Record<string, unknown>[], cols: {structure: string | null, label: string | null,
  score: string | null}, prefix: string): ProposalInput[] {
  const width = String(rows.length).length < 2 ? 2 : String(rows.length).length;
  return rows.map((r, i) => {
    const s = cols.score ? r[cols.score] : null;
    const label = cols.label && r[cols.label] !== null && r[cols.label] !== undefined && String(r[cols.label]).trim() ?
      String(r[cols.label]).trim() : `${prefix}-${String(i + 1).padStart(width, '0')}`;
    return {label, structure: cols.structure ? String(r[cols.structure] ?? '').trim() : '',
      score: s === null || s === undefined || s === '' || Number.isNaN(Number(s)) ? null : Number(s)};
  });
}

/** A design set's code from its name and date: `KIN-A-2026-10-03`, with a counter when taken. */
export function designSetCode(series: string, date: Date, taken: string[]): string {
  const base = `${series}-${date.toISOString().slice(0, 10)}`;
  let code = base;
  for (let n = 2; taken.includes(code); n++)
    code = `${base}-${n}`;
  return code;
}

export interface DesignSetSummary {
  proposed: number;
  chosen: number;
  /** The best-ranked label, if any is scored. */
  best: string | null;
}

export function summarizeDesignSet(proposals: {label: string, rank?: number | null, chosen?: boolean | null}[]): DesignSetSummary {
  const ranked = proposals.filter((p) => p.rank).sort((a, b) => a.rank! - b.rank!);
  return {proposed: proposals.length, chosen: proposals.filter((p) => p.chosen).length, best: ranked[0]?.label ?? null};
}

/** "12 proposed · 3 chosen for synthesis · best KIN-A-D07" */
export function describeDesignSet(s: DesignSetSummary): string {
  return [`${s.proposed} proposed`, `${s.chosen} chosen for synthesis`, ...(s.best ? [`best ${s.best}`] : [])].join(' · ');
}

/** Contributors and approvers design and choose, as they select runs. Not on an ended program. */
export function checkDesignAccess(a: {isApprover: boolean, isContributor: boolean, programStatus?: string | null}): string[] {
  const reasons: string[] = [];
  if (!a.isApprover && !a.isContributor)
    reasons.push('Only the program\'s contributors and approvers propose and choose designs');
  if (['terminated', 'completed', 'transitioned'].includes(a.programStatus ?? ''))
    reasons.push(`The program is ${a.programStatus}`);
  return reasons;
}

export interface SelectionSetInput {
  name: string;
  campaign?: string | null;
  rationale: string;
  compounds: string[];
}

/** Problems with a hit list, as sentences. [campaigns] and [compounds] are the program's own. */
export function validateSelectionSet(s: Partial<SelectionSetInput>, campaigns: string[], compounds: string[]): string[] {
  const errors: string[] = [];
  if (!s.name?.trim())
    errors.push('a selection set needs a name');
  if (s.campaign && !campaigns.includes(s.campaign))
    errors.push(`the program has no campaign ${s.campaign}`);
  if (!s.rationale?.trim())
    errors.push('a selection set needs a rationale: why these compounds');
  const list = s.compounds ?? [];
  if (!list.length)
    errors.push('a selection set holds at least one compound');
  const unknown = list.filter((c) => !compounds.includes(c));
  if (unknown.length)
    errors.push(`not compounds of the program: ${unknown.join(', ')}`);
  if (new Set(list).size !== list.length)
    errors.push('a compound is listed twice');
  return errors;
}

/** The compounds as stored (`{items: [...]}`: a json column holds a map). */
export function selectionCompounds(raw: unknown): string[] {
  const items = Array.isArray(raw) ? raw : Array.isArray((raw as any)?.items) ? (raw as any).items : [];
  return (items as unknown[]).filter((x) => typeof x === 'string' && x.trim()).map((x) => (x as string).trim());
}
