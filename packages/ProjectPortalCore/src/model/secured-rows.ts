/** Who reads the records that belong to one portfolio or one user. Table security is one grant for
 * the whole table, so approval_request, sync_run, follow and compound, readable by All users, let
 * anyone read every portfolio's requests and rationales, sync history, compound structures and
 * everyone's follows. The server cannot change a deployed table's security mode, so the records
 * moved:
 *
 * - approval requests to `approval` and sync runs to `sync_log`: row-secured, each row linked into
 *   a space. A request goes into its program's workspace folder (the portfolio's groups and the
 *   program's roles view it), a model-card request into `_Approvals` and a sync run into
 *   `_Sync history` under the portfolio's root (before 0.3.0: the root itself). The
 *   server lets whoever may view a space view the rows linked into it, View only (the e5adf17c01
 *   pattern), and nobody else; the row's author reaches it too.
 * - follows to `following`: row-secured with no row grants, so each user reaches only their own.
 * - a compound's name, structure, aliases and batches to its program_compound links, which are
 *   secured through the program. `compound` keeps the identity (registry and code).
 *
 * Pure: the caller reads the rows, the folders and their links. */
import {APPROVALS_FOLDER, SYNC_HISTORY_FOLDER} from './workspace-plan';

/** Columns the server writes on every row; a copy leaves them out, except the id. */
export const SYSTEM_COLUMNS = ['version', 'created_on', 'updated_on', 'author_id', 'is_deleted', 'deleted'];

/** [row] as values for an insert into the table that replaces its own: same id, no system columns,
 * no empty values. */
export function copyValues(row: Record<string, unknown>): Record<string, unknown> {
  const out: Record<string, unknown> = {};
  for (const [k, v] of Object.entries(row)) {
    if (!SYSTEM_COLUMNS.includes(k) && v !== null && v !== undefined)
      out[k] = v;
  }
  return out;
}

/** The ids of [old] not yet copied into the new table, which holds [copied]. */
export function toCopy<T extends {id: string}>(old: T[], copied: Set<string>): T[] {
  return old.filter((r) => !copied.has(r.id));
}

/** A row to place: its id and the space it belongs in (null: nowhere to put it). */
export interface Placement {
  id: string;
  label: string;
  spaceId: string | null;
}

/** The space a request belongs in: its program's folder when the program has one, else the
 * portfolio's `_Approvals` folder ([recordsId]). A model card has no program and goes there. */
export function requestSpace(programSpaceId: string | null | undefined, recordsId: string | null): string | null {
  return programSpaceId || recordsId || null;
}

export interface LinkPlan {
  link: Placement[];
  linked: number;
  /** Rows with nowhere to go: they stay with their author (and admins). */
  unplaced: Placement[];
}

/** What to link: [links] holds each existing space's linked ids; a space missing from it is gone
 * and counts as nowhere. */
export function planLinks(rows: Placement[], links: Map<string, Set<string>>): LinkPlan {
  const plan: LinkPlan = {link: [], linked: 0, unplaced: []};
  for (const r of rows) {
    if (!r.spaceId || !links.has(r.spaceId))
      plan.unplaced.push(r);
    else if (links.get(r.spaceId)!.has(r.id))
      plan.linked++;
    else
      plan.link.push(r);
  }
  return plan;
}

/** [rows] of the portfolios in [mine]; every row when [mine] is null (an admin). The server already
 * keeps other portfolios' rows away; this is the app's second check. */
export function inPortfolios<T extends {portfolio_id: string}>(rows: T[], mine: Set<string> | null): T[] {
  return mine === null ? rows : rows.filter((r) => mine.has(r.portfolio_id));
}

/** The portfolios a user belongs to: those whose viewers, contributors or approvers group is one of
 * [groupIds]. */
export function portfoliosOfGroups(portfolios: {id: string, viewers?: string | null, contributors?: string | null,
  approvers?: string | null}[], groupIds: Set<string>): Set<string> {
  return new Set(portfolios.filter((p) => [p.viewers, p.contributors, p.approvers].some((g) => g && groupIds.has(g)))
    .map((p) => p.id));
}

// ---------------------------------------------------------------------------------------------
// Compounds under their programs

/** The compound's values a program_compound link carries: [compound column, link column]. */
export const LINK_VALUES: [string, string][] = [
  ['code', 'compound_code'], ['name', 'compound_name'], ['smiles', 'smiles'], ['aliases', 'aliases'],
  ['batches', 'batches'], ['details', 'details'],
];

/** The values that leaked through `compound`: emptied there once each link carries them. */
export const COMPOUND_CONTENT = ['name', 'smiles', 'aliases', 'batches', 'details'];

/** A compound's values as link columns. Missing values are left out. */
export function linkValuesOf(compound: Record<string, any>): Record<string, unknown> {
  const out: Record<string, unknown> = {};
  for (const [c, l] of LINK_VALUES) {
    if (compound[c] !== null && compound[c] !== undefined)
      out[l] = compound[c];
  }
  return out;
}

/** The compound a link names, from the link's own values; [compound] (the identity row, when the
 * reader may see it) fills what an older link lacks. */
export function compoundOfLink(link: Record<string, any>, compound?: Record<string, any> | null): {id: string, code: string,
  name: string | null, smiles: string | null, aliases: string | null, batches: string | null, details: any} {
  const pick = (c: string, l: string) => link[l] ?? compound?.[c] ?? null;
  return {id: link.compound_id, code: String(pick('code', 'compound_code') ?? link.compound_id), name: pick('name', 'compound_name'),
    smiles: pick('smiles', 'smiles'), aliases: pick('aliases', 'aliases'), batches: pick('batches', 'batches'),
    details: pick('details', 'details')};
}

/** Whether a compound row still holds content that every user could once read. */
export const hasContent = (compound: Record<string, any>): boolean =>
  COMPOUND_CONTENT.some((k) => compound[k] !== null && compound[k] !== undefined && compound[k] !== '');

/** One report of a move. */
export interface MoveReport {
  portfolio: string;
  copied: Record<string, number>;
  already: Record<string, number>;
  linked: Record<string, number>;
  unplaced: string[];
  linkValues: number;
  emptied: number;
  errors: string[];
}

export function describeMove(r: MoveReport, dryRun = false): string[] {
  const would = dryRun ? 'would ' : '';
  const lines: string[] = [];
  for (const t of Object.keys(r.copied))
    lines.push(`${r.portfolio} ${t}: ${would}copy ${r.copied[t]}, ${r.already[t] ?? 0} copied before; ${would}link ${r.linked[t] ?? 0}`);
  if (r.linkValues || r.emptied)
    lines.push(`${r.portfolio} compounds: ${would}copy values onto ${r.linkValues} program link(s), ${would}empty ${r.emptied} compound row(s)`);
  for (const u of r.unplaced)
    lines.push(`${r.portfolio} ${u}: nowhere to link it (no workspace): its author and admins read it`);
  lines.push(...r.errors.map((e) => `${r.portfolio}: ${e}`));
  return lines;
}

// ---------------------------------------------------------------------------------------------
// Records out of the portfolio's root

/** The folder under the root a record linked into the root belongs in, by its entity name
 * (`portal.sync_log:<id>`): sync runs in `_Sync history`, requests in `_Approvals`; null for
 * anything else (a program folder, the root's files). */
export function recordFolderFor(entityName: string): string | null {
  const type = String(entityName ?? '').split(':')[0].toLowerCase();
  return type === 'portal.sync_log' ? SYNC_HISTORY_FOLDER : type === 'portal.approval' ? APPROVALS_FOLDER : null;
}

/** One record to take out of the root: linked into [folder] first, then unlinked from the root. */
export interface RecordMove {
  id: string;
  folder: string;
  /** Already linked into [folder]: only the root's link goes. */
  linked: boolean;
}

/** What taking the records out of the root does: [rootChildren] are the root's children with
 * links; [inFolder] holds, by folder name, the ids already linked into that folder. */
export function planRecordMoves(rootChildren: {id: string, name: string}[], inFolder: Map<string, Set<string>>): RecordMove[] {
  const moves: RecordMove[] = [];
  for (const c of rootChildren) {
    const folder = recordFolderFor(c.name);
    if (folder)
      moves.push({id: c.id, folder, linked: inFolder.get(folder)?.has(c.id) ?? false});
  }
  return moves;
}

/** `EX-SMD: would move 249 record(s) out of the root: 249 into _Sync history`; or that none is left. */
export function describeRecordMoves(code: string, moves: RecordMove[], dryRun: boolean, errors: string[] = []): string[] {
  if (!moves.length && !errors.length)
    return [`${code}: no records in the root`];
  const by = new Map<string, number>();
  for (const m of moves)
    by.set(m.folder, (by.get(m.folder) ?? 0) + 1);
  const into = [...by].map(([f, n]) => `${n} into ${f}`).join(', ');
  return [...(moves.length ? [`${code}: ${dryRun ? 'would move' : 'moved'} ${moves.length - (dryRun ? 0 : errors.length)} record(s) ` +
    `out of the root: ${into}`] : []), ...errors.map((e) => `${code}: ${e}`)];
}
