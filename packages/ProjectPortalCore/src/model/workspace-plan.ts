/** What a portfolio's workspace tree looks like: root → one folder per program → its workstreams.
 * Pure data and pure checks, so the naming rules hold wherever they are read. Compounds, plates,
 * results and predictions are never folders (FR-26); the deepest level is a workstream. */
import {normalizeSpaceName} from './config';

export interface SpaceRef {
  id: string;
  name: string;
  friendlyName?: string | null;
}

export interface WorkspaceProgram {
  id: string;
  code: string;
  name: string;
  spaceId?: string | null;
  stage?: string | null;
  status?: string | null;
}

/** A program's folder: its code first, so folders sort and read the way the registry lists them. */
export function programFolderName(p: Pick<WorkspaceProgram, 'code' | 'name'>): string {
  const name = (p.name ?? '').trim();
  return name === '' || name === p.code ? p.code : `${p.code} · ${name}`;
}

/** Whether [space] is the folder called [name]: the server may hold either spelling. */
export function isSpaceNamed(space: SpaceRef, name: string): boolean {
  const key = normalizeSpaceName(name).toLowerCase();
  return normalizeSpaceName(space.friendlyName ?? '').toLowerCase() === key ||
    normalizeSpaceName(space.name).toLowerCase() === key;
}

/** How the builder finds a program's folder among the root's children. The recorded space id
 * wins, so a program renamed in its registry keeps its folder; the name is the fallback for a
 * folder made before the id was recorded, or by hand. */
export function findProgramFolder(p: WorkspaceProgram, children: SpaceRef[]):
  {space: SpaceRef, by: 'id' | 'name'} | null {
  if (p.spaceId) {
    const byId = children.find((c) => c.id === p.spaceId);
    if (byId)
      return {space: byId, by: 'id'};
  }
  const byName = children.find((c) => isSpaceNamed(c, programFolderName(p)));
  return byName ? {space: byName, by: 'name'} : null;
}

/** Display attributes written on a program's folder. Space properties are not queryable
 * (PLAT-04), so they repeat the `portal` row rather than replace it. */
export function programAttributes(portfolioCode: string, p: WorkspaceProgram): Record<string, string> {
  return {
    portal_portfolio: portfolioCode,
    portal_program: p.code,
    portal_stage: p.stage ?? '',
    portal_status: p.status ?? '',
  };
}

export function parseWorkstreams(stored: string | null | undefined): string[] {
  return (stored ?? '').split('|').map((w) => w.trim()).filter((w) => w.length > 0);
}

/** A stage folder path as declared: segments under the program's folder, '/'-separated. */
export function splitFolderPath(path: string): string[] {
  return String(path ?? '').split('/').map((s) => s.trim()).filter((s) => s.length > 0);
}

export const MAX_FOLDER_DEPTH = 3;

/** Problems with one stage's folder paths, as sentences. */
export function validateStageFolders(stageCode: string, folders: string[]): string[] {
  const errors: string[] = [];
  const seen = new Set<string>();
  for (const f of folders) {
    const parts = splitFolderPath(f);
    if (parts.length === 0)
      errors.push(`stage ${stageCode}: folder "${f}" is empty`);
    else if (parts.length > MAX_FOLDER_DEPTH)
      errors.push(`stage ${stageCode}: folder "${f}" is deeper than ${MAX_FOLDER_DEPTH} levels`);
    else if (f.includes('|'))
      errors.push(`stage ${stageCode}: folder "${f}" may not contain '|'`);
    else if (parts.some((p) => normalizeSpaceName(p) === ''))
      errors.push(`stage ${stageCode}: folder "${f}" has a segment with no letters or digits`);
    const key = parts.map((p) => normalizeSpaceName(p).toLowerCase()).join('/');
    if (seen.has(key))
      errors.push(`stage ${stageCode}: folder "${f}" is declared twice`);
    seen.add(key);
  }
  return errors;
}

export interface StageFolders {
  code: string;
  ordinal: number;
  folders: string[];
}

/** The folders a program at [ordinal] has earned: those of every stage up to and including its
 * own (PORT-35). An unstaged program has none. A stage reached later only adds; a program that
 * recycles keeps what it had, since the builder never removes a folder. */
export function stageFolderPaths(stages: StageFolders[], ordinal: number | null): string[][] {
  if (ordinal === null)
    return [];
  const out: string[][] = [];
  const seen = new Set<string>();
  for (const s of [...stages].sort((a, b) => a.ordinal - b.ordinal)) {
    if (s.ordinal > ordinal)
      break;
    for (const f of s.folders) {
      const parts = splitFolderPath(f);
      const key = parts.map((p) => normalizeSpaceName(p).toLowerCase()).join('/');
      if (parts.length && !seen.has(key)) {
        seen.add(key);
        out.push(parts);
      }
    }
  }
  return out;
}

export interface Orphan {
  /** The folder's space id: what parking moves. */
  id: string;
  path: string;
  kind: 'unclaimed' | 'park';
  reason: string;
}

const PARK_STATUSES = ['terminated', 'transitioned', 'missing_from_feed', 'completed'];

/** The folder under a portfolio's root that parked folders go into (FR-22), one subfolder per day. */
export const PARKED_FOLDER = '_Parked';

/** Whether [space] is the parking folder. The server may hold it as `Parked`. */
export const isParkingFolder = (space: SpaceRef): boolean => isSpaceNamed(space, PARKED_FOLDER);

/** What a prune-or-park pass would look at (FR-22): root folders no program of the portfolio
 * claims, and folders of programs that have ended or left their feed. Nothing below a program
 * folder is judged: people put their own folders there. The parking folder is never listed. */
export function findOrphans(root: string, children: SpaceRef[], programs: WorkspaceProgram[]): Orphan[] {
  const claimed = new Map<string, WorkspaceProgram>();
  for (const p of programs) {
    const found = findProgramFolder(p, children);
    if (found)
      claimed.set(found.space.id, p);
  }
  const out: Orphan[] = [];
  for (const c of [...children].sort((a, b) => (a.friendlyName ?? a.name).localeCompare(b.friendlyName ?? b.name))) {
    // parked folders are already dealt with
    if (isParkingFolder(c))
      continue;
    const label = `${root}/${c.friendlyName ?? c.name}`;
    const p = claimed.get(c.id);
    if (!p)
      out.push({id: c.id, path: label, kind: 'unclaimed', reason: 'no program of the portfolio claims it'});
    else if (PARK_STATUSES.includes(p.status ?? ''))
      out.push({id: c.id, path: label, kind: 'park', reason: `${p.code} is ${p.status}`});
  }
  return out;
}
