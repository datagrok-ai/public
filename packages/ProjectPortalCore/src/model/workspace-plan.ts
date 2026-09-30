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
