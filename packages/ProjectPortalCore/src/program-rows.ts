/** Program rows made visible to the portfolio's groups by linking each row into its workspace
 * folder. The rules and why are in ./model/program-rows. The workspace builder links a program when
 * it makes or reuses its folder; `linkProgramRows` does it for a portfolio's existing folders. */
import * as grok from 'datagrok-api/grok';
import {portalDb} from './generated/db';
import {portfolioByCode} from './portal';
import {describeRowLinks, planRowLinks, ProgramFolderRef} from './model/program-rows';
import {filterValue, ROW_LIMIT} from './rows';

const idOf = (x: any): string => typeof x === 'string' ? x : x?.id;

/** The entity ids linked into each of [spaceIds]; a space that cannot be read is left out. */
export async function folderLinks(spaceIds: string[]): Promise<Map<string, Set<string>>> {
  const out = new Map<string, Set<string>>();
  for (const id of new Set(spaceIds)) {
    const space = await grok.dapi.projects.find(id).catch(() => null);
    if (space)
      out.set(id, new Set((space.links ?? []).map(idOf)));
  }
  return out;
}

/** Links the program row into its folder, unless it is there. True when it linked. The server
 * refuses a second link to the same folder, hence the check. */
export async function linkProgramRow(programId: string, folderId: string): Promise<boolean> {
  const links = await folderLinks([folderId]);
  if (links.get(folderId)?.has(programId))
    return false;
  await grok.dapi.spaces.id(folderId).addEntity(programId, true);
  return true;
}

/** Links every program of portfolio [code] that has a workspace folder into it. Idempotent; with
 * [dryRun] only reports. Creates no folder: building the workspace does that, and links too. */
export async function linkProgramRows(code: string, dryRun = false): Promise<string[]> {
  const portfolio = await portfolioByCode(code);
  const rows = await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)}`, limit: ROW_LIMIT});
  const programs: ProgramFolderRef[] = rows.map((p) => ({id: p.id, code: p.code, spaceId: p.space_id ?? null}))
    .sort((a, b) => a.code.localeCompare(b.code));
  const plan = planRowLinks(programs, await folderLinks(programs.map((p) => p.spaceId).filter((x): x is string => !!x)));
  if (!dryRun) {
    for (const p of [...plan.link]) {
      try {
        await grok.dapi.spaces.id(p.spaceId!).addEntity(p.id, true);
      } catch (e: any) {
        plan.link.splice(plan.link.indexOf(p), 1);
        plan.skipped.push({program: p, reason: `not linked: ${e?.message ?? e}`});
      }
    }
  }
  return describeRowLinks(plan, dryRun);
}
