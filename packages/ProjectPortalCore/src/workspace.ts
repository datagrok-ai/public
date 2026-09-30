/** The workspace builder: a Spaces tree per portfolio, root → program → workstream.
 *
 * Additive: it creates what is missing and never renames, moves or deletes. A program renamed in
 * its registry keeps its folder (the folder's id is recorded on the program row) and the report
 * names the difference. Learned building an earlier workspace generator:
 * - The server stores a space name in normalized form ("Molecule A" → "MoleculeA"), so
 *   `subspaceExists(display name)` misses a folder that is there. Each parent's children are
 *   listed once and matched in both forms.
 * - Creation is slow (~200 ms a folder). Two builds against one root race and duplicate folders,
 *   so a build is refused while another for the same root runs in this session. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {normalizeSpaceName} from './model/config';
import {findProgramFolder, isSpaceNamed, parseWorkstreams, programAttributes, programFolderName, SpaceRef,
  WorkspaceProgram} from './model/workspace-plan';
import {portalDb, PortfolioRow} from './generated/db';
import {portfolioByCode} from './portal';

export interface WorkspaceEntry {
  path: string;
  action: 'create' | 'reuse';
  status: string;
}

export interface WorkspaceReport {
  portfolio: string;
  root: string;
  dryRun: boolean;
  entries: WorkspaceEntry[];
  created: number;
  reused: number;
  failed: number;
  /** Things the builder saw but, being additive, did not change. */
  notes: string[];
}

const building = new Set<string>();

const asRef = (e: DG.Entity): SpaceRef => ({id: e.id, name: e.name, friendlyName: (e as any).friendlyName ?? null});

async function subspaces(parentId: string): Promise<DG.Project[]> {
  return await grok.dapi.spaces.id(parentId).children.ofTypes('Project', false).list({pageSize: 1000}) as DG.Project[];
}

async function findRoot(name: string): Promise<DG.Project | null> {
  const found = await grok.dapi.spaces.filter(`name = "${name}"`).list({pageSize: 100});
  return found.find((s) => isSpaceNamed(asRef(s), name)) ?? null;
}

/** View to viewers, Edit to contributors; subspaces inherit. Unlike program rows, spaces take
 * `permissions.grant`, so the workspace is isolated per portfolio today. */
async function shareRoot(root: DG.Project, portfolio: PortfolioRow, notes: string[]): Promise<void> {
  // a space nobody has been granted comes back without the view/edit keys
  const current = await grok.dapi.permissions.get(root) as {view?: DG.Group[], edit?: DG.Group[]};
  const targets: [string | undefined, boolean][] = [[portfolio.viewers, false], [portfolio.contributors, true]];
  for (const [groupId, edit] of targets) {
    if (!groupId)
      continue;
    const held = ((edit ? current.edit : current.view) ?? []).some((g) => g.id === groupId);
    if (held)
      continue;
    const group = await grok.dapi.groups.find(groupId);
    await grok.dapi.permissions.grant(root, group, edit);
    notes.push(`shared with ${group.friendlyName ?? group.name} (${edit ? 'edit' : 'view'})`);
  }
}

async function writeAttributes(space: DG.Project, attributes: Record<string, string>): Promise<boolean> {
  // they come back namespaced: portal_stage is read as default:portal_stage
  const current = await space.getProperties().catch(() => ({} as Record<string, any>));
  const read = (k: string) => String(current[k] ?? current[`default:${k}`] ?? '');
  if (Object.entries(attributes).every(([k, v]) => read(k) === v))
    return false;
  await space.setProperties(attributes);
  return true;
}

/** Builds (or with [dryRun] previews) the workspace of one portfolio. */
export async function buildPortfolioWorkspace(code: string, dryRun = false): Promise<WorkspaceReport> {
  const portfolio = await portfolioByCode(code);
  const rootName = portfolio.workspace_root;
  if (!rootName)
    throw new Error(`portfolio ${code} has no workspaceRoot`);
  const lock = normalizeSpaceName(rootName).toLowerCase();
  if (building.has(lock))
    throw new Error(`a workspace build under "${rootName}" is already running; wait for it to finish`);
  building.add(lock);
  try {
    return await build(portfolio, rootName, dryRun);
  } finally {
    building.delete(lock);
  }
}

async function build(portfolio: PortfolioRow, rootName: string, dryRun: boolean): Promise<WorkspaceReport> {
  const report: WorkspaceReport = {portfolio: portfolio.code, root: rootName, dryRun, entries: [],
    created: 0, reused: 0, failed: 0, notes: []};
  const record = (path: string, action: 'create' | 'reuse', status = 'ok') => {
    report.entries.push({path, action, status});
    if (status !== 'ok')
      report.failed++;
    else if (action === 'create')
      report.created++;
    else
      report.reused++;
  };
  const workstreams = parseWorkstreams(portfolio.workstreams);

  let root = await findRoot(rootName);
  record(rootName, root ? 'reuse' : 'create');
  if (!root && !dryRun)
    root = await grok.dapi.spaces.createRootSpace(rootName);
  if (root && !dryRun)
    await shareRoot(root, portfolio, report.notes);

  const children = root ? await subspaces(root.id) : [];
  const stages = new Map((await portalDb.stages.query({filter: `stage_model_id = "${portfolio.stage_model_id}"`}))
    .map((s) => [s.id, s.code]));
  const programs: WorkspaceProgram[] = (await portalDb.programs.query({filter: `portfolio_id = "${portfolio.id}"`}))
    .map((p) => ({id: p.id, code: p.code, name: p.name, spaceId: p.space_id ?? null,
      stage: p.stage_id ? stages.get(p.stage_id) ?? null : null, status: p.status ?? null}))
    .sort((a, b) => a.code.localeCompare(b.code));

  let attributesWritten = 0;
  for (const p of programs) {
    const folderName = programFolderName(p);
    const path = `${rootName}/${folderName}`;
    try {
      const found = findProgramFolder(p, children.map(asRef));
      let folder = found ? children.find((c) => c.id === found.space.id)! : null;
      record(path, folder ? 'reuse' : 'create');
      if (found && found.by === 'id' && !isSpaceNamed(found.space, folderName)) {
        report.notes.push(`${p.code}: folder is "${found.space.friendlyName ?? found.space.name}", ` +
          `the registry now says "${folderName}"; not renamed`);
      }
      if (p.spaceId && found?.by !== 'id')
        report.notes.push(`${p.code}: its recorded folder is gone from ${rootName}`);
      if (dryRun) {
        const existing = folder ? await subspaces(folder.id) : [];
        for (const w of workstreams)
          record(`${path}/${w}`, existing.some((c) => isSpaceNamed(asRef(c), w)) ? 'reuse' : 'create');
        continue;
      }
      if (!folder) {
        folder = await grok.dapi.spaces.id(root!.id).addSubspace(folderName);
        children.push(folder);
      }
      if (folder.id !== p.spaceId)
        await portalDb.programs.update(p.id, {space_id: folder.id});
      if (await writeAttributes(folder, programAttributes(portfolio.code, p)))
        attributesWritten++;
      const existing = workstreams.length > 0 ? await subspaces(folder.id) : [];
      for (const w of workstreams) {
        const there = existing.some((c) => isSpaceNamed(asRef(c), w));
        if (!there)
          existing.push(await grok.dapi.spaces.id(folder.id).addSubspace(w));
        record(`${path}/${w}`, there ? 'reuse' : 'create');
      }
    } catch (e: any) {
      record(path, 'create', `error: ${e?.message ?? e}`);
    }
  }
  if (attributesWritten > 0)
    report.notes.push(`attributes written on ${attributesWritten} program folder(s)`);
  return report;
}

export function summarizeWorkspace(r: WorkspaceReport): string {
  const head = `${r.portfolio} workspace "${r.root}"${r.dryRun ? ' (preview)' : ''}: ` +
    `${r.dryRun ? 'would create' : 'created'} ${r.created}, reused ${r.reused}, failed ${r.failed}`;
  const lines = [head, ...r.notes.map((n) => `  ${n}`)];
  if (r.dryRun)
    lines.push(...r.entries.filter((e) => e.action === 'create').map((e) => `  + ${e.path}`));
  lines.push(...r.entries.filter((e) => e.status !== 'ok').map((e) => `  ! ${e.path}: ${e.status}`));
  return lines.join('\n');
}
