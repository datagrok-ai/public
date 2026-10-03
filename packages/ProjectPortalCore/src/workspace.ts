/** The workspace builder: a Spaces tree per portfolio, root → program → workstream, and with
 * `flowFolders` the unit-operation folders of each program's process flows.
 *
 * Additive: it creates what is missing and never renames, moves or deletes. A program renamed in
 * its registry keeps its folder (the folder's id is recorded on the program row) and the report
 * names the difference. Learned building an earlier workspace generator:
 * - The server stores a space name in normalized form ("Program A" → "ProgramA"), so
 *   `subspaceExists(display name)` misses a folder that is there. Each parent's children are
 *   listed once and matched in both forms.
 * - Creation is slow (~200 ms a folder). Two builds against one root race and duplicate folders,
 *   so a build is refused while another for the same root runs in this session. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {normalizeSpaceName} from './model/config';
import {findOrphans, findProgramFolder, isSpaceNamed, noNewFoldersReason, Orphan, PARKED_FOLDER, parseWorkstreams, programAttributes, programFolderName,
  SpaceRef, stageFolderPaths, StageFolders, WorkspaceProgram} from './model/workspace-plan';
import {portalDb, PortfolioRow} from './generated/db';
import {portfolioByCode} from './portal';
import {syncFolderGrants} from './program-roles';
import {filterValue, ROW_LIMIT} from './rows';
import {FlowFoldersConfig, flowFolderPaths, flowFoldersDue, flowsForProgram, ProgramFlow, TableFlow} from './model/flow-folders';

export interface WorkspaceEntry {
  path: string;
  /** skip: the folder is missing and is not made, for [reason]. */
  action: 'create' | 'reuse' | 'skip';
  status: string;
  reason?: string;
}

export interface WorkspaceReport {
  portfolio: string;
  root: string;
  dryRun: boolean;
  entries: WorkspaceEntry[];
  created: number;
  reused: number;
  failed: number;
  /** Folders missing but not made: those of programs that have ended. */
  skipped: number;
  /** Things the builder saw but, being additive, did not change. */
  notes: string[];
  /** Folders a prune-or-park pass would look at; the builder never touches them (parkOrphans
   * moves them, on request). Full builds only. */
  orphans: Orphan[];
}

const building = new Set<string>();

/** Runs [fn] while no build, park or unpark under the same root runs in this session. */
export async function withRootLock<T>(rootName: string, fn: () => Promise<T>): Promise<T> {
  const lock = normalizeSpaceName(rootName).toLowerCase();
  if (building.has(lock))
    throw new Error(`a workspace build or park under "${rootName}" is already running; wait for it to finish`);
  building.add(lock);
  try {
    return await fn();
  } finally {
    building.delete(lock);
  }
}

const asRef = (e: DG.Entity): SpaceRef => ({id: e.id, name: e.name, friendlyName: (e as any).friendlyName ?? null});

/** The ids among `ids` whose workspace folder still exists for the caller. A recorded space_id outlives a folder
 *  deleted by hand, so whether a program has a workspace is a question for the server, not for the row. */
export async function existingSpaces(ids: (string | null | undefined)[]): Promise<Set<string>> {
  const unique = [...new Set(ids.filter((x): x is string => !!x))];
  const found = new Set<string>();
  for (let i = 0; i < unique.length; i += 200)
    for (const e of await grok.dapi.getEntities(unique.slice(i, i + 200))) found.add(e.id);
  return found;
}

export async function subspaces(parentId: string): Promise<DG.Project[]> {
  return await grok.dapi.spaces.id(parentId).children.ofTypes('Project', false).list({pageSize: 1000}) as DG.Project[];
}

export async function findRoot(name: string): Promise<DG.Project | null> {
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

/** Builds (or with [dryRun] previews) the workspace of one portfolio; with [onlyProgramId], of one
 * of its programs, which is what a gate decision asks for after the program advances. */
export async function buildPortfolioWorkspace(code: string, dryRun = false, onlyProgramId?: string): Promise<WorkspaceReport> {
  const portfolio = await portfolioByCode(code);
  const rootName = portfolio.workspace_root;
  if (!rootName)
    throw new Error(`portfolio ${code} has no workspaceRoot`);
  return withRootLock(rootName, () => build(portfolio, rootName, dryRun, onlyProgramId));
}

async function build(portfolio: PortfolioRow, rootName: string, dryRun: boolean, onlyProgramId?: string):
  Promise<WorkspaceReport> {
  const report: WorkspaceReport = {portfolio: portfolio.code, root: rootName, dryRun, entries: [],
    created: 0, reused: 0, failed: 0, skipped: 0, notes: [], orphans: []};
  const record: RecordEntry = (path, action, status = 'ok', reason) => {
    report.entries.push({path, action, status, ...(reason ? {reason} : {})});
    if (status !== 'ok')
      report.failed++;
    else if (action === 'skip')
      report.skipped++;
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
  const stageRows = await portalDb.stages.query({filter: `stage_model_id = ${filterValue(portfolio.stage_model_id)}`, limit: ROW_LIMIT});
  const stages = new Map(stageRows.map((s) => [s.id, s.code]));
  const ordinals = new Map(stageRows.map((s) => [s.id, s.ordinal]));
  const stageFolders: StageFolders[] = stageRows.map((s) => ({code: s.code, ordinal: s.ordinal,
    folders: (s.folders ?? '').split('|').filter((f) => f.length > 0)}));
  const rows = await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)}`, limit: ROW_LIMIT});
  const all: (WorkspaceProgram & {ordinal: number | null})[] = rows
    .map((p) => ({id: p.id, code: p.code, name: p.name, spaceId: p.space_id ?? null,
      stage: p.stage_id ? stages.get(p.stage_id) ?? null : null, status: p.status ?? null,
      ordinal: p.stage_id ? ordinals.get(p.stage_id) ?? null : null}))
    .sort((a, b) => a.code.localeCompare(b.code));
  const programs = onlyProgramId ? all.filter((p) => p.id === onlyProgramId) : all;
  // programs with roles of their own get their folder shared to match
  const withRoles = dryRun ? new Set<string>() : new Set((await portalDb.programRoles.query({limit: ROW_LIMIT}).catch(() => []))
    .map((r) => r.program_id));
  // a parked program folder stays parked: the builder neither rebuilds it nor makes a new one
  const parked = new Set((await portalDb.parkedFolders.query({filter: `portfolio_id = ${filterValue(portfolio.id)}`, limit: ROW_LIMIT}).catch(() => []))
    .filter((r) => r.status === 'parked').map((r) => r.space_id));

  // unit-operation folders, when the portfolio asks for them
  const flowConfig = (portfolio.flow_folders ?? null) as FlowFoldersConfig | null;
  const flowsOf = new Map<string, ProgramFlow[]>();
  if (flowConfig) {
    const table = await tableFlows(programs.map((p) => p.id));
    const from = flowConfig.fromStage ? stageRows.find((s) => s.code === flowConfig.fromStage)?.ordinal ?? null : null;
    for (const p of programs) {
      const row = rows.find((r) => r.id === p.id)!;
      if (flowFoldersDue(from, p.ordinal))
        flowsOf.set(p.id, flowsForProgram(row.modality, table.get(p.id) ?? [], flowConfig));
    }
    const withFlows = programs.filter((p) => flowsOf.get(p.id)?.length);
    report.notes.push(`flow folders: ${withFlows.length ? withFlows.map((p) => `${p.code} (${flowsOf.get(p.id)!
      .map((f) => f.code).join(', ')})`).join(', ') : 'no program has a flow'}` +
      (programs.length > withFlows.length ? `; ${programs.length - withFlows.length} program(s) without` +
        `${flowConfig.fromStage ? ` (or before ${flowConfig.fromStage})` : ''}` : ''));
  }

  let attributesWritten = 0;
  for (const p of programs) {
    const folderName = programFolderName(p);
    const path = `${rootName}/${folderName}`;
    if (p.spaceId && parked.has(p.spaceId)) {
      report.notes.push(`${p.code}: folder parked under ${PARKED_FOLDER}; not rebuilt (unpark it to use it)`);
      continue;
    }
    // an ended program keeps the folders it has and gets no new ones
    const ended = noNewFoldersReason(p);
    try {
      const found = findProgramFolder(p, children.map(asRef));
      let folder = found ? children.find((c) => c.id === found.space.id)! : null;
      if (!folder && ended) {
        record(path, 'skip', 'ok', ended);
        continue;
      }
      record(path, folder ? 'reuse' : 'create');
      if (found && found.by === 'id' && !isSpaceNamed(found.space, folderName)) {
        report.notes.push(`${p.code}: folder is "${found.space.friendlyName ?? found.space.name}", ` +
          `the registry now says "${folderName}"; not renamed`);
      }
      if (p.spaceId && found?.by !== 'id')
        report.notes.push(`${p.code}: its recorded folder is gone from ${rootName}`);
      const extra = [...stageFolderPaths(stageFolders, p.ordinal),
        ...flowFolderPaths(flowsOf.get(p.id) ?? [], flowConfig?.placement)];
      if (dryRun) {
        const existing = folder ? await subspaces(folder.id) : [];
        for (const w of workstreams) {
          const there = existing.some((c) => isSpaceNamed(asRef(c), w));
          record(`${path}/${w}`, there ? 'reuse' : ended ? 'skip' : 'create', 'ok', there ? undefined : ended ?? undefined);
        }
        const counted = new Set(workstreams.map((w) => normalizeSpaceName(`${path}/${w}`).toLowerCase()));
        for (const parts of extra)
          await ensurePath(folder, parts, path, true, counted, record, ended);
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
      if (withRoles.has(p.id)) {
        const row = rows.find((r) => r.id === p.id)!;
        for (const n of await syncFolderGrants(row, portfolio, folder).catch((e) => [`roles not applied: ${e?.message ?? e}`]))
          report.notes.push(`${p.code}: ${n}`);
      }
      const existing = workstreams.length > 0 ? await subspaces(folder.id) : [];
      for (const w of workstreams) {
        const there = existing.some((c) => isSpaceNamed(asRef(c), w));
        if (!there && ended) {
          record(`${path}/${w}`, 'skip', 'ok', ended);
          continue;
        }
        if (!there)
          existing.push(await grok.dapi.spaces.id(folder.id).addSubspace(w));
        record(`${path}/${w}`, there ? 'reuse' : 'create');
      }
      const counted = new Set(workstreams.map((w) => normalizeSpaceName(`${path}/${w}`).toLowerCase()));
      for (const parts of extra)
        await ensurePath(folder, parts, path, false, counted, record, ended);
    } catch (e: any) {
      record(path, 'create', `error: ${e?.message ?? e}`);
    }
  }
  if (attributesWritten > 0)
    report.notes.push(`attributes written on ${attributesWritten} program folder(s)`);
  if (!onlyProgramId)
    report.orphans = findOrphans(rootName, children.map(asRef), all);
  return report;
}

/** The program_flow rows of [programIds], each joined to its flow and steps. A loader of process
 * flows writes them; the builder only reads. */
async function tableFlows(programIds: string[]): Promise<Map<string, TableFlow[]>> {
  const out = new Map<string, TableFlow[]>();
  if (!programIds.length)
    return out;
  const links = (await Promise.all(chunks(programIds, 200).map((ids) => portalDb.programFlows.query({
    filter: `program_id in (${ids.map(filterValue).join(', ')})`, limit: ROW_LIMIT})))).flat();
  const flowIds = [...new Set(links.map((l) => l.process_flow_id))];
  if (!flowIds.length)
    return out;
  const inFlows = flowIds.map(filterValue).join(', ');
  const [flows, steps] = await Promise.all([
    portalDb.processFlows.query({filter: `id in (${inFlows})`, limit: ROW_LIMIT}),
    portalDb.processSteps.query({filter: `process_flow_id in (${inFlows})`, limit: ROW_LIMIT}),
  ]);
  for (const l of links) {
    const f = flows.find((x) => x.id === l.process_flow_id);
    if (!f)
      continue;
    out.set(l.program_id, [...out.get(l.program_id) ?? [], {code: f.code, processType: f.process_type, rowSource: l.source,
      steps: steps.filter((s) => s.process_flow_id === f.id).map((s) => ({step: s.step_number, folderName: s.folder_name}))}]);
  }
  return out;
}

const chunks = <T>(xs: T[], n: number): T[][] => Array.from({length: Math.ceil(xs.length / n)}, (_, i) => xs.slice(i * n, i * n + n));

type RecordEntry = (path: string, action: 'create' | 'reuse' | 'skip', status?: string, reason?: string) => void;

/** Makes sure [parts] exists under [folder], one segment at a time. A segment already there under
 * either spelling is reused, so a stage folder under a workstream (DMPK/In vivo PK) lands in the
 * workstream; [counted] keeps a segment shared by two paths, or a workstream, from counting twice.
 * With [skip] (why no new folders), a missing segment and those below it are listed as skipped. */
async function ensurePath(folder: DG.Project | null, parts: string[], base: string, dryRun: boolean,
  counted: Set<string>, record: RecordEntry, skip: string | null = null): Promise<void> {
  let parent = folder;
  let path = base;
  for (const part of parts) {
    path = `${path}/${part}`;
    const found = (parent ? await subspaces(parent.id) : []).find((c) => isSpaceNamed(asRef(c), part)) ?? null;
    const next = found ?? (dryRun || !parent || skip ? null : await grok.dapi.spaces.id(parent.id).addSubspace(part));
    const key = normalizeSpaceName(path).toLowerCase();
    if (!counted.has(key)) {
      counted.add(key);
      if (found || !skip)
        record(path, found ? 'reuse' : 'create');
      else
        record(path, 'skip', 'ok', skip);
    }
    parent = next;
  }
}

/** Root folders no program claims, and folders of programs that have ended: a report only. */
export async function workspaceOrphans(code: string): Promise<WorkspaceReport> {
  const portfolio = await portfolioByCode(code);
  if (!portfolio.workspace_root)
    throw new Error(`portfolio ${code} has no workspaceRoot`);
  const report: WorkspaceReport = {portfolio: code, root: portfolio.workspace_root, dryRun: true, entries: [],
    created: 0, reused: 0, failed: 0, skipped: 0, notes: [], orphans: []};
  const root = await findRoot(portfolio.workspace_root);
  if (!root) {
    report.notes.push('the workspace is not built yet');
    return report;
  }
  const programs: WorkspaceProgram[] = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)}`, limit: ROW_LIMIT}))
    .map((p) => ({id: p.id, code: p.code, name: p.name, spaceId: p.space_id ?? null, status: p.status ?? null}));
  report.orphans = findOrphans(portfolio.workspace_root, (await subspaces(root.id)).map(asRef), programs);
  return report;
}

export function summarizeWorkspace(r: WorkspaceReport): string {
  if (r.entries.length === 0 && r.created + r.reused + r.failed === 0) {
    const lines = [`${r.portfolio} workspace "${r.root}": ${r.orphans.length} folder(s) for prune-or-park to look at; none touched` +
      (r.orphans.length ? ' (parkOrphans moves them)' : ''),
      ...r.notes.map((n) => `  ${n}`)];
    lines.push(...r.orphans.map((o) => `  ? ${o.path}: ${o.kind === 'park' ? 'park candidate, ' : ''}${o.reason}`));
    return lines.join('\n');
  }
  const head = `${r.portfolio} workspace "${r.root}"${r.dryRun ? ' (preview)' : ''}: ` +
    `${r.dryRun ? 'would create' : 'created'} ${r.created}, reused ${r.reused}, failed ${r.failed}` +
    (r.skipped ? `, skipped ${r.skipped}` : '');
  const lines = [head, ...r.notes.map((n) => `  ${n}`)];
  if (r.dryRun) {
    lines.push(...r.entries.filter((e) => e.action === 'create').map((e) => `  + ${e.path}`));
    lines.push(...r.entries.filter((e) => e.action === 'skip').map((e) => `  - ${e.path}: skipped, ${e.reason}`));
  }
  lines.push(...r.entries.filter((e) => e.status !== 'ok').map((e) => `  ! ${e.path}: ${e.status}`));
  if (r.orphans.length) {
    lines.push(`  Not touched (parking is a separate, explicit step: parkOrphans):`);
    lines.push(...r.orphans.map((o) => `  ? ${o.path}: ${o.kind === 'park' ? 'park candidate, ' : ''}${o.reason}`));
  }
  return lines.join('\n');
}

/** The folder at [parts] under [folderId], segment by segment under either spelling; null when a
 * segment is missing (the stage has not been built into the workspace yet). Never creates. */
export async function findSubfolder(folderId: string, parts: string[]): Promise<DG.Project | null> {
  let current: DG.Project | null = null;
  let parentId = folderId;
  for (const part of parts) {
    current = (await subspaces(parentId)).find((c) => isSpaceNamed(asRef(c), part)) ?? null;
    if (!current)
      return null;
    parentId = current.id;
  }
  return current;
}
