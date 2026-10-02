/** Saving a model run to a program: the call is saved, annotated, linked into the program's
 * workspace, and recorded as a `model_run` row a gate criterion can select.
 *
 * `model_run.func_call` is a soft reference to `Core.func_calls`, and the server refuses to write
 * it today (`column p.func_id does not exist`, the soft-ref alias bug; see the README). Every run
 * therefore also carries `run_id`, the call's id as plain text. The insert tries the soft
 * reference first and, when the server refuses it, writes the row without it, so the soft
 * reference starts being written on its own once the platform is fixed. Rows written without it
 * are migration debt: `func_call` is backfilled from `run_id` then. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {CriterionArtifactRow, ModelRow, ModelRunInsert, ModelRunRow, portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {portfolioByCode} from './portal';
import {programAccess} from './program-roles';
import {checkRunSelection, compareRuns, ComparisonRow, matchesRunFilter, ParamLabel, RunFilter, RunValues,
  SelectionEntry, selectionHistory} from './model/run-selection';
import {ProgramStatus} from './model/config';
import {findSubfolder} from './workspace';
import {RUN_OPTIONS, RunPurpose, runFolderPath, runOptions} from './model/model-runs';
import {programFolderName} from './model/workspace-plan';

export interface SaveModelRunOptions {
  /** Unit operation code (`unit_operation` catalogue). */
  unitOperation?: string;
  /** Run type code (`run_type` catalogue). */
  runType?: string;
  purpose?: RunPurpose;
  title?: string;
  description?: string;
  /** Folder under the program's folder, '/'-separated; default: the current stage's first folder. */
  folder?: string;
  /** Link the run into the program's workspace (default true). */
  link?: boolean;
  isReference?: boolean;
}

export interface SavedModelRun {
  call: DG.FuncCall;
  /** The saved call's id; the same as `modelRun.run_id`. */
  runId: string;
  modelRun: ModelRunRow;
  /** `func_call` when the soft reference was written, `run_id` when only the plain id could be. */
  recordedBy: 'func_call' | 'run_id';
  /** Why the soft reference was refused, when it was. */
  softRefError?: string;
  /** The space the run is linked into, with its path under the workspace root. */
  folder: {id: string, path: string} | null;
  /** Tables uploaded for the call's dataframe inputs and outputs. */
  tableIds: string[];
  notes: string[];
}

async function programByCode(portfolio: PortfolioRow, code: string): Promise<ProgramRow> {
  const p = (await portalDb.programs.query({filter: `portfolio_id = "${portfolio.id}" and code = "${code}"`}))[0];
  if (!p)
    throw new Error(`${portfolio.code} has no program ${code}, or it is not shared with you`);
  return p;
}

async function catalogueId(table: 'unitOperations' | 'runTypes', code: string | undefined, what: string):
  Promise<string | null> {
  if (!code)
    return null;
  const row = await portalDb[table].getByKey({code});
  if (!row)
    throw new Error(`no ${what} ${code} in the catalogue`);
  return row.id;
}

/** Uploads the call's dataframe inputs and outputs, as a saved call holds only references to
 * them, and shares them with the portfolio's viewers. A FuncCall has no ACL of its own. */
async function uploadTables(call: DG.FuncCall, portfolio: PortfolioRow, notes: string[]): Promise<string[]> {
  const ids: string[] = [];
  const viewers = portfolio.viewers ? await grok.dapi.groups.find(portfolio.viewers).catch(() => null) : null;
  const params = [...call.inputParams.values(), ...call.outputParams.values()] as DG.FuncCallParam[];
  for (const p of params) {
    if (p.property.propertyType !== DG.TYPE.DATA_FRAME || !(p.value instanceof DG.DataFrame))
      continue;
    const df = p.value as DG.DataFrame;
    await grok.dapi.tables.uploadDataFrame(df);
    const info = df.getTableInfo();
    ids.push(info.id);
    if (viewers) {
      await grok.dapi.permissions.grant(info, viewers, false)
        .catch((e: any) => notes.push(`table ${p.name} not shared with ${viewers.friendlyName}: ${e?.message ?? e}`));
    }
  }
  return ids;
}

/** The space to link the run into: the asked-for or current-stage folder under the program's
 * folder when it exists, else the program folder; null when the program has no workspace. */
async function targetFolder(portfolio: PortfolioRow, program: ProgramRow, stageFolders: string | null,
  asked: string | undefined, notes: string[]): Promise<{id: string, path: string} | null> {
  if (!program.space_id) {
    notes.push(`${program.code} has no workspace folder yet (buildWorkspace ${portfolio.code}); the run is not linked`);
    return null;
  }
  const base = `${portfolio.workspace_root ?? portfolio.code}/${programFolderName(program)}`;
  const parts = runFolderPath(stageFolders, asked);
  if (parts.length === 0)
    return {id: program.space_id, path: base};
  const sub = await findSubfolder(program.space_id, parts);
  if (sub)
    return {id: sub.id, path: `${base}/${parts.join('/')}`};
  notes.push(`${parts.join('/')} is not in ${program.code}'s workspace; linked to the program folder`);
  return {id: program.space_id, path: base};
}

/** Inserts the `model_run` row: with the soft reference first, by the plain run id when the
 * server refuses it. An existing row for the same run is returned as it is. */
async function recordRun(values: Omit<ModelRunInsert, 'func_call'> & {run_id: string}):
  Promise<{row: ModelRunRow, recordedBy: 'func_call' | 'run_id', softRefError?: string}> {
  const existing = (await portalDb.modelRuns.query({filter: `run_id = "${values.run_id}"`}))[0];
  if (existing)
    return {row: existing, recordedBy: existing.func_call ? 'func_call' : 'run_id'};
  try {
    const [row] = await portalDb.modelRuns.insert({...values, func_call: values.run_id} as any);
    return {row: await portalDb.modelRuns.get(row.id), recordedBy: 'func_call'};
  } catch (e: any) {
    const softRefError = String(e?.message ?? e);
    const [row] = await portalDb.modelRuns.insert(values as any);
    grok.log.warning(`Project Portal: model run ${values.run_id} recorded by run id only; ` +
      `the soft reference to Core.func_calls was refused: ${softRefError}`);
    return {row: await portalDb.modelRuns.get(row.id), recordedBy: 'run_id', softRefError};
  }
}

/** Saves [call], which has run, to a program of a portfolio: annotates it (portfolio, program,
 * stage at run, unit operation, run type, purpose), uploads its tables, saves it, links it into
 * the program's workspace and records it as a `model_run`. */
export async function saveModelRun(call: DG.FuncCall, portfolioCode: string, programCode: string,
  options: SaveModelRunOptions = {}): Promise<SavedModelRun> {
  const notes: string[] = [];
  const portfolio = await portfolioByCode(portfolioCode);
  const program = await programByCode(portfolio, programCode);
  const known = call.id ? (await portalDb.modelRuns.query({filter: `run_id = "${call.id}"`}))[0] : undefined;
  if (known) {
    if (known.program_id !== program.id)
      throw new Error(`run ${call.id} is already recorded for another program`);
    return {call, runId: call.id, modelRun: known, recordedBy: known.func_call ? 'func_call' : 'run_id',
      folder: null, tableIds: [], notes: ['already recorded; nothing changed']};
  }
  const stage = program.stage_id ? await portalDb.stages.get(program.stage_id) : null;
  const unitOperationId = await catalogueId('unitOperations', options.unitOperation, 'unit operation');
  const runTypeId = await catalogueId('runTypes', options.runType, 'run type');
  const func = call.func;
  const model = await portalDb.models.getByKey({code: func.nqName});
  if (!model)
    notes.push(`${func.nqName} is not in the model catalogue; the run has no model_id`);

  const title = options.title ?? `${func.friendlyName ?? func.name} · ${program.code}`;
  const tags = Array.isArray(call.options['tags']) ? call.options['tags'] as string[] : [];
  const annotations = runOptions({portfolio: portfolio.code, program: program.code, stage: stage?.code,
    unitOperation: options.unitOperation, runType: options.runType, purpose: options.purpose, title,
    description: options.description}, tags);
  for (const [k, v] of Object.entries(annotations))
    call.options[k] = v;

  const tableIds = await uploadTables(call, portfolio, notes);
  // the uploaded tables, so whoever removes the run can remove them too
  if (tableIds.length)
    call.options[RUN_OPTIONS.tables] = tableIds;
  const saved = await grok.dapi.functions.calls.save(call);

  const folder = options.link === false ? null :
    await targetFolder(portfolio, program, stage?.folders ?? null, options.folder, notes);
  if (folder)
    await grok.dapi.spaces.id(folder.id).addEntity(saved.id, true);

  const started = saved.started ?? call.started;
  const recorded = await recordRun({
    run_id: saved.id,
    program_id: program.id,
    model_id: model?.id,
    model_name: func.nqName,
    model_version: (func as any).package?.version ?? undefined,
    title,
    stage_id: stage?.id,
    unit_operation_id: unitOperationId ?? undefined,
    run_type_id: runTypeId ?? undefined,
    purpose: options.purpose,
    is_reference: options.isReference ?? false,
    run_author: DG.User.current().id,
    run_started: started ? dayjs(started as any).toISOString() : undefined,
    description: options.description,
  });
  if (recorded.recordedBy === 'run_id')
    notes.push('recorded by run id: the server refused the soft reference to the call (migration debt)');
  return {call: saved, runId: saved.id, modelRun: recorded.row, recordedBy: recorded.recordedBy,
    softRefError: recorded.softRefError, folder, tableIds, notes};
}

/** The runs recorded for a program, newest first. */
export async function modelRunsOf(programId: string): Promise<ModelRunRow[]> {
  const runs = await portalDb.modelRuns.query({filter: `program_id = "${programId}"`});
  return runs.sort((a, b) => b.created_on.valueOf() - a.created_on.valueOf());
}

/** A saved run as the program page lists it: the row and the names it refers to. */
export interface SavedRunItem {
  run: ModelRunRow;
  /** The catalogued model; null when the run's function is not catalogued. */
  model: ModelRow | null;
  stage: string | null;
  unitOperation: string | null;
  runType: string | null;
  /** The catalogue codes, which a criterion's run filter names. */
  unitOperationCode?: string | null;
  runTypeCode?: string | null;
}

/** The program's runs, newest first, with the model, stage, unit operation and run type by name.
 * Each catalogue table is read once, not once per run. */
export async function savedRunsOf(programId: string): Promise<SavedRunItem[]> {
  const runs = await modelRunsOf(programId);
  if (runs.length === 0)
    return [];
  const ids = (pick: (r: ModelRunRow) => string | undefined) => [...new Set(runs.map(pick).filter((x): x is string => !!x))];
  const byIds = async <T extends {id: string}>(table: {query(o: any): Promise<T[]>}, list: string[]) =>
    new Map((list.length ? await table.query({filter: list.map((id) => `id = "${id}"`).join(' or ')}) : [])
      .map((x) => [x.id, x] as const));
  const [models, stages, ops, types] = await Promise.all([
    byIds(portalDb.models, ids((r) => r.model_id)), byIds(portalDb.stages, ids((r) => r.stage_id)),
    byIds(portalDb.unitOperations, ids((r) => r.unit_operation_id)), byIds(portalDb.runTypes, ids((r) => r.run_type_id))]);
  return runs.map((run) => ({run, model: run.model_id ? models.get(run.model_id) ?? null : null,
    stage: run.stage_id ? stages.get(run.stage_id)?.code ?? null : null,
    unitOperation: run.unit_operation_id ? ops.get(run.unit_operation_id)?.name ?? null : null,
    runType: run.run_type_id ? types.get(run.run_type_id)?.name ?? null : null,
    unitOperationCode: run.unit_operation_id ? ops.get(run.unit_operation_id)?.code ?? null : null,
    runTypeCode: run.run_type_id ? types.get(run.run_type_id)?.code ?? null : null}));
}

/** Selects a run for a program's criterion: the run becomes the `selected` artifact, a run
 * selected before is `superseded`, and the criterion is met unless it was waived. Selection
 * instead of approval: the rationale says why this run over the others. */
export async function selectRunForCriterion(programId: string, criterionId: string, modelRunId: string,
  rationale: string): Promise<CriterionArtifactRow> {
  const run = await portalDb.modelRuns.get(modelRunId);
  if (run.program_id !== programId)
    throw new Error('the run belongs to another program');
  // asked again here, not only by the gate panel
  const program = await portalDb.programs.get(programId);
  const portfolio = await portalDb.portfolios.get(program.portfolio_id);
  const access = await programAccess(portfolio, programId);
  const refused = checkRunSelection({isApprover: access.approver, isContributor: access.contributor,
    programStatus: (program.status ?? 'active') as ProgramStatus, rationale});
  if (refused.length)
    throw new Error(refused.join('; '));
  let status = (await portalDb.criterionStatuses.query(
    {filter: `program_id = "${programId}" and gate_criterion_id = "${criterionId}"`}))[0];
  if (!status)
    status = await portalDb.criterionStatuses.get((await portalDb.criterionStatuses.insert(
      {program_id: programId, gate_criterion_id: criterionId, status: 'in_progress'}))[0].id);
  const artifacts = await portalDb.criterionArtifacts.query({filter: `criterion_status_id = "${status.id}"`});
  for (const a of artifacts.filter((a) => a.role === 'selected' && a.model_run_id !== modelRunId))
    await portalDb.criterionArtifacts.update(a.id, {role: 'superseded'});
  const now = dayjs();
  const fields = {role: 'selected' as const, selected_by: DG.User.current().id, selected_on: now, rationale};
  let id = artifacts.find((a) => a.model_run_id === modelRunId)?.id;
  if (id)
    await portalDb.criterionArtifacts.update(id, fields);
  else {
    id = (await portalDb.criterionArtifacts.insert({program_id: programId, criterion_status_id: status.id,
      model_run_id: modelRunId, ...fields}))[0].id;
  }
  if (status.status !== 'waived') {
    await portalDb.criterionStatuses.update(status.id, {status: 'met',
      evidence: `run "${run.title ?? run.run_id}" (${run.run_id ?? run.func_call}) selected: ${rationale}`});
  }
  return await portalDb.criterionArtifacts.get(id);
}

/** What the gate panel and the choose-run dialog show for a criterion that takes a run. */
export interface CriterionRuns {
  /** The program's runs the criterion's filter admits, newest first. */
  candidates: SavedRunItem[];
  /** Runs of the program the filter leaves out. */
  excluded: number;
  /** The model_run id of the run selected now, if any. */
  selectedRunId: string | null;
  selected: SelectionEntry | null;
  superseded: SelectionEntry[];
  /** Why the current user may not select (rationale aside); empty when they may. */
  refused: string[];
}

const itemCodes = (i: SavedRunItem) => ({model: i.run.model_name, unitOperation: i.unitOperationCode ?? null,
  runType: i.runTypeCode ?? null});

/** The candidates, the current choice and the history of choices of one criterion of a program. */
export async function criterionRuns(programId: string, criterionId: string): Promise<CriterionRuns> {
  const [program, criterion, items] = await Promise.all([portalDb.programs.get(programId),
    portalDb.gateCriterions.get(criterionId), savedRunsOf(programId)]);
  const portfolio = await portalDb.portfolios.get(program.portfolio_id);
  const filter = (criterion.run_filter as RunFilter | undefined) ?? null;
  const candidates = items.filter((i) => matchesRunFilter(itemCodes(i), filter));
  const status = (await portalDb.criterionStatuses.query(
    {filter: `program_id = "${programId}" and gate_criterion_id = "${criterionId}"`}))[0];
  const artifacts = status ? await portalDb.criterionArtifacts.query({filter: `criterion_status_id = "${status.id}"`}) : [];
  const titles = new Map(items.map((i) => [i.run.id, i.run.title ?? i.run.run_id ?? i.run.model_name]));
  const users = new Map<string, string>();
  for (const id of new Set(artifacts.map((a) => a.selected_by).filter((x): x is string => !!x)))
    users.set(id, (await grok.dapi.users.find(id).catch(() => null))?.friendlyName ?? '');
  const entries: SelectionEntry[] = artifacts.map((a) => ({runTitle: titles.get(a.model_run_id) ?? a.model_run_id,
    role: a.role, selectedOn: a.selected_on ? dayjs(a.selected_on as any).toISOString() : null,
    selectedBy: a.selected_by ? users.get(a.selected_by) || null : null, rationale: a.rationale ?? null}));
  const history = selectionHistory(entries);
  const access = await programAccess(portfolio, programId);
  return {candidates, excluded: items.length - candidates.length,
    selectedRunId: artifacts.find((a) => a.role === 'selected')?.model_run_id ?? null, ...history,
    refused: checkRunSelection({isApprover: access.approver, isContributor: access.contributor,
      programStatus: (program.status ?? 'active') as ProgramStatus})};
}

/** The runs side by side: their inputs and outputs, read back from the saved calls, labelled from
 * the function's parameters. A run whose call cannot be loaded shows empty values. */
export async function runComparison(items: SavedRunItem[]): Promise<ComparisonRow[]> {
  const calls = await Promise.all(items.map((i) => i.run.run_id ?
    grok.dapi.functions.calls.include('inputs, outputs').find(i.run.run_id).catch(() => null) : Promise.resolve(null)));
  // a saved call holds a dataframe parameter as the uploaded table's id, not the table
  const value = (p: DG.Property, v: unknown) => p.propertyType === DG.TYPE.DATA_FRAME && !(v instanceof DG.DataFrame) ?
    (v ? 'saved table' : null) : v;
  const read = (c: DG.FuncCall | null, which: 'inputs' | 'outputs') => !c ? {} :
    Object.fromEntries(c.func[which].map((p) => [p.name, value(p, c[which][p.name])]));
  const values: RunValues[] = items.map((i, k) => ({title: i.run.title ?? i.run.run_id ?? '',
    inputs: read(calls[k], 'inputs'), outputs: read(calls[k], 'outputs')}));
  const func = calls.find((c) => c)?.func ?? null;
  const label = (p: DG.Property): ParamLabel => ({name: p.name, caption: p.caption || null,
    units: p.units || p.options?.['units'] || null});
  return compareRuns(values, func ? {inputs: func.inputs.map(label), outputs: func.outputs.map(label)} : {});
}

export function summarizeSavedRun(r: SavedModelRun): string {
  return [`run ${r.runId} saved: model_run ${r.modelRun.id} (by ${r.recordedBy})` +
    (r.folder ? `, linked into ${r.folder.path}` : ''), ...r.notes.map((n) => `  ${n}`)].join('\n');
}
