/** Provenance against the platform (F-28): finds the program(s) an object belongs to, reads their
 * records (decisions, evidence, runs and the tables their saved calls read and wrote, run
 * selections, documents, snapshots, criteria, the program dataset), and walks the graph. The
 * graph and the walk live in ./model/provenance. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {portalDb, ProgramCompoundRow} from './generated/db';
import {stagesOf} from './portal';
import {datasetsOf} from './program-data';
import {MAIN_DATASET} from './model/program-data';
import type {RecordedSource} from './model/evidence-tables';
import {
  Provenance, provenanceGraph, ProvenanceRecords, ProvenanceStart, PROVENANCE_STARTS, startKey, traverse,
} from './model/provenance';
import {newestFirstBy} from './model/time';
import {filterValue, inFilter, queryIn, ROW_LIMIT} from './rows';

const iso = (d: any) => d ? dayjs(d).toISOString() : '';

/** Runs whose saved calls are read for their tables; more are shown without them. */
const MAX_CALLS = 40;
/** Programs read for a compound held by many. */
const MAX_PROGRAMS = 10;

/** The programs [kind] [id] belongs to; a compound may be in several. */
async function programsOf(kind: ProvenanceStart, id: string): Promise<string[]> {
  const one = async (p: Promise<{program_id: string} | undefined>) => {
    const row = await p.catch(() => undefined);
    return row ? [row.program_id] : [];
  };
  switch (kind) {
  case 'program':
    return [id];
  case 'decision':
    return one(portalDb.stageTransitions.get(id));
  case 'evidence':
    return one(portalDb.evidences.get(id));
  case 'snapshot':
    return one(portalDb.programSnapshots.get(id));
  case 'document':
    return one(portalDb.documents.get(id));
  case 'model_run':
    return one(portalDb.modelRuns.query({filter: `id = ${filterValue(id)} or run_id = ${filterValue(id)}`, limit: ROW_LIMIT}).then((r) => r[0]));
  case 'compound': {
    // through the links, which are secured by program: only programs the reader sees
    const links = await compoundLinks(id);
    return [...new Set(links.map((l) => l.program_id))].slice(0, MAX_PROGRAMS);
  }
  }
}

const UUID = /^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i;

/** The program links of compound [idOrCode] (and [and]): by the code the link carries, or by the
 * compound row's id for a link from before 0.3.0 when the reader may read the compound table. */
async function compoundLinks(idOrCode: string, and?: string): Promise<ProgramCompoundRow[]> {
  const also = and ? ` and ${and}` : '';
  const key = UUID.test(idOrCode) ? `compound_id = ${filterValue(idOrCode)}` : `compound_code = ${filterValue(idOrCode)}`;
  const byCode = await portalDb.programCompounds.query({filter: `${key}${also}`, limit: ROW_LIMIT}).catch(() => []);
  if (byCode.length)
    return byCode;
  const c = (await portalDb.compounds.query({filter: `id = ${filterValue(idOrCode)} or code = ${filterValue(idOrCode)}`, limit: ROW_LIMIT}).catch(() => []))[0];
  return c ? portalDb.programCompounds.query({filter: `compound_id = ${filterValue(c.id)}${also}`, limit: ROW_LIMIT}).catch(() => []) : [];
}

/** The tables a saved call read and wrote: its dataframe parameters, held as the uploaded tables' ids. */
async function callTables(runId: string): Promise<{inputs: {id: string, name: string}[], outputs: {id: string, name: string}[]}> {
  const call = await grok.dapi.functions.calls.include('inputs, outputs').find(runId).catch(() => null);
  if (!call)
    return {inputs: [], outputs: []};
  const of = (which: 'inputs' | 'outputs') => call.func[which]
    .filter((p) => p.propertyType === DG.TYPE.DATA_FRAME)
    .map((p) => ({value: call[which][p.name], name: p.caption || p.name}))
    .map(({value, name}) => ({id: value instanceof DG.DataFrame ? value.getTableInfo()?.id ?? '' : String(value ?? ''), name}))
    .filter((t) => t.id && t.id !== 'null');
  return {inputs: of('inputs'), outputs: of('outputs')};
}

/** Everything the walk needs about [programIds]. [compound] is the compound asked about, if any. */
export async function provenanceRecords(programIds: string[], compound?: string): Promise<ProvenanceRecords> {
  const r: ProvenanceRecords = {programs: [], decisions: [], evidence: [], runs: [], selections: [], documents: [],
    snapshots: [], criteria: [], datasets: [], compounds: []};
  if (!programIds.length)
    return r;
  const by = inFilter('program_id', programIds)!;
  const programs = await queryIn(portalDb.programs, 'id', programIds);
  const [transitions, evidence, runs, artifacts, statuses, documents, snapshots] = await Promise.all([
    portalDb.stageTransitions.query({filter: by, limit: ROW_LIMIT}), portalDb.evidences.query({filter: by, limit: ROW_LIMIT}),
    portalDb.modelRuns.query({filter: by, limit: ROW_LIMIT}), portalDb.criterionArtifacts.query({filter: by, limit: ROW_LIMIT}),
    portalDb.criterionStatuses.query({filter: by, limit: ROW_LIMIT}), portalDb.documents.query({filter: by, limit: ROW_LIMIT}),
    portalDb.programSnapshots.query({filter: by, limit: ROW_LIMIT}),
  ]);
  const predictions = await portalDb.predictions.query({filter: by, limit: ROW_LIMIT}).catch(() => []);
  const stages = (await Promise.all([...new Set(programs.map((p) => p.stage_model_id))].map((m) => stagesOf(m)))).flat();
  const stageCode = (id?: string | null) => stages.find((s) => s.id === id)?.code ?? null;
  const criteria = await queryIn(portalDb.gateCriterions, 'stage_id', stages.map((s) => s.id));
  const modelIds = [...new Set(runs.map((x) => x.model_id).filter((x): x is string => !!x))];
  const models = await queryIn(portalDb.models, 'id', modelIds);
  const portfolios = await queryIn(portalDb.portfolios, 'id', programs.map((p) => p.portfolio_id));

  r.programs = programs.map((p) => ({id: p.id, code: p.code, name: p.name}));
  r.datasets = programs.flatMap((p) => {
    const pf = portfolios.find((x) => x.id === p.portfolio_id);
    const main = pf ? datasetsOf(pf).find((d) => d.code === MAIN_DATASET) : undefined;
    return main ? [{programId: p.id, name: main.name}] : [];
  });
  r.criteria = criteria.map((c) => ({id: c.id, code: c.code, name: c.name, stage: stageCode(c.stage_id) ?? undefined}));
  r.decisions = transitions.map((t) => ({id: t.id, programId: t.program_id, decision: t.decision,
    on: iso(t.decided_on ?? t.created_on), from: stageCode(t.from_stage_id), to: stageCode(t.to_stage_id)}));
  r.evidence = evidence.map((e) => ({id: e.id, decisionId: e.stage_transition_id ?? null, kind: e.kind, target: e.target,
    name: e.name, context: (e.context as Record<string, any> | undefined) ?? null}));
  const sorted = [...runs].sort(newestFirstBy((x) => x.created_on));
  r.runs = await Promise.all(sorted.map(async (x, i) => {
    const m = models.find((y) => y.id === x.model_id);
    const tables = x.run_id && i < MAX_CALLS ? await callTables(x.run_id) : {inputs: [], outputs: []};
    return {id: x.id, programId: x.program_id, runId: x.run_id ?? null, title: x.title ?? x.run_id ?? x.model_name,
      model: m ? {code: m.code, name: m.name} : {code: x.model_name, name: x.model_name}, ...tables};
  }));
  r.selections = artifacts.filter((a) => a.role === 'selected').flatMap((a) => {
    const s = statuses.find((x) => x.id === a.criterion_status_id);
    return s ? [{programId: a.program_id, criterionId: s.gate_criterion_id, runId: a.model_run_id}] : [];
  });
  r.documents = documents.map((d) => ({id: d.id, programId: d.program_id, title: d.title, version: d.doc_version ?? null,
    status: d.status, criterionId: d.gate_criterion_id ?? null, path: d.path,
    generatedFrom: ((d.generated as {sources?: RecordedSource[]} | undefined)?.sources ?? []).map((x) => ({kind: x.kind, ref: x.ref}))}));
  const scored = new Map<string, {runId: string, predicts: string, count: number}>();
  for (const p of predictions.filter((x) => x.model_run_id)) {
    const k = `${p.model_run_id}|${p.predicts}`;
    scored.set(k, {runId: p.model_run_id!, predicts: p.predicts, count: (scored.get(k)?.count ?? 0) + 1});
  }
  r.predictions = [...scored.values()];
  r.snapshots = snapshots.map((s) => ({id: s.id, programId: s.program_id, label: s.label, takenOn: iso(s.taken_on)}));
  if (compound) {
    const links = await compoundLinks(compound, by);
    r.compounds = links.map((l) => ({id: l.compound_id, code: l.compound_code ?? compound, programId: l.program_id, role: l.role}));
  }
  return r;
}

/** The provenance of one object: what it rests on and what rests on it, [depth] steps each way. */
export async function provenanceOf(kind: string, id: string, depth = 4): Promise<Provenance & {records: ProvenanceRecords}> {
  if (!(PROVENANCE_STARTS as string[]).includes(kind))
    throw new Error(`provenance starts from one of ${PROVENANCE_STARTS.join(', ')}, not "${kind}"`);
  const k = kind as ProvenanceStart;
  const programIds = await programsOf(k, id);
  if (!programIds.length)
    throw new Error(`no ${kind.replace('_', ' ')} ${id}, or it is not shared with you`);
  const records = await provenanceRecords(programIds, k === 'compound' ? id : undefined);
  const g = provenanceGraph(records);
  const start = startKey(g, records, k, id);
  if (!start)
    throw new Error(`no ${kind.replace('_', ' ')} ${id} among its program's records`);
  return {...traverse(g, start, depth), records};
}
