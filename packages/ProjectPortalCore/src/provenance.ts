/** Provenance against the platform (F-28): finds the program(s) an object belongs to, reads their
 * records (decisions, evidence, runs and the tables their saved calls read and wrote, run
 * selections, documents, snapshots, criteria, the program dataset), and walks the graph. The
 * graph and the walk live in ./model/provenance. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {portalDb} from './generated/db';
import {stagesOf} from './portal';
import {datasetsOf} from './program-data';
import {MAIN_DATASET} from './model/program-data';
import {
  Provenance, provenanceGraph, ProvenanceRecords, ProvenanceStart, PROVENANCE_STARTS, startKey, traverse,
} from './model/provenance';

const inList = (ids: string[]) => ids.map((id) => `"${id}"`).join(', ');
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
    return one(portalDb.modelRuns.query({filter: `id = "${id}" or run_id = "${id}"`}).then((r) => r[0]));
  case 'compound': {
    const c = (await portalDb.compounds.query({filter: `id = "${id}" or code = "${id}"`}).catch(() => []))[0];
    if (!c)
      return [];
    const links = await portalDb.programCompounds.query({filter: `compound_id = "${c.id}"`});
    return [...new Set(links.map((l) => l.program_id))].slice(0, MAX_PROGRAMS);
  }
  }
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
  const by = `program_id in (${inList(programIds)})`;
  const programs = await portalDb.programs.query({filter: `id in (${inList(programIds)})`});
  const [transitions, evidence, runs, artifacts, statuses, documents, snapshots] = await Promise.all([
    portalDb.stageTransitions.query({filter: by}), portalDb.evidences.query({filter: by}),
    portalDb.modelRuns.query({filter: by}), portalDb.criterionArtifacts.query({filter: by}),
    portalDb.criterionStatuses.query({filter: by}), portalDb.documents.query({filter: by}),
    portalDb.programSnapshots.query({filter: by}),
  ]);
  const predictions = await portalDb.predictions.query({filter: by}).catch(() => []);
  const stages = (await Promise.all([...new Set(programs.map((p) => p.stage_model_id))].map((m) => stagesOf(m)))).flat();
  const stageCode = (id?: string | null) => stages.find((s) => s.id === id)?.code ?? null;
  const criteria = stages.length ? await portalDb.gateCriterions.query({filter: `stage_id in (${inList(stages.map((s) => s.id))})`}) : [];
  const modelIds = [...new Set(runs.map((x) => x.model_id).filter((x): x is string => !!x))];
  const models = modelIds.length ? await portalDb.models.query({filter: `id in (${inList(modelIds)})`}) : [];
  const portfolios = await portalDb.portfolios.query({filter: `id in (${inList([...new Set(programs.map((p) => p.portfolio_id))])})`});

  r.programs = programs.map((p) => ({id: p.id, code: p.code, name: p.name}));
  r.datasets = programs.flatMap((p) => {
    const main = datasetsOf(portfolios.find((x) => x.id === p.portfolio_id)!).find((d) => d.code === MAIN_DATASET);
    return main ? [{programId: p.id, name: main.name}] : [];
  });
  r.criteria = criteria.map((c) => ({id: c.id, code: c.code, name: c.name, stage: stageCode(c.stage_id) ?? undefined}));
  r.decisions = transitions.map((t) => ({id: t.id, programId: t.program_id, decision: t.decision,
    on: iso(t.decided_on ?? t.created_on), from: stageCode(t.from_stage_id), to: stageCode(t.to_stage_id)}));
  r.evidence = evidence.map((e) => ({id: e.id, decisionId: e.stage_transition_id ?? null, kind: e.kind, target: e.target,
    name: e.name, context: (e.context as Record<string, any> | undefined) ?? null}));
  const sorted = [...runs].sort((a, b) => String(b.created_on).localeCompare(String(a.created_on)));
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
    status: d.status, criterionId: d.gate_criterion_id ?? null, path: d.path}));
  const scored = new Map<string, {runId: string, predicts: string, count: number}>();
  for (const p of predictions.filter((x) => x.model_run_id)) {
    const k = `${p.model_run_id}|${p.predicts}`;
    scored.set(k, {runId: p.model_run_id!, predicts: p.predicts, count: (scored.get(k)?.count ?? 0) + 1});
  }
  r.predictions = [...scored.values()];
  r.snapshots = snapshots.map((s) => ({id: s.id, programId: s.program_id, label: s.label, takenOn: iso(s.taken_on)}));
  if (compound) {
    const c = (await portalDb.compounds.query({filter: `id = "${compound}" or code = "${compound}"`}).catch(() => []))[0];
    if (c) {
      const links = await portalDb.programCompounds.query({filter: `compound_id = "${c.id}" and ${by}`});
      r.compounds = links.map((l) => ({id: c.id, code: c.code, programId: l.program_id, role: l.role}));
    }
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
