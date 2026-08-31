/* eslint-disable camelcase */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import type {BuilderType} from '../escher_src/src/Builder';
import type {CobraModelData, CobraReaction, CobraMetabolite, CobraGene, ReactionBounds} from '../escher_src/src/ts/types';
import {MetabolicGraphView} from './view';
import type {MetabolicAIContext} from './view';
import {WorkerCobraSolver} from './cobra';
import {sampleFluxAveragesWasm, samplesToDataFrame, loadCampaigns, saveAnalysisState,
  analysisExists, readAnalysis, loadStateProxy, getCurrentAnalysis} from './utils';
import {clearTimeCourseSlider, isTimeCourseActive, applyStepDistribution, runTimeCourseSampling} from './timeCourse';
import type {TimeCourseSamplingParams} from './timeCourse';

/** Reaction bounds as the AI passes them: a [lower, upper] pair or an object in either naming style. */
export type ReactionBoundsInput = [number, number] |
  {lower_bound: number, upper_bound: number} | {lowerBound: number, upperBound: number};
export type BoundsMap = {[reactionId: string]: ReactionBoundsInput};
export type ObjectivesMap = {[reactionId: string]: string}; // 'maximize' | 'minimize', validated at runtime
export type FluxMap = {[reactionId: string]: number | string}; // string values are coerced

const DEFAULT_FIND_LIMIT = 15;
const DEFAULT_SAMPLES = 10000;
const MAX_SAMPLES = 100000;

function samplingParamError(samples: number, bins: number, thinning: number): string | null {
  if (samples < 10 || samples > MAX_SAMPLES)
    return `samples must be between 10 and ${MAX_SAMPLES}`;
  if (bins < 2 || bins > 500)
    return 'bins must be between 2 and 500';
  if (thinning < 1 || thinning > 1000)
    return 'thinning must be between 1 and 1000';
  return null;
}

/** The injected view is either the Dart wrapper ({@link DG.View}, unwrapped via `jsView`)
 * or the {@link MetabolicGraphView} itself. */
function mgCtx(view: DG.ViewBase): MetabolicAIContext {
  const v = (view as DG.View)?.jsView ?? view;
  const ctx = v instanceof MetabolicGraphView ? v.aiContext() : null;
  if (!ctx)
    throw new Error('The current view is not a Metabolic Graph view (or its map is still loading)');
  return ctx;
}

function requireModel(builder: BuilderType): CobraModelData {
  const model = builder.model_data;
  if (!model?.reactions?.length)
    throw new Error('No COBRA model is loaded in this view — only the map drawing is shown');
  return model;
}

function requireMap(builder: BuilderType) {
  if (!builder.map)
    throw new Error('No map is loaded in this view');
  return builder.map;
}

const norm = (s: string | null | undefined) => (s ?? '').toLowerCase();

/** Case-insensitive lookup by id, bigg_id or name in a model entity list. */
function findEntity<T extends {id?: string, bigg_id?: string, name?: string | null}>(
  list: T[] | undefined, id: string): T | null {
  const q = norm(id);
  return list?.find((e) => norm(e.id) === q || norm(e.bigg_id) === q || norm(e.name) === q) ?? null;
}

const findReaction = (model: CobraModelData, id: string): CobraReaction | null => findEntity(model.reactions, id);
const findMetabolite = (model: CobraModelData, id: string): CobraMetabolite | null => findEntity(model.metabolites, id);
const findGene = (model: CobraModelData, id: string): CobraGene | null => findEntity(model.genes, id);

/** Map key of the reaction drawn on the map matching any of [ids], or null. */
function mapReactionKey(builder: BuilderType, ...ids: (string | null | undefined)[]): string | null {
  const reactions = builder.map?.reactions;
  if (!reactions)
    return null;
  const qs = ids.filter(Boolean).map((s) => norm(s));
  for (const k of Object.keys(reactions)) {
    const r = reactions[k];
    if (qs.some((q) => norm(r.bigg_id) === q || norm(r.name) === q))
      return k;
  }
  return null;
}

/** Map key of the metabolite node matching any of [ids] (primary nodes preferred), or null. */
function mapNodeKey(builder: BuilderType, ...ids: (string | null | undefined)[]): string | null {
  const nodes = builder.map?.nodes;
  if (!nodes)
    return null;
  const qs = ids.filter(Boolean).map((s) => norm(s));
  let fallback: string | null = null;
  for (const k of Object.keys(nodes)) {
    const n = nodes[k];
    if (n.node_type !== 'metabolite')
      continue;
    if (qs.some((q) => norm(n.bigg_id) === q || norm(n.name) === q)) {
      if (n.node_is_primary)
        return k;
      fallback ??= k;
    }
  }
  return fallback;
}

function reactionData(builder: BuilderType): {[id: string]: number | string | null} | null {
  return builder.settings.get('reaction_data') as {[id: string]: number | string | null} | null;
}

function reactionFlux(builder: BuilderType, r: {id?: string, bigg_id?: string, name?: string | null}): number | null {
  const data = reactionData(builder);
  const v = (r.id ? data?.[r.id] : null) ?? (r.bigg_id ? data?.[r.bigg_id] : null) ?? (r.name ? data?.[r.name] : null);
  return typeof v === 'number' ? v : null;
}

function toBounds(value: ReactionBoundsInput | unknown): ReactionBounds | null {
  if (Array.isArray(value) && value.length === 2 && value.every((x) => typeof x === 'number'))
    return {lower_bound: value[0], upper_bound: value[1]};
  const o = value as {lower_bound?: unknown, upper_bound?: unknown, lowerBound?: unknown, upperBound?: unknown} | null;
  const lb = o?.lower_bound ?? o?.lowerBound;
  const ub = o?.upper_bound ?? o?.upperBound;
  return typeof lb === 'number' && typeof ub === 'number' ? {lower_bound: lb, upper_bound: ub} : null;
}

export function getMetabolicState(view: DG.ViewBase) {
  const {builder} = mgCtx(view);
  const model = builder.model_data;
  const map = builder.map;
  const data = reactionData(builder);
  const numeric = Object.values(data ?? {}).filter((x): x is number => typeof x === 'number');
  const selected = map ?
    [...new Set(Object.values(map.getSelectedNodes()).map((n) => n?.bigg_id).filter(Boolean))] : [];
  return {
    model: model ? {
      id: model.id ?? null,
      name: model.name ?? null,
      reactions: model.reactions?.length ?? 0,
      metabolites: model.metabolites?.length ?? 0,
      genes: model.genes?.length ?? 0,
    } : null,
    map: map ? {
      reactions: Object.keys(map.reactions ?? {}).length,
      metabolites: Object.values(map.nodes ?? {}).filter((n) => n?.node_type === 'metabolite').length,
    } : null,
    objectives: (model?.reactions ?? []).filter((r) => !!r.objective_coefficient)
      .map((r) => ({reactionId: r.id, direction: r.objective_coefficient! > 0 ? 'maximize' : 'minimize'})),
    reactionData: data ? {
      reactions: Object.keys(data).length,
      ...(numeric.length ? {min: Math.min(...numeric), max: Math.max(...numeric)} : {}),
    } : null,
    samplingHistograms: (builder.reaction_sampling_distribution?.data?.size ?? 0) > 0,
    selectedMetabolites: selected,
    timeCourseActive: isTimeCourseActive(),
    currentAnalysis: getCurrentAnalysis(builder),
  };
}

/** The properties {@link findMetabolicEntities} reads off a match — all three cobra entity
 * kinds fit it structurally, and the map-only branch fabricates minimal literals of it. */
type EntityBrief = {
  id?: string, bigg_id?: string, name?: string | null,
  lower_bound?: number, upper_bound?: number, objective_coefficient?: number,
  formula?: string, compartment?: string,
};

export function findMetabolicEntities(view: DG.ViewBase, query: string, kind?: string, limit?: number) {
  const {builder} = mgCtx(view);
  const wanted = kind ? norm(kind) : null;
  if (wanted && !['reaction', 'metabolite', 'gene'].includes(wanted))
    return {success: false, error: `Unknown kind '${kind}' — use reaction, metabolite or gene`};
  const q = norm(query);
  const hit = (...hs: (string | null | undefined)[]) => hs.some((h) => norm(h).includes(q));

  type Row = {kind: 'reaction' | 'metabolite' | 'gene', entity: EntityBrief};
  const rows: Row[] = [];
  const model = builder.model_data;
  if (model?.reactions?.length) {
    if (!wanted || wanted === 'reaction') {
      for (const r of model.reactions)
        if (hit(r.id, r.bigg_id, r.name)) rows.push({kind: 'reaction', entity: r});
    }
    if (!wanted || wanted === 'metabolite') {
      for (const m of model.metabolites ?? [])
        if (hit(m.id, m.bigg_id, m.name)) rows.push({kind: 'metabolite', entity: m});
    }
    if (!wanted || wanted === 'gene') {
      for (const g of model.genes ?? [])
        if (hit(g.id, g.bigg_id, g.name)) rows.push({kind: 'gene', entity: g});
    }
  } else if (builder.map) {
    // map-only view (no model): search what is drawn
    if (!wanted || wanted === 'reaction') {
      for (const k of Object.keys(builder.map.reactions ?? {})) {
        const r = builder.map.reactions[k];
        if (hit(r.bigg_id, r.name)) rows.push({kind: 'reaction', entity: {id: r.bigg_id, name: r.name}});
      }
    }
    if (!wanted || wanted === 'metabolite') {
      const seen = new Set<string>();
      for (const k of Object.keys(builder.map.nodes ?? {})) {
        const n = builder.map.nodes[k];
        if (n.node_type === 'metabolite' && !seen.has(n.bigg_id) && hit(n.bigg_id, n.name)) {
          seen.add(n.bigg_id);
          rows.push({kind: 'metabolite', entity: {id: n.bigg_id, name: n.name}});
        }
      }
    }
  }

  const max = limit && limit > 0 ? limit : DEFAULT_FIND_LIMIT;
  return {
    total: rows.length,
    matches: rows.slice(0, max).map(({kind: k, entity: e}) => {
      const flux = k === 'reaction' && model ? reactionFlux(builder, e) : null;
      return {
        kind: k,
        id: e.id ?? e.bigg_id,
        ...(e.name && e.name !== (e.id ?? e.bigg_id) ? {name: e.name} : {}),
        ...(k === 'reaction' && typeof e.lower_bound === 'number' ?
          {lowerBound: e.lower_bound, upperBound: e.upper_bound} : {}),
        ...(k === 'reaction' && e.objective_coefficient ?
          {objective: e.objective_coefficient > 0 ? 'maximize' : 'minimize'} : {}),
        ...(flux != null ? {flux} : {}),
        ...(k === 'metabolite' && e.formula ? {formula: e.formula} : {}),
        ...(k === 'metabolite' && e.compartment ? {compartment: e.compartment} : {}),
      };
    }),
  };
}

export function getMetabolicEntityDetails(view: DG.ViewBase, entityId: string) {
  const {builder} = mgCtx(view);
  const model = builder.model_data;
  if (!model?.reactions?.length) {
    // map-only view (no model): report what the drawing knows
    const rk = mapReactionKey(builder, entityId);
    if (rk != null) {
      const mr = builder.map!.reactions[rk];
      return {kind: 'reaction', id: mr.bigg_id, name: mr.name, geneRule: mr.gene_reaction_rule || undefined,
        metabolites: Object.fromEntries((mr.metabolites ?? []).map((m) => [m.bigg_id, m.coefficient])),
        onMap: true, note: 'No COBRA model is loaded — only map information is available'};
    }
    const nk = mapNodeKey(builder, entityId);
    if (nk != null) {
      const mn = builder.map!.nodes[nk];
      return {kind: 'metabolite', id: mn.bigg_id, name: mn.name, onMap: true,
        note: 'No COBRA model is loaded — only map information is available'};
    }
    return {success: false, error: `No reaction or metabolite '${entityId}' on the map — call findMetabolicEntities`};
  }

  const r = findReaction(model, entityId);
  if (r) {
    return {
      kind: 'reaction',
      id: r.id,
      name: r.name,
      metabolites: r.metabolites ?? {}, // metabolite id -> stoichiometric coefficient (negative = consumed)
      lowerBound: r.lower_bound ?? null,
      upperBound: r.upper_bound ?? null,
      reversibility: r.reversibility ?? null,
      ...(r.subsystem ? {subsystem: r.subsystem} : {}),
      ...(r.gene_reaction_rule ? {geneRule: r.gene_reaction_rule} : {}),
      ...(r.objective_coefficient ? {objective: r.objective_coefficient > 0 ? 'maximize' : 'minimize'} : {}),
      ...(reactionFlux(builder, r) != null ? {flux: reactionFlux(builder, r)} : {}),
      onMap: mapReactionKey(builder, r.id, r.bigg_id, r.name) != null,
    };
  }

  const m = findMetabolite(model, entityId);
  if (m) {
    const mid: string = m.id ?? m.bigg_id;
    const qm = norm(mid);
    const reactions = model.reactions
      .filter((rx) => Object.keys(rx.metabolites ?? {}).some((k) => norm(k) === qm))
      .map((rx) => ({id: rx.id, coefficient: Object.entries(rx.metabolites)
        .find(([k]) => norm(k) === qm)![1]}));
    return {
      kind: 'metabolite',
      id: mid,
      name: m.name,
      ...(m.formula ? {formula: m.formula} : {}),
      ...(m.compartment ? {compartment: m.compartment} : {}),
      reactions, // negative coefficient = consumed by the reaction, positive = produced
      onMap: mapNodeKey(builder, mid, m.bigg_id, m.name) != null,
    };
  }

  const g = findGene(model, entityId);
  if (g) {
    const gid: string = g.id ?? g.bigg_id;
    return {
      kind: 'gene',
      id: gid,
      name: g.name,
      reactions: model.reactions
        .filter((rx) => (rx.gene_reaction_rule ?? '').includes(gid) ||
          rx.genes?.some((gg) => gg.bigg_id === gid || gg.id === gid))
        .map((rx) => rx.id),
    };
  }

  return {success: false, error: `No reaction, metabolite or gene '${entityId}' — call findMetabolicEntities to search`};
}

export function setReactionBounds(view: DG.ViewBase, bounds: BoundsMap) {
  const {builder} = mgCtx(view);
  const model = requireModel(builder);
  const entries = Object.entries(bounds ?? {});
  if (!entries.length)
    return {success: false, error: 'Pass bounds as a map of reaction id to [lowerBound, upperBound]'};
  const parsed: {[id: string]: ReactionBounds} = {};
  const unknown: string[] = [];
  for (const [id, value] of entries) {
    const b = toBounds(value);
    if (!b)
      return {success: false, error: `Invalid bounds for '${id}' — pass [lowerBound, upperBound] numbers`};
    if (b.lower_bound > b.upper_bound)
      return {success: false, error: `Lower bound above upper bound for '${id}'`};
    const r = findReaction(model, id);
    if (!r) {
      unknown.push(id);
      continue;
    }
    parsed[r.id] = b;
  }
  if (!Object.keys(parsed).length)
    return {success: false, error: `No matching reactions: ${unknown.join(', ')} — call findMetabolicEntities for ids`};
  clearTimeCourseSlider();
  builder.set_reaction_bounds(parsed);
  return {success: true, applied: Object.keys(parsed), ...(unknown.length ? {unknownReactions: unknown} : {})};
}

export async function runFba(view: DG.ViewBase, objectives?: ObjectivesMap) {
  const {builder} = mgCtx(view);
  const model = requireModel(builder);
  if (objectives && Object.keys(objectives).length) {
    const resolved: {r: CobraReaction, coef: number}[] = [];
    for (const [id, dir] of Object.entries(objectives)) {
      const r = findReaction(model, id);
      if (!r)
        return {success: false, error: `Reaction '${id}' not found — call findMetabolicEntities for ids`};
      const d = norm(dir);
      if (d !== 'maximize' && d !== 'minimize')
        return {success: false, error: `Objective for '${id}' must be 'maximize' or 'minimize'`};
      resolved.push({r, coef: d === 'maximize' ? 1 : -1});
    }
    for (const r of model.reactions) {
      if (r.objective_coefficient)
        delete r.objective_coefficient;
    }
    for (const {r, coef} of resolved)
      r.objective_coefficient = coef;
  }
  const objectiveReactions = model.reactions.filter((r) => !!r.objective_coefficient);
  if (!objectiveReactions.length)
    return {success: false, error: 'No objective set — pass objectives, a map of reaction id to maximize or minimize'};
  clearTimeCourseSlider();
  const solution = await WorkerCobraSolver.run_optimization(model);
  if (!solution)
    return {success: false, error: 'FBA failed to run'};
  const fluxes: Record<string, number> = {};
  solution.reactionNames.forEach((name, i) => fluxes[name] = solution.fluxes[i]);
  builder.set_reaction_data(fluxes,
    'FBA Result: ' + objectiveReactions.map((r) =>
      `${r.objective_coefficient! > 0 ? 'Maximize' : 'Minimize'} ${r.id}`).join(', '));
  const objectiveFluxes: Record<string, number> = {};
  let objectiveValue = 0;
  for (const r of objectiveReactions) {
    objectiveFluxes[r.id] = fluxes[r.id] ?? 0;
    objectiveValue += (r.objective_coefficient ?? 0) * (fluxes[r.id] ?? 0);
  }
  return {success: true, objectiveValue, objectiveFluxes, note: 'The map is now colored by the FBA fluxes'};
}

export async function sampleReactionFluxes(view: DG.ViewBase, samples?: number, bins?: number,
  thinning?: number, addDataFrame?: boolean) {
  const {builder} = mgCtx(view);
  const model = requireModel(builder);
  const nSamples = samples ?? DEFAULT_SAMPLES;
  const paramError = samplingParamError(nSamples, bins ?? 20, thinning ?? 10);
  if (paramError)
    return {success: false, error: paramError};
  clearTimeCourseSlider();
  const res = await sampleFluxAveragesWasm(model, bins ?? 20, nSamples, thinning ?? 10);
  if (!res)
    return {success: false, error: 'Invalid reaction bounds in the model — check reaction lower and upper bounds'};
  builder.set_reaction_data(res.reactionData, 'Sampling Histogram');
  applyStepDistribution(builder, res.distribution);
  if (addDataFrame && res.results) {
    const table = samplesToDataFrame(model, res.results, nSamples);
    table.name = 'Sampler results';
    grok.shell.addTableView(table);
  }
  const topAvgs = Object.entries(res.reactionData).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1])).slice(0, 10);
  return {
    success: true,
    samples: nSamples,
    reactions: Object.keys(res.reactionData).length,
    topAbsoluteAverageFluxes: Object.fromEntries(topAvgs),
    note: 'The map is colored by average sampled flux and reaction tooltips show the flux histograms',
  };
}

export async function runTimeCourse(view: DG.ViewBase, startBounds: BoundsMap, endBounds: BoundsMap,
  steps?: number, samples?: number) {
  const {builder} = mgCtx(view);
  const model = requireModel(builder);

  const parse = (o: BoundsMap, label: string): Map<string, ReactionBounds> | {error: string} => {
    const entries = Object.entries(o ?? {});
    if (!entries.length)
      return {error: `${label} is empty — pass a map of reaction id to [lowerBound, upperBound]`};
    const out = new Map<string, ReactionBounds>();
    for (const [id, value] of entries) {
      const b = toBounds(value);
      if (!b)
        return {error: `${label}: invalid bounds for '${id}' — pass [lowerBound, upperBound] numbers`};
      if (b.lower_bound > b.upper_bound)
        return {error: `${label}: lower bound above upper bound for '${id}'`};
      const r = findReaction(model, id);
      if (!r)
        return {error: `${label}: reaction '${id}' not found — call findMetabolicEntities for ids`};
      out.set(r.id, b);
    }
    return out;
  };

  const start = parse(startBounds, 'startBounds');
  if (!(start instanceof Map))
    return {success: false, ...start};
  const end = parse(endBounds, 'endBounds');
  if (!(end instanceof Map))
    return {success: false, ...end};
  const mismatch = [...start.keys()].filter((k) => !end.has(k)).concat([...end.keys()].filter((k) => !start.has(k)));
  if (mismatch.length)
    return {success: false, error: `startBounds and endBounds must list the same reactions. Mismatch: ${mismatch.join(', ')}`};
  const n = steps ?? 10;
  if (n < 2 || n > 200)
    return {success: false, error: 'steps must be between 2 and 200'};
  const nSamples = samples ?? DEFAULT_SAMPLES;
  const paramError = samplingParamError(nSamples, 20, 10);
  if (paramError)
    return {success: false, error: paramError};

  const params: TimeCourseSamplingParams = {
    samples: nSamples, thinning: 10, bins: 20,
    addDf: false, runInPython: false, usePythonFBA: false,
  };
  const res = await runTimeCourseSampling(model, builder, start, end, n, params);
  if (!res.succeeded)
    return {success: false, error: 'Every step had infeasible bounds — nothing was sampled'};
  return {
    success: true,
    steps: n,
    reactions: [...start.keys()],
    ...(res.failedSteps.length ? {infeasibleSteps: res.failedSteps.map((i) => i + 1)} : {}),
    note: 'A time-course slider was added under the map and a table view with average flux per step was opened',
  };
}

export function setReactionFluxes(view: DG.ViewBase, values?: FluxMap, source?: string) {
  const {builder} = mgCtx(view);
  requireMap(builder);
  clearTimeCourseSlider();
  const entries = Object.entries(values ?? {});
  if (!entries.length) {
    builder.set_reaction_data(null, '');
    return {success: true, cleared: true};
  }
  const data: {[id: string]: number} = {};
  for (const [id, v] of entries) {
    const num = typeof v === 'number' ? v : Number(v);
    if (!Number.isFinite(num))
      return {success: false, error: `Value for '${id}' is not a number`};
    data[id] = num;
  }
  builder.set_reaction_data(data, source ?? 'AI-provided data');
  const model = builder.model_data;
  const known = new Set<string>();
  if (model?.reactions?.length) {
    for (const r of model.reactions)
      [r.id, r.bigg_id, r.name].forEach((x) => x && known.add(norm(x)));
  } else {
    for (const k of Object.keys(builder.map.reactions ?? {})) {
      const r = builder.map.reactions[k];
      [r.bigg_id, r.name].forEach((x) => x && known.add(norm(x)));
    }
  }
  const unmatched = known.size ? Object.keys(data).filter((id) => !known.has(norm(id))) : [];
  return {
    success: true,
    applied: Object.keys(data).length,
    ...(unmatched.length ? {unmatchedReactions: unmatched} : {}),
    note: 'The map is colored by the provided values',
  };
}

export function showMetabolicEntity(view: DG.ViewBase, entityId: string) {
  const {builder} = mgCtx(view);
  const map = requireMap(builder);
  const model = builder.model_data;
  const mr = model?.reactions?.length ? findReaction(model, entityId) : null;
  const rk = mapReactionKey(builder, entityId, mr?.id, mr?.bigg_id, mr?.name);
  if (rk != null) {
    map.zoom_to_reaction(rk);
    map.highlight_reaction(rk);
    return {success: true, kind: 'reaction', id: map.reactions[rk].bigg_id};
  }
  const mm = model?.metabolites?.length ? findMetabolite(model, entityId) : null;
  const nk = mapNodeKey(builder, entityId, mm?.id, mm?.bigg_id, mm?.name);
  if (nk != null) {
    map.zoom_to_node(nk);
    map.highlight_node(nk);
    return {success: true, kind: 'metabolite', id: map.nodes[nk].bigg_id};
  }
  if (mr || mm)
    return {success: false, error: `'${entityId}' is in the model but not drawn on the map`};
  return {success: false, error: `No reaction or metabolite '${entityId}' on the map — call findMetabolicEntities`};
}

export function highlightShortestPath(view: DG.ViewBase, fromMetabolite: string, toMetabolite: string, k?: number) {
  const {builder} = mgCtx(view);
  const map = requireMap(builder);
  const model = builder.model_data;
  const resolve = (id: string) => {
    const m = model?.metabolites?.length ? findMetabolite(model, id) : null;
    return mapNodeKey(builder, id, m?.id, m?.bigg_id, m?.name);
  };
  const fromKey = resolve(fromMetabolite);
  if (fromKey == null)
    return {success: false, error: `Metabolite '${fromMetabolite}' not found on the map — call findMetabolicEntities`};
  const toKey = resolve(toMetabolite);
  if (toKey == null)
    return {success: false, error: `Metabolite '${toMetabolite}' not found on the map — call findMetabolicEntities`};
  if (fromKey === toKey)
    return {success: false, error: 'Pick two different metabolites'};
  map.select_none();
  // findAndHighlightShortest orders the endpoints by selectionOrder
  map.nodes[fromKey].selectionOrder = 1;
  map.nodes[toKey].selectionOrder = 2;
  map.select_nodes([fromKey, toKey]);
  map.selectionChanged?.();
  const kth = Math.max(1, Math.round(k ?? 1));
  // a string argument sets k explicitly instead of incrementing it
  map.findAndHighlightShortest(String(kth));
  if (map.last_action?.function !== 'findAndHighlightShortest') {
    return {success: false,
      error: `No path found from '${fromMetabolite}' to '${toMetabolite}' — direction matters, try swapping them`};
  }
  return {success: true, from: map.nodes[fromKey].bigg_id, to: map.nodes[toKey].bigg_id,
    k: kth, note: 'The path is highlighted on the map'};
}

export async function listMetabolicAnalyses(view: DG.ViewBase) {
  mgCtx(view);
  const names = await loadCampaigns();
  return {analyses: names ?? []};
}

const ANALYSIS_NAME_RE = /^[\w\- ]+$/;

export async function saveMetabolicAnalysis(view: DG.ViewBase, name: string,
  description?: string, overwrite?: boolean) {
  const {builder} = mgCtx(view);
  requireMap(builder);
  const clean = (name ?? '').trim();
  if (!clean)
    return {success: false, error: 'Name is required'};
  if (!ANALYSIS_NAME_RE.test(clean))
    return {success: false, error: 'Name may only contain letters, digits, spaces, dashes and underscores'};
  if (!overwrite && await analysisExists(clean))
    return {success: false, error: `Analysis '${clean}' already exists — pass overwrite true to replace it`};
  await saveAnalysisState(builder, clean, description);
  return {success: true, name: clean};
}

export async function loadMetabolicAnalysis(view: DG.ViewBase, name: string) {
  const {builder} = mgCtx(view);
  const clean = (name ?? '').trim();
  if (!ANALYSIS_NAME_RE.test(clean) || !await analysisExists(clean))
    return {success: false, error: `No analysis '${name}' — call listMetabolicAnalyses for names`};
  loadStateProxy(builder, await readAnalysis(clean), clean);
  return {success: true, name: clean, note: 'The saved analysis (map, model and flux data) replaced the current one'};
}
