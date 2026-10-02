/** The model catalogue as configuration declares it: unit operations, run types, and models with
 * their cards (F-12, FR-50). A portfolio file carries its catalogue under `catalogue`; the loader
 * merges the catalogues of every portfolio it loads. Pure data and pure checks. */

import {ModelApplicability, validateApplicability} from './recommend';

export const APPROVAL_STATES = ['draft', 'in_review', 'approved', 'rejected', 'retired'] as const;
export type ApprovalState = typeof APPROVAL_STATES[number];

export const MODEL_STATUSES = ['active', 'deprecated', 'experimental'] as const;
export type ModelStatus = typeof MODEL_STATUSES[number];

export interface UnitOperationConfig {
  code: string;
  name: string;
  ordinal?: number;
  description?: string;
}

export interface RunTypeConfig {
  code: string;
  name: string;
  definition?: string;
  ordinal?: number;
  /** False for values on another axis (input provenance, a run's role). */
  isExecutionMode?: boolean;
}

/** One input or output on a card. For an input, `min` and `max` are the validity range: the
 * values the model was built for, which may be narrower than what the function accepts. */
export interface CardParam {
  name: string;
  caption?: string;
  units?: string;
  min?: number;
  max?: number;
  description?: string;
}

export interface ModelCardConfig {
  /** The model version the card describes. */
  version: string;
  /** The scientific question the model answers (PORT-41). */
  question?: string;
  intendedUse?: string;
  /** Read from the function's signature when absent. */
  inputs?: CardParam[];
  outputs?: CardParam[];
  /** The applicability domain in words. */
  applicability?: string;
  performance?: string;
  lineage?: string;
  /** The version this one replaces. */
  supersedes?: string;
  limitations?: string;
  /** The person or team who answers for the model. */
  owner?: string;
  approval?: ApprovalState;
  /** ISO date. */
  approvedOn?: string;
  /** Where the model applies, by key (task, stage, modality, work-unit kind) and the dataset
   * columns it predicts; what model recommendation reads. Absent: any, within its portfolios. */
  appliesTo?: ModelApplicability;
}

export interface ModelConfig {
  /** The function's nqName (`Package:Name`): the stable handle, since function ids change every release. */
  code: string;
  name: string;
  package?: string;
  category?: string;
  description?: string;
  status?: ModelStatus;
  /** Unit operations the model applies to; the first is the default when a run is saved. */
  unitOperations?: string[];
  /** Run type a run is saved with by default. */
  defaultRunType?: string;
  card?: ModelCardConfig;
}

export interface ModelCatalogueConfig {
  unitOperations?: UnitOperationConfig[];
  runTypes?: RunTypeConfig[];
  models?: ModelConfig[];
}

/** A catalogue merged from several sources: each model knows the portfolios that declare it. */
export interface MergedCatalogue {
  unitOperations: (UnitOperationConfig & {portfolios?: string[]})[];
  runTypes: (RunTypeConfig & {portfolios?: string[]})[];
  models: (ModelConfig & {portfolios: string[]})[];
}

const CODE = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;
const NQ_NAME = /^[A-Za-z0-9_]+:[A-Za-z0-9_]+$/;

/** Every problem with one catalogue, as sentences. Codes it names but does not declare must be
 * declared somewhere in the set ([known]); the set check passes them in. */
export function validateCatalogue(c: ModelCatalogueConfig,
  known: {unitOperations?: Set<string>, runTypes?: Set<string>} = {}): string[] {
  const errors: string[] = [];
  const codes = (items: {code: string}[] | undefined, what: string): Set<string> => {
    const seen = new Set<string>();
    for (const x of items ?? []) {
      if (!x.code || !CODE.test(x.code))
        errors.push(`${what} "${x.code ?? ''}": a code may hold only letters, digits, '_', '.', '-'`);
      else if (seen.has(x.code))
        errors.push(`${what} ${x.code} is declared twice`);
      if (!(x as any).name)
        errors.push(`${what} ${x.code}: name is required`);
      seen.add(x.code);
    }
    return seen;
  };
  const ops = codes(c.unitOperations, 'unit operation');
  const types = codes(c.runTypes, 'run type');
  const opKnown = (code: string) => ops.has(code) || !!known.unitOperations?.has(code);
  const typeKnown = (code: string) => types.has(code) || !!known.runTypes?.has(code);
  const models = new Set<string>();
  for (const m of c.models ?? []) {
    const where = `model ${m.code ?? '?'}`;
    if (!m.code || !NQ_NAME.test(m.code))
      errors.push(`${where}: code must be the function's nqName, Package:Name`);
    else if (models.has(m.code))
      errors.push(`${where} is declared twice`);
    models.add(m.code);
    if (!m.name)
      errors.push(`${where}: name is required`);
    if (m.status && !(MODEL_STATUSES as readonly string[]).includes(m.status))
      errors.push(`${where}: status "${m.status}" is not one of ${MODEL_STATUSES.join(', ')}`);
    for (const op of m.unitOperations ?? []) {
      if (!opKnown(op))
        errors.push(`${where}: unit operation ${op} is not declared`);
    }
    if (m.defaultRunType && !typeKnown(m.defaultRunType))
      errors.push(`${where}: run type ${m.defaultRunType} is not declared`);
    if (m.card)
      errors.push(...validateCard(m.card).map((e) => `${where}: ${e}`));
  }
  return errors;
}

export function validateCard(card: ModelCardConfig): string[] {
  const errors: string[] = [];
  if (!card.version)
    errors.push('card: version is required');
  if (card.approval && !(APPROVAL_STATES as readonly string[]).includes(card.approval))
    errors.push(`card: approval "${card.approval}" is not one of ${APPROVAL_STATES.join(', ')}`);
  if (card.approvedOn && Number.isNaN(Date.parse(card.approvedOn)))
    errors.push(`card: approvedOn "${card.approvedOn}" is not a date`);
  if (card.supersedes && card.supersedes === card.version)
    errors.push('card: a version cannot supersede itself');
  errors.push(...validateApplicability(card.appliesTo).map((e) => `card: ${e}`));
  for (const [kind, params] of [['input', card.inputs], ['output', card.outputs]] as const) {
    const names = new Set<string>();
    for (const p of params ?? []) {
      if (!p.name)
        errors.push(`card: an ${kind} has no name`);
      else if (names.has(p.name))
        errors.push(`card: ${kind} ${p.name} is listed twice`);
      names.add(p.name);
      if (p.min !== undefined && p.max !== undefined && p.min > p.max)
        errors.push(`card: ${kind} ${p.name} has min ${p.min} above max ${p.max}`);
    }
  }
  return errors;
}

/** Merges the catalogues of several portfolios. A unit operation, run type or model declared by
 * two portfolios must be declared the same way; a model then lists both portfolios. */
export function mergeCatalogues(sources: {portfolio: string | null, catalogue: ModelCatalogueConfig}[]):
  {catalogue: MergedCatalogue, errors: string[]} {
  const errors: string[] = [];
  const known = {
    unitOperations: new Set(sources.flatMap((s) => (s.catalogue.unitOperations ?? []).map((x) => x.code))),
    runTypes: new Set(sources.flatMap((s) => (s.catalogue.runTypes ?? []).map((x) => x.code))),
  };
  const merge = <T extends {code: string}>(what: string, pick: (c: ModelCatalogueConfig) => T[] | undefined) => {
    const out = new Map<string, {item: T, from: string}>();
    for (const s of sources) {
      for (const item of pick(s.catalogue) ?? []) {
        const seen = out.get(item.code);
        if (seen && JSON.stringify(seen.item) !== JSON.stringify(item))
          errors.push(`${what} ${item.code} is declared differently by ${seen.from} and ${s.portfolio ?? 'a catalogue'}`);
        else if (!seen)
          out.set(item.code, {item, from: s.portfolio ?? 'a catalogue'});
      }
    }
    return out;
  };
  for (const s of sources) {
    errors.push(...validateCatalogue(s.catalogue, known).map((e) => s.portfolio ? `${s.portfolio}: ${e}` : e));
  }
  const declaring = (pick: (c: ModelCatalogueConfig) => {code: string}[] | undefined, code: string) =>
    sources.filter((s) => s.portfolio && (pick(s.catalogue) ?? []).some((x) => x.code === code)).map((s) => s.portfolio!);
  const unitOperations = [...merge('unit operation', (c) => c.unitOperations).values()]
    .map(({item}) => ({...item, portfolios: declaring((c) => c.unitOperations, item.code)}));
  const runTypes = [...merge('run type', (c) => c.runTypes).values()]
    .map(({item}) => ({...item, portfolios: declaring((c) => c.runTypes, item.code)}));
  const models = [...merge('model', (c) => c.models).values()].map(({item}) => ({...item,
    portfolios: sources.filter((s) => s.portfolio && (s.catalogue.models ?? []).some((m) => m.code === item.code))
      .map((s) => s.portfolio!)}));
  return {catalogue: {unitOperations, runTypes, models}, errors};
}

/** The union of two '|'-separated code lists, the first list's order first. Loading one portfolio's
 * catalogue adds that portfolio to a model and never drops another's. */
export function unionCodes(existing: string | null | undefined, added: string[]): string {
  const out = (existing ?? '').split('|').map((s) => s.trim()).filter((s) => s.length > 0);
  for (const a of added)
    if (!out.includes(a))
      out.push(a);
  return out.join('|');
}

export const splitCodes = (s: string | null | undefined): string[] =>
  (s ?? '').split('|').map((x) => x.trim()).filter((x) => x.length > 0);

/** A function parameter as the platform describes it, reduced to what a card shows. */
export interface FuncParamInfo {
  name: string;
  caption?: string | null;
  units?: string | null;
  min?: number | null;
  max?: number | null;
  description?: string | null;
}

/** The card's inputs or outputs: those the card lists, in its order, else the function's own. A
 * card entry fills in what it leaves out from the function's parameter of the same name. */
export function cardParams(declared: CardParam[] | null | undefined, func: FuncParamInfo[]): CardParam[] {
  const clean = (p: FuncParamInfo): CardParam => {
    const out: CardParam = {name: p.name};
    if (p.caption) out.caption = p.caption;
    if (p.units) out.units = p.units;
    if (p.min !== null && p.min !== undefined) out.min = p.min;
    if (p.max !== null && p.max !== undefined) out.max = p.max;
    if (p.description) out.description = p.description;
    return out;
  };
  if (!declared?.length)
    return func.map(clean);
  return declared.map((d) => ({...clean(func.find((f) => f.name === d.name) ?? {name: d.name}), ...d}));
}

/** A json column holds an object, not an array: card parameters are stored as `{params: [...]}`. */
export const storedParams = (params: CardParam[] | null): {params: CardParam[]} | null => params ? {params} : null;
export const readParams = (stored: unknown): CardParam[] | null =>
  Array.isArray((stored as any)?.params) ? (stored as any).params : Array.isArray(stored) ? stored as CardParam[] : null;

/** Card inputs and outputs the function does not have: a card written for another version. */
export function cardMismatches(card: {inputs?: CardParam[] | null, outputs?: CardParam[] | null},
  func: {inputs: string[], outputs: string[]}): string[] {
  const out: string[] = [];
  for (const p of card.inputs ?? [])
    if (!func.inputs.includes(p.name))
      out.push(`input ${p.name} is not an input of the function`);
  for (const p of card.outputs ?? [])
    if (!func.outputs.includes(p.name))
      out.push(`output ${p.name} is not an output of the function`);
  return out;
}

/** Input values outside the card's validity range, as sentences: a run there is an extrapolation. */
export function outsideValidity(inputs: CardParam[], values: Record<string, unknown>): string[] {
  const out: string[] = [];
  for (const p of inputs) {
    const v = values[p.name];
    if (typeof v !== 'number' || Number.isNaN(v))
      continue;
    const label = p.caption ?? p.name;
    if (p.min !== undefined && v < p.min)
      out.push(`${label} ${v} is below the validity range (${p.min}–${p.max ?? '…'}${p.units ? ` ${p.units}` : ''})`);
    else if (p.max !== undefined && v > p.max)
      out.push(`${label} ${v} is above the validity range (${p.min ?? '…'}–${p.max}${p.units ? ` ${p.units}` : ''})`);
  }
  return out;
}
