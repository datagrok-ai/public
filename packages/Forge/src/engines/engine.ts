import * as DG from 'datagrok-api/dg';
import type {ModelTask} from '../generated/db';
import type {MissingValuesSettings} from '../preparation/missing-values';
import {isRecord} from '../preparation/preparation-options';

export const ENGINE_ROLES = ['train', 'apply', 'isApplicable', 'isInteractive', 'visualize'] as const;
export type EngineRole = typeof ENGINE_ROLES[number];
export type EngineKind = 'function' | 'script';
export type Hyperparameters = {[name: string]: number | string | boolean};
const TRAIN_DATA_INPUTS = ['df', 'predictColumn'] as const;
export const NUMBER_TYPES: string[] = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.NUM, DG.TYPE.QNUM];
/** Methods that run on the server although their functions do not say so (`meta.mlserver`) yet. */
export const SERVER_ENGINES = ['Chemprop'];

export interface Engine {
  name: string;
  namespace: string;
  kind: EngineKind;
  functions: Partial<Record<EngineRole, DG.Func>>;
  isLiveUpdate: boolean;
}

export function isComplete(engine: Engine): boolean {
  const {train, apply, isApplicable} = engine.functions;
  return train !== undefined && apply !== undefined && isApplicable !== undefined;
}

/** The training data leaves the browser: a script in a language other than JavaScript, a `train` with
 * `meta.mlserver: true`, or one of {@link SERVER_ENGINES}. */
export function isServerEngine(engine: Engine): boolean {
  const train = engine.functions.train;
  const isServerScript = train instanceof DG.Script && train.language !== 'javascript';
  return isServerScript || train?.options['mlserver'] === 'true' || SERVER_ENGINES.includes(engine.name);
}

/** The model file embeds training rows (`meta.mlhasrows: true` on `train`), such as SVM's support vectors. */
export function hasTrainingRows(engine: Engine): boolean {
  return engine.functions.train?.options['mlhasrows'] === 'true';
}

/** The method takes rows with missing feature values as they are (`meta.mlnativenans: true` on `train`). */
export function handlesMissingValues(engine: Engine): boolean {
  return engine.functions.train?.options['mlnativenans'] === 'true';
}

/** The method returns class probabilities (`meta.mlprobabilities` on `train`: 'true', or the name of the bool train
 * input that turns them on). */
export function givesProbabilities(engine: Engine): boolean {
  return engine.functions.train?.options['mlprobabilities'] === 'true' || probabilityInput(engine) !== undefined;
}

/** The bool `train` input that turns the class probabilities on, as `meta.mlprobabilities` names it. */
export function probabilityInput(engine: Engine): string | undefined {
  const text: unknown = engine.functions.train?.options['mlprobabilities'];
  if (typeof text !== 'string' || text === 'true')
    return undefined;
  const input = hyperparametersOf(engine).find((p) => p.name === text);
  return input?.propertyType === DG.TYPE.BOOL ? text : undefined;
}

/** [hyperparameters] with the method's probability input set to [on]; unchanged for a method without one. */
export function withProbabilityInput(engine: Engine, hyperparameters: Hyperparameters, on: boolean): Hyperparameters {
  const name = probabilityInput(engine);
  return name === undefined ? hyperparameters : {...hyperparameters, [name]: on};
}

/** The method accepts concurrent train and apply calls (`meta.mlparallel: true` on `train`): a search runs its fits
 * in parallel. */
export function isParallel(engine: Engine): boolean {
  return engine.functions.train?.options['mlparallel'] === 'true';
}

/** The method's search space, as its `train` declares it (`meta.mlspace`, JSON text); undefined without one. */
export function searchSpaceText(engine: Engine): string | undefined {
  const text: unknown = engine.functions.train?.options['mlspace'];
  return typeof text === 'string' && text !== '' ? text : undefined;
}

/** A `meta.mlenabled` condition: every clause holds. A clause on another setting holds when its value is one of a list
 * or above `gt` / below `lt`; `task` holds when the method trains one of the tasks listed. */
export type EnabledCondition = {[name: string]: (string | number | boolean)[] | {gt?: number; lt?: number}};
/** When each setting of a method is used: any of its conditions holds; a setting without one is always used. */
export type EnabledRules = {[setting: string]: EnabledCondition[]};

/** The rules of the method's `meta.mlenabled` on `train` (JSON text `{"<setting>": <condition> | [<condition>, ...]}`):
 * empty without one, undefined when it is malformed. */
export function enabledRules(engine: Engine): EnabledRules | undefined {
  const text: unknown = engine.functions.train?.options['mlenabled'];
  if (typeof text !== 'string' || text === '')
    return {};
  let parsed: unknown;
  try {
    parsed = JSON.parse(text);
  } catch {
    return undefined;
  }
  if (!isRecord(parsed))
    return undefined;
  const rules: EnabledRules = {};
  for (const [setting, value] of Object.entries(parsed)) {
    const conditions = Array.isArray(value) ? value : [value];
    if (conditions.length === 0 || !conditions.every(isCondition))
      return undefined;
    rules[setting] = conditions;
  }
  return rules;
}

/** Whether the method uses [setting] with the settings [values] when it trains a [task] (unknown: any task); a clause
 * on a [free] setting, whose value a search chooses, may hold. */
export function isEnabled(rules: EnabledRules, setting: string, values: Hyperparameters,
  task: ModelTask | undefined, free?: ReadonlySet<string>): boolean {
  const holds = ([name, clause]: [string, EnabledCondition[string]]) => {
    if (name !== 'task' && free?.has(name))
      return true;
    if (Array.isArray(clause)) {
      const value = name === 'task' ? task : values[name];
      return value === undefined ? name === 'task' : clause.some((v) => String(v) === String(value));
    }
    const value = Number(values[name]);
    return (clause.gt === undefined || value > clause.gt) && (clause.lt === undefined || value < clause.lt);
  };
  const conditions = rules[setting];
  return conditions === undefined || conditions.some((condition) => Object.entries(condition).every(holds));
}

export function isChoice(value: unknown): value is string | number | boolean {
  return typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean';
}

function isCondition(value: unknown): value is EnabledCondition {
  return isRecord(value) && Object.values(value).every((clause) => Array.isArray(clause) ? clause.every(isChoice) :
    isRecord(clause) && Object.keys(clause).length > 0 &&
      Object.entries(clause).every(([op, n]) => (op === 'gt' || op === 'lt') && typeof n === 'number'));
}

/** The missing-values handling when none is chosen: `keep` for a method that handles missing values, else `skip`. */
export function defaultMissingValues(engine: Engine): MissingValuesSettings {
  return {mode: handlesMissingValues(engine) ? 'keep' : 'skip'};
}

export function hyperparametersOf(engine: Engine): DG.Property[] {
  const inputs = engine.functions.train?.inputs ?? [];
  return inputs.filter((p) => !TRAIN_DATA_INPUTS.some((name) => name === p.name));
}

/** The settings the form shows: {@link hyperparametersOf} without the probability input, which Predict probability
 * sets. */
export function visibleHyperparametersOf(engine: Engine): DG.Property[] {
  const hidden = probabilityInput(engine);
  return hyperparametersOf(engine).filter((p) => p.name !== hidden);
}

/** [values] without the method's probability input, as {@link visibleHyperparametersOf}. */
export function visibleValuesOf(engine: Engine, values: Hyperparameters): Hyperparameters {
  const hidden = probabilityInput(engine);
  return Object.fromEntries(Object.entries(values).filter(([name]) => name !== hidden));
}

export function defaultHyperparameters(engine: Engine): Hyperparameters {
  return defaultValuesOf(hyperparametersOf(engine));
}

export function defaultValuesOf(props: DG.Property[]): Hyperparameters {
  const values: Hyperparameters = {};
  for (const p of props) {
    const value = initialValueOf(p);
    if (value !== undefined)
      values[p.name] = value;
  }
  return values;
}

// A header default (`= 20`) is kept as the input's initial value, in text; `defaultValue` stays empty.
// The server keeps a string default in single quotes (`'RBF'`), an annotation in double quotes.
function initialValueOf(p: DG.Property): number | string | boolean | undefined {
  const text: unknown = p.initialValue;
  if (typeof text !== 'string' || text === '')
    return undefined;
  if (p.propertyType === DG.TYPE.BOOL)
    return text.toLowerCase() === 'true';
  if (p.propertyType === DG.TYPE.STRING)
    return text.replace(/^(["'])(.*)\1$/, '$2');
  if (!NUMBER_TYPES.includes(p.propertyType))
    return undefined;
  const value = Number(text);
  return Number.isNaN(value) ? undefined : value;
}

export function rolesOf(engine: Engine): EngineRole[] {
  return ENGINE_ROLES.filter((role) => engine.functions[role] !== undefined);
}
