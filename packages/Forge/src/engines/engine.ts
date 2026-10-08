import * as DG from 'datagrok-api/dg';

export const ENGINE_ROLES = ['train', 'apply', 'isApplicable', 'isInteractive', 'visualize'] as const;
export type EngineRole = typeof ENGINE_ROLES[number];
export type EngineKind = 'function' | 'script';
export type Hyperparameters = {[name: string]: number | string | boolean};
const TRAIN_DATA_INPUTS = ['df', 'predictColumn'] as const;
const NUMBER_TYPES: string[] = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.NUM, DG.TYPE.QNUM];
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

export function hyperparametersOf(engine: Engine): DG.Property[] {
  const inputs = engine.functions.train?.inputs ?? [];
  return inputs.filter((p) => !TRAIN_DATA_INPUTS.some((name) => name === p.name));
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
