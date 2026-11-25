import * as DG from 'datagrok-api/dg';

export type InputValue = number | string | boolean;

const NUMBER_TYPES: string[] = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.NUM, DG.TYPE.QNUM];

const cache = new Map<string, {[name: string]: InputValue}>();

/** A header default is kept as the input's initial value, in text; a string one in single quotes (`'RBF'`) */
function initialValueOf(p: DG.Property): InputValue | undefined {
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

/** The header defaults of `f`'s inputs that have one */
export function defaultsOf(f: DG.Func): {[name: string]: InputValue} {
  const values: {[name: string]: InputValue} = {};
  for (const p of f.inputs) {
    const value = initialValueOf(p);
    if (value !== undefined)
      values[p.name] = value;
  }
  return values;
}

/** `values` with every null or undefined input of the Eda function `name` set to its header default:
 * `grok.functions.call` passes an input left out as null */
export function withDefaults<T extends {[name: string]: unknown}>(name: string, values: T): T {
  let defaults = cache.get(name);
  if (defaults === undefined) {
    defaults = defaultsOf(DG.Func.find({package: 'Eda', name})[0]);
    cache.set(name, defaults);
  }
  const filled: {[name: string]: unknown} = {...values};
  for (const [input, value] of Object.entries(defaults))
    filled[input] ??= value;
  return filled as T;
}
