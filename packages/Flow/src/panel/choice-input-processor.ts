/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChoiceFuncRef, parseChoiceFuncRef, isChoiceQueryRef, literalChoiceItems} from '../utils/choice-refs';

export {parseChoiceFuncRef, isChoiceQueryRef};
export type {ChoiceFuncRef};

/** Resolve by function name; ties broken by the truncated package prefix (`Chem`.endsWith(`hem`)). */
export function resolveChoiceFunc(ref: ChoiceFuncRef): DG.Func | null {
  const candidates = DG.Func.find({name: ref.funcName});
  if (candidates.length === 0) return null;
  if (candidates.length === 1) return candidates[0];
  const prefix = ref.packagePrefix.toLowerCase();
  const byPackage = candidates.filter((f) => {
    try {
      return (f.package?.name ?? '').toLowerCase().endsWith(prefix);
    } catch {
      return false;
    }
  });
  return byPackage[0] ?? candidates[0];
}

/** A function reference as the platform reads it in a header — `Pkg:func(args)` with the function's
 *  real nqName — or null when the reference resolves to nothing or can't be spelled there: the
 *  platform splits a `choices:` value on commas, so a reference with several arguments would break. */
export function headerChoicesRef(raw: string): string | null {
  const ref = parseChoiceFuncRef(raw);
  if (!ref || ref.args.length > 1) return null;
  const func = resolveChoiceFunc(ref);
  if (!func) return null;
  return `${func.nqName}(${ref.args.map((a) => `"${a}"`).join('')})`;
}

/** What a typed choices string is: the canonical text to store, or why it can't be used. A literal
 *  list is stored as `a, b`; a function reference (`Pkg:func()`, also typed without the parentheses)
 *  must name a function that exists, takes the arguments given and returns something to choose from. */
export function checkChoices(raw: string): {value: string} | {error: string} {
  const s = (raw ?? '').trim();
  if (s === '' || isChoiceQueryRef(s)) return {value: s}; // a query reference resolves on its adopted connection
  const bare = /^(\w+):(\w+)$/.exec(s);
  const ref = parseChoiceFuncRef(bare ? `${s}()` : s);
  if (ref) {
    const func = resolveChoiceFunc(ref);
    const inPackage = func != null && ((): boolean => {
      try {
        return func.nqName.split(':')[0].toLowerCase().endsWith(ref.packagePrefix.toLowerCase());
      } catch {
        return true;
      }
    })();
    if (!func || !inPackage) {
      // `Pkg:value` with no parentheses may just be a literal value with a colon in it
      if (bare) return {value: literalChoiceItems(s).join(', ')};
      return {error: `No function ${ref.packagePrefix}:${ref.funcName} is registered`};
    }
    const required = func.inputs.filter((p) => !p.nullable && p.defaultValue == null).length;
    if (ref.args.length < required || ref.args.length > func.inputs.length)
      return {error: `${func.nqName} takes ${func.inputs.length} argument(s), ${ref.args.length} given`};
    if (func.outputs.length === 0)
      return {error: `${func.nqName} returns nothing to choose from`};
    if (ref.args.length > 1)
      return {error: 'A choices function can take one argument at most: the run dialog splits choices on commas'};
    return {value: `${func.nqName}(${ref.args.map((a) => `"${a}"`).join('')})`};
  }
  return {value: literalChoiceItems(s).join(', ')};
}

export function choiceValuesFrom(result: unknown): string[] {
  if (Array.isArray(result))
    return result.map((x) => String(x)).filter((x) => x.length > 0);
  if (result instanceof DG.DataFrame && result.columns.length > 0 && result.rowCount > 0)
    return Array.from(result.columns.byIndex(0).categories).filter((x) => x.length > 0);
  if (result instanceof DG.Column)
    return Array.from(result.categories).filter((x) => x.length > 0);
  return [];
}

/** Resolve a choices REFERENCE (mangled func call or `query("…")` string) into literal
 *  items; null when `raw` isn't one or can't be resolved. A query reference runs on
 *  `connection` — without one it is unresolvable. */
export async function resolveChoicesRef(raw: string, connection: DG.DataConnection | null): Promise<string[] | null> {
  raw = raw.replaceAll('\\"', '"'); // the platform may hand the reference JSON-escaped
  let choiceString = raw.toLowerCase();
  // property tends to strip the ends...
  if (choiceString.startsWith('uery("') && choiceString.endsWith('"'))
    choiceString = `q${choiceString})`;
  if (choiceString.startsWith('query("') && choiceString.endsWith('")')) {
    if (!connection) return null;
    const query =
        '--name: choicesQuery\n--output: dataframe out\n' + choiceString.substring(7, choiceString.length - 2);
    const chq = connection.query('choicesQuery', query);
    return choiceValuesFrom(await chq.apply({}));
  }
  const ref = parseChoiceFuncRef(raw);
  if (!ref) return null;
  const choicesFunc = resolveChoiceFunc(ref);
  if (!choicesFunc) {
    console.warn(`Flow: choices function "${ref.funcName}" not found`);
    return null;
  }
  return choiceValuesFrom(await choicesFunc.apply(ref.args));
}

/** Items for an input NODE's stored choices reference (`properties['choices']`);
 *  query references resolve through the connection captured at adoption
 *  (`properties['choicesConnection']`). Null → not resolvable. */
export async function loadChoicesRefItems(raw: string, connectionId: string): Promise<string[] | null> {
  let connection: DG.DataConnection | null = null;
  if (connectionId) {
    try {
      connection = await grok.dapi.connections.find(connectionId);
    } catch {
      // connection gone — a func reference still resolves
    }
  }
  return resolveChoicesRef(raw, connection);
}

export async function processChoiceInput(input: DG.ChoiceInput<any>, func: DG.Func, inputProperty: DG.Property) {
  ui.setUpdateIndicator(input.input, true, 'loadeing');

  try {
    if (!inputProperty.choices || inputProperty.choices.length !== 1)
      return;

    const prevInputValue = input.value;
    const raw = String(inputProperty.choices[0]).replaceAll('\\"', '"');
    const connection = func instanceof DG.DataQuery ? func.connection ?? null : null;
    const items = await resolveChoicesRef(raw, connection);
    if (items) {
      // Always replace the item list, empty included — the single "choice" the input starts with IS the mangled reference string.
      input.items = items;
      if (prevInputValue && items.includes(prevInputValue) && input.value !== prevInputValue)
        input.value = prevInputValue;
    }
  } catch (e) {
    console.error('error loading options for input ' + input.caption, e);
  } finally {
    ui.setUpdateIndicator(input.input, false);
  }
}
