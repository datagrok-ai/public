/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChoiceFuncRef, parseChoiceFuncRef, isChoiceQueryRef} from '../utils/choice-refs';

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
