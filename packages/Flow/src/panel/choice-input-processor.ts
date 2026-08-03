/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChoiceFuncRef, parseChoiceFuncRef, isChoiceQueryRef} from '../utils/choice-refs';

export {parseChoiceFuncRef, isChoiceQueryRef};
export type {ChoiceFuncRef};

/** Resolve a parsed reference to a single function. Matches on the function
 *  name; when several packages expose that name, the truncated package prefix
 *  breaks the tie (a real package name ends with it — `Chem`.endsWith(`hem`)). */
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

/** Coerce whatever a choices function returned into a string list: a plain
 *  list, or the first column of a dataframe. */
export function choiceValuesFrom(result: unknown): string[] {
  if (Array.isArray(result))
    return result.map((x) => String(x)).filter((x) => x.length > 0);
  if (result instanceof DG.DataFrame && result.columns.length > 0 && result.rowCount > 0)
    return Array.from(result.columns.byIndex(0).categories).filter((x) => x.length > 0);
  if (result instanceof DG.Column)
    return Array.from(result.categories).filter((x) => x.length > 0);
  return [];
}

export async function processChoiceInput(input: DG.ChoiceInput<any>, func: DG.Func, inputProperty: DG.Property) {
  ui.setUpdateIndicator(input.input, true, 'loadeing');

  try {
    if (!inputProperty.choices || inputProperty.choices.length !== 1)
      return;

    const prevInputValue = input.value;
    const raw = String(inputProperty.choices[0]).replaceAll('\\"', '"');
    // Always replaces the item list, empty included: the single "choice" the
    // input starts with IS the mangled reference string, so bailing out on an
    // empty result would leave `hem:getMpoProfileNames(` selectable.
    const applyItems = (items: string[]): void => {
      input.items = items;
      if (prevInputValue && items.includes(prevInputValue) && input.value !== prevInputValue)
        input.value = prevInputValue;
    };

    let choiceString = raw.toLowerCase();
    if (func instanceof DG.DataQuery && func.connection) {
      // property tends to strip the ends...
      if (choiceString.startsWith('uery("') && choiceString.endsWith('"'))
        choiceString = `q${choiceString})`;
      if (choiceString.startsWith('query("') && choiceString.endsWith('")')) {
        const query =
            '--name: choicesQuery\n--output: dataframe out\n' + choiceString.substring(7, choiceString.length - 2);
        const chq = func.connection.query('choicesQuery', query);
        const res = await chq.apply({});
        applyItems(choiceValuesFrom(res));
        return;
      }
    }

    // `choices: Pkg:funcName(…)` — resolve the function and use what it returns.
    const ref = parseChoiceFuncRef(raw);
    if (ref) {
      const choicesFunc = resolveChoiceFunc(ref);
      if (!choicesFunc) {
        console.warn(`Flow: choices function "${ref.funcName}" not found for input ${input.caption}`);
        return;
      }
      applyItems(choiceValuesFrom(await choicesFunc.apply(ref.args)));
    }
  } catch (e) {
    console.error('error loading options for input ' + input.caption, e);
  } finally {
    ui.setUpdateIndicator(input.input, false);
  }
}
