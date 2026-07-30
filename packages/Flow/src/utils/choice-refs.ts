/** Parsing for `choices:` values that are a REFERENCE (a function call or a
 *  SQL query) rather than a literal list.
 *
 *  Core mangles a single-entry `choices` string on its way to
 *  `Property.choices`: the **first character and the trailing `)` are
 *  stripped**. So `Chem:getMpoProfileNames()` arrives as
 *  `hem:getMpoProfileNames(` and `query("select …")` as `uery("select …"`.
 *  Everything here is written against that damage, and is pure — no DG, no
 *  backend — so both the property panel (which resolves the reference) and
 *  `FuncNode` (which must NOT seed a reference as if it were a value) can use
 *  it. */

/** A parsed `Pkg:funcName(arg, …)` reference. */
export interface ChoiceFuncRef {
  /** The (first-character-truncated) package prefix, e.g. `hem` — a hint only. */
  packagePrefix: string;
  funcName: string;
  /** Literal arguments, unquoted. */
  args: string[];
}

/** Parse a mangled function-call choices string, or null when it isn't one. */
export function parseChoiceFuncRef(raw: string): ChoiceFuncRef | null {
  const s = (raw ?? '').trim();
  const colon = s.indexOf(':');
  const open = s.indexOf('(');
  if (colon <= 0 || open <= colon + 1) return null;

  const packagePrefix = s.slice(0, colon);
  const funcName = s.slice(colon + 1, open).trim();
  if (!/^\w+$/.test(funcName)) return null;

  // Everything after '(' is the argument list; the closing ')' may or may not
  // have survived the truncation.
  let argsPart = s.slice(open + 1).trim();
  if (argsPart.endsWith(')')) argsPart = argsPart.slice(0, -1).trim();
  const args = argsPart.length === 0 ? [] : argsPart.split(',')
    .map((a) => a.trim().replace(/^['"]|['"]$/g, ''))
    .filter((a) => a.length > 0);

  return {packagePrefix, funcName, args};
}

/** Whether a choices value is a SQL-query reference (`query("…")`, first char
 *  stripped). */
export function isChoiceQueryRef(raw: string): boolean {
  return /^q?uery\("/i.test((raw ?? '').trim());
}

/** Whether a `choices` list is a literal set of selectable values, as opposed to
 *  a single reference that has to be resolved first. Used to decide whether the
 *  first choice is a legitimate default to seed a node with — seeding
 *  `hem:getMpoProfileNames(` as a value would be worse than seeding nothing. */
export function isLiteralChoiceList(choices: readonly unknown[]): boolean {
  if (choices.length === 0) return false;
  if (choices.length > 1) return true;
  const only = String(choices[0]);
  return !isChoiceQueryRef(only) && parseChoiceFuncRef(only) === null;
}
