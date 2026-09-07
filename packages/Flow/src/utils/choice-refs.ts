/** Parsing for `choices:` values that are a REFERENCE (function call / SQL query). Core mangles a
 *  single-entry choices string: the FIRST character and the trailing `)` are stripped —
 *  `Chem:getMpoProfileNames()` arrives as `hem:getMpoProfileNames(` — and everything here is
 *  written against that damage. */

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

  // The closing ')' may or may not have survived the truncation.
  let argsPart = s.slice(open + 1).trim();
  if (argsPart.endsWith(')')) argsPart = argsPart.slice(0, -1).trim();
  const args = argsPart.length === 0 ? [] : argsPart.split(',')
    .map((a) => a.trim().replace(/^['"]|['"]$/g, ''))
    .filter((a) => a.length > 0);

  return {packagePrefix, funcName, args};
}

/** A SQL-query reference (`query("…")`, first char stripped). The quote may arrive
 *  JSON-escaped from the platform (`uery(\"SELECT …\"`) — both forms count. */
export function isChoiceQueryRef(raw: string): boolean {
  return /^q?uery\(\\?"/i.test((raw ?? '').trim());
}

/** A stored `choices` string that is a single REFERENCE (func call / query) rather
 *  than a literal comma list — resolved live by the value editors, kept out of
 *  emitted headers (quotes/parens/commas would corrupt the options block). */
export function isChoicesRefString(s: string): boolean {
  const t = (s ?? '').trim();
  return t !== '' && (isChoiceQueryRef(t) || parseChoiceFuncRef(t) !== null);
}

/** Literal selectable values vs a single reference to resolve first — seeding
 *  `hem:getMpoProfileNames(` as a value would be worse than seeding nothing. */
export function isLiteralChoiceList(choices: readonly unknown[]): boolean {
  if (choices.length === 0) return false;
  if (choices.length > 1) return true;
  const only = String(choices[0]);
  return !isChoiceQueryRef(only) && parseChoiceFuncRef(only) === null;
}
