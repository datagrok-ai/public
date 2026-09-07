/* The expression and validator rules behind FuncCallForm's W4 surfaces (`visible:` / `enabled:` /
   `validator:` options, named-validator verdicts, categoryGroups): evaluation over the
   feature-detected raw `grok_ScriptSync` global — GrokScript itself stays in Dart — plus the
   regex-literal parser, verdict mapping, message plaining and the categoryGroups plan expander.
   Internal to `src/dg/funcs` — func-form.ts is the only importer; not exported through
   `src/dg/index.ts`. */

const api = globalThis as {
  grok_ScriptSync?: (script: string, variables?: Record<string, unknown>) => unknown,
};

function hasScriptSync(): boolean {
  return typeof api.grok_ScriptSync === 'function';
}

/** One `console.warn` per key until {@link clear} — a bad expression re-fails on every keystroke
 * and must not spam the console (warns, never errors: the e2e gate counts errors). Each form owns
 * an instance: two live forms over one call warn — and re-arm — independently. */
export class WarnOnce {
  private readonly warned = new Set<string>();

  warn(key: string, message: string): void {
    if (this.warned.has(key))
      return;
    this.warned.add(key);
    console.warn(`u2 func-form: ${key} — ${message}`);
  }

  clear(): void {
    this.warned.clear();
  }
}

/** The `visible:`/`enabled:` mapping (ib:339-341): anything but `false` is true; a script
 * failure — or a missing seam — keeps the previous state (the ib:328-331 counterpart). */
export function evalRule(expr: unknown, vars: Record<string, unknown>, previous: boolean,
  key: string, warnings: WarnOnce): boolean {
  if (!hasScriptSync()) {
    warnings.warn(key, 'grok_ScriptSync is not available — keeping the previous state');
    return previous;
  }
  try {
    return api.grok_ScriptSync!(String(expr), vars) !== false;
  } catch (e) {
    warnings.warn(key, 'expression failed: ' + (e instanceof Error ? e.message : String(e)));
    return previous;
  }
}

/** ib:230-237 verbatim: false → the expression text as the message; a string → that
 * message; a throw → the Dart onError wording; anything else → valid. */
export function evalValidatorExpression(expr: string, vars: Record<string, unknown>,
  key: string, warnings: WarnOnce): string | null {
  if (!hasScriptSync()) {
    warnings.warn(key, 'grok_ScriptSync is not available — skipping validation');
    return null;
  }
  let result: unknown;
  try {
    result = api.grok_ScriptSync!(expr, vars);
  } catch {
    return `Error during validation: "${expr}"`;
  }
  if (result === false)
    return expr;
  return typeof result === 'string' ? result : null;
}

/** `Validators.tryParseRegexLiteral` (vld:135-142) + `regex()` (vld:125-131): `/pattern/flags`,
 * flags i/m honored, non-strings fail with the same message a non-match gets. */
export function regexValidator(s: string): ((v: unknown) => string | null) | null {
  if (s.length < 2 || !s.startsWith('/'))
    return null;
  const end = s.lastIndexOf('/');
  if (end <= 0)
    return null;
  const pattern = s.slice(1, end);
  const flags = s.slice(end + 1);
  const re = new RegExp(pattern, flags.replace(/[^im]/g, ''));
  const label = `/${pattern}/${flags}`;
  return (v) => typeof v === 'string' && re.test(v) ? null : `Value doesn't match ${label}`;
}

/** Strips the `#{…}` command markup and `**` emphasis the Dart builtin messages carry
 * (vld:43-51) — u2 renders plain validity text, not ValidatorHelper markup. */
export function plainMessage(s: string): string {
  return s.replace(/#\{[^}]*\}/g, '').replace(/\*\*/g, '').replace(/ {2,}/g, ' ').trim();
}

export type PlanItem = {isHeader: boolean, name: string, level: number};

/** fpe:379-398 verbatim: `{header: value}` recurses, lists/strings are leaf categories. */
export function expandCategoryGroups(node: unknown, depth: number): PlanItem[] {
  const result: PlanItem[] = [];
  if (node !== null && typeof node === 'object' && !Array.isArray(node)) {
    for (const key of Object.keys(node as Record<string, unknown>)) {
      result.push({isHeader: true, name: key, level: depth});
      result.push(...expandCategoryGroups((node as Record<string, unknown>)[key], depth + 1));
    }
  }
  else if (Array.isArray(node)) {
    for (const cat of node)
      if (typeof cat === 'string')
        result.push({isHeader: false, name: cat, level: depth});
  }
  else if (typeof node === 'string')
    result.push({isHeader: false, name: node, level: depth});
  return result;
}
