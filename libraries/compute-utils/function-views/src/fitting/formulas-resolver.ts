// Compile a Datagrok bound-formula expression into a fast (ctx) → number
// closure. The compiled body destructures only the referenced names from
// the supplied context object — V8 turns those into monomorphic property
// loads, dramatically faster than the previous per-call `new Function` +
// `with(this)` form which was both reparsed every call and deoptimized.
//
// No internal cache: callers that need "compile once" hold the returned
// closure (bounds-checker precompiles each formula bound at fit setup).

export type CompiledFormula = (context: Record<string, any>) => number | null;

// JS-style identifier-shaped tokens. Doesn't handle Unicode identifiers;
// Datagrok variable names are ASCII in practice. Swap to a unicode-aware
// regex (\p{ID_Start}/\p{ID_Continue}) if that ever changes.
const ID_RX = /\b[A-Za-z_$][A-Za-z0-9_$]*\b/g;

// Compile a formula expression once for a given known set of context names.
// `contextNames` is the universe of names the caller knows the formula MAY
// reference. Names not in this list resolve via the function's lexical /
// global scope (e.g. `Math`, `Number`). Extra names not actually referenced
// by the formula are filtered out by intersecting with regex-extracted
// tokens — the regex over-approximates but the intersect is the safety net.
export function compileFormula(
  formula: string,
  contextNames: readonly string[],
): CompiledFormula {
  const ids = new Set(formula.match(ID_RX) ?? []);
  const used = contextNames.filter((n) => ids.has(n));
  let fn: ((ctx: Record<string, any>) => any) | null = null;
  try {
    const destruct = used.length ? `const {${used.join(', ')}} = ctx;\n` : '';
    fn = new Function('ctx', `${destruct}return (${formula})`) as any;
  } catch {
    return () => null;
  }
  return (ctx) => {
    try {
      const v = fn!(ctx);
      return typeof v === 'number' && !isNaN(v) ? v : null;
    } catch {
      return null;
    }
  };
}

// Legacy entry point — UI / one-shot path. Derives `contextNames` from the
// provided context object's keys and recompiles per call. Hot callers should
// precompile via `compileFormula(formula, knownNames)` once and reuse the
// returned closure across calls.
export function runFormula(formula: string, context: Record<string, any>): number | null {
  return compileFormula(formula, Object.keys(context))(context);
}
