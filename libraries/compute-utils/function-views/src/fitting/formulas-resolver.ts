
function evalInScope(formula: string, context: Record<string, any>) {
  return new Function(`with (this) { return (${formula}); }`).call(context);
}

export function compileFormula(formula: string) {
  return new Function(`with (this) { return (${formula}); }`);
}

function getNumberOrNull(res: any) {
  return typeof res === 'number' && !isNaN(res) ? res : null;
}

export function runFormula(formula: string, context: Record<string, any>): number | null {
  try {
    const res = evalInScope(formula, context);
    return getNumberOrNull(res);
  } catch {
    return null;
  }
}
