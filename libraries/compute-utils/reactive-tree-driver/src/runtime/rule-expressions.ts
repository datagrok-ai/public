import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as jsonLogic from 'json-logic-js';
import {isClientAtLeast} from '@datagrok-libraries/utils/src/format-version-utils';
import {IControllerBase} from '../RuntimeControllers';
import {RuleExpr, RuleTargets} from '../config/PipelineConfiguration';

export type RuleContext = Record<string, any> & {
  all: Record<string, any[]>,
};

let opsRegistered = false;

// the context of the expression being evaluated, read by the script ops
let activeCtx: RuleContext = {all: {}};
let scriptSupport: boolean | undefined;

// GrokScript with a variables map is js-api 1.28+; the same method ignores the map on
// older clients, so the version decides, not the method's presence
function hasScriptSupport(): boolean {
  if (scriptSupport === undefined) {
    scriptSupport = typeof (grok.functions as any).scriptSync === 'function' && isClientAtLeast('1.28.0');
    if (!scriptSupport)
      console.warn('RTD: GrokScript expressions need platform 1.28.0 or later, expression checks are skipped');
  }
  return scriptSupport;
}

function scriptVariables(ctx: RuleContext): Record<string, any> {
  return Object.fromEntries(Object.entries(ctx).filter(([key]) => key !== 'all' && key !== 'literals'));
}

/** Evaluates a GrokScript expression over the rule context; `undefined` when unsupported or failing. */
export function runScript(expr: string, ctx: RuleContext): any {
  if (!hasScriptSupport())
    return undefined;
  try {
    return grok.functions.scriptSync(expr, scriptVariables(ctx));
  } catch {
    return undefined;
  }
}

/** The platform's `validator:` expression contract: `null` when valid, else the message. */
export function scriptVerdict(expr: string, ctx: RuleContext): string | null {
  if (!hasScriptSupport())
    return null;
  let result: unknown;
  try {
    result = grok.functions.scriptSync(expr, scriptVariables(ctx));
  } catch {
    return `Error during validation: "${expr}"`;
  }
  if (result === false)
    return expr;
  return typeof result === 'string' ? result : null;
}

/** Platform column kinds accepted by `type:` annotations, plus a column type or a semantic type. */
export function columnIs(col: any, kind: string): boolean {
  if (!(col instanceof DG.Column))
    return false;
  switch (kind) {
  case 'numerical': return col.isNumerical || col.type === DG.TYPE.DATE_TIME;
  case 'numerical_no_datetime': return col.isNumerical && col.type !== DG.TYPE.DATE_TIME;
  case 'categorical': return col.isCategorical;
  case 'datetime': return col.type === DG.TYPE.DATE_TIME;
  case 'categorical_or_datetime': return col.isCategorical || col.type === DG.TYPE.DATE_TIME;
  default: return col.type === kind || col.semType === kind;
  }
}

function columnsOf(df: any, kind?: string): DG.Column[] {
  if (!(df instanceof DG.DataFrame))
    return [];
  const cols = df.columns.toList();
  return kind == null ? cols : cols.filter((col) => columnIs(col, kind));
}

function registerOps() {
  if (opsRegistered)
    return;
  opsRegistered = true;
  jsonLogic.add_operation('columns', (df: any, kind?: string) => columnsOf(df, kind).map((col) => col.name));
  jsonLogic.add_operation('columnsMissing', (df: any, spec: any[]) => {
    const missing: string[] = [];
    for (const entry of spec ?? []) {
      const [name, kind] = Array.isArray(entry) ? entry : [entry];
      if (!columnsOf(df, kind).some((col) => col.name === name))
        missing.push(kind == null ? name : `${name} (${kind})`);
    }
    return missing;
  });
  jsonLogic.add_operation('columnIs', columnIs);
  jsonLogic.add_operation('nulls', (col: any) => col instanceof DG.Column ? col.stats.missingValueCount : 0);
  jsonLogic.add_operation('column', (df: any, name: string) =>
    df instanceof DG.DataFrame ? df.col(name)?.toList() ?? [] : []);
  jsonLogic.add_operation('row', (df: any, keyColumn: string, key: any) => {
    const col = df instanceof DG.DataFrame ? df.col(keyColumn) : null;
    if (!col || key == null)
      return null;
    for (let i = 0; i < col.length; i++) {
      const value = col.get(i);
      if (value === key || String(value) === String(key))
        return Object.fromEntries(df.columns.names().map((name: string) => [name, df.get(name, i)]));
    }
    return null;
  });
  jsonLogic.add_operation('regex', (val: any, pattern: string, flags?: string) =>
    typeof val === 'string' && new RegExp(pattern, flags ?? '').test(val));
  jsonLogic.add_operation('script', (expr: string) => runScript(expr, activeCtx));
  jsonLogic.add_operation('scriptVerdict', (expr: string) => scriptVerdict(expr, activeCtx));
  jsonLogic.add_operation('len', (val: any) => {
    if (val == null)
      return 0;
    if (val instanceof DG.DataFrame)
      return val.rowCount;
    return val.length ?? 0;
  });
}

export function ruleTargets(targets: RuleTargets): string[] {
  return Array.isArray(targets) ? targets : [targets];
}

const LITERAL = 'literal';

function isLiteral(node: any): boolean {
  return node != null && typeof node === 'object' && !Array.isArray(node) &&
    Object.keys(node).length === 1 && LITERAL in node;
}

// `{literal: x}` shields x from JSON Logic, which would otherwise read any
// one-key object as an operation; the value is moved into the data context
function extractLiterals(expr: any, literals: any[]): any {
  if (Array.isArray(expr))
    return expr.map((item) => extractLiterals(item, literals));
  if (expr == null || typeof expr !== 'object')
    return expr;
  if (isLiteral(expr)) {
    literals.push(expr[LITERAL]);
    return {var: `literals.${literals.length - 1}`};
  }
  return Object.fromEntries(Object.entries(expr).map(([key, value]) => [key, extractLiterals(value, literals)]));
}

export function evaluate(expr: RuleExpr, ctx: RuleContext): any {
  registerOps();
  const literals: any[] = [];
  const logic = extractLiterals(expr, literals);
  const previous = activeCtx;
  activeCtx = ctx;
  try {
    return jsonLogic.apply(logic, literals.length ? {...ctx, literals} : ctx);
  } finally {
    activeCtx = previous;
  }
}

export function isOn(when: RuleExpr | undefined, ctx: RuleContext): boolean {
  return when == null || jsonLogic.truthy(evaluate(when, ctx));
}

export function buildRuleContext(controller: IControllerBase): RuleContext {
  const ctx: RuleContext = {all: {}};
  for (const name of controller.getMatchedInputs()) {
    const all = controller.getAll(name) ?? [];
    ctx[name] = all[0];
    ctx.all[name] = all;
  }
  return ctx;
}

const scopedOps = new Set(['map', 'filter', 'reduce', 'all', 'some', 'none']);

/** Root input aliases referenced by `var`, `missing` and `missing_some`, with the `all.` prefix stripped. */
export function usedAliases(expr: RuleExpr | undefined): string[] {
  const aliases = new Set<string>();
  const addPath = (path: any) => {
    if (typeof path !== 'string' || !path)
      return;
    const segments = path.split('.');
    const root = segments[0] === 'all' ? segments[1] : segments[0];
    if (root)
      aliases.add(root);
  };
  const visit = (node: any) => {
    if (Array.isArray(node)) {
      node.forEach(visit);
      return;
    }
    if (node == null || typeof node !== 'object' || isLiteral(node))
      return;
    const keys = Object.keys(node);
    if (keys.length !== 1) {
      Object.values(node).forEach(visit);
      return;
    }
    const op = keys[0];
    const values: any[] = Array.isArray(node[op]) ? node[op] : [node[op]];
    if (op === 'var')
      addPath(values[0]);
    else if (op === 'missing')
      values.forEach(addPath);
    else if (op === 'missing_some')
      (Array.isArray(values[1]) ? values[1] : [values[1]]).forEach(addPath);
    // the second argument of an array operation runs over the element, not the rule context
    if (scopedOps.has(op)) {
      visit(values[0]);
      if (op === 'reduce')
        visit(values[2]);
      return;
    }
    values.forEach(visit);
  };
  visit(expr);
  return [...aliases];
}
