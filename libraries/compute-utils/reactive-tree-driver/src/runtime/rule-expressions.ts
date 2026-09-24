import * as DG from 'datagrok-api/dg';
import * as jsonLogic from 'json-logic-js';
import {IControllerBase} from '../RuntimeControllers';
import {RuleExpr, RuleTargets} from '../config/PipelineConfiguration';

export type RuleContext = Record<string, any> & {
  all: Record<string, any[]>,
};

let opsRegistered = false;

function registerOps() {
  if (opsRegistered)
    return;
  opsRegistered = true;
  jsonLogic.add_operation('columns', (df: any, kind?: string) => {
    if (!(df instanceof DG.DataFrame))
      return [];
    const cols = df.columns.toList();
    const filtered = kind == null ? cols :
      kind === 'numerical' ? cols.filter((col) => col.isNumerical) :
        kind === 'categorical' ? cols.filter((col) => col.isCategorical) :
          cols.filter((col) => col.type === kind || col.semType === kind);
    return filtered.map((col) => col.name);
  });
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
  return jsonLogic.apply(logic, literals.length ? {...ctx, literals} : ctx);
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

/** Root input aliases referenced by `var`, `missing` and `missing_some`, with the `all.` prefix stripped. */
export function usedAliases(expr: RuleExpr): string[] {
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
    values.forEach(visit);
  };
  visit(expr);
  return [...aliases];
}
