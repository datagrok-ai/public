import {isGroup, isRef, isSpan, walk} from './model.js';
import type {FilterCondition, FilterGroup, FilterKind, FilterNode, FilterProblem, FilterScalar} from './model.js';
import {checkLocks, property} from './schema.js';
import type {FilterProperty, FilterSchema, FilterTemplate} from './schema.js';
import {operators} from './operators.js';
import type {FilterOperator} from './operators.js';
import {kindOf} from './kinds.js';
import {isSpanText} from '../span.js';

export type FilterTarget = 'domain' | 'dataframe';

const DIGITS = /^-?\d+$/;

function badScalar(v: FilterScalar, kind: FilterKind, prop: FilterProperty, op: FilterOperator): string | null {
  const bad = (expected: string) => `Expected ${expected} for "${prop.name}"`;
  switch (kind) {
    case 'int':
      if (typeof v !== 'number' || !Number.isInteger(v))
        return bad('an integer');
      break;
    case 'float':
      if (typeof v !== 'number' || !Number.isFinite(v))
        return bad('a number');
      break;
    case 'bigint':
      if (!(typeof v === 'number' && Number.isInteger(v)) && !(typeof v === 'string' && DIGITS.test(v)))
        return bad('an integer');
      break;
    case 'datetime':
      if (v instanceof Date ? Number.isNaN(v.getTime()) : isSpan(v) ? !isSpanText(v.span) :
        typeof v !== 'string' || Number.isNaN(Date.parse(v)))
        return bad('a date or a time span');
      break;
    case 'bool':
      if (typeof v !== 'boolean')
        return bad('true or false');
      break;
    case 'ref':
      if (!isRef(v) && typeof v !== 'string')
        return bad('a reference');
      break;
    default:
      if (typeof v !== 'string' && typeof v !== 'number')
        return bad('a string');
  }
  if (typeof v === 'number') {
    if (prop.min !== undefined && v < prop.min)
      return `"${prop.name}" is at least ${prop.min}`;
    if (prop.max !== undefined && v > prop.max)
      return `"${prop.name}" is at most ${prop.max}`;
  }
  if (prop.choices && (op.id === '=' || op.id === 'in') && !prop.choices.includes(String(v)))
    return `"${v}" is not one of the choices of "${prop.name}"`;
  return null;
}

function checkCondition(c: FilterCondition, schema: FilterSchema, target: FilterTarget | undefined,
  problems: FilterProblem[]): void {
  const problem = (code: FilterProblem['code'], message: string) => problems.push({nodeId: c.id, code, message});
  const prop = property(schema, c.property);
  if (!prop) {
    problem('unknown-property', `Unknown property "${c.property}"`);
    return;
  }
  if (c.property.includes('.'))
    return;
  const op = operators.for(prop).find((o) => o.id === c.operator);
  if (!op) {
    problem('operator-not-applicable', `Operator "${c.operator}" is not applicable to "${c.property}"`);
    return;
  }
  const v = c.value;
  const empty = v === undefined || v === null || v === '';
  if (op.arity === 1 && empty || op.arity === 2 && !(Array.isArray(v) && v.length === 2) ||
    op.arity === 'n' && !(Array.isArray(v) && v.length > 0)) {
    problem('missing-value', op.arity === 2 ? `"${c.property}" needs two bounds` :
      op.arity === 'n' ? `"${c.property}" needs at least one value` : `"${c.property}" needs a value`);
    return;
  }
  if (op.arity !== 0) {
    const kind = kindOf(prop);
    for (const item of Array.isArray(v) ? v : [v as FilterScalar]) {
      const message = badScalar(item, kind, prop, op);
      if (message) {
        problem('invalid-value', message);
        return;
      }
    }
  }
  const expressible = target === 'domain' ? op.domain !== undefined :
    target === 'dataframe' ? op.mask !== undefined || op.bitset !== undefined : true;
  if (!expressible)
    problem('not-expressible', `"${c.property} ${c.operator}" cannot be expressed for ${target}`);
}

export function validate(root: FilterGroup, schema: FilterSchema, target?: FilterTarget,
  template?: FilterTemplate): FilterProblem[] {
  const problems: FilterProblem[] = [];
  walk(root, (node) => {
    if (!isGroup(node))
      checkCondition(node, schema, target, problems);
  });
  if (template)
    problems.push(...checkLocks(template, root));
  return problems;
}

/** The tree without the nodes `validate` reports, and without the groups that leaves empty:
 * what can be evaluated while the rest is still being edited. The same root when nothing goes. */
export function pruneInvalid(root: FilterGroup, schema: FilterSchema, target?: FilterTarget,
  template?: FilterTemplate): FilterGroup {
  const bad = new Set(validate(root, schema, target, template).map((p) => p.nodeId));
  const prune = (group: FilterGroup): FilterGroup => {
    const nodes: FilterNode[] = [];
    for (const n of group.nodes) {
      if (bad.has(n.id))
        continue;
      const kept = isGroup(n) ? prune(n) : n;
      if (!isGroup(kept) || kept.nodes.length > 0)
        nodes.push(kept);
    }
    const same = nodes.length === group.nodes.length && nodes.every((n, i) => n === group.nodes[i]);
    return same ? group : {...group, nodes};
  };
  return prune(root);
}
