import {isRef, isSpan} from './model.js';
import type {FilterKind, FilterValue} from './model.js';
import type {FilterProperty} from './schema.js';
import {resolveSpan} from '../span.js';

/** The kind behind `propertyType ?? type`: `double/float/num/qnum → float`,
 * `list → string_list`, `ref → ref`, anything unknown → string; an explicit `kind` wins. */
export function kindOf(prop: FilterProperty): FilterKind {
  if (prop.kind)
    return prop.kind;
  if (prop.ref)
    return 'ref';
  switch (prop.propertyType ?? prop.type) {
    case 'int': return 'int';
    case 'bigint': return 'bigint';
    case 'double': case 'float': case 'num': case 'qnum': return 'float';
    case 'bool': return 'bool';
    case 'datetime': return 'datetime';
    case 'list': case 'string_list': return 'string_list';
    default: return 'string';
  }
}

function scalarOut(v: unknown, now: Date): unknown {
  if (v instanceof Date)
    return v.toISOString();
  if (isSpan(v))
    return resolveSpan(v.span, now).toISOString();
  if (isRef(v))
    return v.id;
  return v;
}

/** A model value as the domain tree carries it: ISO dates, spans resolved against `now`,
 * ref ids, lists element-wise. */
export function domainValue(value: FilterValue | undefined, ctx: {now: Date}): unknown {
  return Array.isArray(value) ? value.map((v) => scalarOut(v, ctx.now)) : scalarOut(value, ctx.now);
}
