import {TYPE} from 'datagrok-api/u2core';
import {KIND, isRef, isSpan} from './model.js';
import type {FilterKind, FilterValue} from './model.js';
import type {FilterProperty} from './schema.js';
import {resolveSpan} from '../span.js';

/** The kind behind `propertyType ?? type`: `double/float/num/qnum → float`,
 * `list → string_list`, `ref → ref`, anything unknown → string; an explicit `kind` wins. */
export function kindOf(prop: FilterProperty): FilterKind {
  if (prop.kind)
    return prop.kind;
  if (prop.ref)
    return KIND.REF;
  switch (prop.propertyType ?? prop.type) {
    case TYPE.INT: return KIND.INT;
    case TYPE.BIG_INT: return KIND.BIG_INT;
    case TYPE.FLOAT: case KIND.FLOAT: case TYPE.NUM: case TYPE.QNUM: return KIND.FLOAT;
    case TYPE.BOOL: return KIND.BOOL;
    case TYPE.DATE_TIME: return KIND.DATE_TIME;
    case TYPE.LIST: case TYPE.STRING_LIST: return KIND.STRING_LIST;
    default: return KIND.STRING;
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
