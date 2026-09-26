/// One row normalizer for `check` and `build` (build-plan.md WO-1): nulls dropped, defaults applied,
/// dates canonical, lists deduplicated and scalars coerced into lists, then every present member
/// shape-checked. Required members are not checked here: a node may still be assembled from
/// several extractors and stubs before `finalize()` decides.
import {TypeSystem, EdgeType, Member, ValueHooks, checkValue} from './types';

export type Row = Record<string, unknown>;

export interface RowProblem {
  key: string;
  code: string;
  message: string;
  /** The path or id the problem is about, for the stale report. */
  target?: string;
}

export interface Normalized {
  row: Row;
  problems: RowProblem[];
  /** Members whose value is the declared default, not an authored one. */
  defaulted: string[];
}

/** Keys of an edge row that are not properties. */
export const EDGE_ROW_KEYS = ['type', 'from', 'to', 'name'];

/** List members whose order and repetitions are the fact, not a set: a signature has four `list<string>` inputs and
 * says so four times, a route's parameters come in the order the template spells them. Every other list is a set. */
const ORDERED_MEMBERS = ['input_types', 'output_types', 'path_params', 'query_params', 'actions'];

export function isOrdered(member: string): boolean {
  return ORDERED_MEMBERS.includes(member);
}

/** Code-point order, the same on every platform and locale. */
export function compare(a: string, b: string): number {
  return a < b ? -1 : a > b ? 1 : 0;
}

export interface NormalizeOptions extends Partial<ValueHooks> {
  /** Apply the declared defaults (real rows); a stub carries only what created it. */
  defaults?: boolean;
}

export function normalizeRow(system: TypeSystem, row: Row, options: NormalizeOptions = {}): Normalized {
  const typeName = String(row.type ?? '');
  const type = system.nodes.get(typeName);
  if (!type || type.abstract) {
    const message = type ? `type '${typeName}' is abstract` : `unknown type '${typeName}'`;
    return {row: {...row}, problems: [{key: 'type', code: 'bad-type', message}], defaulted: []};
  }
  const members: Record<string, Member> = {...type.members};
  for (const m of system.buildFields.node) members[m.name] = m;
  return normalizeMembers(system, row, members, ['type'], typeName, options);
}

export function normalizeEdgeRow(system: TypeSystem, edge: EdgeType, row: Row, options: NormalizeOptions = {}): Normalized {
  const members: Record<string, Member> = {...edge.properties};
  for (const m of system.buildFields.edge) members[m.name] = m;
  return normalizeMembers(system, row, members, EDGE_ROW_KEYS, edge.name, options);
}

function normalizeMembers(system: TypeSystem, row: Row, members: Record<string, Member>, passthrough: string[], typeName: string,
  options: NormalizeOptions): Normalized {
  const out: Row = {};
  const problems: RowProblem[] = [];
  for (const key of passthrough)
    if (row[key] !== undefined && row[key] !== null) out[key] = row[key];
  for (const [key, raw] of Object.entries(row)) {
    if (passthrough.includes(key) || raw === null || raw === undefined) continue;
    const member = members[key];
    if (!member) {
      problems.push({key, code: 'unknown-key', message: `unknown key '${key}' for type ${typeName}`});
      continue;
    }
    const value = normalizeValue(member, raw);
    const problem = checkValue(member, value, {provenance: system.provenance, path: options.path, ref: options.ref});
    if (problem) {
      const code = member.kind === 'ref' ? 'unresolved-ref' : member.scalar === 'Path' ? 'missing-path' : 'bad-value';
      problems.push({key, code, message: `${key}: ${problem}`, ...(member.scalar === 'Path' ? {target: String(raw)} : {})});
      continue;
    }
    out[key] = value;
  }
  const defaulted: string[] = [];
  if (options.defaults !== false)
    for (const m of Object.values(members))
      if (m.default !== undefined && out[m.name] === undefined) {
        out[m.name] = m.default;
        defaulted.push(m.name);
      }
  return {row: out, problems, defaulted};
}

function normalizeValue(member: Member, raw: unknown): unknown {
  if (!member.list) return scalar(member, raw);
  const items = (Array.isArray(raw) ? raw : [raw]).filter((v) => v !== null && v !== undefined).map((v) => scalar(member, v));
  if (isOrdered(member.name)) return items;
  const seen = new Set<string>();
  return items.filter((v) => {
    const key = typeof v === 'object' ? JSON.stringify(v) : `${typeof v}:${String(v)}`;
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  });
}

/** A YAML date becomes its ISO text, date-only when it has no time of day; a folded YAML scalar loses its trailing newline. */
function scalar(member: Member, value: unknown): unknown {
  if (value instanceof Date && !Number.isNaN(value.getTime())) {
    const iso = value.toISOString();
    return iso.endsWith('T00:00:00.000Z') ? iso.slice(0, 10) : iso;
  }
  return typeof value === 'string' && (member.scalar === 'string' || member.scalar === 'Text') ? value.trim() : value;
}
