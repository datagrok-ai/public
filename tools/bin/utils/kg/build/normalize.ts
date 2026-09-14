/// One row normalizer for `check` and `build` (build-plan.md WO-1): nulls dropped, defaults applied,
/// dates canonical, lists deduplicated and scalars coerced into lists, then every present member
/// shape-checked. Required members are not checked here: a node may still be assembled from
/// several extractors and stubs before `finalize()` decides.
import {TypeSystem, EdgeType, Member, ValueHooks, checkValue} from '../types';

export type Row = Record<string, unknown>;

export interface RowProblem {
  key: string;
  code: string;
  message: string;
}

export interface Normalized {
  row: Row;
  problems: RowProblem[];
  /** Members whose value is the declared default, not an authored one. */
  defaulted: string[];
}

/** Keys of an edge row that are not properties. */
export const EDGE_ROW_KEYS = ['type', 'from', 'to', 'name'];

export function normalizeRow(system: TypeSystem, row: Row, hooks: Partial<ValueHooks> = {}): Normalized {
  const typeName = String(row.type ?? '');
  const type = system.nodes.get(typeName);
  if (!type || type.abstract) {
    const message = type ? `type '${typeName}' is abstract` : `unknown type '${typeName}'`;
    return {row: {...row}, problems: [{key: 'type', code: 'bad-type', message}], defaulted: []};
  }
  const members: Record<string, Member> = {...type.members};
  for (const m of system.buildFields.node) members[m.name] = m;
  return normalizeMembers(system, row, members, ['type'], typeName, hooks);
}

export function normalizeEdgeRow(system: TypeSystem, edge: EdgeType, row: Row, hooks: Partial<ValueHooks> = {}): Normalized {
  const members: Record<string, Member> = {...edge.properties};
  for (const m of system.buildFields.edge) members[m.name] = m;
  return normalizeMembers(system, row, members, EDGE_ROW_KEYS, edge.name, hooks);
}

function normalizeMembers(system: TypeSystem, row: Row, members: Record<string, Member>, passthrough: string[], typeName: string,
  hooks: Partial<ValueHooks>): Normalized {
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
    const problem = checkValue(member, value, {provenance: system.provenance, path: hooks.path, ref: hooks.ref});
    if (problem) {
      const code = member.kind === 'ref' ? 'unresolved-ref' : member.scalar === 'Path' ? 'missing-path' : 'bad-value';
      problems.push({key, code, message: `${key}: ${problem}`});
      continue;
    }
    out[key] = value;
  }
  const defaulted: string[] = [];
  for (const m of Object.values(members))
    if (m.default !== undefined && out[m.name] === undefined) {
      out[m.name] = m.default;
      defaulted.push(m.name);
    }
  return {row: out, problems, defaulted};
}

function normalizeValue(member: Member, raw: unknown): unknown {
  if (!member.list) return member.scalar === 'Date' ? canonicalDate(raw) : raw;
  const items = (Array.isArray(raw) ? raw : [raw]).filter((v) => v !== null && v !== undefined)
    .map((v) => member.scalar === 'Date' ? canonicalDate(v) : v);
  const seen = new Set<string>();
  return items.filter((v) => {
    const key = typeof v === 'object' ? JSON.stringify(v) : `${typeof v}:${String(v)}`;
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  });
}

/** A YAML date becomes its ISO text; a date without a time of day stays date-only. */
function canonicalDate(value: unknown): unknown {
  if (!(value instanceof Date) || Number.isNaN(value.getTime())) return value;
  const iso = value.toISOString();
  return iso.endsWith('T00:00:00.000Z') ? iso.slice(0, 10) : iso;
}
