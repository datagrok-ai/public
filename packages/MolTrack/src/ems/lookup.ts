import type * as DG from 'datagrok-api/dg';

export const MATCH_RULES = ['ALL_LAYERS', 'STEREO_INSENSITIVE_LAYERS', 'TAUTOMER_INSENSITIVE_LAYERS'] as const;

export type MatchRule = typeof MATCH_RULES[number];

export type HashColumn = 'hash_all_layers' | 'hash_stereo_insensitive' | 'hash_tautomer_insensitive';

export const HASH_COLUMNS: Record<MatchRule, HashColumn> = {
  ALL_LAYERS: 'hash_all_layers',
  STEREO_INSENSITIVE_LAYERS: 'hash_stereo_insensitive',
  TAUTOMER_INSENSITIVE_LAYERS: 'hash_tautomer_insensitive',
};

// the server answers 100 rows to a read without a limit
const CHUNK = 500;
const PAGE = 1000;

type Keyed = {molregno: number} & Partial<Record<HashColumn, unknown>>;

export async function lookupFirstByHash<T extends Keyed>(hashes: string[], rule: MatchRule,
  query: (spec: DG.DomainQuerySpec<HashColumn | 'molregno', never>) => Promise<T[]>): Promise<Map<string, T>> {
  const column = HASH_COLUMNS[rule];
  const found = new Map<string, T>();
  for (let at = 0; at < hashes.length; at += CHUNK) {
    const chunk = {property: column, operator: '=' as const, value: hashes.slice(at, at + CHUNK)};
    for (let after = 0; ;) {
      const rows = await query({filter: [chunk, {property: 'molregno', operator: '>', value: after}],
        columns: [column, 'molregno'], sort: 'molregno', limit: PAGE});
      for (const row of rows) {
        const key = String(row[column]);
        if (!found.has(key))
          found.set(key, row);
      }
      if (rows.length < PAGE)
        break;
      after = rows[rows.length - 1].molregno;
    }
  }
  return found;
}
