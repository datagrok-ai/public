// The server answers 100 rows when a read names no limit, so every read names one and pages.
const CHUNK = 500;
const PAGE = 1000;

export interface LookupQuery {
  filter: {property: string; operator: '=' | '>'; value: string[] | number}[];
  columns: string[];
  sort: string;
  limit: number;
}

export async function lookupFirstByHash<T extends {molregno: number}>(
  hashes: string[], column: string, query: (spec: LookupQuery) => Promise<T[]>): Promise<Map<string, T>> {
  const found = new Map<string, T>();
  for (let at = 0; at < hashes.length; at += CHUNK) {
    const chunk = {property: column, operator: '=' as const, value: hashes.slice(at, at + CHUNK)};
    // paged by key: an offset shifts when a row is deleted between two pages
    for (let after = 0; ;) {
      const rows = await query({filter: [chunk, {property: 'molregno', operator: '>', value: after}],
        columns: [column, 'molregno'], sort: 'molregno', limit: PAGE});
      for (const row of rows) {
        const key = String((row as Record<string, unknown>)[column]);
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
