// Pure (no platform imports) so a headless test can run it. The server answers 100 rows when a read names no limit,
// so every read here names one and pages until a short page: under a stereo- or tautomer-insensitive rule several
// compounds share a hash, and the rows can outnumber the hashes asked for.
export const LOOKUP_CHUNK = 500;
export const LOOKUP_PAGE = 1000;

export interface LookupQuery {
  filter: {property: string; operator: '='; value: string[]};
  columns: string[];
  sort: string;
  limit: number;
  offset: number;
}

// the first (lowest molregno) row of every hash found among `hashes`
export async function lookupFirstByHash<T extends {molregno: number}>(
  hashes: string[], column: string,
  query: (spec: LookupQuery) => Promise<T[]>,
  chunk: number = LOOKUP_CHUNK, page: number = LOOKUP_PAGE): Promise<Map<string, T>> {
  const found = new Map<string, T>();
  for (let at = 0; at < hashes.length; at += chunk) {
    const filter = {property: column, operator: '=' as const, value: hashes.slice(at, at + chunk)};
    for (let offset = 0; ; offset += page) {
      const rows = await query({filter, columns: [column, 'molregno'], sort: 'molregno', limit: page, offset});
      for (const row of rows) {
        const key = String((row as Record<string, unknown>)[column]);
        if (!found.has(key))
          found.set(key, row);
      }
      if (rows.length < page)
        break;
    }
  }
  return found;
}
