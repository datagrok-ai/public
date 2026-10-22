/* Reads bigger than one request: the server answers a JSON query with at most MAX_ROWS rows, and a long `=` list
   goes a chunk at a time. */
export class DomainQueries {
  /** The row cap of a JSON query. */
  static readonly MAX_ROWS = 10000;
  /** Values per `=` list, rows per page. */
  static readonly CHUNK = 500;

  /** Every row of a paged `query`, which must sort by a unique column for the pages to line up. */
  static async all<T>(query: (page: {limit: number, offset: number}) => Promise<T[]>,
    size = DomainQueries.CHUNK): Promise<T[]> {
    const rows: T[] = [];
    for (let offset = 0; ; offset += size) {
      const page = await query({limit: size, offset});
      rows.push(...page);
      if (page.length < size)
        return rows;
    }
  }

  /** `run` over `values` a chunk at a time, the results in order. */
  static async chunked<V, T>(values: readonly V[], run: (chunk: V[]) => Promise<T[]>,
    size = DomainQueries.CHUNK): Promise<T[]> {
    const results: T[] = [];
    for (let at = 0; at < values.length; at += size)
      results.push(...await run(values.slice(at, at + size)));
    return results;
  }
}
