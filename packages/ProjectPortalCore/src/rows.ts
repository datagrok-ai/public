/** Reading the portal's tables. A `query({filter})` without a limit gets the server's default of
 * 100 rows, silently: the sync history, a busy program's evidence or the follows of a portfolio
 * outgrow it. Every list read passes ROW_LIMIT, the client's cap. */
export const ROW_LIMIT = 10000;

/** [v] as a quoted value in a filter string: a `"` or `\` in it cannot end the value early. */
export const filterValue = (v: string): string => `"${String(v ?? '').replace(/[\\"]/g, '\\$&')}"`;

/** `column in (...)` over the distinct non-empty [values], or null when there are none: the server
 * cannot parse an empty list (`Cannot parse filter: "id in ()"`). */
export function inFilter(column: string, values: Iterable<string | null | undefined>): string | null {
  const vs = [...new Set([...values].filter((v): v is string => v != null && v !== ''))];
  return vs.length ? `${column} in (${vs.map(filterValue).join(', ')})` : null;
}

/** Values per `in (...)`: keeps a filter short for a list of hundreds. */
const IN_BATCH = 100;

/** The rows of [table] whose [column] is one of [values], and [and] when given; none, without asking
 * the server, when there are no values. A long list is read in batches. */
export async function queryIn<T>(table: {query(q: {filter: string, limit: number}): Promise<T[]>}, column: string,
  values: Iterable<string | null | undefined>, and?: string): Promise<T[]> {
  const vs = [...new Set([...values].filter((v): v is string => v != null && v !== ''))];
  const out: T[] = [];
  for (let i = 0; i < vs.length; i += IN_BATCH) {
    const f = inFilter(column, vs.slice(i, i + IN_BATCH))!;
    out.push(...await table.query({filter: and ? `${and} and ${f}` : f, limit: ROW_LIMIT}));
  }
  return out;
}
