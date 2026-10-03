/** Reading the portal's tables. A `query({filter})` without a limit gets the server's default of
 * 100 rows, silently: the sync history, a busy program's evidence or the follows of a portfolio
 * outgrow it. Every list read passes ROW_LIMIT, the client's cap. */
export const ROW_LIMIT = 10000;

/** [v] as a quoted value in a filter string: a `"` or `\` in it cannot end the value early. */
export const filterValue = (v: string): string => `"${String(v ?? '').replace(/[\\"]/g, '\\$&')}"`;
