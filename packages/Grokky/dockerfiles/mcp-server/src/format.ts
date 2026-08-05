// Tool-result shaping. A tool result is not a one-off cost: it stays in the transcript and is
// re-read on every later hop of the turn, so an unpaged pretty-printed list taxes the whole turn
// (LATENCY.md Tier 4). Results are therefore compact JSON, array results are paged, and everything
// is capped — a naive `list` on a large instance used to inject hundreds of KB mid-turn.

export type ToolResult = {content: {type: 'text'; text: string}[]};

const DEFAULT_LIMIT = 50;
const MAX_CHARS = 20_000;
const MAX_RAW_CHARS = 40_000;

const wrap = (text: string): ToolResult => ({content: [{type: 'text', text}]});

function cap(text: string, max: number): string {
  return text.length <= max ? text :
    `${text.slice(0, max)}\n…truncated at ${max} chars — narrow the filter or use limit/offset`;
}

export interface FormatOpts {
  paged?: boolean;
  raw?: boolean;
  args?: Record<string, any>;
}

export function formatResult(data: unknown, opts: FormatOpts = {}): ToolResult {
  if (typeof data === 'string')
    return wrap(cap(data, opts.raw ? MAX_RAW_CHARS : MAX_CHARS));

  let payload = data;
  if (opts.paged && Array.isArray(data)) {
    const offset = Number(opts.args?.offset) || 0;
    const limit = Number(opts.args?.limit) || DEFAULT_LIMIT;
    const items = data.slice(offset, offset + limit);
    payload = items.length < data.length ?
      {total: data.length, showing: `${offset + 1}-${offset + items.length}`,
        note: 'pass offset/limit for more', items} :
      items;
  }
  return wrap(cap(JSON.stringify(payload), opts.raw ? MAX_RAW_CHARS : MAX_CHARS));
}

/** Errors come back as a result, not a throw: the model can read and correct them in-place,
 * whereas a protocol-level error ends the tool call with nothing actionable in it. */
export function formatError(message: string, detail?: unknown): ToolResult {
  return wrap(JSON.stringify(detail === undefined ? {error: message} : {error: message, ...(detail as object)}));
}
