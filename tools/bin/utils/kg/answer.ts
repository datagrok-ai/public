/// What an operation answers (conventions.md §11.1, build-plan.md WO-7): sections of rows over a resolved target,
/// led by the caveats of the sources the manifest could not read. Data only, shaped once: the CLI prints it
/// (`print.ts`), the browser and an MCP server serialize it as it is.

export interface Section {
  title: string;
  rows: Record<string, unknown>[];
  /** Rows before paging, when the section was paged: it is counted whole first, so a header can say `50 of 394`. */
  total?: number;
  /** The whole header line of a section that came back empty, saying why it did. */
  empty?: string;
}

export interface Answer {
  op: string;
  target: Record<string, unknown> | null;
  /** One line per source the manifest does not report as `ok`; the Dart clause leads. */
  notes?: string[];
  sections: Section[];
}

/** One row of `explain`'s `edges` section: an edge type in one direction, what vouches for it and where it leads. */
export type EdgeGroup = {
  /** The `edges/` folder of the type, `reference` for a reference property. */
  group: string;
  edge: string;
  direction: 'in' | 'out';
  count: number;
  derived_by: string[];
  /** The lowest and the highest confidence among the edges; none when no edge carries one. */
  confidence: [number, number] | null;
  /** The evidence of the first edge that has any. */
  evidence: string[];
  /** The first neighbours, at most `HOP_TARGETS` of `count`. */
  targets: {id: string, name?: string, type?: string}[];
};

/** Every section is counted whole and paged after, so the header can tell a page from the total. */
export function section(title: string, rows: Record<string, unknown>[], limit: number, empty?: string): Section {
  const paged: Section = {title, rows: rows.slice(0, limit), total: rows.length};
  if (!rows.length && empty) paged.empty = empty;
  return paged;
}

/** What a manifest says about its sources: their status, and what each contributes (`provides`). */
export type Coverage = {sources?: Record<string, string>, provides?: Record<string, string>};

export interface Caveat {
  /** The manifest source, `dart` for the Dart clause. */
  source: string;
  status: string;
  note: string;
}

/** What the manifest says about the Dart pass; every op and every report leads with it unless the pass is `ok`. */
export function coverageNote(sources: Record<string, string> | undefined): string | undefined {
  const status = sources?.dart;
  if (status === undefined || status === 'missing') return 'Dart coverage unknown (the dart extractor did not run)';
  return /^ok\b/.test(status) ? undefined : 'Dart coverage partial (some markers did not resolve)';
}

/** One caveat per source the manifest does not report as `ok`, Dart first: what an answer here cannot be based on. */
export function sourceCaveats(coverage: Coverage | undefined): Caveat[] {
  const sources = coverage?.sources;
  const caveats: Caveat[] = [];
  const dart = coverageNote(sources);
  if (dart) caveats.push({source: 'dart', status: sources?.dart ?? 'missing', note: dart});
  for (const [name, status] of Object.entries(sources ?? {}).sort(([a], [b]) => a < b ? -1 : 1)) {
    if (name === 'dart' || /^ok\b/.test(status)) continue;
    const means = coverage?.provides?.[name] ?? `what ${name} contributes`;
    const reach = status === 'missing' ? 'no' : status === 'stale' ? 'possibly out-of-date' : 'incomplete';
    caveats.push({source: name, status, note: `${name} ${status}: ${reach} coverage of ${means}`});
  }
  return caveats;
}

/** The caveats as the lines an answer carries. */
export function caveatNotes(coverage: Coverage | undefined): string[] {
  return sourceCaveats(coverage).map((c) => c.note);
}
