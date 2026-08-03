import * as DG from 'datagrok-api/dg';

/**
 * Single source of truth for the organisms the package supports for g:Profiler
 * enrichment and UniProt subcellular-location narrowing. `code` is the g:Profiler
 * organism id; `display` is the user-facing dropdown label. Moved here (a leaf
 * module) from `analysis/enrichment.ts` so parsers, the annotate dialog, and
 * enrichment can all share it without a circular import.
 */
export const ORGANISM_LIST = [
  {display: 'Homo sapiens (Human)', code: 'hsapiens'},
  {display: 'Mus musculus (Mouse)', code: 'mmusculus'},
  {display: 'Rattus norvegicus (Rat)', code: 'rnorvegicus'},
  {display: 'Saccharomyces cerevisiae (Yeast)', code: 'scerevisiae'},
  {display: 'Escherichia coli (K12)', code: 'ecoli'},
  {display: 'Danio rerio (Zebrafish)', code: 'drerio'},
  {display: 'Drosophila melanogaster (Fruit fly)', code: 'dmelanogaster'},
  {display: 'Arabidopsis thaliana', code: 'athaliana'},
  {display: 'Caenorhabditis elegans', code: 'celegans'},
] as const;

/** Scientific-name prefix (lowercased, parenthetical stripped) → g:Profiler code,
 * derived once from ORGANISM_LIST. e.g. 'homo sapiens' → 'hsapiens'. */
const SCIENTIFIC_TO_CODE: ReadonlyArray<{scientific: string; code: string}> =
  ORGANISM_LIST.map((o) => ({
    scientific: o.display.split('(')[0].trim().toLowerCase(),
    code: o.code,
  }));

/** Display label → code and code → display, for dialog seeding. */
export function organismDisplayForCode(code: string | undefined): string | undefined {
  return ORGANISM_LIST.find((o) => o.code === code)?.display;
}
export function organismCodeForDisplay(display: string | undefined): string | undefined {
  return ORGANISM_LIST.find((o) => o.display === display)?.code;
}

/**
 * Maps a raw UniProt/Spectronaut organism-name string to a g:Profiler code, or
 * undefined if it matches none of the supported organisms. Matches on the
 * scientific-name prefix so strain-qualified values resolve too — e.g.
 * "Escherichia coli (strain K12)" → 'ecoli', "Rattus norvegicus" → 'rnorvegicus'.
 */
export function resolveOrganismCode(name: string | null | undefined): string | undefined {
  if (!name) return undefined;
  const lower = name.toLowerCase();
  for (const {scientific, code} of SCIENTIFIC_TO_CODE)
    if (lower.includes(scientific)) return code;
  return undefined;
}

/** Column-name hints for the per-protein organism column that some vendors emit
 * (Spectronaut's `PG.Organisms`; generic "Organism"/"Species"). */
const ORGANISM_COLUMN_HINTS = ['pg.organisms', 'organism', 'organisms', 'species'];

function findOrganismColumn(df: DG.DataFrame): DG.Column | null {
  for (const col of df.columns.toList()) {
    const n = col.name.toLowerCase();
    if (ORGANISM_COLUMN_HINTS.some((h) => n === h || n.includes(h))) return col;
  }
  return null;
}

/**
 * Infers the experiment's organism from an organism column in the data, when the
 * data is unambiguously a single supported species. Returns the g:Profiler code,
 * or undefined when there is no organism column, nothing resolves, or MORE THAN
 * ONE supported organism is present (e.g. a HYE multi-species mix) — in which
 * case the caller should leave the choice to the user rather than guess.
 */
export function detectOrganismCode(df: DG.DataFrame): string | undefined {
  const col = findOrganismColumn(df);
  if (!col) return undefined;
  const codes = new Set<string>();
  for (let i = 0; i < col.length; i++) {
    if (col.isNone(i)) continue;
    const code = resolveOrganismCode(col.get(i) as string);
    if (code) codes.add(code);
    if (codes.size > 1) return undefined; // ambiguous — don't guess
  }
  return codes.size === 1 ? [...codes][0] : undefined;
}
