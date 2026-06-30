import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {findColumn} from './column-detection';
import {SEMTYPE} from './proteomics-types';
import {runWithConcurrency} from '../analysis/subcellular-location';

// ---------------------------------------------------------------------------
// LOCKED CLIENT CONTRACT — verbatim port of CKomics_tool2.py
//   • improve_gene_labels_with_ensrnog_marking  (lines 895-1059)
//   • extract_readable_description              (lines 1062-1104)
//   • get_ensembl_annotations request shape     (lines 756-855)
// D-08: inline `*` (grouped) / `†` (predicted/reclassified) markers; raw ID
//        kept in `Source ID` column.
// D-09: Ensembl REST POST /lookup/id batched ≤1000 IDs per request, cached
//        cross-session via grok.dapi.userDataStorage with `__schema_v` invalidation.
// D-10: duplicate gene names disambiguated with `(Source ID)` suffix; one
//        grok.shell.warning emitted per import with affected-row count.
// ---------------------------------------------------------------------------

/** userDataStorage key for the cross-session id → EnsemblEntry cache.
 * Exported so callers (and tests) can peek the cache directly. */
export const STORE_GENE_LABELS = 'proteomics-gene-labels';

/** Bump only if the locked CK-omics resolver contract changes. Mismatch in
 * userDataStorage → stale cache discarded. */
export const SCHEMA_V_GENE_LABELS = '14-r1-1';
const SCHEMA_KEY = '__schema_v';

/** Ensembl rest.ensembl.org /lookup/id hard cap (Pitfall 3). */
const ENS_CHUNK = 1000;

/** Matches the Phase 13 UniProt cap. Ensembl rate budget is 55k req/hour
 * (≈15 req/s); FETCH_CONCURRENCY = 6 keeps us well under that bound. */
const FETCH_CONCURRENCY = 6;

/** Predicted-pattern prefix list (CK-omics line 898-909, in this exact order). */
export const PREDICTED_PREFIXES = [
  'ENSRNOG', 'ENSMUSG', 'LOC', 'ENSG', 'ENSDARG',
  'ENSRNO', 'MGP_', 'RGD', 'AABR',
] as const;

/** Ensembl-eligible test per CK-omics line 933: only IDs starting with 'ENS'
 * OR 'MGP_' are sent to /lookup/id. LOC / RGD / AABR are NCBI / RGD /
 * Affymetrix respectively — they keep the raw ID + `†` marker only. */
export function isEnsemblEligible(id: string): boolean {
  return id.startsWith('ENS') || id.startsWith('MGP_');
}

/** True iff `id` starts with any locked predicted-pattern prefix. */
export function isPredicted(id: string): boolean {
  for (const p of PREDICTED_PREFIXES)
    if (id.startsWith(p)) return true;
  return false;
}

export type SpeciesCode =
  'homo_sapiens' | 'mus_musculus' | 'rattus_norvegicus' | 'danio_rerio';

/** Prefix → Ensembl species code. ENSRNO falls back to rattus_norvegicus.
 * Non-Ensembl-eligible prefixes (LOC / RGD / AABR) return null so callers
 * skip them at the POST stage. */
export function detectSpecies(id: string): SpeciesCode | null {
  if (id.startsWith('ENSG')) return 'homo_sapiens';
  if (id.startsWith('ENSMUSG')) return 'mus_musculus';
  if (id.startsWith('MGP_')) return 'mus_musculus';
  if (id.startsWith('ENSRNOG')) return 'rattus_norvegicus';
  if (id.startsWith('ENSRNO')) return 'rattus_norvegicus';
  if (id.startsWith('ENSDARG')) return 'danio_rerio';
  return null;
}

export interface EnsemblEntry {
  display_name?: string;
  external_name?: string;
  description?: string;
  species?: string;
  biotype?: string;
  object_type?: string;
}

/** Verbatim port of CKomics extract_readable_description (lines 1062-1104).
 * Returns null when the cleaned description is empty, still says
 * "uncharacterized", or is too short to be useful. */
export function extractReadableDescription(description: string | null | undefined): string | null {
  if (description == null || description === '') return null;
  let s = String(description);

  // Strip [organism] suffix.
  s = s.split('[')[0].trim();

  // Strip the various "Predicted to ..." preambles, case-insensitive, in
  // order — `enable` / `be involved in` / `be located in` / `be part of` /
  // `be` first because the bare `Predicted to ` would otherwise consume them.
  s = s.replace(/^Predicted to enable\s+/i, '');
  s = s.replace(/^Predicted to be involved in\s+/i, '');
  s = s.replace(/^Predicted to be located in\s+/i, '');
  s = s.replace(/^Predicted to be part of\s+/i, '');
  s = s.replace(/^Predicted to be\s+/i, '');
  s = s.replace(/^Predicted to\s+/i, '');

  // First sentence only.
  const sentences = s.split('.');
  if (sentences.length > 0) s = sentences[0].trim();

  // Strip species suffixes.
  s = s.replace(/_RAT$/i, '');
  s = s.replace(/_MOUSE$/i, '');
  s = s.replace(/_HUMAN$/i, '');

  // Capitalize first letter if lowercase.
  if (s.length > 0 && s[0] === s[0].toLowerCase() && s[0] !== s[0].toUpperCase())
    s = s[0].toUpperCase() + s.slice(1);

  // Truncate runaway descriptions (CK-omics line 1097-1098).
  if (s.length > 60) s = s.slice(0, 57) + '...';

  if (s === '' || s.length < 3) return null;
  if (s.toLowerCase().includes('uncharacterized')) return null;

  return s;
}

/** CK-omics line 1011 — `improved_name + '*'? + '†'?`. The dagger is appended
 * UNCONDITIONALLY when the resolver acted on a predicted row (CK-omics
 * `final_label = improved_name + ('*' if has_asterisk else '') + '†'`). */
export function applyMarkerRules(improvedName: string,
  hadAsterisk: boolean, wasReclassified: boolean): string {
  return improvedName + (hadAsterisk ? '*' : '') + (wasReclassified ? '†' : '');
}

function chunk<T>(arr: T[], size: number): T[][] {
  const out: T[][] = [];
  for (let i = 0; i < arr.length; i += size) out.push(arr.slice(i, i + size));
  return out;
}

async function loadCache(): Promise<Record<string, EnsemblEntry | string>> {
  try {
    const raw = (await grok.dapi.userDataStorage.get(STORE_GENE_LABELS)) ?? {};
    if ((raw as any)[SCHEMA_KEY] !== SCHEMA_V_GENE_LABELS)
      return {[SCHEMA_KEY]: SCHEMA_V_GENE_LABELS};
    return raw as Record<string, EnsemblEntry | string>;
  } catch {
    return {[SCHEMA_KEY]: SCHEMA_V_GENE_LABELS};
  }
}

async function flushCache(cache: Record<string, EnsemblEntry | string>,
  fetched: Record<string, EnsemblEntry>): Promise<void> {
  try {
    await grok.dapi.userDataStorage.put(STORE_GENE_LABELS,
      {...cache, ...fetched, [SCHEMA_KEY]: SCHEMA_V_GENE_LABELS});
  } catch (e: any) {
    console.warn(`Gene-label cache write failed: ${e?.message ?? e}`);
  }
}

/** One POST to https://rest.ensembl.org/lookup/id with body `{ids: [...]}`.
 * Honours a single 429 retry-after backoff per chunk. Warn-and-continue on
 * any other error; the resolver must not block import on Ensembl flakes. */
export async function lookupEnsemblBatch(ids: string[]): Promise<Map<string, EnsemblEntry>> {
  const out = new Map<string, EnsemblEntry>();
  if (ids.length === 0) return out;

  const post = async (): Promise<Response> => grok.dapi.fetchProxy(
    'https://rest.ensembl.org/lookup/id',
    {
      method: 'POST',
      headers: {'Content-Type': 'application/json', 'Accept': 'application/json'},
      body: JSON.stringify({ids}),
    },
  );

  let resp: Response;
  try {
    resp = await post();
  } catch (e: any) {
    console.warn(`Ensembl lookup/id batch failed (network): ${e?.message ?? e}`);
    return out;
  }

  if (resp.status === 429) {
    // Single retry per CK-omics line 822-825 contract; Retry-After header in
    // seconds (fractional allowed). Cap the back-off at 30s to keep import
    // responsive — the resolver degrades gracefully on continued failure.
    const retryAfter = parseFloat(resp.headers.get('Retry-After') ?? '1');
    const waitMs = Math.min(30_000, Math.max(0, retryAfter * 1000));
    await new Promise((res) => setTimeout(res, waitMs));
    try {
      resp = await post();
    } catch (e: any) {
      console.warn(`Ensembl lookup/id retry failed: ${e?.message ?? e}`);
      return out;
    }
  }

  if (!resp.ok) {
    console.warn(`Ensembl lookup/id returned status ${resp.status}; continuing`);
    return out;
  }

  let body: any;
  try {
    body = await resp.json();
  } catch (e: any) {
    console.warn(`Ensembl lookup/id JSON parse failed: ${e?.message ?? e}`);
    return out;
  }

  if (body && typeof body === 'object') {
    for (const [id, entry] of Object.entries(body)) {
      if (entry && typeof entry === 'object') {
        const e = entry as Record<string, unknown>;
        // Type-guard every field read (threat T-14-01-T1).
        const safe: EnsemblEntry = {};
        if (typeof e.display_name === 'string') safe.display_name = e.display_name;
        if (typeof e.external_name === 'string') safe.external_name = e.external_name;
        if (typeof e.description === 'string') safe.description = e.description;
        if (typeof e.species === 'string') safe.species = e.species;
        if (typeof e.biotype === 'string') safe.biotype = e.biotype;
        if (typeof e.object_type === 'string') safe.object_type = e.object_type;
        out.set(id, safe);
      }
    }
  }
  return out;
}

/** Three-level fallback per CK-omics lines 800-804:
 *   external_name → display_name → description.split('[')[0].strip() → raw ID.
 * Acceptance gate per CK-omics lines 977-988: reject the candidate if it equals
 * the raw ID OR starts with a predicted prefix. */
function pickBestName(rawId: string, entry: EnsemblEntry | undefined): string | null {
  if (!entry) return null;
  const candidates: Array<string | undefined> = [
    entry.external_name,
    entry.display_name,
    entry.description ? entry.description.split('[')[0].trim() : undefined,
  ];
  for (const c of candidates) {
    if (c && c !== rawId && !isPredicted(c)) return c;
  }
  return null;
}

export type GeneLabelProgress = (done: number, total: number,
  phase: 'lookup' | 'apply') => void;

/** Public entry — runs after every parser produces its DataFrame. Always
 * creates `Display Name` and `Source ID` columns even when no predicted IDs
 * are present (Pitfall 9: downstream Volcano label bindings depend on the
 * `Display Name` invariant). The resolver degrades gracefully on Ensembl
 * outages — analyst sees raw IDs + a single warning toast.
 */
export async function resolveGeneLabels(df: DG.DataFrame,
  progress?: GeneLabelProgress): Promise<void> {
  const geneCol = findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene names', 'gene symbol']);
  const n = df.rowCount;
  const displayArr: string[] = new Array(n).fill('');
  const sourceArr: string[] = new Array(n).fill('');

  // Pre-fill from Gene name (when present); the resolver overrides specific
  // rows below. Rows without a Gene name source still get '' for both, so the
  // columns are always rectangular.
  if (geneCol != null) {
    for (let i = 0; i < n; i++) {
      const raw = geneCol.get(i);
      const s = (raw == null) ? '' : String(raw);
      displayArr[i] = s;
    }
  }

  // Identify rows whose Gene name (stripped of pre-existing * and †) matches a
  // predicted-pattern prefix. CK-omics line 916-925: also walks FirstGene /
  // FirstUniProt columns. Datagrok's parsers normalize to a single Gene name
  // column, so checking that column is sufficient — but if the proteomics
  // primary protein id column happens to also be a predicted Ensembl ID, we
  // honour the same logic on it.
  const idCol = findColumn(df, SEMTYPE.PROTEIN_ID,
    ['protein ids', 'primary protein id', 'protein id', 'majority protein ids', 'accession', 'uniprot']);

  /** Per-row work item carrying the raw ID we want to resolve. */
  type Predicted = {row: number; clean: string; hadAsterisk: boolean; source: 'gene' | 'protein'};
  const predicted: Predicted[] = [];

  const recordIfPredicted = (row: number, raw: string, source: 'gene' | 'protein'): boolean => {
    if (!raw) return false;
    const hadAsterisk = raw.includes('*');
    const clean = raw.replace(/[*†]/g, '');
    if (!isPredicted(clean)) return false;
    predicted.push({row, clean, hadAsterisk, source});
    return true;
  };

  for (let i = 0; i < n; i++) {
    const g = geneCol != null ? geneCol.get(i) : null;
    const gs = g == null ? '' : String(g);
    if (recordIfPredicted(i, gs, 'gene')) continue;
    // Fall through to protein id only when gene name isn't a predicted ID.
    if (idCol != null) {
      const p = idCol.get(i);
      const ps = p == null ? '' : String(p).split(';')[0];
      recordIfPredicted(i, ps, 'protein');
    }
  }

  // Always ensure the two new columns exist (Pitfall 9). Remove first if a
  // re-run is in progress so we re-init cleanly.
  for (const name of ['Display Name', 'Source ID']) {
    const existing = df.col(name);
    if (existing != null) df.columns.remove(name);
  }
  // Fast path: no predicted IDs → write columns and exit, but always create
  // them. Display Name = raw gene name (or ''), Source ID = ''.
  if (predicted.length === 0) {
    const dCol = df.columns.addNewString('Display Name');
    dCol.init((i) => displayArr[i]);
    dCol.semType = SEMTYPE.DISPLAY_NAME;
    const sCol = df.columns.addNewString('Source ID');
    sCol.init((i) => sourceArr[i]);
    sCol.semType = SEMTYPE.SOURCE_ID;
    return;
  }

  // Collect unique predicted IDs and the Ensembl-eligible subset.
  const uniqueIds = [...new Set(predicted.map((p) => p.clean))];
  const ensIds = uniqueIds.filter((id) => isEnsemblEligible(id));

  const cache = await loadCache();
  const fetched: Record<string, EnsemblEntry> = {};
  const lookupResult = new Map<string, EnsemblEntry>();

  // Seed from cache.
  for (const id of ensIds) {
    const c = cache[id];
    if (c && typeof c === 'object') lookupResult.set(id, c as EnsemblEntry);
  }
  const misses = ensIds.filter((id) => !lookupResult.has(id));

  let unresolvedAny = false;
  if (misses.length > 0) {
    const chunks = chunk(misses, ENS_CHUNK);
    const total = chunks.length;
    let done = 0;
    try {
      await runWithConcurrency(chunks, FETCH_CONCURRENCY, async (group) => {
        const got = await lookupEnsemblBatch(group);
        for (const [id, entry] of got) {
          lookupResult.set(id, entry);
          fetched[id] = entry;
        }
        // Anything in `group` that didn't come back is unresolved by Ensembl —
        // we still mark the row reclassified, but Display Name falls back to
        // the raw ID. Track this so the warning toast can name it.
        for (const id of group) if (!got.has(id)) unresolvedAny = true;
        done++;
        progress?.(done, total, 'lookup');
      });
    } finally {
      await flushCache(cache, fetched);
    }
  }

  // Apply the marker / fallback rules per predicted row.
  // CK-omics fallback chain when Ensembl yielded nothing: try the
  // ProteinDescriptions column (CK-omics line 997-1005). Datagrok parsers
  // surface this as "Protein Descriptions" (Spectronaut) or "Protein names"
  // (MaxQuant); we accept either by name.
  const descCol = df.col('Protein Descriptions') ?? df.col('Protein names')
    ?? df.col('Description') ?? df.col('protein_descriptions');

  let reclassified = 0;
  for (const p of predicted) {
    const entry = lookupResult.get(p.clean);
    let improved: string | null = pickBestName(p.clean, entry);

    if (improved == null && entry != null) {
      const cleaned = extractReadableDescription(entry.description ?? null);
      if (cleaned != null && cleaned !== p.clean && !isPredicted(cleaned))
        improved = cleaned;
    }

    // ProteinDescriptions fallback.
    if (improved == null && descCol != null) {
      const dRaw = descCol.get(p.row);
      const cleaned = extractReadableDescription(dRaw == null ? null : String(dRaw));
      if (cleaned != null && cleaned !== p.clean && !isPredicted(cleaned) && cleaned.length > 3)
        improved = cleaned;
    }

    // If still nothing, keep raw ID — marker still appended (CK-omics line 1011).
    const finalName = applyMarkerRules(improved ?? p.clean, p.hadAsterisk, true);
    displayArr[p.row] = finalName;
    sourceArr[p.row] = p.clean;
    reclassified++;
  }
  progress?.(predicted.length, predicted.length, 'apply');

  // D-10 duplicate disambiguation. Group rows by Display Name; if a name maps
  // to >1 distinct Source ID, append `(Source ID)` to each of those rows.
  const byName = new Map<string, Set<string>>();
  const rowsByName = new Map<string, number[]>();
  for (let i = 0; i < n; i++) {
    const name = displayArr[i];
    const src = sourceArr[i];
    if (!name || !src) continue;
    if (!byName.has(name)) byName.set(name, new Set());
    byName.get(name)!.add(src);
    if (!rowsByName.has(name)) rowsByName.set(name, []);
    rowsByName.get(name)!.push(i);
  }
  let duplicateRows = 0;
  for (const [name, srcs] of byName) {
    if (srcs.size <= 1) continue;
    for (const r of rowsByName.get(name) ?? []) {
      const sid = sourceArr[r];
      if (!sid) continue;
      displayArr[r] = `${name} (${sid})`;
      duplicateRows++;
    }
  }
  if (duplicateRows > 0) {
    grok.shell.warning(`${duplicateRows} duplicate gene names disambiguated with source IDs.`);
  }

  if (unresolvedAny) {
    const unresolvedCount = predicted.filter((p) =>
      isEnsemblEligible(p.clean) && !lookupResult.has(p.clean)).length;
    grok.shell.warning(
      `Ensembl gene-label resolution unavailable. Showing raw IDs (${unresolvedCount} unresolved). Re-import to retry.`,
    );
  }

  const dCol = df.columns.addNewString('Display Name');
  dCol.init((i) => displayArr[i]);
  dCol.semType = SEMTYPE.DISPLAY_NAME;
  const sCol = df.columns.addNewString('Source ID');
  sCol.init((i) => sourceArr[i]);
  sCol.semType = SEMTYPE.SOURCE_ID;

  // Silence unused-warning for `reclassified` — kept for future telemetry.
  void reclassified;
}
