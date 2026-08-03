import * as grok from 'datagrok-api/grok';

/**
 * Shared UniProt subcellular-location service.
 *
 * The classifier, keyword map, and hex palette are ported VERBATIM from the
 * locked CK-omics contract (~/Downloads/ck/CKomics_tool2.py:1384-1436
 * `parse_subcellular_location` / `get_location_colors` and
 * Subcellular_Location_Classification_README.txt). Any divergence is a parity
 * regression — do not re-derive. Pure functions do no I/O.
 */

// --- Verbatim CK-omics 11-category keyword map (dict insertion order matters:
//     it is the tie-break for equal-position matches). ---

export const LOCATION_KEYWORDS: Record<string, string[]> = {
  'Nucleus': ['nucleus', 'nuclear', 'nucleoplasm', 'chromatin', 'nucleolus',
    'nuclear envelope', 'nuclear membrane', 'nuclear pore', 'nuclear speck',
    'cajal body', 'nuclear body'],
  'Cytoplasm': ['cytoplasm', 'cytosol', 'cytoplasmic', 'cytoskeleton',
    'microtubule', 'actin', 'intermediate filament', 'sarcomere', 'myofibril',
    'z disc', 'z disk', 'spindle', 'centriole', 'centrosome', 'cilium',
    'flagellum', 'stress fiber'],
  'Mitochondria': ['mitochondria', 'mitochondrial', 'mitochondrion',
    'mitochondrial matrix', 'mitochondrial membrane',
    'mitochondrial inner membrane', 'mitochondrial outer membrane',
    'mitochondrial intermembrane space'],
  'ER': ['endoplasmic reticulum', 'sarcoplasmic reticulum',
    'rough endoplasmic reticulum', 'smooth endoplasmic reticulum',
    'er-golgi intermediate compartment', 'ergic'],
  'Golgi': ['golgi', 'trans-golgi', 'golgi apparatus', 'golgi membrane',
    'cis-golgi', 'trans-golgi network', 'tgn'],
  'Plasma Membrane': ['plasma membrane', 'cell membrane', 'cell surface',
    'tight junction', 'gap junction', 'adherens junction', 'desmosome',
    'focal adhesion', 'cell junction', 'synapse', 'postsynaptic',
    'presynaptic', 'neuromuscular junction'],
  'Lysosome': ['lysosome', 'lysosomal', 'vacuole', 'lysosomal membrane',
    'autophagosome', 'phagosome'],
  'Peroxisome': ['peroxisome', 'peroxisomal', 'glyoxysome'],
  'Ribosome': ['ribosome', 'ribosomal', 'ribosomal protein', 'polysome'],
  'Extracellular': ['extracellular', 'secreted', 'extracellular space',
    'extracellular matrix', 'basement membrane', 'basal lamina', 'collagen',
    'ecm'],
  'Vesicles': ['vesicle', 'endosome', 'transport vesicle', 'secretory vesicle',
    'early endosome', 'late endosome', 'recycling endosome',
    'clathrin-coated vesicle', 'coated vesicle'],
};

/** Locked hex palette (README §COLOR SCHEME) as ARGB ints — same 0xFF… form as
 * volcano.ts:60-62. 11 categories + Unknown. */
export const LOCATION_COLORS: Record<string, number> = {
  'Nucleus': 0xFF000000 | parseInt('FF6B6B', 16),
  'Cytoplasm': 0xFF000000 | parseInt('ECDC44', 16),
  'Mitochondria': 0xFF000000 | parseInt('45B7D1', 16),
  'ER': 0xFF000000 | parseInt('96CEB4', 16),
  'Golgi': 0xFF000000 | parseInt('FAAFFE', 16),
  'Plasma Membrane': 0xFF000000 | parseInt('DDA0DD', 16),
  'Lysosome': 0xFF000000 | parseInt('F39C12', 16),
  'Peroxisome': 0xFF000000 | parseInt('E17055', 16),
  'Ribosome': 0xFF000000 | parseInt('A29BFE', 16),
  'Extracellular': 0xFF000000 | parseInt('FD79A8', 16),
  'Vesicles': 0xFF000000 | parseInt('FDCB6E', 16),
  'Unknown': 0xFF000000 | parseInt('CCCCCC', 16),
};

/** Escapes regex metacharacters (mirrors Python `re.escape` faithfully — the
 * current keywords contain none, but stay defensive). */
function escapeRegex(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

/** First location keyword (by earliest character position across ALL
 * categories) in the lowercased RAW text. Ties keep the first-iterated
 * category (strict `<`), so keyword-map insertion order is the tie-break.
 * The raw string is NOT pre-stripped (Pitfall 2). */
function findFirstLocation(text: string): string | null {
  if (!text) return null;
  const lower = text.toLowerCase();
  let earliestMatch: string | null = null;
  let earliestPos = Infinity;
  for (const category of Object.keys(LOCATION_KEYWORDS)) {
    for (const keyword of LOCATION_KEYWORDS[category]) {
      const m = new RegExp('\\b' + escapeRegex(keyword) + '\\b', 'i').exec(lower);
      if (m && m.index < earliestPos) {
        earliestPos = m.index;
        earliestMatch = category;
      }
    }
  }
  return earliestMatch;
}

/** Verbatim port of CK-omics `parse_subcellular_location`: subcellular field is
 * scanned fully first, GO cellular-component is fallback only, else 'Unknown'. */
export function parseSubcellularLocation(
  subcell: string | null | undefined,
  go: string | null | undefined,
): string {
  const subcellText = subcell == null ? '' : subcell;
  const goText = go == null ? '' : go;
  let location = findFirstLocation(subcellText);
  if (location) return location;
  location = findFirstLocation(goText);
  if (location) return location;
  return 'Unknown';
}

/** Result of parsing one UniProt stream TSV body. */
export interface StreamTsvResult {
  locByAcc: Map<string, string>;
  geneByAcc: Map<string, string>;
  /** Whether the winning location for an accession came from a reviewed entry. */
  reviewedByAcc: Map<string, boolean>;
}

/**
 * Pure positional parse of a UniProt stream TSV (fields order:
 * accession, cc_subcellular_location, go_c, reviewed, gene_primary). The first
 * line is the display-name header and is skipped; columns are read by POSITION
 * (`parts[0..4]`), never by header name (Pitfall 1). Within the text, the
 * priority merge is: first occurrence sets; a reviewed non-Unknown overwrites;
 * an existing Unknown is overwritten by any non-Unknown.
 */
export function mergeStreamTsv(text: string): StreamTsvResult {
  const locByAcc = new Map<string, string>();
  const geneByAcc = new Map<string, string>();
  const reviewedByAcc = new Map<string, boolean>();
  const lines = text.trim().split('\n');
  for (let li = 1; li < lines.length; li++) {
    const line = lines[li];
    if (!line) continue;
    const p = line.split('\t');
    const acc = (p[0] ?? '').trim();
    if (!acc) continue;
    const subcell = p[1] ?? '';
    const goc = p[2] ?? '';
    const reviewed = (p[3] ?? '').toLowerCase() === 'reviewed';
    const gene = (p[4] ?? '').trim();
    const cat = parseSubcellularLocation(subcell, goc);

    if (!locByAcc.has(acc)) {
      locByAcc.set(acc, cat);
      reviewedByAcc.set(acc, reviewed);
    } else if (reviewed && cat !== 'Unknown') {
      locByAcc.set(acc, cat);
      reviewedByAcc.set(acc, true);
    } else if (locByAcc.get(acc) === 'Unknown' && cat !== 'Unknown') {
      locByAcc.set(acc, cat);
      reviewedByAcc.set(acc, reviewed);
    }
    if (gene && (!geneByAcc.has(acc) || reviewed))
      geneByAcc.set(acc, gene);
  }
  return {locByAcc, geneByAcc, reviewedByAcc};
}

// --- Fetch + cache + D-03 fallback (Task 2) ---

/** Organism g:Profiler code → NCBI (species-level) taxonomy id. Covers all 9
 * organisms in enrichment's `ORGANISM_LIST`. Species-level ids are used so the
 * UniProt `taxonomy_id:` filter (hierarchical — matches the taxon and its
 * descendants) stays inclusive of strains (e.g. E. coli K-12 under 562).
 * Only the gene-fallback pass narrows by this; accession queries do not (A1:
 * accessions are globally unique, so filtering them risks excluding valid hits). */
export const ORGANISM_TAXONOMY: Record<string, number> = {
  hsapiens: 9606, mmusculus: 10090, rnorvegicus: 10116, scerevisiae: 4932,
  ecoli: 562, drerio: 7955, dmelanogaster: 7227, athaliana: 3702, celegans: 6239,
};

/** userDataStorage key for the cross-session accession → category cache.
 * Exported so the volcanoOptions OK handler can peek the cache to decide
 * which pre-OK toast wording (cold vs warm) to show — see package.ts
 * volcanoOptions. */
export const STORE = 'proteomics-subcell-loc';
/** Bump only if the locked keyword map ever changes (D-04: it is a contract,
 * so this is effectively fixed). Mismatch → stale cache discarded. */
const SCHEMA_V = '13-04-1';
const SCHEMA_KEY = '__schema_v';
const ACC_CHUNK = 100;
const GENE_CHUNK = 20;

/** Progress callback shape used by getSubcellularLocations and forwarded by
 * ensureLocationColumn / recomputeVolcano / volcanoOptions. The dialog's
 * DG.TaskBarProgressIndicator.update is the call site that consumes this.
 *
 * - 'fetch-acc'    — Pass 1 chunked accession-stream UniProt queries.
 * - 'fetch-gene'   — Pass 2 reviewed-by-gene fallback queries.
 * - 'init-column'  — bulk-init of the LOCATION_COL column on the DataFrame
 *                    (in volcano.ts) AFTER the network work completes. */
export type ProgressCb = (done: number, total: number,
  phase: 'fetch-acc' | 'fetch-gene' | 'init-column') => void;

/** Maximum in-flight UniProt stream chunks. 6 is a compromise: 4× the
 * sequential baseline on a real ~80-chunk Spectronaut Candidates file while
 * staying well under rest.uniprot.org's documented limits. Increase only
 * with measured evidence; UniProt's stream endpoint is shared infrastructure. */
const FETCH_CONCURRENCY = 6;

/** Wall-clock interval between incremental cache flushes during a fetch.
 * A single setInterval-driven flush from the orchestrator (not from inside
 * workers) sidesteps the per-chunk-counter race that a previous draft hit
 * when two workers reached the flush check between `await` points. */
const CACHE_FLUSH_INTERVAL_MS = 5000;

const VALID_CATEGORIES = new Set<string>([...Object.keys(LOCATION_KEYWORDS), 'Unknown']);

const STREAM_FIELDS = 'accession,cc_subcellular_location,go_c,reviewed,gene_primary';
const GENE_FIELDS = 'accession,gene_primary,cc_subcellular_location,go_c';

function chunk<T>(arr: T[], size: number): T[][] {
  const out: T[][] = [];
  for (let i = 0; i < arr.length; i += size) out.push(arr.slice(i, i + size));
  return out;
}

/** Index-based worker pool. Maintains a shared `nextIndex` counter; launches
 * `limit` async workers; each worker pulls indices and awaits `work(items[i], i)`,
 * storing results at `results[i]`; resolves with the full results array when
 * every worker drains. Used by getSubcellularLocations Pass 1/2 to cap
 * concurrent UniProt fetchProxy calls at FETCH_CONCURRENCY.
 *
 * Exported only so the test suite can verify the concurrency cap directly
 * with synthetic deferred-promise tasks (no network). */
export async function runWithConcurrency<T, R>(
  items: T[], limit: number, work: (item: T, idx: number) => Promise<R>,
): Promise<R[]> {
  const results: R[] = new Array(items.length);
  let nextIndex = 0;
  const workerCount = Math.min(Math.max(1, limit), items.length || 1);
  const workers: Promise<void>[] = [];
  for (let w = 0; w < workerCount; w++) {
    workers.push((async () => {
      while (true) {
        const i = nextIndex++;
        if (i >= items.length) return;
        results[i] = await work(items[i], i);
      }
    })());
  }
  await Promise.all(workers);
  return results;
}

function taxonomyClause(organism?: string): string {
  const tax = organism ? ORGANISM_TAXONOMY[organism] : undefined;
  return tax ? ` AND (taxonomy_id:${tax})` : '';
}

async function loadCache(): Promise<Record<string, string>> {
  try {
    const raw = (await grok.dapi.userDataStorage.get(STORE)) ?? {};
    if (raw[SCHEMA_KEY] !== SCHEMA_V)
      return {[SCHEMA_KEY]: SCHEMA_V}; // stale → discard
    return raw as Record<string, string>;
  } catch {
    return {[SCHEMA_KEY]: SCHEMA_V};
  }
}

/**
 * Resolves accession → subcellular category for every requested accession.
 * Cache-first (cross-session userDataStorage map, `__schema_v`-keyed), then a
 * chunked UniProt stream fetch via grok.dapi.fetchProxy (no raw fetch — CORS),
 * then the D-03 reviewed-by-gene fallback for Unknown+gene accessions. Never
 * throws out of the fetch loop: a failed batch warns and continues.
 */
export async function getSubcellularLocations(
  accessions: string[],
  organism?: string,
  progress?: ProgressCb,
): Promise<Map<string, string>> {
  const unique = [...new Set(accessions.filter((a) => a && a.length > 0))];
  const cache = await loadCache();

  const resolved = new Map<string, string>();
  const misses: string[] = [];
  for (const acc of unique) {
    const c = cache[acc];
    if (c !== undefined && VALID_CATEGORIES.has(c)) resolved.set(acc, c);
    else misses.push(acc);
  }

  const geneByAcc = new Map<string, string>();
  const fetched: Record<string, string> = {};
  // Applied ONLY to the Pass 2 gene fallback below — gene symbols are shared
  // across species, so without this a rat gene resolves to the human entry.
  const taxClause = taxonomyClause(organism);

  // Single timer-driven incremental flush. Workers DO NOT touch the cache —
  // they only mutate `resolved`/`fetched`/`geneByAcc`. The timer reads the
  // current `fetched` map at its tick and writes through, so a session
  // interrupted mid-fetch keeps every accession that finished. The previous
  // draft folded the flush into the per-chunk worker body and hit a documented
  // double-flush race when two workers reached the check between `await`
  // points — see plan 13-08 BLOCKER FIX notes.
  const flushCache = async (): Promise<void> => {
    if (Object.keys(fetched).length === 0) return;
    try {
      await grok.dapi.userDataStorage.put(STORE,
        {...cache, ...fetched, [SCHEMA_KEY]: SCHEMA_V});
    } catch (e: any) {
      console.warn(`Subcellular-location cache write failed: ${e?.message ?? e}`);
    }
  };
  const flushTimer = setInterval(() => { flushCache().catch(() => {}); },
    CACHE_FLUSH_INTERVAL_MS);

  try {
    // Pass 1: chunked accession stream queries, run with bounded concurrency.
    const accChunks = chunk(misses, ACC_CHUNK);
    const totalAcc = accChunks.length;
    let doneAcc = 0;
    await runWithConcurrency(accChunks, FETCH_CONCURRENCY, async (group) => {
      // Pass 1 is deliberately NOT taxonomy-filtered: accessions are globally
      // unique, so a wrong/mis-set organism must never exclude a valid protein.
      const q = group.map((a) => `accession:${a}`).join(' OR ');
      const url = 'https://rest.uniprot.org/uniprotkb/stream' +
        `?query=(${encodeURIComponent(q)})` +
        `&fields=${STREAM_FIELDS}&format=tsv`;
      try {
        const resp = await grok.dapi.fetchProxy(url);
        if (!resp.ok) {
          console.warn(`UniProt stream returned status ${resp.status}; continuing`);
          return;
        }
        const {locByAcc, geneByAcc: g} = mergeStreamTsv(await resp.text());
        for (const [acc, loc] of locByAcc) {
          resolved.set(acc, loc);
          fetched[acc] = loc;
        }
        for (const [acc, gene] of g) geneByAcc.set(acc, gene);
      } catch (e: any) {
        console.warn(`UniProt stream batch failed: ${e?.message ?? e}; continuing`);
      }
      // Increment AFTER the chunk completes — bounded concurrency means
      // worker order can interleave, so the callback signal is non-decreasing
      // (not strictly monotonic). Tests assert the non-decreasing invariant.
      doneAcc++;
      progress?.(doneAcc, totalAcc, 'fetch-acc');
    });

    // Pass 2 (D-03): reviewed-by-gene fallback for Unknown accessions with a gene.
    const geneToUnknownAccs = new Map<string, string[]>();
    for (const acc of misses) {
      const loc = resolved.get(acc) ?? 'Unknown';
      const gene = geneByAcc.get(acc);
      if (loc === 'Unknown' && gene) {
        if (!geneToUnknownAccs.has(gene)) geneToUnknownAccs.set(gene, []);
        geneToUnknownAccs.get(gene)!.push(acc);
      }
    }
    const genes = [...geneToUnknownAccs.keys()];
    const geneChunks = chunk(genes, GENE_CHUNK);
    const totalGene = geneChunks.length;
    let doneGene = 0;
    await runWithConcurrency(geneChunks, FETCH_CONCURRENCY, async (group) => {
      const q = group.map((g) => `gene_exact:${g}`).join(' OR ');
      const url = 'https://rest.uniprot.org/uniprotkb/stream' +
        `?query=(${encodeURIComponent(q)}) AND (reviewed:true)${encodeURIComponent(taxClause)}` +
        `&fields=${GENE_FIELDS}&format=tsv`;
      try {
        const resp = await grok.dapi.fetchProxy(url);
        if (!resp.ok) {
          console.warn(`UniProt gene fallback returned status ${resp.status}; continuing`);
          return;
        }
        // Positional: accession, gene_primary, cc_subcellular_location, go_c.
        const lines = (await resp.text()).trim().split('\n');
        for (let li = 1; li < lines.length; li++) {
          const parts = lines[li].split('\t');
          const gene = (parts[1] ?? '').trim();
          const loc = parseSubcellularLocation(parts[2] ?? '', parts[3] ?? '');
          if (!gene || loc === 'Unknown') continue;
          for (const acc of geneToUnknownAccs.get(gene) ?? []) {
            resolved.set(acc, loc);
            fetched[acc] = loc;
          }
        }
      } catch (e: any) {
        console.warn(`UniProt gene fallback batch failed: ${e?.message ?? e}; continuing`);
      }
      doneGene++;
      progress?.(doneGene, totalGene, 'fetch-gene');
    });

    // Anything still unresolved → Unknown (explicit, also cached so we don't
    // re-query it every session).
    for (const acc of unique) {
      if (!resolved.has(acc)) resolved.set(acc, 'Unknown');
      if (fetched[acc] === undefined) fetched[acc] = resolved.get(acc)!;
    }
  } finally {
    clearInterval(flushTimer);
    // Single final write-through. Schema invariant preserved on every flush —
    // tests assert this.
    await flushCache();
  }

  const out = new Map<string, string>();
  for (const acc of unique) out.set(acc, resolved.get(acc) ?? 'Unknown');
  return out;
}

// --- Cached per-accession UniProt JSON (closes the folded cache-uniprot todo;
//     uniprot-panel.ts delegates here so re-clicking a protein doesn't refetch). ---

const _entryCache = new Map<string, unknown>();

/** Cached single-accession UniProt JSON fetch. Mirrors the original panel
 * discipline: grok.dapi.fetchProxy, warn + null on !resp.ok, warn + null on
 * throw. Successful and definitive-miss (non-OK) results are cached for the
 * session; transient throws are not. */
export async function fetchUniProtEntry(accession: string): Promise<unknown | null> {
  if (_entryCache.has(accession)) return _entryCache.get(accession) ?? null;
  const url = `https://rest.uniprot.org/uniprotkb/${encodeURIComponent(accession)}.json` +
    '?fields=accession,protein_name,gene_names,organism_name,cc_function,go';
  try {
    const resp = await grok.dapi.fetchProxy(url);
    if (!resp.ok) {
      console.warn(`UniProt fetch for ${accession} returned status ${resp.status}`);
      _entryCache.set(accession, null);
      return null;
    }
    const json = await resp.json();
    _entryCache.set(accession, json);
    return json;
  } catch (e: any) {
    console.warn(`UniProt fetch for ${accession} failed:`, e?.message ?? e);
    return null;
  }
}
