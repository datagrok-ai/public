/*
 * Stage 1 — RCSB Data API client (vendored).
 *
 * Single GraphQL query (ENTRY_INFO). Body lifted verbatim from
 * packages/BiostructureViewer/src/utils/rcsb-gql-adapter.ts:40-121 (BSV's
 * `RcsbGraphQLAdapter` is not callable cross-package — not registered as a
 * platform function — so we vendor what we need rather than wait for it to
 * be exposed).
 *
 * Drops the polymer/batch query branches from the vendored copy. v1.5 will
 * re-vendor them when the target-name Search-API leg lands.
 *
 * Blueprint reference - Phase 1, section 1 reuse inventory (BSV vendoring), Q11.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';

import {_package} from './package';
import type {PipelineOptions} from './orchestrator-types';
import {isLocalId, localQcRows} from './local-pdb-store';

const GRAPHQL_ENDPOINT = 'https://data.rcsb.org/graphql';
const SEARCH_ENDPOINT = 'https://search.rcsb.org/rcsbsearch/v2/query';

// Vendored from packages/BiostructureViewer/src/utils/rcsb-gql-adapter.ts:40-121.
// Stripped of polymer organism/host fields v1 doesn't consume; otherwise verbatim.
const ENTRY_INFO_QUERY = `
  query GetEntryInfo($entryId: String!) {
    entry(entry_id: $entryId) {
      rcsb_id
      struct { title }
      rcsb_entry_info {
        resolution_combined
        experimental_method
      }
      nonpolymer_entities {
        rcsb_nonpolymer_entity_container_identifiers {
          auth_asym_ids
        }
        nonpolymer_comp {
          chem_comp { id name formula }
          rcsb_chem_comp_descriptor { SMILES SMILES_stereo }
        }
      }
    }
  }
`;

interface NonpolymerEntity {
  rcsb_nonpolymer_entity_container_identifiers?: {auth_asym_ids?: string[]};
  nonpolymer_comp?: {
    chem_comp?: {id: string; name?: string; formula?: string};
    rcsb_chem_comp_descriptor?: {SMILES?: string; SMILES_stereo?: string};
  };
}

interface EntryInfo {
  rcsb_id: string;
  struct?: {title?: string};
  rcsb_entry_info?: {resolution_combined?: number[]; experimental_method?: string};
  nonpolymer_entities?: NonpolymerEntity[];
}

/** GraphQL POST via fetchProxy (CORS-safe). Throws on HTTP / GraphQL errors. */
async function executeEntryInfo(pdbId: string): Promise<EntryInfo | null> {
  const resp = await grok.dapi.fetchProxy(GRAPHQL_ENDPOINT, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({query: ENTRY_INFO_QUERY, variables: {entryId: pdbId.toUpperCase()}}),
  });
  if (!resp.ok)
    throw new Error(`RCSB Data API HTTP ${resp.status} for ${pdbId}`);
  const json = await resp.json();
  if (json.errors)
    throw new Error(`RCSB Data API errors for ${pdbId}: ${JSON.stringify(json.errors)}`);
  return json.data?.entry ?? null;
}

/** Load the bundled cofactor deny-list as a Set of upper-cased CCD IDs. */
async function loadCofactorDenylist(): Promise<Set<string>> {
  const df = await _package.files.readCsv('cofactor-denylist.csv');
  const idCol = df.col('ccd_id');
  const out = new Set<string>();
  if (!idCol) return out;
  for (let i = 0; i < df.rowCount; i++) {
    const v = String(idCol.get(i) ?? '').trim().toUpperCase();
    if (v) out.add(v);
  }
  return out;
}

/** Validate SMILES + extract exactmw via WASM RDKit; returns null on parse failure. */
function validateAndWeigh(rdkit: {get_mol: (s: string) => any}, smiles: string): {mw: number | null} | null {
  let mol: any = null;
  try {
    mol = rdkit.get_mol(smiles);
    if (!mol || !mol.is_valid?.()) return null;
    let mw: number | null = null;
    try {
      const d = JSON.parse(mol.get_descriptors());
      if (typeof d.exactmw === 'number') mw = d.exactmw;
    } catch { /* keep mw null */ }
    return {mw};
  } catch {
    return null;
  } finally {
    mol?.delete?.();
  }
}

/**
 * Stage 1 — enrich a list of PDB IDs into a pdb_qc DataFrame.
 *
 * For each PDB:
 *   - Fetch ENTRY_INFO via RCSB Data API GraphQL.
 *   - Apply QC: resolution <= maxResolution; X-ray-only if requireXray.
 *   - For each nonpolymer ligand:
 *     - Reject if in cofactor deny-list.
 *     - Validate SMILES via WASM RDKit; reject if unparseable.
 *     - Compute exactmw; reject if < minLigandMw.
 *     - Emit one row per (pdb, ligand) — multi-ligand entries yield multiple rows.
 *
 * Schema (per Appendix B):
 *   pdb_id, resolution, experimental_method, r_free, ligand_comp_id, ligand_chain,
 *   ligand_formula_weight, ligand_smiles (semType=Molecule), ligand_rscc
 *
 * r_free and ligand_rscc are nullable in v1 — the ENTRY_INFO query doesn't fetch
 * them. v1.5 can add them via a second query if Stage 5a's enrichment relies on
 * the values.
 */
/**
 * Per-PDB record explaining why one of the input IDs did NOT make it into
 * the accepted pdb_qc DataFrame. The wizard's Step 1 panel renders this
 * as a "Dropped" section so users see explicitly what their filters
 * removed and at what RCSB-reported value.
 */
export interface DroppedPdb {
  pdb_id: string;
  /** Plain-English reason (e.g. "Resolution 2.60 Å > max 2.50 Å"). */
  reason: string;
  /** Best-effort metadata from RCSB before the drop (for the table). */
  resolution: number | null;
  method: string;
}

/** Result envelope returned by enrichPdbList: the accepted rows as a
 *  DataFrame plus a list of dropped PDBs with reasons. */
export interface EnrichResult {
  accepted: DG.DataFrame;
  dropped: DroppedPdb[];
}

export async function enrichPdbList(
  pdbIds: string[], options: PipelineOptions,
): Promise<EnrichResult> {
  const cofactors = await loadCofactorDenylist();
  const rdkit = await getRdKitModule();

  const rows: Record<string, any>[] = [];
  const droppedList: DroppedPdb[] = [];
  let dropped = 0;
  let droppedResolution = 0;
  let droppedXray = 0;
  let droppedNoLigand = 0;
  let droppedCofactor = 0;
  let droppedSmiles = 0;
  let droppedMw = 0;

  // Fetch all PDB entries in parallel — the previous sequential `for...await`
  // serialized 5 EGFR demo PDBs into 5 sequential RCSB GraphQL round-trips
  // (~5s total instead of ~1s). Each fetch settles independently so one
  // 404 / network error doesn't take the batch down with it; failures are
  // tracked alongside successes and surfaced via the same warning toast.
  // User-dropped local structures have no RCSB entry — exclude them from the
  // GraphQL fetch; their QC rows are synthesized below via localQcRows().
  const upper = pdbIds.map((s) => s.toUpperCase()).filter((id) => !isLocalId(id));
  const entries: Array<{pdbId: string; entry: EntryInfo | null; error?: Error}> =
    await Promise.all(upper.map(async (pdbId) => {
      try {
        return {pdbId, entry: await executeEntryInfo(pdbId)};
      } catch (e) {
        return {pdbId, entry: null, error: e as Error};
      }
    }));

  for (const {pdbId, entry, error} of entries) {
    if (error) {
      grok.shell.warning(`${pdbId}: enrichment failed (${error.message}). Skipped.`);
      droppedList.push({pdb_id: pdbId, reason: `RCSB fetch failed: ${error.message}`,
        resolution: null, method: ''});
      dropped++;
      continue;
    }
    if (!entry) {
      droppedList.push({pdb_id: pdbId, reason: 'RCSB returned no entry record',
        resolution: null, method: ''});
      dropped++; continue;
    }

    const resolution = entry.rcsb_entry_info?.resolution_combined?.[0];
    const method = entry.rcsb_entry_info?.experimental_method ?? '';
    if (typeof resolution === 'number' && resolution > options.maxResolution) {
      droppedList.push({pdb_id: pdbId,
        reason: `Resolution ${resolution.toFixed(2)} Å > max ${options.maxResolution.toFixed(2)} Å`,
        resolution, method});
      droppedResolution++; dropped++; continue;
    }
    // Classify the experimental method from the RCSB record. The string
    // varies in case across entries — substring match keeps the rules robust:
    //  - "X-ray" / "X-RAY DIFFRACTION" → X-ray
    //  - "Solution NMR" / "Solid-state NMR" / "NMR" → NMR
    //  - "Electron microscopy" / "Cryo-EM" → cryo-EM
    //  - "Computational" / "Predicted" / "AlphaFold" → predicted model
    const m = method.toLowerCase();
    const isXray = m.includes('x-ray');
    const isNmr = m.includes('nmr');
    const isCryoEm = m.includes('electron microscopy') || m.includes('cryo-em');
    const isAlphaFold = m.includes('alphafold') || m.includes('predicted') ||
                        m.includes('computational') || m.includes('computed');
    const allowed =
      (isXray && options.allowXray) ||
      (isNmr && options.allowNmr) ||
      (isCryoEm && options.allowCryoEm) ||
      (isAlphaFold && options.allowAlphaFold);
    if (!allowed) {
      const methodLabel = isXray ? 'X-ray' : isNmr ? 'NMR' :
        isCryoEm ? 'Cryo-EM' : isAlphaFold ? 'AI predicted' : (method || 'unknown');
      droppedList.push({pdb_id: pdbId,
        reason: `Method "${methodLabel}" not enabled in the filter`,
        resolution: resolution ?? null, method});
      droppedXray++; dropped++; continue;
    }

    let ligandsEmitted = 0;
    for (const ne of entry.nonpolymer_entities ?? []) {
      const compId = (ne.nonpolymer_comp?.chem_comp?.id ?? '').toUpperCase();
      if (!compId) continue;
      if (cofactors.has(compId)) { droppedCofactor++; continue; }
      const smiles = ne.nonpolymer_comp?.rcsb_chem_comp_descriptor?.SMILES_stereo ??
                     ne.nonpolymer_comp?.rcsb_chem_comp_descriptor?.SMILES;
      if (!smiles) { droppedSmiles++; continue; }
      const validated = validateAndWeigh(rdkit, smiles);
      if (!validated) { droppedSmiles++; continue; }
      if (validated.mw !== null && validated.mw < options.minLigandMw) {
        droppedMw++; continue;
      }

      const chains = ne.rcsb_nonpolymer_entity_container_identifiers?.auth_asym_ids ?? [];
      rows.push({
        pdb_id: pdbId,
        resolution: resolution ?? null,
        experimental_method: method,
        r_free: null,
        ligand_comp_id: compId,
        ligand_chain: chains[0] ?? '',
        ligand_formula_weight: validated.mw,
        ligand_smiles: smiles,
        ligand_rscc: null,
      });
      ligandsEmitted++;
    }
    if (ligandsEmitted === 0) {
      droppedList.push({pdb_id: pdbId,
        reason: 'No usable ligand (cofactors / unparseable SMILES / below MW cutoff)',
        resolution: resolution ?? null, method});
      droppedNoLigand++; dropped++;
    }
  }

  if (dropped > 0) {
    const parts: string[] = [];
    if (droppedResolution) parts.push(`resolution>${options.maxResolution}A: ${droppedResolution}`);
    if (droppedXray) parts.push(`non-X-ray: ${droppedXray}`);
    if (droppedNoLigand) parts.push(`no usable ligand: ${droppedNoLigand}`);
    if (droppedCofactor) parts.push(`only cofactors: ${droppedCofactor}`);
    if (droppedSmiles) parts.push(`bad SMILES: ${droppedSmiles}`);
    if (droppedMw) parts.push(`MW<${options.minLigandMw}: ${droppedMw}`);
    const breakdown = parts.length ? ` (${parts.join('; ')})` : '';
    grok.shell.info(`Stage 1: dropped ${dropped} entries${breakdown}.`);
  }

  // Append synthesized rows for user-dropped local structures (auto-detected
  // ligand; resolution from REMARK 2; no RCSB metadata or filters applied).
  rows.push(...localQcRows(pdbIds));

  // Build with EXPLICIT column types (not DG.DataFrame.fromObjects). fromObjects
  // infers each column's type from its values, and columns that are null for
  // every row (r_free, ligand_rscc — always null in v1) make that inference do
  // a `.first` on an empty set and throw "Bad state: No element". Constructing
  // the columns directly keeps the schema stable whether rows is empty or not.
  const numCol = (k: string): (number | null)[] =>
    rows.map((r) => (typeof r[k] === 'number' ? r[k] as number : null));
  const strCol = (k: string): string[] => rows.map((r) => String(r[k] ?? ''));
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('pdb_id', strCol('pdb_id')),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'resolution', numCol('resolution')),
    DG.Column.fromStrings('experimental_method', strCol('experimental_method')),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'r_free', numCol('r_free')),
    DG.Column.fromStrings('ligand_comp_id', strCol('ligand_comp_id')),
    DG.Column.fromStrings('ligand_chain', strCol('ligand_chain')),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'ligand_formula_weight',
      numCol('ligand_formula_weight')),
    DG.Column.fromStrings('ligand_smiles', strCol('ligand_smiles')),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'ligand_rscc', numCol('ligand_rscc')),
  ]);
  const smilesCol = df.col('ligand_smiles');
  if (smilesCol) smilesCol.semType = DG.SEMTYPE.MOLECULE;
  return {accepted: df, dropped: droppedList};
}

/**
 * A ligand-similarity hit: a PDB id and its RCSB similarity score (0..1). */
export interface LigandHit {
  id: string;
  score: number;
}

/** Lightweight per-PDB metadata for the ligand-search selection dialog. */
export interface PdbSummary {
  id: string;
  title: string;
  resolution: number | null;
  method: string;
  ligands: {compId: string; name: string; smiles: string}[];
}

/**
 * Search-by-ligand: find PDB entries whose bound ligand is similar to a query
 * SMILES, via the RCSB Search API's `chemical` service (Tanimoto
 * fingerprint-similarity). Returns experimental-structure hits ordered by
 * similarity (best first), de-duped, capped at `maxResults`.
 *
 * Empty result is a normal outcome (RCSB replies HTTP 204 No Content). An
 * invalid SMILES surfaces as an HTTP 400 from the API.
 */
export async function searchPdbsByLigand(
  smiles: string, maxResults = 25,
): Promise<LigandHit[]> {
  const query = {
    query: {
      type: 'terminal',
      service: 'chemical',
      parameters: {
        value: smiles,
        type: 'descriptor',
        descriptor_type: 'SMILES',
        match_type: 'fingerprint-similarity',
      },
    },
    return_type: 'entry',
    request_options: {
      paginate: {start: 0, rows: maxResults},
      results_content_type: ['experimental'],
    },
  };
  const resp = await grok.dapi.fetchProxy(SEARCH_ENDPOINT, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(query),
  });
  if (resp.status === 204) return []; // RCSB returns 204 when nothing matches
  if (!resp.ok) {
    let detail = '';
    try { detail = (await resp.text()).slice(0, 200); } catch { /* ignore */ }
    throw new Error(`RCSB Search API HTTP ${resp.status}${detail ? `: ${detail}` : ''}`);
  }
  const json = await resp.json();
  const seen = new Set<string>();
  const hits: LigandHit[] = [];
  for (const r of json.result_set ?? []) {
    const id = String(r.identifier ?? '').toUpperCase();
    if (!/^[0-9][A-Z0-9]{3}$/.test(id) || seen.has(id)) continue;
    seen.add(id);
    hits.push({id, score: typeof r.score === 'number' ? r.score : 0});
  }
  return hits;
}

/**
 * Fetch lightweight metadata (title, resolution, method, bound ligands) for a
 * set of PDB ids — used to populate the ligand-search selection dialog. Each
 * entry is fetched independently; ids that fail to resolve are omitted.
 */
/** Morgan fingerprint (bitstring) for a SMILES, or null if it won't parse. */
function morganFp(rdkit: {get_mol: (s: string) => any}, smiles: string): string | null {
  let mol: any = null;
  try {
    mol = rdkit.get_mol(smiles);
    if (!mol || !mol.is_valid?.()) return null;
    return typeof mol.get_morgan_fp === 'function' ? mol.get_morgan_fp() : null;
  } catch {
    return null;
  } finally {
    mol?.delete?.();
  }
}

/** Tanimoto coefficient between two equal-length fingerprint bitstrings. */
function tanimoto(a: string, b: string): number {
  const n = Math.min(a.length, b.length);
  let inter = 0;
  let union = 0;
  for (let i = 0; i < n; i++) {
    const x = a.charCodeAt(i) === 49; // '1'
    const y = b.charCodeAt(i) === 49;
    if (x && y) inter++;
    if (x || y) union++;
  }
  return union ? inter / union : 0;
}

/**
 * Fetch lightweight metadata (title, resolution, method, bound ligand) for a
 * set of PDB ids — used to populate the ligand-search selection dialog.
 *
 * When `querySmiles` is given, only the single bound ligand most chemically
 * similar to the query (Morgan/Tanimoto) is kept per entry — i.e. the ligand
 * that actually drove the RCSB hit — so ions, buffers and cryoprotectants
 * (IOD, MPD, …) that the cofactor deny-list misses don't clutter the row.
 * Without a query, cofactors are filtered by the deny-list as a fallback.
 */
export async function fetchPdbSummaries(
  ids: string[], querySmiles?: string,
): Promise<PdbSummary[]> {
  const cofactors = await loadCofactorDenylist();
  const rdkit = await getRdKitModule();
  const queryFp = querySmiles ? morganFp(rdkit, querySmiles) : null;

  const out = await Promise.all(ids.map(async (id): Promise<PdbSummary | null> => {
    try {
      const e = await executeEntryInfo(id);
      if (!e) return null;
      const allLigands = (e.nonpolymer_entities ?? [])
        .map((ne) => ({
          compId: (ne.nonpolymer_comp?.chem_comp?.id ?? '').toUpperCase(),
          name: ne.nonpolymer_comp?.chem_comp?.name ?? '',
          smiles: ne.nonpolymer_comp?.rcsb_chem_comp_descriptor?.SMILES_stereo ??
                  ne.nonpolymer_comp?.rcsb_chem_comp_descriptor?.SMILES ?? '',
        }))
        .filter((l) => l.compId);

      let ligands = allLigands;
      if (queryFp) {
        let best: typeof allLigands[number] | null = null;
        let bestSim = -1;
        for (const l of allLigands) {
          if (!l.smiles) continue;
          const fp = morganFp(rdkit, l.smiles);
          if (!fp) continue;
          const sim = tanimoto(queryFp, fp);
          if (sim > bestSim) { bestSim = sim; best = l; }
        }
        // Keep just the matched ligand; fall back to the deny-list filter if
        // none of the ligands produced a fingerprint.
        ligands = best ? [best] : allLigands.filter((l) => !cofactors.has(l.compId));
      } else {
        ligands = allLigands.filter((l) => !cofactors.has(l.compId));
      }

      return {
        id: id.toUpperCase(),
        title: e.struct?.title ?? '',
        resolution: e.rcsb_entry_info?.resolution_combined?.[0] ?? null,
        method: e.rcsb_entry_info?.experimental_method ?? '',
        ligands,
      };
    } catch {
      return null;
    }
  }));
  return out.filter((s): s is PdbSummary => s !== null);
}
