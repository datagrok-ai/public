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

const GRAPHQL_ENDPOINT = 'https://data.rcsb.org/graphql';

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
export async function enrichPdbList(
  pdbIds: string[], options: PipelineOptions,
): Promise<DG.DataFrame> {
  const cofactors = await loadCofactorDenylist();
  const rdkit = await getRdKitModule();

  const rows: Record<string, any>[] = [];
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
  const upper = pdbIds.map((s) => s.toUpperCase());
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
      dropped++;
      continue;
    }
    if (!entry) { dropped++; continue; }

    const resolution = entry.rcsb_entry_info?.resolution_combined?.[0];
    const method = entry.rcsb_entry_info?.experimental_method ?? '';
    if (typeof resolution === 'number' && resolution > options.maxResolution) {
      droppedResolution++; dropped++; continue;
    }
    // RCSB returns 'X-ray' (sentence-cased, hyphenated) for X-ray diffraction
    // entries via this field. Substring match also covers legacy 'X-RAY DIFFRACTION'.
    const isXray = method.toLowerCase().includes('x-ray');
    if (options.requireXray && !isXray) { droppedXray++; dropped++; continue; }

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
    if (ligandsEmitted === 0) { droppedNoLigand++; dropped++; }
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

  // Always build the DataFrame from a schema-stable column set so empty results
  // still have the expected columns downstream.
  const df = rows.length
    ? DG.DataFrame.fromObjects(rows)!
    : DG.DataFrame.fromColumns([
      DG.Column.fromStrings('pdb_id', []),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'resolution', []),
      DG.Column.fromStrings('experimental_method', []),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'r_free', []),
      DG.Column.fromStrings('ligand_comp_id', []),
      DG.Column.fromStrings('ligand_chain', []),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'ligand_formula_weight', []),
      DG.Column.fromStrings('ligand_smiles', []),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'ligand_rscc', []),
    ]);
  const smilesCol = df.col('ligand_smiles');
  if (smilesCol) smilesCol.semType = DG.SEMTYPE.MOLECULE;
  return df;
}
