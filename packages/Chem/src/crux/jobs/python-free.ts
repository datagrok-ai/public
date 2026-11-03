/**
 * The Crux worker's jobs for Chem's Python and server code in Crux mode (python-free.ts sends them, crux.worker.ts runs
 * them; crux-core decision chem-python-free): chem-wasm's replacements of Curate, Mutate, Generate Conformers, Murcko
 * Scaffolds (scripts/stand-alone/murcko_scaffold.py), Filter by Catalogs (filter_catalogs.py) and Map Identifiers'
 * chemistry (the chem-chem container's conversions), and the column work of Butina clustering (the script's Morgan
 * fingerprints, the neighbours within the distance cutoff). The Gasteiger panel's one molecule is drawn on the main
 * thread's module; Calculate logP is the descriptor job's MolLogP.
 */
import type {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {WorkerJob, throwIfTrap} from '../crux-api';

/** butina_cluster.py's fingerprint: `GetMorganFingerprintAsBitVect(mol, 2, 1024)`. */
const BUTINA_FP_BITS = 1024;
const BUTINA_FP_DETAILS = JSON.stringify({radius: 2, nBits: BUTINA_FP_BITS});
/** A fingerprint's 32-bit words. */
export const BUTINA_FP_WORDS = BUTINA_FP_BITS / 32;

/** The rows' Butina neighbours: each row's count (itself included) and, row after row, their indices. */
export interface ButinaNeighbours {counts: Int32Array, neighbours: Int32Array}

/**
 * butina_cluster.py's fingerprints of `molecules`, BUTINA_FP_WORDS words each: a molecule read as the script reads it
 * (`MolFromSmiles` / `MolFromMolBlock`, RDKit's defaults: RDKit.js's `get_mol`), no bits where it reads none (the
 * script's `Chem.Mol()`) or the cell is empty.
 */
function butinaFingerprints(module: RDModule, molecules: (string | null)[]): Uint32Array {
  const fps = new Uint32Array(molecules.length * BUTINA_FP_WORDS);
  const bytes = new Uint8Array(fps.buffer);
  for (let i = 0; i < molecules.length; i++) {
    let mol: RDMol | null = null;
    try {
      mol = molecules[i] ? module.get_mol(molecules[i]!) : null;
      if (mol?.is_valid())
        bytes.set(mol.get_morgan_fp_as_uint8array(BUTINA_FP_DETAILS), i * BUTINA_FP_WORDS * 4);
    } catch (e) {
      throwIfTrap(e);
    } finally {
      mol?.delete();
    }
  }
  return fps;
}

const popcount = (v: number): number => {
  v -= (v >>> 1) & 0x55555555;
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return (((v + (v >>> 4)) & 0x0F0F0F0F) * 0x01010101) >>> 24;
};

/**
 * Each of `rows`' neighbours within `cutoff` (Tanimoto distance, `1 - similarity`), in index order and itself included,
 * as RDKit's `Butina.ClusterData(dists, n, cutoff, isDistData=True)` lists them from butina_cluster.py's distances
 * (`1 - BulkTanimotoSimilarity`; two fingerprints without bits: similarity 0, RDKit's; a row is its own neighbour).
 */
function butinaNeighbours(fps: Uint32Array, rows: number[], cutoff: number): ButinaNeighbours {
  const n = fps.length / BUTINA_FP_WORDS;
  const ones = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    let c = 0;
    for (let w = i * BUTINA_FP_WORDS; w < (i + 1) * BUTINA_FP_WORDS; w++)
      c += popcount(fps[w]);
    ones[i] = c;
  }
  const counts = new Int32Array(rows.length);
  const neighbours: number[] = [];
  rows.forEach((i, k) => {
    const start = neighbours.length;
    for (let j = 0; j < n; j++) {
      if (j === i) {
        neighbours.push(j);
        continue;
      }
      const total = ones[i] + ones[j];
      let common = 0;
      if (total > 0) {
        for (let w = 0; w < BUTINA_FP_WORDS; w++)
          common += popcount(fps[i * BUTINA_FP_WORDS + w] & fps[j * BUTINA_FP_WORDS + w]);
      }
      // Chem's tanimotoSimilarity (common / (total - common)) but for two empty fingerprints: RDKit's 0
      const similarity = total === 0 ? 0 : common / (total - common);
      if (1 - similarity <= cutoff)
        neighbours.push(j);
    }
    counts[k] = neighbours.length - start;
  });
  return {counts, neighbours: Int32Array.from(neighbours)};
}

/** The jobs by name (the Crux worker's `op`). */
export const SCRIPT_JOBS: {[op: string]: WorkerJob} = {
  // scripts/curate.py's `curated_molecule` per molecule ('' where it writes '', null where it raises)
  chemCurate: ({crux}, molecules: (string | null)[], options: string): (string | null)[] =>
    JSON.parse(crux.chemCurate(molecules, options)),
  // scripts/mutate.py's `mutations` (a throw where the script raises)
  chemMutate: ({crux}, molecule: string, options: string): string[] => JSON.parse(crux.chemMutate(molecule, options)),
  // scripts/generate_conformers.py's table: {columns, rows} (a throw where it raises or the deadline passes)
  chemConformers: ({crux}, molecule: string, options: string): unknown =>
    JSON.parse(crux.chemConformers(molecule, options)),
  butinaFingerprints: ({module}, molecules: (string | null)[]): Uint32Array => butinaFingerprints(module, molecules),
  butinaNeighbours: (_ctx, rows: number[], fps: Uint32Array, cutoff: number): ButinaNeighbours =>
    butinaNeighbours(fps, rows, cutoff),
  // murcko_scaffold.py's `scaffolds` per molecule ('' without rings, null where the molecule does not read)
  chemMurckoScaffolds: ({crux}, molecules: (string | null)[]): (string | null)[] =>
    JSON.parse(crux.chemMurckoScaffolds(molecules)),
  // filter_catalogs.py's flag per molecule (1 where the catalog matches); a catalog it has no table of throws
  chemFilterCatalog: ({crux}, molecules: (string | null)[], catalog: string): Uint8Array =>
    crux.chemFilterCatalog(molecules, catalog),
  // the container's molecules_to_canonical: a SMILES per molecule, '' where it does not read
  chemCanonicalSmiles: ({crux}, molecules: (string | null)[]): string[] =>
    JSON.parse(crux.chemCanonicalSmiles(molecules)),
  // the container's molecules_to_inchi up to the InChI library: an input per molecule, null where it writes ''
  chemInchiInputs: ({crux}, molecules: (string | null)[]): (Uint8Array | null)[] =>
    crux.chemInchiInputs(molecules),
};
