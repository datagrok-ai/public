/*
 * Small PDB-block helpers used by the orchestrator's Preview button and
 * (eventually) by the reference-overlay path in renderer.ts.
 *
 *  fetchPdbBlock(id)            -> raw PDB text from RCSB
 *  applyPdbTransform(block, T)  -> rewrites ATOM/HETATM coords with a 4x4 transform (row-major,
 *                                  identity is a no-op)
 *  concatPdbStructures(blocks)  -> stacks N PDBs into one renderable block with globally
 *                                  unique atom serials; ATOM/HETATM/TER preserved, CONECT
 *                                  dropped (CONECT references the original serials and would
 *                                  be invalid after renumbering)
 *
 * Blueprint reference - section 1 reuse inventory ("transform.ts"), Phase 6.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {resolveFamily} from './family-map';

/** PDB column widths from the v3.3 spec. */
const ATOM_PREFIXES = ['ATOM', 'HETATM'];
const COORD_RECORD_PREFIXES = [...ATOM_PREFIXES, 'TER'];

/**
 * 3-letter residue names of HETATM groups that should NEVER appear in the Mol*
 * concat block — waters, ions, cryoprotectants, common buffers, cofactors. Kept
 * in sync with `files/cofactor-denylist.csv` (the CSV is the doc source, this
 * set is what TS-side code uses for filtering).
 *
 * Stripping them at PDB level (vs. setting alpha=0 in Mol*) is necessary because
 * Mol*'s default partition lumps EVERY HETATM that isn't water/ion/non-standard
 * into a single "Ligand" component. GOL/MES/EDO would land next to the drug
 * ligand and we can't selectively hide them post-load via setOptions.
 */
const COFACTOR_DENYLIST = new Set([
  'HOH', 'DOD',
  'SO4', 'PO4',
  'CL', 'NA', 'K', 'CA', 'MG', 'ZN', 'FE', 'MN', 'NI', 'CU',
  'EDO', 'GOL', 'PEG', 'PG4', 'PGE',
  'DMS', 'ACT', 'TRS', 'MES', 'IMD',
  'HEM', 'HEC',
  'NAG', 'BMA', 'MAN', 'FUC', 'GAL',
  'FAD', 'NAD', 'NAP', 'ADP', 'ATP', 'GDP', 'GTP', 'COA', 'FMN', 'SAM', 'SAH',
]);

/** Strip HETATM lines whose residue name (cols 18-20) is in the cofactor deny-list. */
export function stripCofactorsFromPdb(pdb: string): string {
  if (!pdb) return pdb;
  const out: string[] = [];
  for (const line of pdb.split(/\r?\n/)) {
    if (line.startsWith('HETATM')) {
      const resName = line.slice(17, 20).trim().toUpperCase();
      if (COFACTOR_DENYLIST.has(resName)) continue;
    }
    out.push(line);
  }
  return out.join('\n');
}

/** Identity 4x4 in row-major flattening — what stubStage2a emits. */
export const IDENTITY_4X4: number[] = [
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1,
];

/**
 * Fetch a raw PDB block from RCSB. Returns '' on any failure (caller decides
 * whether to bubble or fall back).
 */
export async function fetchPdbBlock(pdbId: string): Promise<string> {
  try {
    const url = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb`;
    const resp = await grok.dapi.fetchProxy(url);
    if (!resp.ok) return '';
    return await resp.text();
  } catch {
    return '';
  }
}

function parseTransformJson(json: string): number[] {
  if (!json) return IDENTITY_4X4;
  try {
    const t = JSON.parse(json);
    if (Array.isArray(t) && t.length === 16) return t.map((v: any) => Number(v));
  } catch { /* fall through */ }
  return IDENTITY_4X4;
}

function isIdentity(t: number[]): boolean {
  for (let i = 0; i < 16; i++) {
    const expected = (i % 5 === 0) ? 1 : 0; // diagonal entries at indices 0/5/10/15
    if (Math.abs(t[i] - expected) > 1e-9) return false;
  }
  return true;
}

/** Apply a row-major 4x4 affine transform to one ATOM/HETATM line; returns the rewritten line. */
function transformAtomLine(line: string, t: number[]): string {
  // PDB columns 31-38 / 39-46 / 47-54 (1-based, inclusive) = x/y/z, each %8.3f.
  const x = parseFloat(line.slice(30, 38));
  const y = parseFloat(line.slice(38, 46));
  const z = parseFloat(line.slice(46, 54));
  if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) return line;
  const nx = t[0] * x + t[1] * y + t[2] * z + t[3];
  const ny = t[4] * x + t[5] * y + t[6] * z + t[7];
  const nz = t[8] * x + t[9] * y + t[10] * z + t[11];
  const f = (v: number) => v.toFixed(3).padStart(8);
  return line.slice(0, 30) + f(nx) + f(ny) + f(nz) + line.slice(54);
}

/**
 * Rewrites the x/y/z columns of every ATOM/HETATM record by applying the row-major
 * 4x4 transform. The transform JSON comes from Stage 2a/2b's `transform_4x4_json` column.
 * Identity transforms short-circuit (no string churn).
 */
export function applyPdbTransform(pdbBlock: string, transformJson: string): string {
  const t = parseTransformJson(transformJson);
  if (isIdentity(t)) return pdbBlock;
  return pdbBlock.split(/\r?\n/).map((line) => {
    return ATOM_PREFIXES.some((p) => line.startsWith(p)) ? transformAtomLine(line, t) : line;
  }).join('\n');
}

/** Rewrite the 5-char serial field (cols 7-11) of an ATOM/HETATM/TER line. */
function renumberSerial(line: string, newSerial: number): string {
  return line.slice(0, 6) + newSerial.toString().padStart(5) + line.slice(11);
}

/** Rewrite the 1-char chainID field (col 22) of an ATOM/HETATM/TER line. */
function rewriteChainId(line: string, newChain: string): string {
  return line.slice(0, 21) + newChain.charAt(0) + line.slice(22);
}

// Pool of single-character chain IDs Mol* will color distinctly. Uppercase letters
// first (standard PDB), then digits, then lowercase — gives ~60 unique chains before
// we'd wrap.
const CHAIN_LETTER_POOL =
  'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz';

/**
 * Stack multiple PDB blocks into a single block with globally unique atom serials
 * AND globally unique chain IDs.
 *
 * Each input block contributes its ATOM/HETATM/TER records (with serials renumbered
 * and chain IDs remapped) surrounded by a tagging REMARK so downstream debugging is
 * obvious. CONECT records are dropped because they reference the original serials.
 * Header records (HEADER/COMPND/etc.) are also dropped from the 2nd+ block — a single
 * concatenated PDB has one header.
 *
 * The chain remap is per-structure: each PDB's original chain IDs are mapped to fresh
 * letters from CHAIN_LETTER_POOL. Mol* defaults to per-chain coloring, so this makes
 * each input structure render in a different color — easy visual differentiation.
 *
 * If labels.length matches blocks.length the labels are used in the per-structure REMARK.
 */
export interface ConcatResult {
  /** Assembled PDB body. Does NOT include trailing END if `terminate: false`. */
  body: string;
  /** Next atom serial number the caller can use for an appended block. */
  nextSerial: number;
}

export function concatPdbStructures(
  blocks: string[],
  labels?: string[],
  options: {terminate?: boolean; startSerial?: number} = {},
): ConcatResult {
  const terminate = options.terminate ?? true;
  const out: string[] = ['REMARK   Consensus Pharmacophore: stacked input structures'];
  let serial = options.startSerial ?? 1;
  let nextChainIdx = 0;
  blocks.forEach((block, i) => {
    if (!block) return;
    const label = labels && labels[i] ? labels[i] : `structure-${i + 1}`;
    out.push(`REMARK   --- ${label} ---`);

    // Each input PDB gets its own chain-id remap. Original chain 'A' in structure 1
    // and original chain 'A' in structure 2 land on DIFFERENT letters in the output.
    const chainMap = new Map<string, string>();
    const remapChain = (originalChain: string): string => {
      if (!chainMap.has(originalChain)) {
        chainMap.set(originalChain,
          CHAIN_LETTER_POOL.charAt(nextChainIdx % CHAIN_LETTER_POOL.length));
        nextChainIdx += 1;
      }
      return chainMap.get(originalChain)!;
    };

    for (const line of block.split(/\r?\n/)) {
      if (!COORD_RECORD_PREFIXES.some((p) => line.startsWith(p))) continue;
      const originalChain = line.length > 21 ? line.charAt(21) : ' ';
      const newChain = remapChain(originalChain);
      const renumbered = renumberSerial(line, serial);
      out.push(rewriteChainId(renumbered, newChain));
      serial += 1;
    }
  });
  if (terminate) out.push('END', '');
  return {body: out.join('\n'), nextSerial: serial};
}

/**
 * Convert a pocket_atoms DataFrame (Stage 3 output) into a HETATM block suitable
 * for overlay in Mol*. By default emits only Cα and ligand-seed atoms (cleaner
 * visual; the pocket cloud otherwise has thousands of side-chain points).
 *
 * All output atoms land on a single distinguished chain ID ('P' by default) so
 * Mol* colors them uniformly and distinctly from the per-structure chains.
 */
export function pocketAtomsToOverlayBlock(
  pocketAtoms: DG.DataFrame,
  options: {chain?: string; resName?: string; cAlphaOnly?: boolean; serialStart?: number} = {},
): string {
  const chain = (options.chain ?? 'P').charAt(0);
  const resName = (options.resName ?? 'POC').padEnd(3).slice(0, 3);
  const cAlphaOnly = options.cAlphaOnly ?? true;
  let serial = options.serialStart ?? 1;

  const cols = {
    pdb: pocketAtoms.col('pdb_id'),
    atomName: pocketAtoms.col('atom_name'),
    element: pocketAtoms.col('element'),
    x: pocketAtoms.col('x_consensus'),
    y: pocketAtoms.col('y_consensus'),
    z: pocketAtoms.col('z_consensus'),
    seed: pocketAtoms.col('ligand_seed'),
  };
  if (!cols.atomName || !cols.x || !cols.y || !cols.z) return '';

  const lines: string[] = ['REMARK   Pocket overlay (Stage 3)'];
  for (let i = 0; i < pocketAtoms.rowCount; i++) {
    const isSeed = cols.seed ? Boolean(cols.seed.get(i)) : false;
    const atomName = String(cols.atomName.get(i) ?? '');
    if (cAlphaOnly && !isSeed && atomName !== 'CA') continue;

    const x = Number(cols.x.get(i));
    const y = Number(cols.y.get(i));
    const z = Number(cols.z.get(i));
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) continue;

    const elem = (String(cols.element?.get(i) ?? atomName.charAt(0))).trim() || 'C';
    const atomNameField = elem.length === 2 ? elem.padEnd(4) : ' ' + elem.padEnd(3);
    const f = (v: number) => v.toFixed(3).padStart(8);
    lines.push(
      'HETATM' + serial.toString().padStart(5) + ' ' +
      atomNameField + ' ' +
      resName + ' ' + chain + (serial % 9999).toString().padStart(4) + '    ' +
      f(x) + f(y) + f(z) +
      '  1.00  1.00          ' + elem.padStart(2),
    );
    serial += 1;
  }
  return lines.join('\n');
}

/**
 * Convert a `ligand_features` DataFrame (Stage 4 output) into a HETATM block
 * suitable for overlay in Mol*. Each feature centroid becomes one HETATM whose
 * element symbol comes from the family-map (Donor→N, Acceptor→O, Aromatic→S,
 * Hydrophobic→C, Positive→NA, Negative→CL, Halogen→BR) — Mol*'s default CPK
 * coloring then gives a distinct color per family without any extra config.
 *
 * Diagnostic rows (skip_reason non-empty) are dropped silently.
 */
export function ligandFeaturesToOverlayBlock(
  features: DG.DataFrame,
  options: {chain?: string; serialStart?: number} = {},
): string {
  const chain = (options.chain ?? 'F').charAt(0);
  let serial = options.serialStart ?? 1;

  const cols = {
    family: features.col('family'),
    x: features.col('x'),
    y: features.col('y'),
    z: features.col('z'),
    skip: features.col('skip_reason'),
  };
  if (!cols.family || !cols.x || !cols.y || !cols.z) return '';

  const lines: string[] = ['REMARK   Ligand features overlay (Stage 4)'];
  for (let i = 0; i < features.rowCount; i++) {
    if (cols.skip && String(cols.skip.get(i) ?? '').trim() !== '') continue;
    const fam = resolveFamily(String(cols.family.get(i) ?? ''));
    const x = Number(cols.x.get(i));
    const y = Number(cols.y.get(i));
    const z = Number(cols.z.get(i));
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) continue;

    const elem = fam.element;
    const atomNameField = elem.length === 2 ? elem.padEnd(4) : ' ' + elem.padEnd(3);
    const f = (v: number) => v.toFixed(3).padStart(8);
    lines.push(
      'HETATM' + serial.toString().padStart(5) + ' ' +
      atomNameField + ' ' +
      fam.resName + ' ' + chain + (serial % 9999).toString().padStart(4) + '    ' +
      f(x) + f(y) + f(z) +
      '  1.00  1.00          ' + elem.padStart(2),
    );
    serial += 1;
  }
  return lines.join('\n');
}
