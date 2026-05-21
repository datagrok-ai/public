/*
 * Stage 5b — render consensus_model to a fake-PDB block.
 * Each consensus feature becomes a HETATM (chain P); B-factor = frequency, occupancy = 1.00.
 * Mol*'s default CPK coloring separates families by element (N, O, S, C, NA, CL, BR).
 *
 * Optional reference-PDB overlay - ATOM/HETATM lines from the reference are renumbered
 * to follow the consensus HETATMs so serials stay unique in the concatenated block.
 *
 * Port of files/_reference/render_consensus_pharmacophore_3d_original.py.
 * Blueprint reference - Phase 5b.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {resolveFamily} from './family-map';

/** Format one fixed-width PDB HETATM record per the PDB v3.3 spec. */
function hetatm(
  serial: number, elem: string, resName: string, chain: string,
  resSeq: number, x: number, y: number, z: number, occ: number, bFactor: number,
): string {
  // Cols 13-16 (atom name): 2-letter elements occupy 13-14; 1-letter elements right-shift one space.
  const atomName = elem.length === 2 ? elem.padEnd(4) : ' ' + elem.padEnd(3);
  const f = (v: number, width: number, dp: number) => v.toFixed(dp).padStart(width);
  return (
    'HETATM' + serial.toString().padStart(5) + ' ' +
    atomName + ' ' +
    resName.padEnd(3) + ' ' + chain + resSeq.toString().padStart(4) + '    ' +
    f(x, 8, 3) + f(y, 8, 3) + f(z, 8, 3) +
    f(occ, 6, 2) + f(bFactor, 6, 2) +
    '          ' + elem.padStart(2)
  );
}

/** Renumber ATOM/HETATM serials in a reference PDB block so they start at `startSerial`. */
function renumberReferenceSerials(refPdbBlock: string, startSerial: number): string {
  let n = startSerial;
  return refPdbBlock.split(/\r?\n/).map((line) => {
    if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      const newSerial = n.toString().padStart(5);
      n += 1;
      return line.slice(0, 6) + newSerial + line.slice(11);
    }
    return line;
  }).join('\n');
}

/** Fetch a reference PDB block; return '' on any failure (caller treats overlay as optional). */
async function fetchReferencePdb(refPdbId: string): Promise<string> {
  try {
    const url = `https://files.rcsb.org/download/${refPdbId.toUpperCase()}.pdb`;
    const resp = await grok.dapi.fetchProxy(url);
    if (!resp.ok) return '';
    const text = await resp.text();
    const keep = ['ATOM', 'HETATM', 'TER', 'SSBOND', 'HELIX', 'SHEET', 'CONECT'];
    return text.split(/\r?\n/).filter((ln) => keep.some((k) => ln.startsWith(k))).join('\n');
  } catch {
    return '';
  }
}

/**
 * Build the consensus pharmacophore PDB block. `refPdbId` is optional; when supplied,
 * the reference structure is fetched, filtered to ATOM/HETATM/connectivity, renumbered,
 * and concatenated *before* the consensus HETATMs so the consensus pops on top in Mol*.
 */
export async function renderConsensusPharmacophore(
  consensus: DG.DataFrame, refPdbId?: string,
): Promise<string> {
  const familyCol = consensus.col('family');
  const xCol = consensus.col('x');
  const yCol = consensus.col('y');
  const zCol = consensus.col('z');
  if (!familyCol || !xCol || !yCol || !zCol)
    throw new Error('consensus_model requires columns: family, x, y, z');
  const freqCol = consensus.col('frequency');

  const lines: string[] = [];
  for (let i = 0; i < consensus.rowCount; i++) {
    const fam = resolveFamily(familyCol.get(i));
    const freq = freqCol ? Number(freqCol.get(i)) : 1.0;
    lines.push(hetatm(
      i + 1, fam.element, fam.resName, 'P', i + 1,
      Number(xCol.get(i)), Number(yCol.get(i)), Number(zCol.get(i)),
      1.0, freq,
    ));
  }
  const consensusBlock = lines.join('\n');
  const nConsensus = lines.length;

  let refBlock = '';
  const cleanedRef = (refPdbId ?? '').trim().toUpperCase();
  const sentinels = new Set(['', 'NONE', 'NA', 'NULL', '-']);
  if (cleanedRef && !sentinels.has(cleanedRef)) {
    const raw = await fetchReferencePdb(cleanedRef);
    if (raw) refBlock = renumberReferenceSerials(raw, nConsensus + 1) + '\n';
  }

  return (
    'REMARK   Consensus pharmacophore -- chain P, B-factor = frequency\n' +
    refBlock + consensusBlock + '\nEND\n'
  );
}
