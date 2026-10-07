/** Residues of a Macromolecule column as HMMER compares them, from any Bio
 * notation: one upper-case letter per non-gap monomer (non-canonical monomers
 * become `X`, HMMER's unknown residue), plus each residue's monomer position
 * in its cell, which is what Bio's annotations index. */
import * as DG from 'datagrok-api/dg';

/** Bio's sequence handler, loaded only for non-FASTA columns (keeps package.js small). */
async function seqHandler(column: DG.Column<string>) {
  const {getSeqHelper} = await import('@datagrok-libraries/bio/src/utils/seq-helper');
  return (await getSeqHelper()).getSeqHandler(column);
}

export interface ColumnResidues {
  residues: string[];
  /** `positions[row][k]`: monomer position of residue `k` in the cell. */
  positions: number[][];
}

/** Plain FASTA (single-letter monomers): residues are the non-gap characters. */
function plain(raw: string): {residues: string; positions: number[]} {
  let residues = '';
  const positions: number[] = [];
  for (let i = 0; i < raw.length; i++) {
    const c = raw[i];
    if (c === '-' || c === '.') continue;
    residues += c.toUpperCase();
    positions.push(i);
  }
  return {residues, positions};
}

export async function columnResidues(column: DG.Column<string>): Promise<ColumnResidues> {
  const n = column.length;
  const residues = new Array<string>(n);
  const positions = new Array<number[]>(n);
  const units = column.meta.units;
  const values = column.toList() as (string | null)[];
  if (!units || units === 'fasta' && !values.some((v) => v?.includes('['))) {
    for (let i = 0; i < n; i++) {
      const row = plain(values[i] ?? '');
      residues[i] = row.residues;
      positions[i] = row.positions;
    }
    return {residues, positions};
  }
  const handler = await seqHandler(column);
  for (let i = 0; i < n; i++) {
    let text = '';
    const at: number[] = [];
    if (values[i]) {
      const split = handler.getSplitted(i);
      for (let p = 0; p < split.length; p++) {
        if (split.isGap(p)) continue;
        const monomer = split.getCanonical(p);
        text += monomer.length === 1 ? monomer.toUpperCase() : 'X';
        at.push(p);
      }
    }
    residues[i] = text;
    positions[i] = at;
  }
  return {residues, positions};
}

/** Residues of one cell (see {@link columnResidues}). */
export async function cellResidues(column: DG.Column<string>, row: number):
  Promise<{residues: string; positions: number[]}> {
  const value = column.get(row) ?? '';
  const units = column.meta.units;
  if (!units || units === 'fasta' && !value.includes('[')) return plain(value);
  const split = (await seqHandler(column)).getSplitted(row);
  let residues = '';
  const positions: number[] = [];
  for (let p = 0; p < split.length; p++) {
    if (split.isGap(p)) continue;
    const monomer = split.getCanonical(p);
    residues += monomer.length === 1 ? monomer.toUpperCase() : 'X';
    positions.push(p);
  }
  return {residues, positions};
}
