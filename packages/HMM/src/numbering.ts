/** ANARCI numbering of a macromolecule column, shaped for Bio's numbering
 * dialog (`meta.role: antibodyNumbering`; see Bio's CLAUDE.md, "Adding a New
 * Antibody Numbering Engine"): `position_names`, `chain_type`,
 * `annotations_json`, `numbering_detail`, `numbering_map`, plus ANARCI's
 * chain, species, germline and E-value columns. */
import * as DG from 'datagrok-api/dg';
import type {AnarciOptions, SequenceResult} from './hmmer/anarci/anarci.ts';
import type {Scheme} from './hmmer/anarci/types.ts';
import {HmmerPool, type Progress} from './pool';
import {columnResidues} from './sequences';
import {type DomainChoice, type NumberingRow, numberingRow} from './numbering-rows';

export {CHAIN_NAMES, domainChoice, regionDefs, SCHEMES, type DomainChoice} from './numbering-rows';

/** The residues Bio's dialog indexes: gaps (`-`, `.`) removed, upper case. */
export function extractSequence(raw: string | null | undefined): string {
  return (raw ?? '').replace(/[-.]/g, '').toUpperCase();
}

/** Run ANARCI on a column in the worker pool; the workers return the frame's
 * rows (numbering indices refer to the gap-free residues, as Bio's numbering
 * dialog expects). */
export async function numberColumn(pool: HmmerPool, column: DG.Column<string>, options: AnarciOptions,
  scheme: Scheme, choice: DomainChoice = 'first', progress?: Progress):
  Promise<{sequences: string[]; rows: NumberingRow[]}> {
  const sequences = (await columnResidues(column)).residues;
  const rows = await pool.anarciRows(sequences.map((s, i) => [String(i), s]), options, scheme, choice, undefined,
    progress);
  return {sequences, rows};
}

/** The numbering DataFrame of `rows` (numberingRow): one row per input sequence. */
export function frameFromRows(rows: NumberingRow[]): DG.DataFrame {
  const n = rows.length;
  const string = (name: keyof NumberingRow) => {
    const col = DG.Column.fromType(DG.COLUMN_TYPE.STRING, name, n);
    for (let i = 0; i < n; i++) {
      const value = rows[i][name] as string;
      if (value) col.set(i, value, false);
    }
    return col;
  };
  // E-values reach 1e-60 and identities are exact ratios, so floats are doubles.
  const double = (name: keyof NumberingRow) => {
    const values = new Float64Array(n).fill(DG.FLOAT_NULL);
    for (let i = 0; i < n; i++) {
      const value = rows[i][name] as number | null;
      if (value !== null) values[i] = value;
    }
    return DG.Column.fromFloat64Array(name, values);
  };
  const domains = DG.Column.fromInt32Array('domains', Int32Array.from(rows, (r) => r.domains));
  const cols = [string('position_names'), string('chain_type'), string('annotations_json'),
    string('numbering_detail'), string('numbering_map'), string('chain'), string('species'), string('v_gene'),
    string('j_gene'), double('v_identity'), double('j_identity'), double('evalue'), double('bitscore'), domains,
    string('error')];
  cols[11].meta.format = 'scientific';
  for (const col of [cols[9], cols[10]]) col.meta.format = '0.00';
  return DG.DataFrame.fromColumns(cols);
}

/** The numbering DataFrame of full ANARCI results (computed here; see numberColumn
 * for the worker path): the N-terminal numbered domain of each sequence, or its
 * first heavy-like or light-like domain (`choice`); `domains` counts all of them. */
export function numberingFrame(scheme: Scheme, results: SequenceResult[], choice: DomainChoice = 'first'):
  DG.DataFrame {
  return frameFromRows(results.map((r) => numberingRow(scheme, r, choice)));
}

/** ANARCI's chain types per scheme (as its command line: TCRs only in IMGT and AHo). */
export function allowedChains(scheme: Scheme): string[] {
  return scheme === 'imgt' || scheme === 'aho' ? ['H', 'K', 'L', 'A', 'B', 'G', 'D'] : ['H', 'K', 'L'];
}

/** Species choices: ANARCI's default (human and mouse preferred), all, or one species. */
// No commas in choices: Datagrok's function headers split choice lists at commas.
export const SPECIES_CHOICES = ['human and mouse', 'all', 'human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig',
  'alpaca', 'cow'];

export function speciesOption(choice: string | null | undefined): string[] | null {
  if (!choice || choice === 'human and mouse' || choice === 'human, mouse') return ['human', 'mouse'];
  if (choice === 'all') return null;
  return [choice];
}
