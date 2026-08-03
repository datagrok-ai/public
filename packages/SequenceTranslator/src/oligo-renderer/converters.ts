/**
 * Column converters for the OligoNucleotide semantic type.
 *
 * - `tagAsOligoNucleotide(col)` — flags an existing HELM column so the
 *   platform's renderer registry picks up the oligo cell renderer.
 * - `convertHelmColumnToOligo(table, helmCol)` — creates a NEW column
 *   with the same HELM values but tagged as OligoNucleotide. Adds it to
 *   the table next to the source column. Original column is left intact.
 * - `combineSenseAntisenseToOligo(table, senseCol, antiCol)` — produces
 *   an OligoNucleotide column from two existing single-strand HELM columns.
 *
 * No autodetection: the user opts in by calling one of these.
 */

import * as DG from 'datagrok-api/dg';

import {looksLikeHelm} from './helm-parser';
import {OLIGO_SEM_TYPE, OLIGO_UNITS} from './types';

/** Set the OligoNucleotide semType + unit + cell-renderer hint on a column. */
export function tagAsOligoNucleotide(col: DG.Column<string>): DG.Column<string> {
  col.semType = OLIGO_SEM_TYPE;
  col.meta.units = OLIGO_UNITS;
  // Make the renderer choice explicit even if the platform's autodiscovery
  // hasn't run yet for this session.
  col.setTag('cell.renderer', OLIGO_SEM_TYPE);
  col.setTag('quality', OLIGO_SEM_TYPE);
  return col;
}

/** Clone a HELM column and tag the clone as OligoNucleotide. */
export function convertHelmColumnToOligo(
  table: DG.DataFrame, helmCol: DG.Column<string>, newName?: string,
): DG.Column<string> {
  const values: string[] = new Array(helmCol.length);
  for (let i = 0; i < helmCol.length; i++)
    values[i] = helmCol.get(i) ?? '';
  const out = DG.Column.fromStrings(table.columns.getUnusedName(newName ?? `${helmCol.name} (oligo)`), values);
  tagAsOligoNucleotide(out);
  table.columns.add(out);
  return out;
}

/** Build an OligoNucleotide column from separate sense + antisense HELM columns.
 * Each cell becomes `RNA1{...}|RNA2{...}` with the original chain bodies.
 * If a row's value isn't HELM, that cell is left empty. */
export function combineSenseAntisenseToOligo(
  table: DG.DataFrame, senseCol: DG.Column<string>, antiCol: DG.Column<string>, newName?: string,
): DG.Column<string> {
  if (senseCol.length !== antiCol.length)
    throw new Error('Sense and antisense columns must have the same length');

  const values: string[] = new Array(senseCol.length);
  for (let i = 0; i < senseCol.length; i++) {
    const s = senseCol.get(i) ?? '';
    const a = antiCol.get(i) ?? '';
    if (!looksLikeHelm(s)) { values[i] = ''; continue; }
    const senseChain = renumberChain(s, 1);
    const antiChain = looksLikeHelm(a) ? renumberChain(a, 2) : '';
    const polymerSection = antiChain ? `${senseChain}|${antiChain}` : senseChain;
    values[i] = `${polymerSection}$$$$`;
  }

  const out = DG.Column.fromStrings(
    table.columns.getUnusedName(newName ?? `${senseCol.name}+${antiCol.name} (oligo)`),
    values,
  );
  tagAsOligoNucleotide(out);
  table.columns.add(out);
  return out;
}

/** Strip the polymer-section index off `RNA1{...}` / `RNA2{...}` and rewrite to `RNA<n>{...}`. */
function renumberChain(helm: string, n: number): string {
  const polymerSection = helm.split('$')[0];
  return polymerSection.replace(/^(RNA|DNA|PEPTIDE|CHEM|BLOB)\d+\{/, `$1${n}{`);
}
