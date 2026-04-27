/**
 * SMILES↔Molecule3D column link helpers for the Chem atom-picker bridge.
 * The link is a reciprocal pair of persistent column tags:
 *   smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL]        = mol3DColName
 *   mol3DCol.tags [CHEM_ATOM_PICKER_LINKED_SMILES_COL] = smilesColName
 * `delete col.tags[KEY]` routes through the Datagrok MapProxy to grok_Map_Delete.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';

/** Mirrors CHEM_ATOM_PICKER_LINKED_COL from Chem's constants.ts.
 *  Duplicated to avoid a cross-package import; keep string values in sync. */
export const CHEM_ATOM_PICKER_LINKED_COL = '%chem-atom-picker-linked-col';

/** Mirrors CHEM_ATOM_PICKER_LINKED_SMILES_COL from Chem's constants.ts. */
export const CHEM_ATOM_PICKER_LINKED_SMILES_COL = '%chem-atom-picker-linked-smiles-col';

/** Mirrors CHEM_INTERACTIVE_SELECTION_EVENT from Chem's constants.ts. */
const CHEM_SELECTION_EVENT = 'chem-interactive-selection-changed';

/** Mirrors ChemTemps.SUBSTRUCT_PROVIDERS from @datagrok-libraries/chem-meta. */
const SUBSTRUCT_PROVIDERS_KEY = 'substruct-providers';

/** Returns the SMILES column name linked to the given Mol3D column, or null.
 *  O(1) via the reciprocal tag on the Mol3D column; falls back to O(n) scan
 *  over forward tags. */
export function findLinkedSmilesColName(
  df: DG.DataFrame, mol3DColName: string): string | null {
  const mol3DCol = df.col(mol3DColName);
  if (mol3DCol) {
    const reciprocal = mol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL];
    if (reciprocal && df.col(reciprocal)) return reciprocal;
  }
  const c = df.columns.toList().find(
    (col) => col.tags[CHEM_ATOM_PICKER_LINKED_COL] === mol3DColName);
  return c?.name ?? null;
}

/** Removes all picker links that point at `mol3DColName` and clears residual highlights. */
export function clearLinksToMol3D(
  df: DG.DataFrame, mol3DColName: string): void {
  for (const c of df.columns.toList()) {
    if (c.tags[CHEM_ATOM_PICKER_LINKED_COL] === mol3DColName) {
      delete c.tags[CHEM_ATOM_PICKER_LINKED_COL];
      clearAtomPickerHighlights(c);
    }
  }
  const mol3DCol = df.col(mol3DColName);
  if (mol3DCol?.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL])
    delete mol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL];
}

/** Writes or clears the reciprocal-tag link between `smilesCol` and `mol3DColName`.
 *  Pass null to unlink (clears both forward and reverse tags + highlights). */
export function setSmilesColLink(
  smilesCol: DG.Column, mol3DColName: string | null): void {
  const df = smilesCol.dataFrame;
  const prevMol3DName = smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL];
  const prevMol3DCol = prevMol3DName && df ? df.col(prevMol3DName) : null;
  if (prevMol3DCol && prevMol3DCol.name !== mol3DColName &&
      prevMol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL])
    delete prevMol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL];

  if (!mol3DColName) {
    if (smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL])
      delete smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL];
    clearAtomPickerHighlights(smilesCol);
    return;
  }

  smilesCol.setTag(CHEM_ATOM_PICKER_LINKED_COL, mol3DColName);
  const mol3DCol = df?.col(mol3DColName);
  if (mol3DCol)
    mol3DCol.setTag(CHEM_ATOM_PICKER_LINKED_SMILES_COL, smilesCol.name);
}

/** Strips picker providers, fires clear-all event, and invalidates the grid. */
export function clearAtomPickerHighlights(smilesCol: DG.Column): void {
  try {
    const providers = (smilesCol.temp?.[SUBSTRUCT_PROVIDERS_KEY] ?? []) as
      Array<{__atomPicker?: boolean}>;
    if (providers.length > 0) {
      smilesCol.temp[SUBSTRUCT_PROVIDERS_KEY] = providers.filter(
        (p) => !p.__atomPicker);
    }
    grok.events.fireCustomEvent(CHEM_SELECTION_EVENT, {
      column: smilesCol, rowIdx: -1, atoms: [],
      persistent: true, clearAll: true,
    });
    const dfName = smilesCol.dataFrame?.name;
    if (dfName) {
      const tv = grok.shell.getTableView(dfName);
      tv?.grid?.invalidate();
    }
  } catch (err: unknown) {
    _package.logger.error(
      `mol3d-link clearAtomPickerHighlights: ${err instanceof Error ? err.message : String(err)}`);
  }
}
