/**
 * Helpers for manual SMILES↔Molecule3D column linking (Chem atom-picker
 * bridge). The link is stored as a **reciprocal pair of persistent tags**,
 * one on each column pointing to its partner:
 *
 *   smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL]        = mol3DColName
 *   mol3DCol.tags [CHEM_ATOM_PICKER_LINKED_SMILES_COL] = smilesColName
 *
 * Two-way tags let either side find its partner in O(1) without iterating
 * the dataframe, and make the relationship self-describing (a glance at
 * either column's tags shows the link).
 *
 * These helpers are used by the link widget in
 * `panels/mol3d-atom-picker-link-panel.ts` to establish or change the
 * link. They take care of:
 *
 * - Reading the link with backward-compatible `col.temp[...]` fallback
 *   for sessions that predate the persistent-tag migration.
 * - Writing/clearing both tags atomically through the Datagrok MapProxy,
 *   which routes raw `delete col.tags[KEY]` to the native `grok_Map_Delete`.
 * - Clearing residual 2D/3D atom-picker highlights on unlink (strips
 *   atom-picker providers, fires the cross-package clear-all event,
 *   invalidates the grid).
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';

/** Persistent column tag used to link a SMILES (Molecule) column to a
 *  specific Molecule3D column for the Chem atom-picker bridge. Mirrors
 *  the same-named constant in Chem's `src/constants.ts` and (eventually)
 *  in `@datagrok-libraries/bio/src/viewers/molecule3d`; duplicated here
 *  to avoid a hard dependency on a bio release that publishes the
 *  constant. Keep the string value in sync across packages. */
export const CHEM_ATOM_PICKER_LINKED_COL = '%chem-atom-picker-linked-col';

/** Reciprocal tag on the Mol3D column pointing back at the linked SMILES
 *  column name. Lets Mol3D→SMILES lookup be O(1) and makes the link
 *  self-describing from either side. Mirror the same-named constant in
 *  Chem's `src/constants.ts`. */
export const CHEM_ATOM_PICKER_LINKED_SMILES_COL = '%chem-atom-picker-linked-smiles-col';

/** Chem cross-package event for atom selection sync — mirrors the
 *  constant in Chem's `src/constants.ts`. Defined here to avoid a
 *  cross-package import. */
const CHEM_SELECTION_EVENT = 'chem-interactive-selection-changed';

/** Key where Chem's cell renderer stores atom-picker providers in
 *  `col.temp`. Mirrors `ChemTemps.SUBSTRUCT_PROVIDERS` from
 *  `@datagrok-libraries/chem-meta`. */
const SUBSTRUCT_PROVIDERS_KEY = 'substruct-providers';

/** Reads the forward-direction link from a SMILES column — the Mol3D
 *  column name the picker targets. Prefers the persistent
 *  `col.tags[...]` location; falls back to `col.temp[...]` for
 *  backward compatibility with older session-only links. */
function getLinkedMol3DColName(smilesCol: DG.Column): string | null {
  const fromTags = smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL];
  if (fromTags) return fromTags;
  const fromTemp = smilesCol.temp?.[CHEM_ATOM_PICKER_LINKED_COL];
  return typeof fromTemp === 'string' ? fromTemp : null;
}

/** Reads the reverse-direction link from a Mol3D column — the SMILES
 *  column name that points back at it. Same tags-first, temp-fallback
 *  strategy as the forward direction. */
function getLinkedSmilesColName(mol3DCol: DG.Column): string | null {
  const fromTags = mol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL];
  if (fromTags) return fromTags;
  const fromTemp = mol3DCol.temp?.[CHEM_ATOM_PICKER_LINKED_SMILES_COL];
  return typeof fromTemp === 'string' ? fromTemp : null;
}

/** Finds the SMILES column currently linked TO the given Mol3D column.
 *  O(1) via the reciprocal tag on the Mol3D side; falls back to scanning
 *  SMILES columns (O(n)) for backward compatibility with old one-way
 *  links written before the reciprocal-tag migration. */
export function findLinkedSmilesColName(
  df: DG.DataFrame, mol3DColName: string): string | null {
  const mol3DCol = df.col(mol3DColName);
  if (mol3DCol) {
    const reciprocal = getLinkedSmilesColName(mol3DCol);
    if (reciprocal && df.col(reciprocal)) return reciprocal;
  }
  // Fallback: one-way-tagged SMILES columns from pre-migration sessions.
  const c = df.columns.toList().find(
    (col) => getLinkedMol3DColName(col) === mol3DColName);
  return c?.name ?? null;
}

/** Removes any existing link that points at `mol3DColName` — strips the
 *  forward tag from every linked SMILES column AND the reciprocal tag
 *  from the Mol3D column itself. `delete col.tags[KEY]` routes through
 *  the Datagrok MapProxy's `deleteProperty` trap which calls the native
 *  `grok_Map_Delete`, so this correctly unsets the persistent tags.
 *  Also wipes residual atom-picker highlights for every unlinked
 *  SMILES column so the 3D viewer repaints. */
export function clearLinksToMol3D(
  df: DG.DataFrame, mol3DColName: string): void {
  for (const c of df.columns.toList()) {
    if (getLinkedMol3DColName(c) === mol3DColName) {
      if (c.tags[CHEM_ATOM_PICKER_LINKED_COL])
        delete c.tags[CHEM_ATOM_PICKER_LINKED_COL];
      if (c.temp?.[CHEM_ATOM_PICKER_LINKED_COL])
        delete c.temp[CHEM_ATOM_PICKER_LINKED_COL];
      clearAtomPickerHighlights(c);
    }
  }
  // Clear the reciprocal tag on the Mol3D side too.
  const mol3DCol = df.col(mol3DColName);
  if (mol3DCol) {
    if (mol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL])
      delete mol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL];
    if (mol3DCol.temp?.[CHEM_ATOM_PICKER_LINKED_SMILES_COL])
      delete mol3DCol.temp[CHEM_ATOM_PICKER_LINKED_SMILES_COL];
  }
}

/** Writes a reciprocal-tag link: forward tag on `smilesCol` pointing at
 *  `mol3DColName`, plus reverse tag on the Mol3D column pointing back at
 *  `smilesCol.name`. If `mol3DColName` is null/empty, clears BOTH tags
 *  (forward on this SMILES col + reverse on whichever Mol3D col it
 *  previously pointed at) and wipes active atom-picker highlights. */
export function setSmilesColLink(
  smilesCol: DG.Column, mol3DColName: string | null): void {
  if (smilesCol.temp?.[CHEM_ATOM_PICKER_LINKED_COL])
    delete smilesCol.temp[CHEM_ATOM_PICKER_LINKED_COL];

  const df = smilesCol.dataFrame;
  // Find the Mol3D column this SMILES col was previously linked to so we
  // can strip its reciprocal tag even when the unlink target is null.
  const prevMol3DName = getLinkedMol3DColName(smilesCol);
  const prevMol3DCol = prevMol3DName && df ? df.col(prevMol3DName) : null;
  if (prevMol3DCol && prevMol3DCol.name !== mol3DColName) {
    if (prevMol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL])
      delete prevMol3DCol.tags[CHEM_ATOM_PICKER_LINKED_SMILES_COL];
    if (prevMol3DCol.temp?.[CHEM_ATOM_PICKER_LINKED_SMILES_COL])
      delete prevMol3DCol.temp[CHEM_ATOM_PICKER_LINKED_SMILES_COL];
  }

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

/** Clears any residual 2D/3D atom-picker highlights for a SMILES column
 *  that just got unlinked. Strips atom-picker providers from the
 *  column's `col.temp['substruct-providers']` list (preserving other
 *  providers), fires the Chem cross-package clear-all event (so Molstar
 *  drops its overpaint and the selection cache is flushed), and
 *  invalidates the grid so 2D cells redraw without highlights. */
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
