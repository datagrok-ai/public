/**
 * Helpers for manual SMILES↔Molecule3D column linking (Chem atom-picker
 * bridge). The link is stored as a persistent tag on the SMILES column:
 *
 *   smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL] = mol3DColName
 *
 * These helpers are used by the link widget in
 * `panels/mol3d-atom-picker-link-panel.ts` to establish or change the
 * link. They take care of:
 *
 * - Reading the link with backward-compatible `col.temp[...]` fallback
 *   for sessions that predate the persistent-tag migration.
 * - Writing/clearing the link through the Datagrok MapProxy, which
 *   routes raw `delete col.tags[KEY]` to the native `grok_Map_Delete`.
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

/** Chem cross-package event for atom selection sync — mirrors the
 *  constant in Chem's `src/constants.ts`. Defined here to avoid a
 *  cross-package import. */
const CHEM_SELECTION_EVENT = 'chem-interactive-selection-changed';

/** Key where Chem's cell renderer stores atom-picker providers in
 *  `col.temp`. Mirrors `ChemTemps.SUBSTRUCT_PROVIDERS` from
 *  `@datagrok-libraries/chem-meta`. */
const SUBSTRUCT_PROVIDERS_KEY = 'substruct-providers';

/** Reads the link on a SMILES column. Prefers the persistent
 *  `col.tags[...]` location; falls back to `col.temp[...]` for
 *  backward compatibility with older session-only links. */
function getLinkedMol3DColName(smilesCol: DG.Column): string | null {
  const fromTags = smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL];
  if (fromTags) return fromTags;
  const fromTemp = smilesCol.temp?.[CHEM_ATOM_PICKER_LINKED_COL];
  return typeof fromTemp === 'string' ? fromTemp : null;
}

/** Finds the SMILES column currently linked TO the given Mol3D column
 *  (i.e., the SMILES column whose tag value equals `mol3DColName`). */
export function findLinkedSmilesColName(
  df: DG.DataFrame, mol3DColName: string): string | null {
  const c = df.columns.toList().find(
    (col) => getLinkedMol3DColName(col) === mol3DColName);
  return c?.name ?? null;
}

/** Removes any existing link that points at `mol3DColName` from any
 *  SMILES column of `df`. `delete col.tags[KEY]` routes through the
 *  Datagrok MapProxy's `deleteProperty` trap which calls the native
 *  `grok_Map_Delete`, so this correctly unsets the persistent tag.
 *  Also wipes residual atom-picker highlights for every unlinked
 *  column so the 3D viewer repaints. */
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
}

/** Writes a link on `smilesCol` pointing at `mol3DColName`, or clears
 *  it if `mol3DColName` is null/empty. On clear, also wipes any active
 *  atom-picker highlights so the 3D viewer drops its overpaint. */
export function setSmilesColLink(
  smilesCol: DG.Column, mol3DColName: string | null): void {
  if (smilesCol.temp?.[CHEM_ATOM_PICKER_LINKED_COL])
    delete smilesCol.temp[CHEM_ATOM_PICKER_LINKED_COL];
  if (!mol3DColName) {
    if (smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL])
      delete smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL];
    clearAtomPickerHighlights(smilesCol);
    return;
  }
  smilesCol.setTag(CHEM_ATOM_PICKER_LINKED_COL, mol3DColName);
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
