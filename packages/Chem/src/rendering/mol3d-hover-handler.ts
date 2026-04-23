/**
 * Reverse 3D→2D hover handler — extracted from `rdkit-cell-renderer.ts`.
 *
 * Fired when BiostructureViewer's Molstar viewer posts `CHEM_MOL3D_HOVER_EVENT`
 * because the user is hovering a ligand atom in a 3D pose. The handler looks
 * up the linked 2D SMILES column (via `CHEM_ATOM_PICKER_LINKED_COL`), reverses
 * the 3D Molstar atom serial back to a 2D atom index using the pre-computed
 * `AtomIndexMapping`, and delegates to the renderer's existing 2D rendering
 * methods — so both directions share the same provider storage.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {mapAtomIndices2Dto3D, AtomIndexMapping} from '../utils/atom-index-mapper';

/** Payload for CHEM_MOL3D_HOVER_EVENT fired by BiostructureViewer's Molstar
 *  viewer (3D→2D reverse bridge). Mirrors `Mol3DHoverEventArgs` in
 *  `molstar-highlight-utils.ts`; duplicated here to avoid a cross-package
 *  import from Chem into BiostructureViewer. */
export interface Mol3DHoverEventArgs {
  mol3DColumnName: string;
  rowIdx: number;
  atom3DSerial: number | null;
  mode: 'preview' | 'paint' | 'erase';
}

/** Minimal topology info needed to render a 2D atom highlight — subset of
 *  the full `CellInteractiveInfo` the renderer keeps per-cell. */
interface CellTopology {
  bondAtoms: Map<number, [number, number]>;
}

/** Callbacks exposed by `RDKitCellRenderer` that the handler needs. Lets us
 *  keep the handler in its own module without a circular import on the
 *  renderer class. Each member mirrors the same-named method on the renderer;
 *  shape kept identical so `this as unknown as Mol3DHoverRendererDeps` is a
 *  safe cast from the class body. */
export interface Mol3DHoverRendererDeps {
  rdKitModule: RDModule;
  _previewFrom3D: boolean;
  _getLinkedMol3DColName(col: DG.Column): string | null;
  _removePreviewAtom(): void;
  _setPreviewAtom(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void;
  _addAtomToRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void;
  _removeAtomFromRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void;
  _getCellAtomPositions(molString: string, cssWidth: number,
    cssHeight: number): CellTopology | null;
  _clearRendersCache(): void;
}

/** Validates a raw event payload against the `Mol3DHoverEventArgs` shape.
 *  Returns the typed payload on success, `null` on any type/shape mismatch
 *  (the subscription uses `args: unknown` so the guard is required). */
export function parseMol3DHoverArgs(args: unknown): Mol3DHoverEventArgs | null {
  if (!args || typeof args !== 'object') return null;
  const {mol3DColumnName, rowIdx, atom3DSerial, mode} = args as Mol3DHoverEventArgs;
  if (typeof mol3DColumnName !== 'string' || typeof rowIdx !== 'number' ||
      rowIdx < 0) return null;
  return {mol3DColumnName, rowIdx, atom3DSerial, mode};
}

/** Inverts the forward 2D→3D atom-index mapping. Mirrors the transform in
 *  BSV's `computeSerials`:
 *    forward:  idx2D → mapped = mapping[idx2D] → serial = pdbSerials[mapped] || mapped+1
 *    reverse:  serial → mapped = (pdbSerials.indexOf(serial) ?? serial-1) → idx2D = mapping.indexOf(mapped)
 *  Returns -1 when the serial is outside the mapped range. */
export function reverseMap3DSerialTo2DIdx(
  mapping: AtomIndexMapping, atom3DSerial: number,
): number {
  const mapped = mapping.pdbSerials ?
    mapping.pdbSerials.indexOf(atom3DSerial) :
    atom3DSerial - 1;
  if (mapped < 0) return -1;
  return mapping.mapping.indexOf(mapped);
}

/**
 * Reverse 3D→2D hover handler. Called when BiostructureViewer's Molstar
 * viewer fires CHEM_MOL3D_HOVER_EVENT — the user is hovering a ligand
 * atom in the 3D pose, and we need to highlight the matching 2D atom.
 *
 * Flow:
 *   1. Resolve the linked 2D SMILES column (the one whose
 *      CHEM_ATOM_PICKER_LINKED_COL tag points to the fired 3D column).
 *   2. Compute (on demand) the 2D↔3D atom-index mapping for the row's
 *      (smiles, pose3D) pair via `mapAtomIndices2Dto3D`.
 *   3. Reverse-look-up the 3D Molstar atom serial → 2D atom index via
 *      `reverseMap3DSerialTo2DIdx`.
 *   4. Route by `mode`: preview / paint / erase — each delegates to the
 *      renderer's existing 2D rendering methods so both directions share
 *      the same provider storage.
 *
 * Events with `atom3DSerial: null` mean "cursor left all ligand atoms" —
 * clear the transient preview (painted atoms are untouched).
 *
 * NOTE: resolves the DataFrame via `grok.shell.tv?.grid?.dataFrame`. In
 * test harnesses that fire multiple events rapidly without a stable
 * focused view, subsequent events may silently bail at the grid-lookup
 * guard — see `atom-picker-3d-hover-tests.ts` for which scenarios can and
 * cannot be automated today.
 */
export function handleMol3DHoverEvent(
  renderer: Mol3DHoverRendererDeps, args: unknown,
): void {
  const parsed = parseMol3DHoverArgs(args);
  if (!parsed) return;
  const {mol3DColumnName, rowIdx, atom3DSerial, mode} = parsed;

  const grid = grok.shell.tv?.grid;
  if (!grid) return;
  const df = grid.dataFrame;
  if (!df) return;

  // Find the 2D SMILES column linked to this 3D pose column.
  const smilesCol = df.columns.toList().find(
    (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE &&
      renderer._getLinkedMol3DColName(c) === mol3DColumnName);
  if (!smilesCol) return;

  // Cursor left an atom → clear the preview only. Painted atoms stay.
  if (atom3DSerial == null) {
    renderer._removePreviewAtom();
    grid.invalidate();
    return;
  }

  const smiles2D = smilesCol.get(rowIdx);
  const mol3DCol = df.col(mol3DColumnName);
  const pose3D = mol3DCol?.get(rowIdx);
  if (!smiles2D || !pose3D) return;

  let mapping: AtomIndexMapping | null;
  try {
    mapping = mapAtomIndices2Dto3D(renderer.rdKitModule, smiles2D, pose3D);
  } catch {
    return;
  }
  if (!mapping) return;

  const idx2D = reverseMap3DSerialTo2DIdx(mapping, atom3DSerial);
  if (idx2D < 0) return;

  // `_getCellAtomPositions` cache key is `molString|WxH`. Using 100x100 here
  // matches what `_removePreviewAtom` already does — we only need
  // `bondAtoms` (topology, dimension-independent), not pixel positions.
  const cellInfo = renderer._getCellAtomPositions(smiles2D, 100, 100);
  const bondAtoms = cellInfo?.bondAtoms ?? new Map<number, [number, number]>();

  if (mode === 'paint')
    renderer._addAtomToRow(smilesCol, rowIdx, idx2D, bondAtoms);
  else if (mode === 'erase')
    renderer._removeAtomFromRow(smilesCol, rowIdx, idx2D, bondAtoms);
  else {
    // Mark the preview as 3D-sourced BEFORE setting it, so the 2D
    // `_onDocumentMouseMove` no-modifier branch does not wipe it on
    // the next DOM mousemove over the 3D canvas.
    renderer._previewFrom3D = true;
    renderer._setPreviewAtom(smilesCol, rowIdx, idx2D, bondAtoms);
  }

  renderer._clearRendersCache();
  grid.invalidate();
}
