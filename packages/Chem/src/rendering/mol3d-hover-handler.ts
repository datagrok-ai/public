/**
 * Reverse 3D→2D hover handler for the interactive atom-highlighting bridge.
 * Extracted from `rdkit-cell-renderer.ts` to keep the renderer file focused.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {mapAtomIndices2Dto3D, AtomIndexMapping} from '../utils/atom-index-mapper';

/** Payload for CHEM_MOL3D_HOVER_EVENT. Mirrors `Mol3DHoverEventArgs` in
 *  `molstar-highlight-utils.ts`; duplicated to avoid a cross-package import. */
export interface Mol3DHoverEventArgs {
  mol3DColumnName: string;
  rowIdx: number;
  atom3DSerial: number | null;
  mode: 'preview' | 'paint' | 'erase';
}

interface CellTopology {
  bondAtoms: Map<number, [number, number]>;
}

/** Subset of RDKitCellRenderer needed by the handler. Shape is identical to
 *  the renderer's private members so `this as unknown as Mol3DHoverRendererDeps`
 *  is a safe cast. Avoids a circular import on the full renderer class. */
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
}

/** Parses and validates a raw CHEM_MOL3D_HOVER_EVENT payload. */
export function parseMol3DHoverArgs(args: unknown): Mol3DHoverEventArgs | null {
  if (!args || typeof args !== 'object') return null;
  const {mol3DColumnName, rowIdx, atom3DSerial, mode} = args as Mol3DHoverEventArgs;
  if (typeof mol3DColumnName !== 'string' || typeof rowIdx !== 'number' ||
      rowIdx < 0) return null;
  return {mol3DColumnName, rowIdx, atom3DSerial, mode};
}

/** Inverts the 2D→3D mapping: serial → mapped idx → 2D atom index.
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
 * Handles CHEM_MOL3D_HOVER_EVENT: reverses the 3D atom serial to a 2D atom
 * index and routes preview/paint/erase to the renderer's shared provider path.
 * `atom3DSerial: null` means "cursor left the atom" — clears the preview only.
 *
 * NOTE: resolves the DataFrame via `grok.shell.tv?.grid?.dataFrame`; events
 * fired without a stable focused view silently bail at the grid-lookup guard.
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

  // 100x100: we only need bondAtoms (topology); pixel positions are dimension-dependent.
  const cellInfo = renderer._getCellAtomPositions(smiles2D, 100, 100);
  const bondAtoms = cellInfo?.bondAtoms ?? new Map<number, [number, number]>();

  if (mode === 'paint')
    renderer._addAtomToRow(smilesCol, rowIdx, idx2D, bondAtoms);
  else if (mode === 'erase')
    renderer._removeAtomFromRow(smilesCol, rowIdx, idx2D, bondAtoms);
  else {
    // Mark as 3D-sourced before setting so 2D mousemoves don't wipe it.
    renderer._previewFrom3D = true;
    renderer._setPreviewAtom(smilesCol, rowIdx, idx2D, bondAtoms);
  }
  grid.invalidate();
}
