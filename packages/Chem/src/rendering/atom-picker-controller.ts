/**
 * Interactive 2D↔3D atom-picker controller for `RDKitCellRenderer`.
 *
 * Holds all hover/preview/paint/erase state and the cross-package selection
 * event plumbing. The renderer keeps slim override shims that delegate to a
 * single instance of this controller; render-machinery (mol fetch, raster
 * cache, atom-position SVG extraction) stays on the renderer and is reached
 * through `AtomPickerRendererDeps`.
 *
 * Extracted from `rdkit-cell-renderer.ts` to keep the renderer file focused.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {addSubstructProvider, ISubstructProvider} from '@datagrok-libraries/chem-meta/src/types';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {
  CHEM_ATOM_SELECTION_EVENT, CHEM_ATOM_PICKER_LINKED_3D_COL_TAG,
  CHEM_MOL3D_SELECTION_EVENT,
} from '@datagrok-libraries/chem-meta/src/types';
import {mapAtomIndices2Dto3D, AtomIndexMapping} from '../utils/atom-index-mapper';
import {findNearestAtom, computeSelectedBonds} from '../utils/chem-atom-picker-utils';
import {handleMol3DHoverEvent, Mol3DHoverRendererDeps} from './mol3d-hover-handler';

/** Shape of the cross-package atom selection event. */
export interface ChemSelectionEvent {
  column: DG.Column;
  rowIdx: number;
  atoms: number[];
  persistent?: boolean;
  clearAll?: boolean;
  mapping3D?: AtomIndexMapping | null;
  mol3DColumnName?: string;
}

/** Substruct provider with atom-picker metadata fields. */
export interface AtomPickerProvider extends ISubstructProvider {
  __atomPicker?: boolean;
  __rowIdx?: number;
  __atoms?: Set<number>;
}

/** Per-cell layout info for atom hit-testing. Cached by "molString|WxH". */
export interface CellInteractiveInfo {
  /** Atom index → center in cell-local CSS pixels. */
  positions: Map<number, {x: number, y: number}>;
  /** Bond index → [atomA, atomB] pair from RDKit SVG class annotations. */
  bondAtoms: Map<number, [number, number]>;
}

/** Subset of `RDKitCellRenderer` the picker controller reads/calls.
 *  Kept narrow so the controller can be unit-tested without the full
 *  renderer surface and mirrors the pattern in `mol3d-hover-handler.ts`. */
export interface AtomPickerRendererDeps {
  readonly rdKitModule: RDModule;
  /** Returns positions/bondAtoms for the cell, populating the cache. */
  _getCellAtomPositions(molString: string, cssWidth: number,
    cssHeight: number): CellInteractiveInfo | null;
}

export class AtomPickerController {
  /** Per-cell layout info cached by "molString|w|h" — atom positions + bond connectivity. */
  atomPositionsCache: DG.LruCache<string, CellInteractiveInfo> =
    new DG.LruCache<string, CellInteractiveInfo>();
  /** Last atom processed by hover — dedup guard so we don't re-fire on every pixel. */
  _lastHoveredAtom: {col: string, rowIdx: number, atomIdx: number, erase?: boolean} | null = null;
  /** Transient preview atom (no modifier hover). Cleared on mouse-leave or shift. */
  _previewAtomIdx: number | null = null;
  /** Row the preview is on. */
  _previewRowIdx: number = -1;
  /** True when the active preview was set by a 3D→2D hover event rather than a 2D hover.
   *  A 3D-sourced preview must survive 2D mousemoves and is cleared only by the
   *  corresponding 3D "cursor left the atom" event (`atom3DSerial === null`). */
  _previewFrom3D: boolean = false;

  /** Registered once on first render; grid rebuilds get re-attached automatically. */
  _globalListenersAttached = false;
  _wiredCanvases: WeakSet<HTMLCanvasElement> = new WeakSet();

  constructor(private readonly deps: AtomPickerRendererDeps) {}

  // -- Provider lookups -----------------------------------------------------

  _getProviders(col: DG.Column): AtomPickerProvider[] {
    return (col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as AtomPickerProvider[];
  }

  _getProviderForRow(col: DG.Column, rowIdx: number): AtomPickerProvider | undefined {
    return this._getProviders(col).find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
  }

  /** Builds a selection event, attaching 3D mapping when a linked Mol3D column exists. */
  _buildSelectionEvent(
    col: DG.Column, rowIdx: number, atoms: number[], persistent: boolean = true,
  ): ChemSelectionEvent {
    const event: ChemSelectionEvent = {column: col, rowIdx, atoms, persistent};
    try {
      const df = col.dataFrame;
      if (df) {
        const linkedColName = this._getLinkedMol3DColName(col);
        const mol3DCol = linkedColName ? df.col(linkedColName) : null;
        if (mol3DCol && rowIdx >= 0) {
          const smiles2D = col.get(rowIdx);
          const pose3D = mol3DCol.get(rowIdx);
          if (smiles2D && pose3D) {
            const mapping = mapAtomIndices2Dto3D(this.deps.rdKitModule, smiles2D, pose3D);
            if (mapping) {
              event.mapping3D = mapping;
              event.mol3DColumnName = mol3DCol.name;
            }
          }
        }
      }
    } catch {
      /* mapping failed */
    }
    return event;
  }

  _isSameAtom(colName: string, rowIdx: number, atomIdx: number, erase?: boolean): boolean {
    const last = this._lastHoveredAtom;
    if (!last) return false;
    return last.col === colName &&
      last.rowIdx === rowIdx &&
      last.atomIdx === atomIdx &&
      (erase === undefined || last.erase === erase);
  }

  /** Returns `true` (same atom, caller should bail) or records it and returns `false`. */
  _trackHoveredAtom(
    colName: string, rowIdx: number, atomIdx: number, erase?: boolean,
  ): boolean {
    if (this._isSameAtom(colName, rowIdx, atomIdx, erase)) return true;
    this._lastHoveredAtom = {col: colName, rowIdx, atomIdx, erase};
    return false;
  }

  /** Returns true when the picker is active for the column (a Mol3D link exists). */
  _isPickerActive(col: DG.Column): boolean {
    return !!this._getLinkedMol3DColName(col);
  }

  /** Returns the linked Mol3D column name, preferring persistent tags over session temp. */
  _getLinkedMol3DColName(col: DG.Column): string | null {
    const fromTags = col.tags[CHEM_ATOM_PICKER_LINKED_3D_COL_TAG];
    if (fromTags) return fromTags;
    const fromTemp = col.temp?.[CHEM_ATOM_PICKER_LINKED_3D_COL_TAG];
    return typeof fromTemp === 'string' ? fromTemp : null;
  }

  // -- Render-time hooks ----------------------------------------------------

  /** Subscribes once to the 3D→2D bridge fired by BiostructureViewer's Molstar viewer. */
  ensureGlobalListeners(): void {
    if (this._globalListenersAttached) return;
    this._globalListenersAttached = true;
    grok.events.onCustomEvent(CHEM_MOL3D_SELECTION_EVENT).subscribe(
      (args: unknown) => this._onMol3DHoverEvent(args));
  }

  /** Prewarms the atom-positions cache for one cell. Called from the
   *  renderer's `onMouseEnter` (per cell the cursor actually enters), not
   *  from `render()`, so scrolling 100s of rows doesn't queue 100s of
   *  SVG-extraction tasks. The work is deferred via `setTimeout(0)` so the
   *  current event handler returns before the heavy SVG generation runs. */
  warmCacheIfNeeded(col: DG.Column, molString: string, cssW: number, cssH: number): void {
    if (!this._isPickerActive(col)) return;
    const dpr = window.devicePixelRatio || 1;
    const cacheKey = molString + '|' + Math.round(cssW * dpr) + 'x' + Math.round(cssH * dpr);
    if (this.atomPositionsCache.has(cacheKey)) return;
    setTimeout(
      () => {this.deps._getCellAtomPositions(molString, cssW, cssH);},
    );
  }

  // -- Hover-paint atom picker ----------------------------------------------
  // No-modifier click shows a transient preview of one atom under the cursor.
  // Shift+click/drag paints atoms persistently; Shift+Ctrl/Cmd erases.
  // Mirrors the 3D Molstar-side modifier scheme so paint/erase keys match.
  // Escape clears all picker selections in the grid's DataFrame.
  // Wired via DG.GridCellRenderer overrides (no document-wide listeners).

  /** Resolves the atom under the cursor. Coordinates are grid-canvas-relative
   *  (e.offsetX/Y) and translated into cell-local space via gridCell.bounds. */
  _hitTestAtomInCell(gridCell: DG.GridCell, e: MouseEvent): {
    cellInfo: CellInteractiveInfo; nearest: number;
  } | null {
    if (gridCell.tableRowIndex == null || gridCell.tableRowIndex < 0) return null;
    const molString: string = gridCell.cell.value;
    if (!molString || DG.chem.Sketcher.isEmptyMolfile(molString)) return null;

    const bounds = gridCell.bounds;
    const pointerCellX = e.offsetX - bounds.left;
    const pointerCellY = e.offsetY - bounds.top;
    if (pointerCellX < 0 || pointerCellY < 0 ||
        pointerCellX > bounds.width || pointerCellY > bounds.height) return null;

    const cellInfo = this.deps._getCellAtomPositions(molString, bounds.width, bounds.height);
    if (!cellInfo || cellInfo.positions.size === 0) return null;

    const nearest = findNearestAtom(cellInfo.positions, pointerCellX, pointerCellY);
    if (nearest === null) return null;

    return {cellInfo, nearest};
  }

  /**
   * No-modifier hover handler — highlights ONLY the single atom under the
   * cursor (preview mode). Moving to a different atom highlights that one
   * instead. Moving off all atoms clears the preview.
   *
   * `onMouseLeave` handles preview-clear when the cursor leaves the cell
   * toward a non-molecule column.
   */
  onMouseDown(gridCell: DG.GridCell, e: MouseEvent): void {
    const col = gridCell.tableColumn;
    if (!col || col.semType !== DG.SEMTYPE.MOLECULE || !this._isPickerActive(col))
      return;

    const hit = this._hitTestAtomInCell(gridCell, e);
    if (hit)
      e.preventDefault();
    const isErase = e.shiftKey && (e.ctrlKey || e.metaKey);
    const isPaint = e.shiftKey && !(e.ctrlKey || e.metaKey);

    // -- Shift: paint. Shift+Ctrl / Shift+Cmd: erase. No-modifier: preview. --
    if (isPaint || isErase) {
      this._previewAtomIdx = null;
      if (!hit) return;
      const rowIdx = gridCell.tableRowIndex!;
      if (this._trackHoveredAtom(col.name, rowIdx, hit.nearest, isErase)) return;
      if (isErase)
        this._removeAtomFromRow(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
      else
        this._addAtomToRow(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
      gridCell.render();
      return;
    }

    // -- No modifier: preview mode --
    if (!hit) {
      // A 3D-sourced preview must survive 2D mousemoves; it is cleared only by
      // the corresponding 3D "cursor left the atom" event (atom3DSerial === null).
      if (this._previewAtomIdx !== null && !this._previewFrom3D)
        this._removePreviewAtom();
      this._lastHoveredAtom = null;
      return;
    }

    const rowIdx = gridCell.tableRowIndex!;
    if (this._trackHoveredAtom(col.name, rowIdx, hit.nearest))
      return;

    this._previewFrom3D = false; // 2D hover takes ownership
    this._setPreviewAtom(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
    gridCell.render();
  }

  /** Clears the 2D-sourced preview when the cursor leaves a molecule cell. */
  onMouseLeave(gridCell: DG.GridCell, _e: MouseEvent): void {
    const col = gridCell.tableColumn;
    if (!col || col.semType !== DG.SEMTYPE.MOLECULE || !this._isPickerActive(col))
      return;
    if (this._previewAtomIdx !== null && !this._previewFrom3D)
      this._removePreviewAtom();
    this._lastHoveredAtom = null;
  }

  // Cast is safe: Mol3DHoverRendererDeps mirrors this controller's members by
  // name, plus `_getCellAtomPositions` / `rdKitModule`
  // exposed via the renderer's shims (the controller's `deps` proxies them).
  _onMol3DHoverEvent(args: unknown): void {
    handleMol3DHoverEvent(this as unknown as Mol3DHoverRendererDeps, args);
  }

  /** Escape when a molecule cell is focused — Datagrok routes the keydown here
   *  at Dart level and it does NOT bubble to grid.root. The grid-root listener
   *  in `_attachGridCanvasShiftListener` covers non-molecule focused cells. */
  onKeyDown(gridCell: DG.GridCell, e: KeyboardEvent): void {
    if (e.key === 'Escape')
      this._clearPickerSelection(gridCell.grid);
  }

  /** Clears all picker providers and broadcasts a persistent clear-all event
   *  so 3D viewers drop their cached overpaint. Called from both Escape paths. */
  _clearPickerSelection(grid: DG.Grid): void {
    const df = grid?.dataFrame;
    if (!df) return;
    const molCols = df.columns.toList().filter(
      (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE && this._isPickerActive(c));
    if (molCols.length === 0) return;
    for (const col of molCols) {
      col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = this._getProviders(col)
        .filter((p) => !p.__atomPicker);
      grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
        column: col, rowIdx: -1, atoms: [], persistent: true, clearAll: true,
      });
    }
    this._lastHoveredAtom = null;
    this._previewAtomIdx = null;
    this._previewFrom3D = false;
    grid.invalidate();
  }

  /** Shift-drag handler from the canvas capture listener. Re-does the grid
   *  hit-test from page coords because Datagrok captures shift-mousemove
   *  at the Dart level and does not forward it to the renderer's onMouseMove. */
  _onShiftDragMouseMove(e: MouseEvent): void {
    const isErase = e.shiftKey && (e.ctrlKey || e.metaKey);
    const isPaint = e.shiftKey && !(e.ctrlKey || e.metaKey);
    if (!isPaint && !isErase) return;

    const grid = grok.shell.tv?.grid;
    if (!grid) return;
    let gridRoot: HTMLElement | null = null;
    try {gridRoot = grid.root;} catch {return;}
    if (!gridRoot) return;

    const gridRect = gridRoot.getBoundingClientRect();
    const localX = e.clientX - gridRect.left;
    const localY = e.clientY - gridRect.top;
    if (localX < 0 || localY < 0 || localX > gridRect.width || localY > gridRect.height) return;

    let gridCell: DG.GridCell | null = null;
    try {gridCell = grid.hitTest(localX, localY);} catch {return;}
    const col = gridCell?.tableColumn;
    if (!gridCell || !col || col.semType !== DG.SEMTYPE.MOLECULE) return;
    if (!this._isPickerActive(col)) return;

    const molString: string = gridCell.cell.value;
    if (!molString || DG.chem.Sketcher.isEmptyMolfile(molString)) return;

    const bounds = gridCell.bounds;
    const pointerCellX = e.clientX - (gridRect.left + bounds.left);
    const pointerCellY = e.clientY - (gridRect.top + bounds.top);
    if (pointerCellX < 0 || pointerCellY < 0 ||
        pointerCellX > bounds.width || pointerCellY > bounds.height) return;

    const cellInfo = this.deps._getCellAtomPositions(molString, bounds.width, bounds.height);
    if (!cellInfo || cellInfo.positions.size === 0) return;

    const nearest = findNearestAtom(cellInfo.positions, pointerCellX, pointerCellY);
    if (nearest === null) return;

    this._previewAtomIdx = null;
    const rowIdx = gridCell.tableRowIndex!;
    if (this._trackHoveredAtom(col.name, rowIdx, nearest, isErase)) return;

    if (isErase)
      this._removeAtomFromRow(col, rowIdx, nearest, cellInfo.bondAtoms);
    else
      this._addAtomToRow(col, rowIdx, nearest, cellInfo.bondAtoms);

    grid.invalidate();
  }

  _getPersistentAtoms(col: DG.Column, rowIdx: number): Set<number> {
    const prov = this._getProviderForRow(col, rowIdx);
    return new Set<number>(prov?.__atoms ?? []);
  }

  _addAtomToRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    const current = this._getPersistentAtoms(col, rowIdx);
    if (current.has(atomIdx)) return; // already selected
    current.add(atomIdx);
    this._updateRowSelection(col, rowIdx, [...current], {add: true}, bondAtoms);
  }

  _removeAtomFromRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    const prov = this._getProviderForRow(col, rowIdx);
    if (!prov?.__atoms?.has(atomIdx)) return;

    // Compute the reduced set on a copy — `prov` is a snapshot from
    // `col.temp`, so in-place `Set.delete` would not survive the next read.
    const reduced = new Set(prov.__atoms);
    reduced.delete(atomIdx);

    if (reduced.size === 0) {
      // Last atom removed — drop the picker provider and broadcast empty.
      col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = this._getProviders(col).filter(
        (p) => !p.__atomPicker || p.__rowIdx !== rowIdx);
      grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT,
        this._buildSelectionEvent(col, rowIdx, []));
    } else {
      // Replace the provider via the canonical write path — this re-builds
      // the closure-captured atomsArr/bonds so 2D renders the reduced set,
      // and pushes the new providers array through the col.temp proxy in
      // one atomic write. (`add: false` clears prior and uses boxed.)
      this._updateRowSelection(col, rowIdx, [...reduced], {add: false},
        bondAtoms, true, true);
    }
  }

  /** Renders persistent atoms + the hovered atom as a non-persistent preview,
   *  but stores ONLY the persistent set in `__atoms` so the preview atom
   *  doesn't leak into the canonical state. (`col.temp` is a Datagrok proxy
   *  that returns fresh snapshots on each read; post-hoc mutation of an
   *  accessed object would be lost — so we let `_updateRowSelection` write
   *  the correct `__atoms` on its first pass via `storeAtoms`.) */
  _setPreviewAtom(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    this._previewAtomIdx = atomIdx;
    this._previewRowIdx = rowIdx;
    const persistent = this._getPersistentAtoms(col, rowIdx);
    const combined = new Set(persistent);
    combined.add(atomIdx);
    this._updateRowSelection(col, rowIdx, [...combined], {add: true},
      bondAtoms, true, false, [...persistent]);
  }

  /** Restores display to the persistent selection, dropping the preview atom. */
  _removePreviewAtom(): void {
    const prevAtom = this._previewAtomIdx;
    const prevRow = this._previewRowIdx;
    this._previewAtomIdx = null;
    this._previewFrom3D = false;
    if (prevAtom === null) return;

    const grid = grok.shell.tv?.grid;
    if (!grid) return;
    const df = grid.dataFrame;
    if (!df) return;
    const molCol = df.columns.toList().find(
      (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE && this._isPickerActive(c));
    if (!molCol) return;

    const persistent = this._getPersistentAtoms(molCol, prevRow);
    if (persistent.size > 0) {
      // Restore just the persistent atoms in 2D and 3D.
      // Fire event as non-persistent so the cache keeps the stable Shift-paint set.
      const cellInfo = this.deps._getCellAtomPositions(
        molCol.get(prevRow), 100, 100);
      const bondAtoms = cellInfo?.bondAtoms ?? new Map();
      this._updateRowSelection(molCol, prevRow, [...persistent], {add: true}, bondAtoms, true, false);
    } else {
      // No persistent atoms — clear everything.
      molCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = this._getProviders(molCol)
        .filter((p) => !p.__atomPicker || p.__rowIdx !== prevRow);
      // Fire event to clear 3D highlights too (non-persistent).
      grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
        column: molCol, rowIdx: prevRow, atoms: [], persistent: false,
      });
    }
    grid.invalidate();
  }


  /**
   * Writes an atom-picker ISubstructProvider for `rowIdx`, replacing the
   * previous one for that row. Other rows' providers are preserved.
   *
   * `storeAtoms` (optional) overrides what's written to `provider.__atoms`.
   * `boxed` controls what's RENDERED (via the closure-captured `atomsArr`);
   * `__atoms` is the canonical persistent set used by callers that read
   * the picker state. They differ only for previews — `_setPreviewAtom`
   * passes a wider `boxed` (persistent + preview atom) and a narrower
   * `storeAtoms` (persistent-only) so the transient atom is shown but not
   * stored. When omitted, `__atoms` defaults to the rendered set.
   */
  _updateRowSelection(col: DG.Column, rowIdx: number, boxed: number[],
    modifiers: {add: boolean}, bondAtoms: Map<number, [number, number]>,
    fire3DEvent: boolean = true, persistent: boolean = true,
    storeAtoms?: number[]): void {
    const prior = this._getProviderForRow(col, rowIdx);
    const current = new Set<number>(prior?.__atoms ?? []);

    // No modifier → replace existing selection with boxed atoms.
    // Shift (add) → keep existing and union with boxed atoms.
    if (!modifiers.add) current.clear();
    for (const a of boxed) current.add(a);

    // Yellow (#FFFF00): avoids collision with Datagrok's green row-selection tint.
    const pickColor: number[] = [1.0, 1.0, 0.0, 1.0];
    const atomsArr = [...current];
    const highlightAtomColors: {[k: number]: number[]} = {};
    for (const a of atomsArr) highlightAtomColors[a] = pickColor;
    const {bondsArr, highlightBondColors} =
      computeSelectedBonds(current, bondAtoms, pickColor);

    const provider: AtomPickerProvider = {
      getSubstruct: (ridx) => {
        if (!this._isPickerActive(col)) return undefined;
        if (ridx !== rowIdx) return undefined;
        return {atoms: atomsArr, bonds: bondsArr, highlightAtomColors, highlightBondColors};
      },
    };
    provider.__atomPicker = true;
    provider.__rowIdx = rowIdx;
    provider.__atoms = storeAtoms !== undefined ? new Set(storeAtoms) : current;

    // Replace only the prior provider for this row; keep picks on other rows.
    col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = this._getProviders(col).filter(
      (p) => !p.__atomPicker || p.__rowIdx !== rowIdx);
    addSubstructProvider(col.temp, provider);

    if (fire3DEvent) {
      grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT,
        this._buildSelectionEvent(col, rowIdx, atomsArr, persistent));
    }
  }

  // -- Mol3DHoverRendererDeps shape -----------------------------------------
  // `handleMol3DHoverEvent` expects the controller-owned hover state plus the
  // three renderer-owned hooks below; we proxy them so `this as unknown as
  // Mol3DHoverRendererDeps` is sound without forcing a wider deps interface.

  get rdKitModule(): RDModule {return this.deps.rdKitModule;}

  _getCellAtomPositions(molString: string, cssWidth: number, cssHeight: number):
      CellInteractiveInfo | null {
    return this.deps._getCellAtomPositions(molString, cssWidth, cssHeight);
  }
}
