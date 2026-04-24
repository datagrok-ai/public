import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {CHEM_ATOM_PICKER_LINKED_COL} from '../constants';
import {GridCellRendererProxy, RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';

/** Test-only structural alias for private hover state on RDKitCellRenderer.
 *  Uses `as unknown as` (not intersection) because TypeScript collapses an
 *  intersection to `never` when any of the named members is already `private`
 *  in the base class. The alias names match the real field/method names exactly
 *  so a rename in the source produces a type error here. */
interface RDKitCellRendererInternals {
  _previewAtomIdx: number | null;
  _previewFrom3D: boolean;
  _lastHoveredAtom: {col: string; rowIdx: number; atomIdx: number; erase?: boolean} | null;
  _trackHoveredAtom(col: string, row: number, atom: number, erase?: boolean): boolean;
  onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void;
  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void;
}

category('renderer-hover', () => {
  let rdkitModule: any;

  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule');
  });

  // ---------------------------------------------------------------------------
  // Test 1 — GridCellRendererProxy must forward all 9 mouse/key events
  // ---------------------------------------------------------------------------

  test('proxy-forwards-onMouseMove', async () => {
    const flags = {
      enter: false, leave: false, down: false, up: false, move: false,
      click: false, dblClick: false, keyDown: false, keyPress: false,
    };

    class StubRenderer extends DG.GridCellRenderer {
      override get name(): string {return 'stub';}
      override get cellType(): string {return 'stub';}
      override onMouseEnter(_gc: DG.GridCell, _e: MouseEvent): void {flags.enter = true;}
      override onMouseLeave(_gc: DG.GridCell, _e: MouseEvent): void {flags.leave = true;}
      override onMouseDown(_gc: DG.GridCell, _e: MouseEvent): void {flags.down = true;}
      override onMouseUp(_gc: DG.GridCell, _e: MouseEvent): void {flags.up = true;}
      override onMouseMove(_gc: DG.GridCell, _e: MouseEvent): void {flags.move = true;}
      override onClick(_gc: DG.GridCell, _e: MouseEvent): void {flags.click = true;}
      override onDoubleClick(_gc: DG.GridCell, _e: MouseEvent): void {flags.dblClick = true;}
      override onKeyDown(_gc: DG.GridCell, _e: KeyboardEvent): void {flags.keyDown = true;}
      override onKeyPress(_gc: DG.GridCell, _e: KeyboardEvent): void {flags.keyPress = true;}
    }

    const stub = new StubRenderer();
    const proxy = new GridCellRendererProxy(stub, 'TestCell');

    const gc = {} as DG.GridCell;
    const me = {} as MouseEvent;
    const ke = {} as KeyboardEvent;

    proxy.onMouseEnter(gc, me);
    proxy.onMouseLeave(gc, me);
    proxy.onMouseDown(gc, me);
    proxy.onMouseUp(gc, me);
    proxy.onMouseMove(gc, me);
    proxy.onClick(gc, me);
    proxy.onDoubleClick(gc, me);
    proxy.onKeyDown(gc, ke);
    proxy.onKeyPress(gc, ke);

    expect(flags.enter, true);
    expect(flags.leave, true);
    expect(flags.down, true);
    expect(flags.up, true);
    expect(flags.move, true);
    expect(flags.click, true);
    expect(flags.dblClick, true);
    expect(flags.keyDown, true);
    expect(flags.keyPress, true);
  });

  // ---------------------------------------------------------------------------
  // Test 2 — onMouseLeave clears 2D preview but keeps 3D-sourced preview
  // ---------------------------------------------------------------------------

  // onMouseLeave reads `gridCell.tableColumn`. We build a minimal fake that
  // has semType=MOLECULE and the picker-activation tag set via col.tags[].
  // When _previewFrom3D=false the method calls _removePreviewAtom(), which
  // unconditionally nulls _previewAtomIdx at its first line (before any
  // grok.shell.tv guard). When _previewFrom3D=true the branch is skipped
  // entirely — preview survives.

  test('onMouseLeave-clears-2d-preview-keeps-3d-preview', async () => {
    const r = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    const col = DG.Column.fromStrings('smiles', ['CCO']);
    col.semType = DG.SEMTYPE.MOLECULE;
    col.tags[CHEM_ATOM_PICKER_LINKED_COL] = 'pose3D';

    const fakeGc = {tableColumn: col} as unknown as DG.GridCell;
    const fakeMe = {} as MouseEvent;

    // 2D-sourced preview: onMouseLeave must null _previewAtomIdx.
    r._previewAtomIdx = 5;
    r._previewFrom3D = false;
    r.onMouseLeave(fakeGc, fakeMe);
    expect(r._previewAtomIdx, null);

    // 3D-sourced preview: onMouseLeave must leave _previewAtomIdx untouched.
    r._previewAtomIdx = 7;
    r._previewFrom3D = true;
    r.onMouseLeave(fakeGc, fakeMe);
    expect(r._previewAtomIdx, 7);
  });

  // ---------------------------------------------------------------------------
  // Test 3 — onMouseMove bails immediately when shiftKey is true
  // ---------------------------------------------------------------------------

  // `if (e.shiftKey) return;` — the very first statement in
  // onMouseMove. No column check, no state mutation. Prove the guard works
  // by leaving _previewAtomIdx and _lastHoveredAtom unchanged after a
  // shiftKey move.

  test('onMouseMove-bails-on-shiftKey', async () => {
    const r = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    r._previewAtomIdx = 42;
    r._lastHoveredAtom = null;

    // The shiftKey guard fires before any gridCell field access, so an empty
    // gridCell object is safe here.
    r.onMouseMove({} as DG.GridCell, {shiftKey: true} as MouseEvent);

    expect(r._previewAtomIdx, 42);
    expect(r._lastHoveredAtom, null);
  });

  // ---------------------------------------------------------------------------
  // Test 4 — _trackHoveredAtom dedup accounts for the erase flag
  // ---------------------------------------------------------------------------

  // Calling with the same (col, row, atom) but different `erase` values must
  // be treated as distinct atoms (returns false = "new"). Calling twice with
  // identical args (including erase) returns true = "same, skip".

  test('trackHoveredAtom-dedup-across-erase-flag', async () => {
    const r = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;
    const track = (col: string, row: number, atom: number, erase?: boolean) =>
      r._trackHoveredAtom(col, row, atom, erase);

    // First call: no prior state → new atom.
    expect(track('m', 0, 3, false), false);
    // Same args: deduped → same.
    expect(track('m', 0, 3, false), true);
    // Same col/row/atom but erase flipped → new.
    expect(track('m', 0, 3, true), false);
    // Repeated with erase=true → same.
    expect(track('m', 0, 3, true), true);
  });
});
