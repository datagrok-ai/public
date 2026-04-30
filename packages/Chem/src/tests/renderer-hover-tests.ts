import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {CHEM_ATOM_PICKER_LINKED_3D_COL_TAG} from '@datagrok-libraries/chem-meta/src/types';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';

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
}

category('renderer hover', () => {
  let rdkitModule: RDModule;

  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule') as RDModule;
  });

  // onMouseLeave reads `gridCell.tableColumn`. We build a minimal fake that
  // has semType=MOLECULE and the picker-activation tag set via col.tags[].
  // When _previewFrom3D=false the method calls _removePreviewAtom(), which
  // unconditionally nulls _previewAtomIdx at its first line (before any
  // grok.shell.tv guard). When _previewFrom3D=true the branch is skipped
  // entirely — preview survives.

  test('onMouseLeave — clears 2D preview, keeps 3D preview', async () => {
    const r = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    const col = DG.Column.fromStrings('smiles', ['CCO']);
    col.semType = DG.SEMTYPE.MOLECULE;
    col.tags[CHEM_ATOM_PICKER_LINKED_3D_COL_TAG] = 'pose3D';

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

  // Calling with the same (col, row, atom) but different `erase` values must
  // be treated as distinct atoms (returns false = "new"). Calling twice with
  // identical args (including erase) returns true = "same, skip".

  test('trackHoveredAtom — dedup across erase flag', async () => {
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
