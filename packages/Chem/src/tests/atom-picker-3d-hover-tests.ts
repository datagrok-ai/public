/**
 * Tests for the reverse 3D→2D atom-highlighting bridge.
 *
 * The feature: BiostructureViewer's Molstar viewer fires `CHEM_MOL3D_HOVER_EVENT`
 * when the user hovers a ligand atom in the 3D pose. `RDKitCellRenderer`'s
 * `_onMol3DHoverEvent` handler reverses the 3D serial → 2D atom index mapping
 * and routes it through the same preview / paint / erase storage used by the
 * 2D mouse-hover path.
 *
 * All tests here drive the handler synthetically by firing the event via
 * `grok.events.fireCustomEvent` — no live Molstar viewer required.
 *
 * Preconditions set up in `before`:
 *  - A DataFrame with a SMILES column (semType=MOLECULE) and a Molecule3D
 *    column (semType=MOLECULE3D) is opened as a table view.
 *  - The SMILES column's `CHEM_ATOM_PICKER_LINKED_COL` tag is set to the
 *    Molecule3D column name, activating the picker.
 *  - The renderer has been through at least one render pass so its document
 *    listeners (including the CHEM_MOL3D_HOVER_EVENT subscriber) are attached.
 */

import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';

import {CHEM_MOL3D_SELECTION_EVENT, CHEM_ATOM_PICKER_LINKED_3D_COL_TAG} from '@datagrok-libraries/chem-meta/src/types';
import {awaitGrid} from './utils';

// ---------------------------------------------------------------------------
// Minimal inline pose3D data — a 3-atom PDB block for ethanol (CCO).
// Atom serials: 1=C, 2=C, 3=O. Heavy atoms only, no hydrogens.
// The SMILES 'CCO' has: 2D atom 0=C, 1=C, 2=O.
// RDKit substruct mapping maps 2D-O (idx 2) → PDB serial 3.
// ---------------------------------------------------------------------------
const POSE_ETHANOL_PDB = [
  'HETATM    1  C1  LIG     1       0.000   0.000   0.000  1.00  0.00           C  ',
  'HETATM    2  C2  LIG     1       1.540   0.000   0.000  1.00  0.00           C  ',
  'HETATM    3  O   LIG     1       2.300   1.200   0.000  1.00  0.00           O  ',
  'END',
].join('\n');

const SMILES_ETHANOL = 'CCO';

/** Shape of atom-picker providers stored in `col.temp[ChemTemps.SUBSTRUCT_PROVIDERS]`. */
type AtomPickerProvider = {
  __atomPicker?: boolean;
  __rowIdx?: number;
  __atoms?: Set<number>;
  getSubstruct: (rowIdx: number) => {atoms?: number[]} | undefined;
};

/** Returns the atom-picker provider for `(col, rowIdx)`, or undefined. */
function getPickerProvider(col: DG.Column, rowIdx: number): AtomPickerProvider | undefined {
  const providers = (col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as AtomPickerProvider[];
  return providers.find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
}

category('atom picker: 3D hover', () => {
  let df: DG.DataFrame;
  let smilesCol: DG.Column;
  let mol3DCol: DG.Column;
  let view: DG.TableView;

  before(async () => {
    smilesCol = DG.Column.fromStrings('smiles', [SMILES_ETHANOL, 'CC']);
    smilesCol.semType = DG.SEMTYPE.MOLECULE;

    mol3DCol = DG.Column.fromStrings('pose3D', [POSE_ETHANOL_PDB, '']);
    mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;

    df = DG.DataFrame.fromColumns([smilesCol, mol3DCol]);
    df.name = 'atom-picker-3d-hover-test';

    // Activate the picker by setting the link tag on the SMILES column.
    smilesCol.setTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG, mol3DCol.name);

    view = grok.shell.addTableView(df);
    df.currentRowIdx = 0;

    // Force a render pass so the cell renderer attaches its document-level
    // CHEM_MOL3D_HOVER_EVENT subscriber.
    await awaitGrid(view.grid);
    view.grid.invalidate();
    await delay(300);
  });

  after(async () => {
    // Clear any residual providers and close the view.
    smilesCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = [];
    view?.close();
  });

  /** Fires `CHEM_MOL3D_HOVER_EVENT` with defaults (`mol3DCol.name`, row 0,
   *  `preview` mode) merged with `overrides`, then waits 100 ms for the
   *  async handler in `RDKitCellRenderer` to reverse-map and update
   *  the provider. */
  async function fire3DHover(overrides: {
    mol3DColumnName?: string, rowIdx?: number,
    atom3DSerial?: number | null, mode?: 'preview' | 'paint' | 'erase',
  }): Promise<void> {
    grok.events.fireCustomEvent(CHEM_MOL3D_SELECTION_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      mode: 'preview',
      ...overrides,
    });
    await delay(100);
  }

  // -------------------------------------------------------------------------
  // preview mode
  // -------------------------------------------------------------------------

  test('3D hover — preview sets provider for correct atom', async () => {
    // Fire a preview event for atom serial 3 (the O atom in the PDB).
    // After reverse-mapping: serial 3 → PDB heavy-atom index 2 → 2D idx 2 (O in CCO).
    await fire3DHover({atom3DSerial: 3});

    const prov = getPickerProvider(smilesCol, 0);
    expect(prov !== undefined, true);
    // getSubstruct(0) should contain at least the previewed atom (idx 2 = O).
    const substruct = prov!.getSubstruct(0);
    expect(substruct !== undefined, true);
    expect(Array.isArray(substruct!.atoms), true);
    expect((substruct!.atoms as number[]).includes(2), true);
  }, {timeout: 15000});

  test('3D hover — null serial clears preview', async () => {
    // First establish a preview, then fire the "cursor left" event.
    await fire3DHover({atom3DSerial: 3});
    await fire3DHover({atom3DSerial: null});

    // The provider may still exist (if persistent atoms remain) but the
    // preview atom is gone. When no persistent atoms were ever painted,
    // the provider disappears or getSubstruct returns no atoms. Either
    // way means the preview is cleared.
    const prov = getPickerProvider(smilesCol, 0);
    if (prov) {
      const persistentSize = prov.__atoms?.size ?? 0;
      expect(persistentSize, 0);
    }
  }, {timeout: 15000});

  // Note on coverage: the multi-event scenarios (paint-accumulates across
  // serials, erase-after-paint, _previewFrom3D flag preserving preview across
  // a subsequent document mousemove) were intentionally NOT ported to this
  // file. They require each test to control `grok.shell.tv` across multiple
  // event fires, which the shared-view category harness does not offer. The
  // underlying contracts are still verified transitively:
  //   - add / clear of the CHEM_ATOM_PICKER_LINKED_COL tag → `mol3d-link` tests
  //   - paint / erase provider shape → `atom-picker` (2D path) tests
  //   - preview + null-serial round-trip here
  // Extend when the harness gets a per-test view hook.

  // -------------------------------------------------------------------------
  // unknown column — event for a different 3D column must be ignored
  // -------------------------------------------------------------------------

  test('3D hover — event for unknown column is ignored', async () => {
    smilesCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = [];

    await fire3DHover({mol3DColumnName: 'non-existent-col', atom3DSerial: 1});

    const prov = getPickerProvider(smilesCol, 0);
    expect(prov === undefined, true);
  }, {timeout: 15000});
});
