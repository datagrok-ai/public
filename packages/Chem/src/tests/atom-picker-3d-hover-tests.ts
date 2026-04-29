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

import {AtomPickerProvider, CHEM_MOL3D_SELECTION_EVENT, CHEM_ATOM_PICKER_LINKED_3D_COL_TAG}
  from '@datagrok-libraries/chem-meta/src/types';
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
    expect(substruct!.atoms!.includes(2), true);
  }, {timeout: 15000});

  // -------------------------------------------------------------------------
  // paint mode — 3D event should add the atom persistently in 2D
  // -------------------------------------------------------------------------

  test('3D hover — paint mode adds atom persistently', async () => {
    smilesCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = [];
    // Paint serial 3 (= 2D O atom, idx 2 in CCO).
    await fire3DHover({atom3DSerial: 3, mode: 'paint'});

    const prov = getPickerProvider(smilesCol, 0);
    expect(prov !== undefined, true);
    // __atoms is the canonical persistent set — paint must have added 2.
    expect(prov!.__atoms?.has(2), true);
  }, {timeout: 15000});

  // -------------------------------------------------------------------------
  // erase mode — 3D event should remove the atom from the persistent set
  // -------------------------------------------------------------------------

  test('3D hover — erase mode removes atom from persistent set', async () => {
    smilesCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = [];
    // First paint atom (serial 3 → 2D atom 2), then erase it via 3D event.
    await fire3DHover({atom3DSerial: 3, mode: 'paint'});
    await fire3DHover({atom3DSerial: 3, mode: 'erase'});

    const prov = getPickerProvider(smilesCol, 0);
    // Last atom removed → provider gone, OR present with __atoms empty.
    if (prov)
      expect(prov.__atoms?.size ?? 0, 0);
  }, {timeout: 15000});

  // -------------------------------------------------------------------------
  // burst preview events — should not accumulate atoms in __atoms
  // -------------------------------------------------------------------------

  // Reproduces the user-reported bug: shift+click in 2D, then 3D hover
  // keeps adding atoms to the persistent set instead of just previewing.
  // After our fix, the persistent set must only contain the originally-
  // painted atom no matter how many preview events fire.
  test('3D hover — burst preview events do not accumulate state', async () => {
    smilesCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = [];
    // Paint atom 0 persistently first (via paint event for serial 1).
    await fire3DHover({atom3DSerial: 1, mode: 'paint'});

    // Burst of preview events on different serials, ending with cursor-leave.
    for (const s of [2, 3, 2, 3, 2, null, 3, null]) {
      grok.events.fireCustomEvent(CHEM_MOL3D_SELECTION_EVENT, {
        mol3DColumnName: mol3DCol.name, rowIdx: 0,
        atom3DSerial: s, mode: 'preview',
      });
    }
    await delay(300);

    // Re-read provider state via fresh col.temp access.
    const prov = getPickerProvider(smilesCol, 0);
    expect(prov !== undefined, true);
    // __atoms must contain ONLY the originally painted atom (idx 0 from
    // serial 1). Preview atoms must not have leaked in.
    const atoms = Array.from(prov!.__atoms ?? []);
    expect(JSON.stringify(atoms), '[0]');
  }, {timeout: 20000});
});
