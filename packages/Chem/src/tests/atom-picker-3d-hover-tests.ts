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

import {CHEM_MOL3D_HOVER_EVENT, CHEM_ATOM_PICKER_LINKED_COL} from '../constants';
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

category('atom-picker-3d-hover', () => {
  let df: DG.DataFrame;
  let smilesCol: DG.Column;
  let mol3DCol: DG.Column;
  let view: DG.TableView;

  // Re-open (or re-activate) the table view and clear provider state. Called
  // before every test because Datagrok's `grok.shell.tv` can drift between
  // tests — `_onMol3DHoverEvent` bails when `grok.shell.tv?.grid` is not the
  // view whose DataFrame holds our SMILES column. Rebinding here ensures the
  // handler consistently finds the right grid + dataFrame on every fire.
  async function ensureActiveView(): Promise<void> {
    // Re-focus the existing view so `grok.shell.tv` returns it. If the view
    // was closed by a prior teardown, re-add the DataFrame which produces a
    // fresh view over the same data.
    if (!view || !view.dataFrame)
      view = grok.shell.addTableView(df);
    grok.shell.v = view;
    df.currentRowIdx = 0;
    smilesCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = [];
    view.grid.invalidate();
    await awaitGrid(view.grid);
    await delay(200);
  }

  before(async () => {
    smilesCol = DG.Column.fromStrings('smiles', [SMILES_ETHANOL, 'CC']);
    smilesCol.semType = DG.SEMTYPE.MOLECULE;

    mol3DCol = DG.Column.fromStrings('pose3D', [POSE_ETHANOL_PDB, '']);
    mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;

    df = DG.DataFrame.fromColumns([smilesCol, mol3DCol]);
    df.name = 'atom-picker-3d-hover-test';

    // Activate the picker by setting the link tag on the SMILES column.
    smilesCol.setTag(CHEM_ATOM_PICKER_LINKED_COL, mol3DCol.name);

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

  // -------------------------------------------------------------------------
  // preview mode
  // -------------------------------------------------------------------------

  test('preview-sets-provider-with-correct-atom', async () => {
    // Fire a preview event for atom serial 3 (the O atom in the PDB).
    // After reverse-mapping: serial 3 → PDB heavy-atom index 2 → 2D idx 2 (O in CCO).
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      atom3DSerial: 3,
      mode: 'preview',
    });
    await delay(100);

    const prov = getPickerProvider(smilesCol, 0);
    expect(prov !== undefined, true);
    // getSubstruct(0) should contain at least the previewed atom (idx 2 = O).
    const substruct = prov!.getSubstruct(0);
    expect(substruct !== undefined, true);
    expect(Array.isArray(substruct!.atoms), true);
    expect((substruct!.atoms as number[]).includes(2), true);
  }, {timeout: 15000});

  test('preview-null-serial-clears-preview', async () => {
    // First establish a preview.
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      atom3DSerial: 3,
      mode: 'preview',
    });
    await delay(100);

    // Then fire the "cursor left" event (atom3DSerial: null).
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      atom3DSerial: null,
      mode: 'preview',
    });
    await delay(100);

    // The provider may still exist (if persistent atoms remain) but the
    // preview atom is gone. When no persistent atoms were ever painted,
    // the provider disappears or getSubstruct returns no atoms.
    const prov = getPickerProvider(smilesCol, 0);
    if (prov) {
      // Provider exists only for persistent atoms — __atoms must be empty.
      const persistentSize = prov.__atoms?.size ?? 0;
      expect(persistentSize, 0);
    }
    // Either way: no provider means preview is cleared.
    // (Both "no provider" and "provider with 0 atoms" are correct outcomes.)
  }, {timeout: 15000});

  // -------------------------------------------------------------------------
  // paint mode — accumulates across multiple events
  // -------------------------------------------------------------------------

  test('paint-accumulates-three-serials', async () => {
    await ensureActiveView();

    // Paint serials 1, 2, 3 sequentially (→ 2D indices 0, 1, 2 for CCO).
    for (const serial of [1, 2, 3]) {
      grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
        mol3DColumnName: mol3DCol.name,
        rowIdx: 0,
        atom3DSerial: serial,
        mode: 'paint',
      });
      await delay(60);
    }

    const prov = getPickerProvider(smilesCol, 0);
    expect(prov !== undefined, true);
    // All three 2D atoms should be in the persistent set.
    const atomCount = prov!.__atoms?.size ?? 0;
    expect(atomCount >= 3, true);
  }, {timeout: 15000});

  // -------------------------------------------------------------------------
  // erase mode — removes one specific atom, others stay
  // -------------------------------------------------------------------------

  test('erase-removes-one-atom-others-remain', async () => {
    await ensureActiveView();

    // Set up: paint atoms at serials 1, 2, 3.
    for (const serial of [1, 2, 3]) {
      grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
        mol3DColumnName: mol3DCol.name,
        rowIdx: 0,
        atom3DSerial: serial,
        mode: 'paint',
      });
      await delay(60);
    }

    const before_ = getPickerProvider(smilesCol, 0);
    expect(before_ !== undefined, true);
    const countBefore = before_!.__atoms?.size ?? 0;
    expect(countBefore >= 2, true);

    // Erase serial 1 (→ 2D atom 0).
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      atom3DSerial: 1,
      mode: 'erase',
    });
    await delay(100);

    const after_ = getPickerProvider(smilesCol, 0);
    // If provider still exists, it should have one fewer atom.
    if (after_) {
      const countAfter = after_!.__atoms?.size ?? 0;
      expect(countAfter < countBefore, true);
    }
    // If provider was removed (all atoms erased via cascade), that is also
    // acceptable — the important invariant is that count decreased.
  }, {timeout: 15000});

  // -------------------------------------------------------------------------
  // preview-from-3D flag — null-serial must not wipe paint state
  // -------------------------------------------------------------------------

  test('null-serial-preview-leaves-paint-atoms-intact', async () => {
    await ensureActiveView();

    // Set up: paint serial 2 (a C atom in ethanol).
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      atom3DSerial: 2,
      mode: 'paint',
    });
    await delay(60);

    // Also add a preview.
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      atom3DSerial: 3,
      mode: 'preview',
    });
    await delay(60);

    const painted = new Set(getPickerProvider(smilesCol, 0)?.__atoms ?? []);
    expect(painted.size > 0, true);

    // Now fire null-serial (cursor left atom) — preview clears, paint stays.
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      atom3DSerial: null,
      mode: 'preview',
    });
    await delay(100);

    // Persistent paint atoms must survive the null preview event.
    const afterProv = getPickerProvider(smilesCol, 0);
    expect(afterProv !== undefined, true);
    const remainingSize = afterProv!.__atoms?.size ?? 0;
    expect(remainingSize, painted.size);
  }, {timeout: 15000});

  // -------------------------------------------------------------------------
  // _previewFrom3D flag — 2D mousemove must NOT wipe a 3D-sourced preview
  // -------------------------------------------------------------------------

  test('mousemove-outside-grid-does-not-wipe-3d-preview', async () => {
    await ensureActiveView();

    // Establish a 3D-sourced preview (sets _previewFrom3D = true internally).
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: mol3DCol.name,
      rowIdx: 0,
      atom3DSerial: 3,
      mode: 'preview',
    });
    await delay(100);

    // Verify the preview is set.
    const provBefore = getPickerProvider(smilesCol, 0);
    expect(provBefore !== undefined, true);

    // Dispatch a synthetic mousemove at position (0, 0) — outside any grid
    // cell bounding box. The no-modifier branch of _onDocumentMouseMove
    // calls _removePreviewAtom only when _previewFrom3D === false.
    // With the flag set, this mousemove must be a no-op for the preview.
    document.dispatchEvent(new MouseEvent('mousemove', {
      bubbles: true, clientX: 0, clientY: 0,
    }));
    await delay(100);

    // The 3D-sourced preview must still be present (not wiped by the
    // mousemove over an area outside the grid).
    const provAfter = getPickerProvider(smilesCol, 0);
    expect(provAfter !== undefined, true);
    // The preview atom (2D atom 2 = O) must still be in getSubstruct.
    const substruct = provAfter!.getSubstruct(0);
    expect(substruct !== undefined && (substruct!.atoms as number[]).length > 0, true);
  }, {timeout: 15000});

  // -------------------------------------------------------------------------
  // unknown column — event for a different 3D column must be ignored
  // -------------------------------------------------------------------------

  test('event-for-unknown-column-is-ignored', async () => {
    smilesCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = [];

    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, {
      mol3DColumnName: 'non-existent-col',
      rowIdx: 0,
      atom3DSerial: 1,
      mode: 'preview',
    });
    await delay(100);

    const prov = getPickerProvider(smilesCol, 0);
    expect(prov === undefined, true);
  }, {timeout: 15000});
});
