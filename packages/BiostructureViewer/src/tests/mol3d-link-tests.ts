/* eslint-disable max-len */
/**
 * Unit tests for `utils/mol3d-link.ts` helpers.
 *
 * These are pure column-tag operations: no Molstar viewer, no RDKit, no live
 * server interaction needed beyond the Datagrok column/DataFrame API.
 * Each test constructs its own DataFrame inline so tests are fully isolated.
 */

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  setSmilesColLink,
  clearLinksToMol3D,
  findLinkedSmilesColName,
  clearAtomPickerHighlights,
} from '../utils/mol3d-link';

import {CHEM_ATOM_PICKER_LINKED_3D_COL_TAG,
  CHEM_ATOM_PICKER_LINKED_SMILES_COL,} from '@datagrok-libraries/chem-meta/src/types';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Key where Chem's cell renderer stores atom-picker providers in col.temp.
 *  Mirrors ChemTemps.SUBSTRUCT_PROVIDERS from @datagrok-libraries/chem-meta. */
const PROVIDERS_KEY = 'substruct-providers';

/** Minimal provider shape used in these tests. */
type TestProvider = {__atomPicker?: boolean; __scaffold?: boolean; __ordinaryMarker?: boolean; getSubstruct: () => undefined};

// Builds a minimal DataFrame with one SMILES col and one Molecule3D col.
function makeDF(): {df: DG.DataFrame; smilesCol: DG.Column; mol3DCol: DG.Column} {
  const smilesCol = DG.Column.fromStrings('smiles', ['CCO', 'CC']);
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  const mol3DCol = DG.Column.fromStrings('pose3D', ['...', '...']);
  mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
  const df = DG.DataFrame.fromColumns([smilesCol, mol3DCol]);
  return {df, smilesCol, mol3DCol};
}

// ---------------------------------------------------------------------------
// setSmilesColLink
// ---------------------------------------------------------------------------

category('Mol3DLink', () => {
  test('setSmilesColLink — writes tag', async () => {
    const {smilesCol, mol3DCol} = makeDF();
    setSmilesColLink(smilesCol, mol3DCol.name);
    const tag = smilesCol.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG);
    expect(tag, mol3DCol.name);
  });

  test('setSmilesColLink — null removes tag', async () => {
    const {smilesCol, mol3DCol} = makeDF();
    setSmilesColLink(smilesCol, mol3DCol.name);
    // Sanity: tag was set.
    expect(smilesCol.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG) !== null, true);

    setSmilesColLink(smilesCol, null);
    const tagAfter = smilesCol.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG);
    // Tag must be absent (null or empty string) after clearing.
    expect(!tagAfter, true);
  });

  test('setSmilesColLink — overwrites previous', async () => {
    const {smilesCol} = makeDF();
    setSmilesColLink(smilesCol, 'pose_A');
    setSmilesColLink(smilesCol, 'pose_B');
    expect(smilesCol.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG), 'pose_B');
  });

  // -------------------------------------------------------------------------
  // Two-way tag pair
  // -------------------------------------------------------------------------

  test('setSmilesColLink — writes reciprocal tag on Mol3D', async () => {
    const {smilesCol, mol3DCol} = makeDF();
    setSmilesColLink(smilesCol, mol3DCol.name);
    // Forward tag on SMILES col.
    expect(smilesCol.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG), mol3DCol.name);
    // Reciprocal tag on Mol3D col.
    expect(mol3DCol.getTag(CHEM_ATOM_PICKER_LINKED_SMILES_COL), smilesCol.name);
  });

  test('setSmilesColLink — null clears both tags', async () => {
    const {smilesCol, mol3DCol} = makeDF();
    setSmilesColLink(smilesCol, mol3DCol.name);
    setSmilesColLink(smilesCol, null);
    expect(!smilesCol.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG), true);
    expect(!mol3DCol.getTag(CHEM_ATOM_PICKER_LINKED_SMILES_COL), true);
  });

  test('setSmilesColLink — overwrite strips old Mol3D reciprocal', async () => {
    // When a SMILES col's link moves from mol3DA to mol3DB, mol3DA should
    // lose its reciprocal tag — otherwise a stale back-link leaks.
    const smilesCol = DG.Column.fromStrings('smiles', ['CCO']);
    smilesCol.semType = DG.SEMTYPE.MOLECULE;
    const mol3DA = DG.Column.fromStrings('pose_A', ['...']);
    mol3DA.semType = DG.SEMTYPE.MOLECULE3D;
    const mol3DB = DG.Column.fromStrings('pose_B', ['...']);
    mol3DB.semType = DG.SEMTYPE.MOLECULE3D;
    DG.DataFrame.fromColumns([smilesCol, mol3DA, mol3DB]);

    setSmilesColLink(smilesCol, mol3DA.name);
    setSmilesColLink(smilesCol, mol3DB.name);

    expect(!mol3DA.getTag(CHEM_ATOM_PICKER_LINKED_SMILES_COL), true);
    expect(mol3DB.getTag(CHEM_ATOM_PICKER_LINKED_SMILES_COL), smilesCol.name);
  });

  // -------------------------------------------------------------------------
  // clearLinksToMol3D
  // -------------------------------------------------------------------------

  test('clearLinksToMol3D — removes link for target col', async () => {
    const {df, smilesCol, mol3DCol} = makeDF();
    setSmilesColLink(smilesCol, mol3DCol.name);

    clearLinksToMol3D(df, mol3DCol.name);

    // Both forward + reciprocal tags must be gone after clearing.
    expect(!smilesCol.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG), true);
    expect(!mol3DCol.getTag(CHEM_ATOM_PICKER_LINKED_SMILES_COL), true);
  });

  test('clearLinksToMol3D — leaves unrelated links intact', async () => {
    // Two SMILES cols; only one is linked to 'pose3D', the other to 'poseX'.
    const smilesA = DG.Column.fromStrings('smilesA', ['CCO']);
    smilesA.semType = DG.SEMTYPE.MOLECULE;
    const smilesB = DG.Column.fromStrings('smilesB', ['CC']);
    smilesB.semType = DG.SEMTYPE.MOLECULE;
    const mol3DA = DG.Column.fromStrings('pose3D', ['...']);
    mol3DA.semType = DG.SEMTYPE.MOLECULE3D;
    const mol3DB = DG.Column.fromStrings('poseX', ['...']);
    mol3DB.semType = DG.SEMTYPE.MOLECULE3D;
    const df = DG.DataFrame.fromColumns([smilesA, smilesB, mol3DA, mol3DB]);

    setSmilesColLink(smilesA, mol3DA.name); // linked to 'pose3D'
    setSmilesColLink(smilesB, mol3DB.name); // linked to 'poseX'

    clearLinksToMol3D(df, mol3DA.name); // clear only 'pose3D' links

    // smilesA tag must be gone.
    expect(!smilesA.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG), true);
    // smilesB tag for 'poseX' must be untouched.
    expect(smilesB.getTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG), mol3DB.name);
  });

  // -------------------------------------------------------------------------
  // findLinkedSmilesColName
  // -------------------------------------------------------------------------

  test('findLinkedSmilesColName — returns linked SMILES col', async () => {
    const {df, smilesCol, mol3DCol} = makeDF();
    setSmilesColLink(smilesCol, mol3DCol.name);

    const found = findLinkedSmilesColName(df, mol3DCol.name);
    expect(found, smilesCol.name);
  });

  test('findLinkedSmilesColName — returns null when no link', async () => {
    const {df, mol3DCol} = makeDF();
    const found = findLinkedSmilesColName(df, mol3DCol.name);
    expect(found === null, true);
  });

  test('findLinkedSmilesColName — returns null for unknown 3D col', async () => {
    const {df, smilesCol, mol3DCol} = makeDF();
    setSmilesColLink(smilesCol, mol3DCol.name);

    const found = findLinkedSmilesColName(df, 'does-not-exist');
    expect(found === null, true);
  });

  // -------------------------------------------------------------------------
  // clearAtomPickerHighlights
  // -------------------------------------------------------------------------

  test('clearAtomPickerHighlights — strips atom-picker providers', async () => {
    const {smilesCol} = makeDF();

    // Simulate two providers: one atom-picker (tagged), one ordinary
    // (marked with `__ordinaryMarker` so we can identify it after the round
    // trip through `col.temp`, which does not preserve object identity).
    const atomPickerProv: TestProvider = {__atomPicker: true, getSubstruct: () => undefined};
    const ordinaryProv: TestProvider = {__ordinaryMarker: true, getSubstruct: () => undefined};
    smilesCol.temp[PROVIDERS_KEY] = [atomPickerProv, ordinaryProv];

    clearAtomPickerHighlights(smilesCol);

    const remaining = (smilesCol.temp[PROVIDERS_KEY] ?? []) as TestProvider[];
    // The ordinary provider must survive.
    expect(remaining.some((p) => p.__ordinaryMarker === true), true);
    // The atom-picker provider must be gone.
    expect(remaining.some((p) => p.__atomPicker === true), false);
  });

  test('clearAtomPickerHighlights — fires clear-all event', async () => {
    const {smilesCol} = makeDF();
    const CHEM_SELECTION_EVENT = 'chem-interactive-selection-changed';

    let receivedClearAll = false;
    const sub = grok.events.onCustomEvent(CHEM_SELECTION_EVENT)
      .subscribe((args: {clearAll?: boolean}) => {
        if (args?.clearAll)
          receivedClearAll = true;
      });

    try {
      clearAtomPickerHighlights(smilesCol);
      // Allow a generous window for RxJS subscriber scheduling to deliver
      // the event across the Datagrok event bus.
      await new Promise<void>((r) => setTimeout(r, 250));
      expect(receivedClearAll, true);
    } finally {
      sub.unsubscribe();
    }
  });

  test('clearAtomPickerHighlights — preserves non-picker providers', async () => {
    const {smilesCol} = makeDF();

    const scaffoldProv: TestProvider = {__scaffold: true, getSubstruct: () => undefined};
    const pickerProv: TestProvider = {__atomPicker: true, getSubstruct: () => undefined};
    smilesCol.temp[PROVIDERS_KEY] = [scaffoldProv, pickerProv];

    clearAtomPickerHighlights(smilesCol);

    const remaining = (smilesCol.temp[PROVIDERS_KEY] ?? []) as TestProvider[];
    expect(remaining.length, 1);
    // `col.temp` round-trips values as data, so check by marker rather than
    // object identity — the surviving provider is the scaffold one.
    expect(remaining[0].__scaffold === true, true);
  });
});
