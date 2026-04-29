import {before, category, test, expect} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {AtomPickerProvider, CHEM_ATOM_PICKER_LINKED_3D_COL_TAG}
  from '@datagrok-libraries/chem-meta/src/types';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';

const PROVIDERS_KEY = ChemTemps.SUBSTRUCT_PROVIDERS;

/** Reads provider state via a fresh `col.temp[KEY]` access — important
 *  because col.temp returns a snapshot proxy on each read. Returns `null`
 *  when no picker provider exists for the row. */
function readProviderFresh(col: DG.Column, rowIdx: number):
    {atoms: number[]; rendered: number[] | null} | null {
  const all = (col.temp[PROVIDERS_KEY] ?? []) as AtomPickerProvider[];
  const p = all.find((x) => x.__atomPicker && x.__rowIdx === rowIdx);
  if (!p) return null;
  return {
    atoms: Array.from(p.__atoms ?? []),
    rendered: p.getSubstruct?.(rowIdx)?.atoms ?? null,
  };
}

function makeLinkedDF(): {smilesCol: DG.Column; mol3DCol: DG.Column} {
  const smilesCol = DG.Column.fromStrings('smiles', ['CCO', 'CCCN']);
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  const mol3DCol = DG.Column.fromStrings('pose3D', ['', '']);
  mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
  DG.DataFrame.fromColumns([smilesCol, mol3DCol]);
  smilesCol.setTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG, mol3DCol.name);
  return {smilesCol, mol3DCol};
}

/** Synthesizes a CCO bondAtoms map (atoms 0-1-2, bonds [0,1] and [1,2]).
 *  The map is what `_addAtomToRow` / `_setPreviewAtom` pass through to
 *  `_updateRowSelection` for bond-highlight computation. The snapshot tests
 *  only check the atom set, so the bond contents don't matter — we just
 *  need a valid Map. Building it inline avoids depending on the renderer's
 *  SVG-based atom-position cache (`_getCellAtomPositions`), which requires
 *  a fully bootstrapped RDKit + DOM layout context. */
function ccoBondAtoms(): Map<number, [number, number]> {
  return new Map<number, [number, number]>([[0, [0, 1]], [1, [1, 2]]]);
}

/** Same as `ccoBondAtoms` but for CCCN (atoms 0-1-2-3, bonds 0-1, 1-2, 2-3). */
function cccnBondAtoms(): Map<number, [number, number]> {
  return new Map<number, [number, number]>([[0, [0, 1]], [1, [1, 2]], [2, [2, 3]]]);
}

category('atom picker', () => {
  let rdkitModule: RDModule;

  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule') as RDModule;
  });

  // Snapshot round-trip tests: after a picker mutation, re-read
  // `col.temp[KEY]` to force a fresh proxy snapshot and verify the change
  // actually persisted. Mutating an object inside a previous snapshot does
  // NOT propagate — see `_setPreviewAtom` / `_updateRowSelection` for the
  // canonical write pattern this guards.

  test('preview atom does NOT leak into __atoms (snapshot round-trip)', async () => {
    const {smilesCol} = makeLinkedDF();
    const renderer = new RDKitCellRenderer(rdkitModule);
    const bondAtoms = ccoBondAtoms();

    // Paint atom 0 (persistent), then set a preview atom 2.
    renderer.picker._addAtomToRow(smilesCol, 0, 0, bondAtoms);
    renderer.picker._setPreviewAtom(smilesCol, 0, 2, bondAtoms);

    // Fresh snapshot read — the bug we caught was that the preview atom
    // leaked into __atoms because the post-_updateRowSelection mutation was
    // performed on a stale snapshot reference.
    const snap = readProviderFresh(smilesCol, 0);
    expect(snap !== null, true);
    expect(JSON.stringify(snap!.atoms), '[0]');
    expect(JSON.stringify(snap!.rendered), '[0,2]');
  });

  test('erase syncs to 2D rendering set immediately (snapshot round-trip)', async () => {
    const {smilesCol} = makeLinkedDF();
    const renderer = new RDKitCellRenderer(rdkitModule);
    const bondAtoms = ccoBondAtoms();

    // Paint [0, 1, 2], then erase atom 1.
    renderer.picker._addAtomToRow(smilesCol, 0, 0, bondAtoms);
    renderer.picker._addAtomToRow(smilesCol, 0, 1, bondAtoms);
    renderer.picker._addAtomToRow(smilesCol, 0, 2, bondAtoms);
    renderer.picker._removeAtomFromRow(smilesCol, 0, 1, bondAtoms);

    // Both __atoms (canonical state) AND getSubstruct().atoms (what 2D
    // renders) must reflect the reduced set on a fresh read. The earlier
    // bug had 3D updating but 2D's getSubstruct closure pointing at the
    // stale unreduced set.
    const snap = readProviderFresh(smilesCol, 0);
    expect(JSON.stringify(snap!.atoms), '[0,2]');
    expect(JSON.stringify(snap!.rendered), '[0,2]');
  });

  test('painting different rows is independent', async () => {
    const {smilesCol} = makeLinkedDF();
    const renderer = new RDKitCellRenderer(rdkitModule);
    const baR0 = ccoBondAtoms();
    const baR1 = cccnBondAtoms();

    renderer.picker._addAtomToRow(smilesCol, 0, 0, baR0);
    renderer.picker._addAtomToRow(smilesCol, 1, 3, baR1);

    // Mutate row 0 only (add atom 2). Row 1 must be untouched.
    renderer.picker._addAtomToRow(smilesCol, 0, 2, baR0);

    const r0 = readProviderFresh(smilesCol, 0);
    const r1 = readProviderFresh(smilesCol, 1);
    expect(JSON.stringify(r0!.atoms), '[0,2]');
    expect(JSON.stringify(r1!.atoms), '[3]');

    // Erase from row 0 — row 1 still untouched.
    renderer.picker._removeAtomFromRow(smilesCol, 0, 0, baR0);
    const r0After = readProviderFresh(smilesCol, 0);
    const r1After = readProviderFresh(smilesCol, 1);
    expect(JSON.stringify(r0After!.atoms), '[2]');
    expect(JSON.stringify(r1After!.atoms), '[3]');
  });
});
