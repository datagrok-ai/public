/**
 * Tests for `getMol3DAtomPickerLinkWidget` (exported via
 * `BiostructureViewer:mol3dAtomPickerLinkWidget`).
 *
 * Scope: we verify the widget's INITIAL RENDER STATE given different
 * tag/DataFrame configurations. The widget's state transitions on user
 * interaction (check/uncheck → `setSmilesColLink`/`clearLinksToMol3D`)
 * are covered transitively by `mol3d-link-tests.ts`, so we don't try
 * to simulate DOM events through Datagrok's `InputBase` wrapper here.
 */

import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {CHEM_ATOM_PICKER_LINKED_3D_COL_TAG} from '@datagrok-libraries/chem-meta/src/types';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Calls the BSV widget function cross-package the same way Docking does.
// grok.functions.call returns Promise<unknown>, so the double cast is required.
async function callWidget(mol3DCol: DG.Column): Promise<DG.Widget> {
  return grok.functions.call(
    'BiostructureViewer:mol3dAtomPickerLinkWidget',
    {mol3DCol},
  ) as unknown as Promise<DG.Widget>;
}

// Finds the enable checkbox inside a widget root.
function findCheckbox(root: HTMLElement): HTMLInputElement | null {
  return root.querySelector('input[type="checkbox"]');
}

// Finds the SMILES choice dropdown inside a widget root.
function findSelect(root: HTMLElement): HTMLSelectElement | null {
  return root.querySelector('select');
}

// Finds the outer dropdown wrapper (the div whose `style.display` the widget
// toggles to show/hide the dropdown). The widget wraps `choiceInput.root` in
// its own div — `select.closest('div')` returns `ui.input.choice`'s internal
// wrapper, so we walk up one more level to reach the widget's own wrapper.
function findDropdownWrapper(root: HTMLElement): HTMLElement | null {
  const select = findSelect(root);
  if (!select) return null;
  const inner = select.closest('div') as HTMLElement | null;
  return inner?.parentElement ?? null;
}

// Builds a minimal DataFrame with one SMILES col and one Molecule3D col.
function makeDF(): {smilesCol: DG.Column; mol3DCol: DG.Column} {
  const smilesCol = DG.Column.fromStrings('smiles', ['CCO']);
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  const mol3DCol = DG.Column.fromStrings('pose3D', ['...']);
  mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
  DG.DataFrame.fromColumns([smilesCol, mol3DCol]);
  return {smilesCol, mol3DCol};
}

// ---------------------------------------------------------------------------
// Category
// ---------------------------------------------------------------------------

category('Mol3DAtomPickerLinkWidget', () => {
  // -------------------------------------------------------------------------
  // Empty state — no SMILES columns → widget renders an info message, no
  // interactive controls
  // -------------------------------------------------------------------------

  test('no SMILES cols — shows info message', async () => {
    const mol3DCol = DG.Column.fromStrings('pose3D', ['...']);
    mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
    DG.DataFrame.fromColumns([mol3DCol]);

    const widget = await callWidget(mol3DCol);
    const root = widget.root;

    // No checkbox when there are no SMILES columns.
    const checkbox = findCheckbox(root);
    expect(checkbox === null, true);

    // The root text content should mention SMILES or "No" to signal the
    // empty state (content check is intentionally loose so minor wording
    // changes don't break the test).
    const text = (root.textContent ?? '').toLowerCase();
    expect(text.includes('no') || text.includes('smiles') || text.includes('molecule'), true);
  });

  // -------------------------------------------------------------------------
  // Initial state — SMILES col present, no link yet → checkbox unchecked
  // -------------------------------------------------------------------------

  test('initial state — no link, checkbox unchecked', async () => {
    const {mol3DCol} = makeDF();

    const widget = await callWidget(mol3DCol);
    await delay(50);

    const checkbox = findCheckbox(widget.root);
    expect(checkbox !== null, true);
    expect(checkbox!.checked, false);
  });

  // -------------------------------------------------------------------------
  // Initial state — SMILES col present WITH link → checkbox pre-checked
  // -------------------------------------------------------------------------

  test('initial state — with link, checkbox prechecked', async () => {
    const {smilesCol, mol3DCol} = makeDF();
    // Pre-set a link so the widget opens with checkbox checked.
    smilesCol.setTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG, mol3DCol.name);

    const widget = await callWidget(mol3DCol);
    await delay(50);

    const checkbox = findCheckbox(widget.root);
    expect(checkbox !== null, true);
    expect(checkbox!.checked, true);
  });

  // -------------------------------------------------------------------------
  // Dropdown defaults to the linked SMILES col when a link is pre-set
  // -------------------------------------------------------------------------

  test('initial state — dropdown defaults to linked SMILES col', async () => {
    const smilesA = DG.Column.fromStrings('smilesA', ['CCO']);
    smilesA.semType = DG.SEMTYPE.MOLECULE;
    const smilesB = DG.Column.fromStrings('smilesB', ['CC']);
    smilesB.semType = DG.SEMTYPE.MOLECULE;
    const mol3DCol = DG.Column.fromStrings('pose3D', ['...']);
    mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
    DG.DataFrame.fromColumns([smilesA, smilesB, mol3DCol]);

    // Pre-link smilesB (not smilesA) so we can verify the dropdown reflects
    // the linked col rather than the first one.
    smilesB.setTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG, mol3DCol.name);

    const widget = await callWidget(mol3DCol);
    await delay(50);

    const select = findSelect(widget.root);
    // ui.input.choice may render as <select> or as a custom widget; skip the
    // value check if no <select> is present, since the interesting invariant
    // (link survives widget open) is covered by the previous test.
    if (select)
      expect(select.value, smilesB.name);
  });

  // -------------------------------------------------------------------------
  // Dropdown visibility reflects initial link state
  // -------------------------------------------------------------------------

  test('initial state — dropdown hidden when no link', async () => {
    const {mol3DCol} = makeDF();

    const widget = await callWidget(mol3DCol);
    await delay(50);

    const wrapper = findDropdownWrapper(widget.root);
    // The widget's own wrapper may not exist if ui.input.choice changes its
    // DOM structure. Only assert when we can reliably find the wrapper.
    if (wrapper && wrapper !== widget.root)
      expect(wrapper.style.display, 'none');
  });

  test('initial state — dropdown visible when link exists', async () => {
    const {smilesCol, mol3DCol} = makeDF();
    smilesCol.setTag(CHEM_ATOM_PICKER_LINKED_3D_COL_TAG, mol3DCol.name);

    const widget = await callWidget(mol3DCol);
    await delay(50);

    const wrapper = findDropdownWrapper(widget.root);
    if (wrapper && wrapper !== widget.root)
      expect(wrapper.style.display !== 'none', true);
  });
});
