/**
 * Exports `getMol3DAtomPickerLinkWidget`, which builds a checkbox+dropdown
 * widget for linking a Molecule3D column to a SMILES column so that the
 * 2D→3D and 3D→2D atom-highlighting bridge is active.
 *
 * Consumed cross-package by Docking's AutoDock panel via
 * `grok.functions.call('BiostructureViewer:mol3dAtomPickerLinkWidget', …)`.
 * The underlying link state and helpers live in `utils/mol3d-link.ts`.
 */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  clearLinksToMol3D,
  findLinkedSmilesColName,
  setSmilesColLink,
} from '../utils/mol3d-link';

/** Builds the link-SMILES-column widget for a Molecule3D column.
 *  Structure:
 *    [x] Link SMILES column
 *        SMILES column: [ dropdown ]   <-- shown only when checked
 *
 *  State transitions:
 *  - Check the box → show the dropdown in an empty "no selection" state.
 *    **No link is written yet.** The user must explicitly pick a SMILES
 *    column via the dropdown — avoiding auto-selection when multiple
 *    molecule columns exist (drizhina review #5).
 *  - Pick a SMILES column → write the link tag.
 *  - Change dropdown (while checked, already linked) → rewrite the link
 *    to the new SMILES column, clearing the previous.
 *  - Uncheck → clear the link (both tag and residual highlights);
 *    dropdown hidden.
 *
 *  Empty state: if the DataFrame has no SMILES columns, the widget
 *  renders a short info message and no controls. */
export function getMol3DAtomPickerLinkWidget(mol3DCol: DG.Column): DG.Widget {
  const df = mol3DCol.dataFrame;
  if (!df)
    return new DG.Widget(ui.divText('No DataFrame.'));

  const smilesCols = df.columns.toList()
    .filter((c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE);
  if (smilesCols.length === 0) {
    return new DG.Widget(ui.divText(
      'No SMILES columns in this DataFrame.',
      {style: {color: 'var(--grey-4)', fontStyle: 'italic', padding: '6px 0'}},
    ));
  }

  const currentLinkedSmiles = findLinkedSmilesColName(df, mol3DCol.name);
  const smilesNames = smilesCols.map((c) => c.name);

  // Dropdown: pre-selects the currently-linked SMILES column if any;
  // otherwise stays empty (`null`) so the user must make an explicit
  // choice — drizhina #5 asked us not to auto-pick when multiple
  // molecule columns are present.
  const choiceInput = ui.input.choice('SMILES column', {
    value: currentLinkedSmiles ?? null,
    items: smilesNames,
    nullable: true,
    onValueChanged: (newName) => {
      if (!enableInput.value) return;
      if (!newName) {
        // Dropdown cleared while checkbox is on: drop the link but keep
        // the checkbox state — user can pick again without re-toggling.
        clearLinksToMol3D(df, mol3DCol.name);
        df.fireValuesChanged();
        return;
      }
      clearLinksToMol3D(df, mol3DCol.name);
      const newCol = df.col(newName);
      if (newCol)
        setSmilesColLink(newCol, mol3DCol.name);
      df.fireValuesChanged();
    },
  });
  const choiceWrap = ui.div([choiceInput.root]);
  choiceWrap.style.display = currentLinkedSmiles ? '' : 'none';

  const enableInput = ui.input.bool('Link SMILES column', {
    value: currentLinkedSmiles !== null,
    onValueChanged: (checked) => {
      if (checked) {
        // Show dropdown; the onValueChanged handler writes the tag on pick.
        choiceWrap.style.display = '';
      } else {
        clearLinksToMol3D(df, mol3DCol.name);
        choiceInput.value = null;
        choiceWrap.style.display = 'none';
        df.fireValuesChanged();
      }
    },
  });
  enableInput.setTooltip(
    'Enable interactive atom highlighting between this 3D pose column and a 2D SMILES column.');

  return new DG.Widget(ui.divV([enableInput.root, choiceWrap]));
}
