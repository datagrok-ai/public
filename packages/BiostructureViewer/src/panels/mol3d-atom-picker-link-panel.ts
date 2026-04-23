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
 *  - Check the box → auto-link to the currently selected dropdown
 *    value (or the first SMILES column if no prior link); dropdown
 *    becomes visible.
 *  - Uncheck → clear the link (both tag and residual highlights);
 *    dropdown hidden.
 *  - Change dropdown (while checked) → rewrite the link to the new
 *    SMILES column, clearing the previous.
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

  // Dropdown reflects the currently-linked SMILES column name, or
  // defaults to the first available SMILES column if no link exists.
  const choiceInput = ui.input.choice('SMILES column', {
    value: currentLinkedSmiles ?? smilesNames[0],
    items: smilesNames,
    onValueChanged: (newName) => {
      // Only write a new link when the enable checkbox is on; otherwise
      // changing the dropdown while disabled is a no-op (prevents
      // accidentally re-linking after an unlink).
      if (!enableInput.value || !newName)
        return;
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
        // Auto-link to the dropdown's current value so the checkbox
        // state always matches a real tag (no "checked but unlinked"
        // limbo state).
        const target = choiceInput.value || smilesNames[0];
        const newCol = df.col(target);
        if (newCol) {
          clearLinksToMol3D(df, mol3DCol.name);
          setSmilesColLink(newCol, mol3DCol.name);
          df.fireValuesChanged();
        }
        choiceWrap.style.display = '';
      } else {
        clearLinksToMol3D(df, mol3DCol.name);
        df.fireValuesChanged();
        choiceWrap.style.display = 'none';
      }
    },
  });
  enableInput.setTooltip(
    'Enable interactive atom highlighting between this 3D pose column and a 2D SMILES column.');

  return new DG.Widget(ui.divV([enableInput.root, choiceWrap]));
}
