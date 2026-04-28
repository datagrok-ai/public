/** Checkbox+dropdown widget for linking a Molecule3D column to a SMILES column.
 *  Consumed cross-package by Docking's AutoDock panel via grok.functions.call. */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  clearLinksToMol3D,
  findLinkedSmilesColName,
  setSmilesColLink,
} from '../utils/mol3d-link';

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

  const choiceInput = ui.input.choice('SMILES column', {
    value: currentLinkedSmiles ?? null,
    items: smilesNames,
    nullable: true,
    tooltipText: 'Enable interactive atom highlighting between this 3D pose column and a 2D SMILES column.',
    onValueChanged: (newName) => {
      if (!newName) {
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
  // choiceWrap.style.display = currentLinkedSmiles ? '' : 'none';

  return new DG.Widget(ui.divV([choiceWrap]));
}
