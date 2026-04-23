import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Tag key mirroring `CHEM_ATOM_PICKER_LINKED_COL` from Chem's constants.ts.
 *  Set on a Molecule (SMILES) column to link it to a specific Molecule3D
 *  column for interactive atom picking. Defined here to avoid a
 *  cross-package import. */
const LINKED_TAG = '%chem-atom-picker-linked-col';

const NONE = '(none)';

/** Context-panel widget for Molecule3D columns. Lets the user explicitly
 *  pick which SMILES (Molecule) column in the same DataFrame is linked
 *  to this docking-poses column. When a link is set, the 2D interactive
 *  atom picker in the Chem cell renderer activates for that SMILES column
 *  and synchronises highlights with the 3D Molstar viewer. */
export function getMol3DColumnPropertyPanel(mol3DCol: DG.Column): DG.Widget {
  const df = mol3DCol.dataFrame;
  if (!df)
    return new DG.Widget(ui.divText('No DataFrame'));

  // Candidates: all Molecule (SMILES) columns in the same DataFrame.
  const candidates = df.columns.toList()
    .filter((c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE)
    .map((c: DG.Column) => c.name);

  // Current link: the SMILES column whose tag points at this Mol3D column.
  const currentSmiles = df.columns.toList()
    .find((c: DG.Column) => c.temp[LINKED_TAG] === mol3DCol.name)?.name ?? NONE;

  const choice = ui.input.choice('Linked SMILES column', {
    value: currentSmiles,
    items: [NONE].concat(candidates),
    onValueChanged: (value) => {
      // Clear any existing tag that points at this Mol3D column
      // (handles switching from one SMILES col to another).
      for (const c of df.columns.toList()) {
        if (c.temp[LINKED_TAG] === mol3DCol.name)
          delete c.temp[LINKED_TAG];
      }
      // Set the tag on the newly chosen SMILES column.
      if (value !== NONE) {
        const newCol = df.col(value);
        if (newCol)
          newCol.temp[LINKED_TAG] = mol3DCol.name;
      }
      df.fireValuesChanged();
    },
  });
  choice.setTooltip('Select the SMILES column to link for interactive atom highlighting');

  return new DG.Widget(ui.form([choice]));
}
