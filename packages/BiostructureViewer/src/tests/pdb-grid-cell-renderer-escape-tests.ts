import {category, expect, test} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import {AtomPickerProvider} from '@datagrok-libraries/chem-meta/src/types';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import {setSmilesColLink} from '../utils/mol3d-link';
import {PdbGridCellRenderer} from '../utils/pdb-grid-cell-renderer';

const PROVIDERS_KEY = ChemTemps.SUBSTRUCT_PROVIDERS;

function plantProvider(col: DG.Column): void {
  const prov: AtomPickerProvider = {__atomPicker: true, getSubstruct: () => undefined};
  const existing = (col.temp[PROVIDERS_KEY] as AtomPickerProvider[] | undefined) ?? [];
  col.temp[PROVIDERS_KEY] = [...existing, prov];
}

function getPickerProviders(col: DG.Column): AtomPickerProvider[] {
  return ((col.temp[PROVIDERS_KEY] ?? []) as AtomPickerProvider[]).filter((p) => p.__atomPicker);
}

function makeLinkedDF(): {smilesCol: DG.Column; mol3DCol: DG.Column} {
  const smilesCol = DG.Column.fromStrings('smiles', ['CCO', 'CC']);
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  const mol3DCol = DG.Column.fromStrings('pose3D', ['', '']);
  mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
  DG.DataFrame.fromColumns([smilesCol, mol3DCol]);
  // Writes both forward (smiles→mol3D) and reciprocal (mol3D→smiles) tags.
  setSmilesColLink(smilesCol, mol3DCol.name);
  return {smilesCol, mol3DCol};
}

category('pdbGridCellRenderer: escape', () => {
  test('Escape — clears providers on linked SMILES column', async () => {
    const {smilesCol, mol3DCol} = makeLinkedDF();
    plantProvider(smilesCol);
    expect(getPickerProviders(smilesCol).length, 1);

    const renderer = new PdbGridCellRenderer();
    renderer.onKeyDown(
      {tableColumn: mol3DCol} as unknown as DG.GridCell,
      new KeyboardEvent('keydown', {key: 'Escape'}),
    );
    expect(getPickerProviders(smilesCol).length, 0);
  });

  test('Escape — targets only the linked SMILES column', async () => {
    const smilesLinked = DG.Column.fromStrings('smiles_linked', ['CCO']);
    smilesLinked.semType = DG.SEMTYPE.MOLECULE;
    const smilesOther = DG.Column.fromStrings('smiles_other', ['CC']);
    smilesOther.semType = DG.SEMTYPE.MOLECULE;
    const mol3DCol = DG.Column.fromStrings('pose3D', ['']);
    mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
    DG.DataFrame.fromColumns([smilesLinked, smilesOther, mol3DCol]);

    setSmilesColLink(smilesLinked, mol3DCol.name);
    plantProvider(smilesLinked);
    plantProvider(smilesOther);

    const renderer = new PdbGridCellRenderer();
    renderer.onKeyDown(
      {tableColumn: mol3DCol} as unknown as DG.GridCell,
      new KeyboardEvent('keydown', {key: 'Escape'}),
    );

    expect(getPickerProviders(smilesLinked).length, 0);
    expect(getPickerProviders(smilesOther).length, 1);
  });
});
