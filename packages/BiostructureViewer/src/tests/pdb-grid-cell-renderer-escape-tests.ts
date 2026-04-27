import {category, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {setSmilesColLink} from '../utils/mol3d-link';
import {PdbGridCellRenderer} from '../utils/pdb-grid-cell-renderer';

const PROVIDERS_KEY = 'substruct-providers';

type PickerProvider = {
  __atomPicker?: boolean;
  getSubstruct: () => undefined;
};

function plantProvider(col: DG.Column): void {
  const prov: PickerProvider = {__atomPicker: true, getSubstruct: () => undefined};
  const existing = (col.temp[PROVIDERS_KEY] as PickerProvider[] | undefined) ?? [];
  col.temp[PROVIDERS_KEY] = [...existing, prov];
}

function getPickerProviders(col: DG.Column): PickerProvider[] {
  return ((col.temp[PROVIDERS_KEY] ?? []) as PickerProvider[]).filter((p) => p.__atomPicker);
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

  test('Escape — fires clear-all selection event', async () => {
    const {smilesCol, mol3DCol} = makeLinkedDF();
    plantProvider(smilesCol);

    const CHEM_SELECTION_EVENT = 'chem-interactive-selection-changed';
    let receivedClearAll = false;
    const sub = grok.events.onCustomEvent(CHEM_SELECTION_EVENT)
      .subscribe((args: {clearAll?: boolean}) => {
        if (args?.clearAll) receivedClearAll = true;
      });

    const renderer = new PdbGridCellRenderer();

    try {
      renderer.onKeyDown(
        {tableColumn: mol3DCol} as unknown as DG.GridCell,
        new KeyboardEvent('keydown', {key: 'Escape'}),
      );
      await new Promise<void>((r) => setTimeout(r, 250));
      expect(receivedClearAll, true);
    } finally {
      sub.unsubscribe();
    }
  }, {timeout: 10000});

  test('non-Escape key is no-op', async () => {
    const {smilesCol, mol3DCol} = makeLinkedDF();
    plantProvider(smilesCol);

    const renderer = new PdbGridCellRenderer();
    renderer.onKeyDown(
      {tableColumn: mol3DCol} as unknown as DG.GridCell,
      new KeyboardEvent('keydown', {key: 'Enter'}),
    );
    expect(getPickerProviders(smilesCol).length, 1);
  });

  test('Escape — no-op without linked SMILES column', async () => {
    const mol3DCol = DG.Column.fromStrings('pose3D', ['']);
    mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
    DG.DataFrame.fromColumns([mol3DCol]);

    const renderer = new PdbGridCellRenderer();
    renderer.onKeyDown(
      {tableColumn: mol3DCol} as unknown as DG.GridCell,
      new KeyboardEvent('keydown', {key: 'Escape'}),
    );
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
