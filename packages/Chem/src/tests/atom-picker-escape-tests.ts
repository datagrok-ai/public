import {before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';

import {CHEM_INTERACTIVE_SELECTION_EVENT, CHEM_ATOM_PICKER_LINKED_COL} from '../constants';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';
import {awaitGrid} from './utils';

const PROVIDERS_KEY = ChemTemps.SUBSTRUCT_PROVIDERS;

type PickerProvider = {
  __atomPicker?: boolean;
  __rowIdx?: number;
  __atoms?: Set<number>;
  getSubstruct: (rowIdx: number) => {atoms?: number[]} | undefined;
};

function plantProvider(col: DG.Column, rowIdx: number, atoms: number[]): void {
  const existing = (col.temp[PROVIDERS_KEY] ?? []) as PickerProvider[];
  const prov: PickerProvider = {
    __atomPicker: true,
    __rowIdx: rowIdx,
    __atoms: new Set(atoms),
    getSubstruct: (r) => r === rowIdx ? {atoms} : undefined,
  };
  col.temp[PROVIDERS_KEY] = [...existing, prov];
}

function getPickerProviders(col: DG.Column): PickerProvider[] {
  return ((col.temp[PROVIDERS_KEY] ?? []) as PickerProvider[]).filter((p) => p.__atomPicker);
}

function makeLinkedDF(): {df: DG.DataFrame; smilesCol: DG.Column; mol3DCol: DG.Column} {
  const smilesCol = DG.Column.fromStrings('smiles', ['CCO', 'CC']);
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  const mol3DCol = DG.Column.fromStrings('pose3D', ['', '']);
  mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
  const df = DG.DataFrame.fromColumns([smilesCol, mol3DCol]);
  smilesCol.setTag(CHEM_ATOM_PICKER_LINKED_COL, mol3DCol.name);
  return {df, smilesCol, mol3DCol};
}

/** Test-only structural alias for private hover state on RDKitCellRenderer.
 *  Uses `as unknown as` (not intersection) because TypeScript collapses an
 *  intersection to `never` when any of the named members is already `private`
 *  in the base class. The alias names match the real field/method names exactly
 *  so a rename in the source produces a type error here. */
interface RDKitCellRendererInternals {
  _previewAtomIdx: number | null;
  _previewFrom3D: boolean;
  _lastHoveredAtom: unknown;
  onKeyDown(gridCell: DG.GridCell, e: KeyboardEvent): void;
}

category('atom picker: escape', () => {
  let rdkitModule: any;

  before(async () => {
    rdkitModule = await grok.functions.call('Chem:getRdKitModule');
  });

  test('Escape — clears atom-picker providers', async () => {
    const {df, smilesCol} = makeLinkedDF();

    plantProvider(smilesCol, 0, [0, 1]);
    plantProvider(smilesCol, 1, [2]);
    expect(getPickerProviders(smilesCol).length, 2);

    const renderer = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    const view = grok.shell.addTableView(df);
    try {
      const fakeGc = {grid: view.grid} as unknown as DG.GridCell;
      renderer.onKeyDown(fakeGc, new KeyboardEvent('keydown', {key: 'Escape'}));
      expect(getPickerProviders(smilesCol).length, 0);
    } finally {
      view.close();
    }
  }, {timeout: 15000});

  test('Escape — preserves non-picker providers', async () => {
    const {df, smilesCol} = makeLinkedDF();

    const scaffoldProv: PickerProvider = {
      __atomPicker: undefined,
      getSubstruct: () => undefined,
    };
    smilesCol.temp[PROVIDERS_KEY] = [scaffoldProv];
    plantProvider(smilesCol, 0, [0]);

    const renderer = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    const view = grok.shell.addTableView(df);
    try {
      renderer.onKeyDown(
        {grid: view.grid} as unknown as DG.GridCell,
        new KeyboardEvent('keydown', {key: 'Escape'}),
      );
      const remaining = (smilesCol.temp[PROVIDERS_KEY] ?? []) as PickerProvider[];
      expect(remaining.some((p) => !p.__atomPicker), true);
      expect(remaining.some((p) => p.__atomPicker === true), false);
    } finally {
      view.close();
    }
  }, {timeout: 15000});

  test('Escape — fires clear-all selection event', async () => {
    const {df, smilesCol} = makeLinkedDF();
    plantProvider(smilesCol, 0, [0, 1, 2]);

    const renderer = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    let received: any = null;
    const sub = grok.events.onCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT)
      .subscribe((args: any) => {received = args;});

    const view = grok.shell.addTableView(df);
    try {
      renderer.onKeyDown(
        {grid: view.grid} as unknown as DG.GridCell,
        new KeyboardEvent('keydown', {key: 'Escape'}),
      );
      await delay(150); // Allow RxJS event delivery.
      expect(received !== null, true);
      expect(received.clearAll, true);
      expect(received.persistent, true);
      expect(Array.isArray(received.atoms) && received.atoms.length === 0, true);
    } finally {
      sub.unsubscribe();
      view.close();
    }
  }, {timeout: 15000});

  test('non-Escape key is no-op', async () => {
    const {df, smilesCol} = makeLinkedDF();
    plantProvider(smilesCol, 0, [0]);

    const renderer = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    const view = grok.shell.addTableView(df);
    try {
      renderer.onKeyDown(
        {grid: view.grid} as unknown as DG.GridCell,
        new KeyboardEvent('keydown', {key: 'Enter'}),
      );
      expect(getPickerProviders(smilesCol).length, 1);
    } finally {
      view.close();
    }
  }, {timeout: 15000});

  test('Escape — no-op without picker-active column', async () => {
    // DataFrame with no MOLECULE semType column — _isPickerActive returns false.
    const col = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [1.0, 2.0]);
    const df = DG.DataFrame.fromColumns([col]);

    const renderer = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    const view = grok.shell.addTableView(df);
    try {
      renderer.onKeyDown(
        {grid: view.grid} as unknown as DG.GridCell,
        new KeyboardEvent('keydown', {key: 'Escape'}),
      );
    } finally {
      view.close();
    }
  }, {timeout: 15000});

  test('Escape — resets preview state', async () => {
    const {df, smilesCol} = makeLinkedDF();
    plantProvider(smilesCol, 0, [0]);

    const r = new RDKitCellRenderer(rdkitModule) as unknown as RDKitCellRendererInternals;

    r._previewAtomIdx = 3;
    r._previewFrom3D = true;
    r._lastHoveredAtom = {col: 'smiles', rowIdx: 0, atomIdx: 3};

    const view = grok.shell.addTableView(df);
    try {
      r.onKeyDown(
        {grid: view.grid} as unknown as DG.GridCell,
        new KeyboardEvent('keydown', {key: 'Escape'}),
      );
      expect(r._previewAtomIdx, null);
      expect(r._previewFrom3D, false);
      expect(r._lastHoveredAtom, null);
    } finally {
      view.close();
    }
  }, {timeout: 15000});

  // Path 2 — grid.root capture-phase listener wired by _attachGridCanvasShiftListener on
  // first render. We wait for a real draw pass then dispatch on grid.root without
  // needing any SMILES rendering.
  test('Escape on grid root — clears from any focused cell', async () => {
    const {df, smilesCol} = makeLinkedDF();
    plantProvider(smilesCol, 0, [0, 1]);

    const view = grok.shell.addTableView(df);
    df.currentRowIdx = 0;

    try {
      await awaitGrid(view.grid);
      await delay(200);

      expect(getPickerProviders(smilesCol).length, 1);

      view.grid.root.dispatchEvent(
        new KeyboardEvent('keydown', {key: 'Escape', bubbles: true, cancelable: true}),
      );

      await delay(150);

      expect(getPickerProviders(smilesCol).length, 0);
    } finally {
      view.close();
    }
  }, {timeout: 30000});
});
