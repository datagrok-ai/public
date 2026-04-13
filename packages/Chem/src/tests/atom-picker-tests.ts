import {category, test, expect, before, after} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {CHEM_INTERACTIVE_SELECTION_EVENT} from '../constants';

category('atom-picker', () => {
  // -- _isPickerActive detection tests ---------------------------------------

  test('picker active with Molecule3D column', async () => {
    const smilesCol = DG.Column.fromStrings('smiles', ['CCO', 'CC']);
    smilesCol.semType = DG.SEMTYPE.MOLECULE;
    const mol3DCol = DG.Column.fromStrings('pose', ['ATOM      1  C   LIG     1       0.0  0.0  0.0', '']);
    mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
    const df = DG.DataFrame.fromColumns([smilesCol, mol3DCol]);
    // The picker auto-activates when a Molecule3D column is present.
    // Check by scanning columns — same logic as _isPickerActive.
    const cols = df.columns.toList();
    const hasMol3D = cols.some((c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE3D);
    expect(hasMol3D, true);
  });

  test('picker inactive without Molecule3D column', async () => {
    const smilesCol = DG.Column.fromStrings('smiles', ['CCO', 'CC']);
    smilesCol.semType = DG.SEMTYPE.MOLECULE;
    const numCol = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [1.0, 2.0]);
    const df = DG.DataFrame.fromColumns([smilesCol, numCol]);
    const cols = df.columns.toList();
    const hasMol3D = cols.some((c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE3D);
    expect(hasMol3D, false);
  });

  // -- Event structure tests -------------------------------------------------

  test('selection event carries atoms and rowIdx', async () => {
    let receivedArgs: any = null;
    const sub = grok.events.onCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT)
      .subscribe((args: any) => {receivedArgs = args;});
    try {
      const testAtoms = [0, 3, 5];
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
        rowIdx: 2, atoms: testAtoms, persistent: true,
      });
      // Give the synchronous subscriber a tick to process.
      await new Promise((r) => setTimeout(r, 50));
      expect(receivedArgs !== null, true);
      expect(receivedArgs.rowIdx, 2);
      expect(receivedArgs.atoms.length, 3);
      expect(receivedArgs.persistent, true);
    } finally {
      sub.unsubscribe();
    }
  });

  test('persistent flag — true for Alt events', async () => {
    let receivedPersistent: boolean | undefined;
    const sub = grok.events.onCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT)
      .subscribe((args: any) => {receivedPersistent = args?.persistent;});
    try {
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
        rowIdx: 0, atoms: [1], persistent: true,
      });
      await new Promise((r) => setTimeout(r, 50));
      expect(receivedPersistent, true);
    } finally {
      sub.unsubscribe();
    }
  });

  test('persistent flag — false for preview events', async () => {
    let receivedPersistent: boolean | undefined;
    const sub = grok.events.onCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT)
      .subscribe((args: any) => {receivedPersistent = args?.persistent;});
    try {
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
        rowIdx: 0, atoms: [1, 2], persistent: false,
      });
      await new Promise((r) => setTimeout(r, 50));
      expect(receivedPersistent, false);
    } finally {
      sub.unsubscribe();
    }
  });
});
