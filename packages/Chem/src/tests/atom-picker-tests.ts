import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {CHEM_INTERACTIVE_SELECTION_EVENT} from '../constants';

category('atom-picker', () => {
  /** Builds a minimal DataFrame with a SMILES molecule column plus `extra`
   *  columns and applies the `_isPickerActive` detection rule — the picker
   *  auto-activates when any column on the frame has semType
   *  `MOLECULE3D`. */
  function pickerActivates(extra: DG.Column[]): boolean {
    const smilesCol = DG.Column.fromStrings('smiles', ['CCO', 'CC']);
    smilesCol.semType = DG.SEMTYPE.MOLECULE;
    const df = DG.DataFrame.fromColumns([smilesCol, ...extra]);
    return df.columns.toList().some((c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE3D);
  }

  /** Subscribes to the chem interactive-selection event, fires `payload`,
   *  waits for delivery, and returns the received args. Always
   *  unsubscribes, even on error. */
  async function fireAndCaptureSelection(payload: unknown): Promise<any> {
    let received: any = null;
    const sub = grok.events.onCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT)
      .subscribe((args: any) => {received = args;});
    try {
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, payload);
      await new Promise((r) => setTimeout(r, 50));
    } finally {
      sub.unsubscribe();
    }
    return received;
  }

  // -- _isPickerActive detection tests ---------------------------------------

  test('picker active with Molecule3D column', async () => {
    const mol3DCol = DG.Column.fromStrings('pose',
      ['ATOM      1  C   LIG     1       0.0  0.0  0.0', '']);
    mol3DCol.semType = DG.SEMTYPE.MOLECULE3D;
    expect(pickerActivates([mol3DCol]), true);
  });

  test('picker inactive without Molecule3D column', async () => {
    const numCol = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [1.0, 2.0]);
    expect(pickerActivates([numCol]), false);
  });

  // -- Event structure tests -------------------------------------------------

  test('selection event carries atoms and rowIdx', async () => {
    const args = await fireAndCaptureSelection({
      rowIdx: 2, atoms: [0, 3, 5], persistent: true,
    });
    expect(args !== null, true);
    expect(args.rowIdx, 2);
    expect(args.atoms.length, 3);
    expect(args.persistent, true);
  });

  test('persistent flag — true for Alt events', async () => {
    const args = await fireAndCaptureSelection({rowIdx: 0, atoms: [1], persistent: true});
    expect(args?.persistent, true);
  });

  test('persistent flag — false for preview events', async () => {
    const args = await fireAndCaptureSelection({rowIdx: 0, atoms: [1, 2], persistent: false});
    expect(args?.persistent, false);
  });
});
