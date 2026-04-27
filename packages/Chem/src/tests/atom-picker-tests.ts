import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';

import {CHEM_INTERACTIVE_SELECTION_EVENT} from '../constants';

category('atom picker', () => {
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
