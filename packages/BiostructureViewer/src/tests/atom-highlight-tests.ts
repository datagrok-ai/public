import * as grok from 'datagrok-api/grok';

import {category, test, testEvent} from '@datagrok-libraries/test/src/test';

import {CHEM_ATOM_SELECTION_EVENT} from '@datagrok-libraries/chem-meta/src/types';

/** Fires the event and resolves once it is delivered: the subscribers that came before have handled it by then. */
async function fire(args: object): Promise<void> {
  await testEvent(grok.events.onCustomEvent(CHEM_ATOM_SELECTION_EVENT), () => {},
    () => grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, args), 5000);
}

category('AtomHighlight', () => {
  test('global cache stores persistent events only', async () => {
    await fire({rowIdx: 999, atoms: [0, 1, 2], persistent: true});

    await fire({rowIdx: 999, atoms: [0, 1, 2, 3, 4, 5], persistent: false});

    await fire({rowIdx: 999, atoms: [], persistent: true});
  });
});
