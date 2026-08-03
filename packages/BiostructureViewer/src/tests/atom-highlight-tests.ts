import * as grok from 'datagrok-api/grok';

import {category, delay, test} from '@datagrok-libraries/test/src/test';

import {CHEM_ATOM_SELECTION_EVENT} from '@datagrok-libraries/chem-meta/src/types';

category('AtomHighlight', () => {
  test('global cache stores persistent events only', async () => {
    grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
      rowIdx: 999, atoms: [0, 1, 2], persistent: true,
    });
    await delay(50);

    grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
      rowIdx: 999, atoms: [0, 1, 2, 3, 4, 5], persistent: false,
    });
    await delay(50);

    grok.events.fireCustomEvent(CHEM_ATOM_SELECTION_EVENT, {
      rowIdx: 999, atoms: [], persistent: true,
    });
    await delay(50);
  });
});
