/** Locks the data-testid convention (ff- namespace, slugified dynamic parts)
 *  so selectors stay stable as the UI grows. */
import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {tid, tidSlug, setTid} from '../utils/test-ids';

category('Flow: test-ids', () => {
  test('tidSlug makes selector-safe tokens', async () => {
    expect(tidSlug('Data Sources'), 'data-sources');
    expect(tidSlug('Chem:addChemProperties'), 'chem-addchemproperties');
    expect(tidSlug('DG Functions/Transform/Join Tables'), 'dg-functions-transform-join-tables');
    expect(tidSlug('   '), 'x', 'empty/blank falls back to a stable token');
    expect(tidSlug('--Already--Slugged--'), 'already-slugged', 'trims leading/trailing separators');
  });

  test('tid composes an ff-namespaced id from all parts', async () => {
    expect(tid('node'), 'ff-node');
    expect(tid('browser-item', 'OpenFile'), 'ff-browser-item-openfile');
    expect(tid('socket-input', 'table'), 'ff-socket-input-table');
    expect(tid('prop-input', 'Param Name'), 'ff-prop-input-param-name');
  });

  test('setTid stamps data-testid on an element', async () => {
    const el = document.createElement('div');
    const ret = setTid(el, 'node-status');
    expect(el.getAttribute('data-testid'), 'ff-node-status');
    expect(ret === el, true, 'returns the same element for chaining');
  });
});
