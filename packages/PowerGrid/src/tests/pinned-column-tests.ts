import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {awaitCheck, category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';

category('pinned-column', () => {
  const csv = `col1, col2, col3, col4
 a11,a12,a13,a14
 a21,a22,a23,a24
 a31,a32,a33,a34
 a41,a42,a43,a44`;

  test('single-first', async () => {
    const view = await _openGrid(csv);
    const pinGCol = view.grid.columns.byIndex(0)!;
    await grok.functions.call('PowerGrid:addPinnedColumn', {gridCol: pinGCol});
  });

  test('single-second', async () => {
    const view = await _openGrid(csv);
    const pinGCol = view.grid.columns.byIndex(1)!;
    await grok.functions.call('PowerGrid:addPinnedColumn', {gridCol: pinGCol});
  });


  test('twice', async () => {
    const view = await _openGrid(csv);
    const pin1GCol = view.grid.columns.byIndex(0)!;
    const pin2GCol = view.grid.columns.byIndex(1)!;
    await grok.functions.call('PowerGrid:addPinnedColumn', {gridCol: pin1GCol});
    await grok.functions.call('PowerGrid:addPinnedColumn', {gridCol: pin2GCol});
  });
}, {owner: 'dkovalyov@datagrok.ai'});

async function _openGrid(csv: string): Promise<DG.TableView> {
  const df = DG.DataFrame.fromCsv(csv);
  const view = grok.shell.addTableView(df);

  await awaitCheck(() => { return $(view.root).find('.d4-grid canvas').length > 0; });
  return view;
}
