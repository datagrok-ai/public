import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

category('GridWithTree', () => {
  let dfList: DG.DataFrame[];
  let vList: DG.ViewBase[];
  let currentView: DG.ViewBase;

  before(async () => {
    dfList = [];
    vList = [];
    currentView = grok.shell.v;
  });

  after(async () => {
    for (const df of dfList) grok.shell.closeTable(df);
    for (const v of vList) v.close();
    grok.shell.v = currentView;
  });

  test('open', async () => {
    // TODO:
  });
});