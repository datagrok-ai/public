import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

export function viewsTests(tests: (ctx: { dfList: DG.DataFrame[], vList: DG.ViewBase[] }) => void): () => void {
  return () => {
    const ctx: {
      dfList: DG.DataFrame[],
      vList: DG.ViewBase[],
      currentView: DG.ViewBase | null
    } = {dfList: [], vList: [], currentView: null};

    before(async () => {
      ctx.dfList = [];
      ctx.vList = [];
      ctx.currentView = grok.shell.v;
    });

    after(async () => {
      // for (const df of ctx.dfList) grok.shell.closeTable(df);
      // for (const v of ctx.vList) v.close();
      grok.shell.v = ctx.currentView!;
    });

    tests(ctx);
  };
}
