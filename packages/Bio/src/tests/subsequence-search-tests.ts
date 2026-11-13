import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, isDialogPresent, test} from '@datagrok-libraries/test/src/test';

category('subsequenceSearch', () => {
  test('noMacromoleculeColumn', async () => {
    const view = grok.shell.addTableView(grok.data.demo.demog(10));
    try {
      const call = DG.Func.find({package: 'Bio', name: 'SubsequenceSearchTopMenu'})[0].prepare({});
      await grok.functions.call('Bio:SearchSubsequenceEditor', {call});
      expect(isDialogPresent('Substructure Search'), false);
    } finally {
      view.close();
    }
  });
});
