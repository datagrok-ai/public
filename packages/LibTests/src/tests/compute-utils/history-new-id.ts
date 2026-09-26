import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

category('ComputeUtils: History new id', async () => {
  test('Saving a loaded run with newId keeps the earlier run', async () => {
    const saved: DG.FuncCall[] = [];
    try {
      const fc = DG.Func.byName('LibTests:TestAdd2').prepare({a: 1, b: 2});
      await fc.call();
      fc.options['title'] = 'first';
      const originalId = fc.id;
      const first = await historyUtils.saveRun(fc, {newId: true});
      saved.push(first);
      expectDeepEqual(first.id, fc.id, {prefix: 'live call carries the saved id'});
      expectDeepEqual(first.id !== originalId, true, {prefix: 'saved id is fresh'});

      const loaded = await historyUtils.loadRun(first.id);
      loaded.inputs['a'] = 5;
      loaded.options['title'] = 'second';
      const second = await historyUtils.saveRun(loaded, {newId: true});
      saved.push(second);
      expectDeepEqual(second.id !== first.id, true, {prefix: 'second save gets a new id'});

      const firstAgain = await historyUtils.loadRun(first.id);
      expectDeepEqual(firstAgain.options['title'], 'first');
      expectDeepEqual(firstAgain.inputs['a'], 1);
      const secondAgain = await historyUtils.loadRun(second.id);
      expectDeepEqual(secondAgain.options['title'], 'second');
      expectDeepEqual(secondAgain.inputs['a'], 5);
    } finally {
      for (const run of saved)
        await grok.dapi.functions.calls.allPackageVersions().delete(run);
    }
  });
});
