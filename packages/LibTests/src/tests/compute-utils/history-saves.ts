import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {deleteRuns} from '../../test-utils';

category('ComputeUtils: History saves', async () => {
  test('A save that fails only to link the run into its space succeeds and shows why', async () => {
    const space = await grok.dapi.spaces.createRootSpace(`Deleted space ${Date.now()}`);
    await grok.dapi.spaces.delete(space);
    const fc = DG.Func.byName('LibTests:TestAdd2').prepare({a: 1, b: 2});
    await fc.call();
    try {
      const saved = await historyUtils.saveRun(fc,
        {newId: true, metadata: {name: 'Unlinked'}, target: {spaceId: space.id}});
      expectDeepEqual(saved.id, fc.id, {prefix: 'the save resolves with the run'});
      expectDeepEqual((await historyUtils.loadRun(fc.id)).options['title'], 'Unlinked', {prefix: 'saved run'});
      await awaitCheck(() => Array.from(document.querySelectorAll('.d4-balloon.error'))
        .some((balloon) => balloon.textContent?.includes('Saved the run, but couldn\'t add it to the space')),
      'a failed link should show an error', 5000);
    } finally {
      await deleteRuns([fc.id]);
    }
  });

  test('A save with only a name keeps the description and tags', async () => {
    const fc = DG.Func.byName('LibTests:TestAdd2').prepare({a: 1, b: 2});
    await fc.call();
    const saved: string[] = [];
    try {
      saved.push((await historyUtils.saveRun(fc, {newId: true,
        metadata: {name: 'Described', description: 'Kept', tags: ['kept']}})).id);
      const loaded = await historyUtils.loadRun(saved[0]);
      saved.push((await historyUtils.saveRun(loaded, {newId: true, metadata: {name: 'Renamed'}})).id);
      expectDeepEqual(historyUtils.getRunMetadata(await historyUtils.loadRun(saved[1])),
        {name: 'Renamed', description: 'Kept', tags: ['kept']});
    } finally {
      await deleteRuns(saved);
    }
  });
});
