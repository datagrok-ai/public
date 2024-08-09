import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {PipelineConfiguration} from '@datagrok-libraries/compute-utils';
import {getProcessedConfig} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import { callHandler } from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import { expectDeepEqual } from '@datagrok-libraries/utils/src/expect';

category('ComputeUtils: Driver state tree persistence', async () => {
  test('Load and save simple config', async () => {
    const config = await callHandler<PipelineConfiguration>('LibTests:MockProvider1', {}).toPromise();
    const pconf = await getProcessedConfig(config);
    const tree = StateTree.fromConfig(pconf);
    await tree.init().toPromise();
    const sc = tree.toSerializedState(true);
    console.log(sc);
    const id = await tree.save().toPromise();
    const loadedTree = await StateTree.load(id, pconf).toPromise();
    await loadedTree.init().toPromise();
    const lc = loadedTree.toSerializedState(true);
    expectDeepEqual(sc, lc);
    console.log(lc);
  });
});
