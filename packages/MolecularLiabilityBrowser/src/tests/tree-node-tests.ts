import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject, expectArray} from '@datagrok-libraries/utils/src/test';
import {getVId, mlbTreeNodeRe} from '../utils/tree-stats';

category('treeNode', () => {

  test('mlbTreeNodeRe1', async () => {
    const ma = mlbTreeNodeRe.exec('TPP000022964|PSMB325|VR000012782|PRJ000000165');
    expectArray(ma.slice(1, 5), ['TPP000022964', 'PSMB325', 'VR000012782', 'PRJ000000165']);
  });

  test('getVId1', async () => {
    expect(getVId('TPP000022964|PSMB325|VR000012782|PRJ000000165'), 'VR000012782');
    expect(getVId('TPP000022964|PSMB325|VR000012782|PRJ000000165_1'), 'VR000012782');
    expect(getVId('TPP000180676|PSMB1039|_VR000035650|PRJ000000165'), '_VR000035650');
    expect(getVId('TPP000180495|PSMB971|_VR000044692|PRJ000000165'), '_VR000044692');
    expect(getVId('TPP000051800|PSMW32|VR000012782|PRJ000000165'), 'VR000012782');
    expect(getVId('TPP000041924|PSMB331|_VR000038164|PRJ000000165'), '_VR000038164');
    expect(getVId('TPP000206916|PSMB1517|VR000012782|PRJ000000165'), 'VR000012782');
  });
});