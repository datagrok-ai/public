import * as DG from 'datagrok-api/dg';
//import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';  await _package.files.readAsText(name);

import {_package} from '../package-test';
import {before, category, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';


category('Viewers', () => {
  let treeHelper: ITreeHelper;

  before(async () => {
    treeHelper = await getTreeHelper();
  });

  const viewers = DG.Func.find({package: 'PhyloTreeViewer', meta: {role: DG.FUNC_TYPES.VIEWER}}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, await (async () => {
        const newickStr: string = await _package.files.readAsText('data/tree95.nwk');
        return treeHelper.newickToDf(newickStr, 'tree95');
      })(), {detectSemanticTypes: true});
    });
  }
});
