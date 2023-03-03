import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

import {category, test, testViewer} from '@datagrok-libraries/utils/src/test';
//import {_package} from './package-test-int';
import {PdbHelper} from '../utils/pdb-helper';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {getPdbHelper} from '../package';
import {_packageName} from './utils';


category('Viewers', () => {

  // -- Viewers tests --
  const viewers = DG.Func.find({package: 'BiostructureViewer', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, await (async () => {
        const ph: IPdbHelper = await getPdbHelper();
        const pdbFn: string = `System:AppData/${_packageName}/samples/1bdq.pdb`;
        const pdbStr: string = await grok.dapi.files.readAsText(pdbFn);
        const df: DG.DataFrame = await ph.pdbToDf(pdbStr, '1bdq');
        return df;
      })(), true);
    });
  }
});
