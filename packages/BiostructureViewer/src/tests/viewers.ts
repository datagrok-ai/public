import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {category, test, testViewer} from '@datagrok-libraries/test/src/test';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {asIViewer, IViewer} from '@datagrok-libraries/bio/src/viewers/viewer';

import {_package} from '../package-test';

category('Viewers', () => {
  const packageName = 'BiostructureViewer';
  // -- Viewers tests --
  const viewers = DG.Func.find({package: packageName, meta: {role: DG.FUNC_TYPES.VIEWER}}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, await (async () => {
        const ph: IPdbHelper = await grok.functions.call('BiostructureViewer:getPdbHelper');
        const pdbFn: string = `System:AppData/${_package.name}/samples/1bdq.pdb`;
        const pdbStr: string = await grok.dapi.files.readAsText(pdbFn);
        const df: DG.DataFrame = await ph.pdbToDf(pdbStr, '1bdq');
        // Put structure data to the tag of the DataFrame object to be displayed by the viewer
        df.setTag(pdbTAGS.PDB, pdbStr);
        return df;
      })(), {
        detectSemanticTypes: true,
        packageName: packageName,
        awaitViewer: async (viewer: DG.Viewer) => {
          const iViewer = asIViewer(viewer);
          if (iViewer)
            await iViewer.awaitRendered();
        }
      });
    }, {timeout: 30000});
  }
});
