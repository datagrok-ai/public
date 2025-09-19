import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import * as api from '../package-api';

import {category, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {asIViewer, IViewer} from '@datagrok-libraries/bio/src/viewers/viewer';

import {_package} from '../package-test';

category('Viewers', () => {
  const packageName = 'BiostructureViewer';
  // -- Viewers tests --
  const viewers = DG.Func.find({package: packageName, tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, await (async () => {
        const ph: IPdbHelper = await api.funcs.getPdbHelper();
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
