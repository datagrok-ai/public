import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import {getDataProviderList} from '@datagrok-libraries/bio/src/utils/data-provider';

import {MolstarViewer} from '../viewers/molstar-viewer';

import {_package} from '../package-test';
import {awaitGrid} from './utils';
import {DebounceIntervals} from '../viewers/molstar-viewer/molstar-viewer';


category('dataProvider', () => {
  const pdbIdCsv: string = `pdb_id
1QBS
1ZP8
2BDJ
1IAN
4UJ1
2BPW`;

  test('getDataProviderList-Molecule3D', async () => {
    // At least two data providers for Molecule3D are registered
    const listRes: DG.Func[] = await getDataProviderList(DG.SEMTYPE.MOLECULE3D);
    expect(listRes.length >= 2, true,
      `Data provider count for '${DG.SEMTYPE.MOLECULE3D}' semantic type.`);
  });

  test('MolstarViewer', async () => {
    const df = DG.DataFrame.fromCsv(pdbIdCsv);
    const tv = grok.shell.addTableView(df);
    const viewer = await df.plot.fromType('Biostructure',
      {
        biostructureIdColumnName: 'pdb_id',
        biostructureDataProvider: 'BiostructureViewer:getBiostructureRcsbMmcif'
      });
    tv.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure with data provider', 0.4);

    // await delay(50); // await for debounce onRebuildViewLigands
    // await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);

    df.currentRowIdx = 0;
    await awaitGrid(tv.grid);
    await delay(DebounceIntervals.currentRow * 2.5);
    await viewer.awaitRendered(5000);
    // @ts-ignore
    expect(viewer.dataEff!.options!.name, '1QBS');

    df.currentRowIdx = 2;
    await delay(DebounceIntervals.currentRow * 2.5);
    await Promise.all([awaitGrid(tv.grid), viewer.awaitRendered(5000)]);
    //@ts-ignore
    expect(viewer.dataEff!.options!.name, '2BDJ');
  }, {timeout: 15000});
});
