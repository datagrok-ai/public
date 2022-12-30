import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';
import {hierarchicalClusteringUI} from '../utils/hierarchical-clustering';
import {_package} from '../package-test';
import {viewsTests} from './utils/views-tests';


category('hierarchicalClustering', viewsTests((ctx: { dfList: DG.DataFrame[], vList: DG.ViewBase[] }) => {
  test('UI', async () => {
    const csv: string = await _package.files.readAsText('data/demog-short.csv');
    const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    dataDf.name = 'testDemogShort';

    const tv: DG.TableView = grok.shell.addTableView(dataDf);

    ctx.vList.push(tv);
    ctx.dfList.push(dataDf);

    hierarchicalClusteringUI(dataDf, ['HEIGHT'], 'euclidean', 'ward');
  });
}));
