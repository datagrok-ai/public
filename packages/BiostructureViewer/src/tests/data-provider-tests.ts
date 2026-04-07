import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {getDataProviderList} from '@datagrok-libraries/bio/src/utils/data-provider';

import {awaitGrid} from './utils';
import {MolstarViewer} from '../viewers/molstar-viewer';
import {DebounceIntervals} from '../viewers/molstar-viewer/molstar-viewer';

import {_package} from '../package-test';


category('dataProvider', () => {
  test('getDataProviderList-Molecule3D', async () => {
    // At least two data providers for Molecule3D are registered
    const listRes: DG.Func[] = await getDataProviderList(DG.SEMTYPE.MOLECULE3D);
    expect(listRes.length >= 2, true,
      `Data provider count for '${DG.SEMTYPE.MOLECULE3D}' semantic type.`);
  });
});
