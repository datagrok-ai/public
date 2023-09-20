import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {TaskBarProgressIndicator} from 'datagrok-api/dg';
import {NglForGridTestApp} from '../apps/ngl-for-grid-test-app';

import {_package} from '../package';

export async function biostructureInGridApp(appName: string, pi: TaskBarProgressIndicator): Promise<void> {
  const _piMsg = pi.description;
  const dfCsv: string = await _package.files.readAsText('pdb_data.csv');

  const df: DG.DataFrame = DG.DataFrame.fromCsv(dfCsv);

  const app = new NglForGridTestApp(appName);
  await app.init({df});
}
