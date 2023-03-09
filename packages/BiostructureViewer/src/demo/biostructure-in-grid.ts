import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {TaskBarProgressIndicator} from 'datagrok-api/dg';
import {NglForGridTestApp} from '../apps/ngl-for-grid-test-app';

export async function biostructureInGridApp(appName: string, pi: TaskBarProgressIndicator): Promise<void> {
  const piMsg = pi.description;
  const dfCsv: string = await _package.files.readAsText('pdb_data.csv');

  const df: DG.DataFrame = DG.DataFrame.fromCsv(dfCsv);

  const app = new NglForGridTestApp(appName);
  await app.init({df});
}
