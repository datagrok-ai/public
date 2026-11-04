import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';

export async function demoBioSimDivLayout() {
  const t = await _package.files.readCsv('samples/peptides-non-natural.csv');
  const tv = grok.shell.addTableView(t);
  await _package.files.readAsText('demo-files/bio_similarity_diversity.layout').then(async (layoutString: string) => {
    const layout = DG.ViewLayout.fromJson(layoutString);
    tv.loadLayout(layout);
  });
}
