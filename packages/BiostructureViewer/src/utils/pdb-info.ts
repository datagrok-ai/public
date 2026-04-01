import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RcsbGraphQLAdapter} from './rcsb-gql-adapter';

export async function pdbInfoWidget(pdbId: string): Promise<DG.Widget> {
  try {
    const info = await RcsbGraphQLAdapter.getEntryInfo(pdbId);
    const map: {[key: string]: any} = {};
    if (info.experimentalMethod)
      map['Method'] = info.experimentalMethod;
    if (info.resolution != null)
      map['Resolution'] = `${info.resolution} \u00C5`;
    if (info.title)
      map['Description'] = info.title;
    if (info.doi)
      map['PDB DOI'] = ui.link(`https://doi.org/${info.doi}`, () => window.open(`https://doi.org/${info.doi}`));
    if (info.classification)
      map['Classification'] = info.classification;
    if (info.organisms)
      map['Organism(s)'] = info.organisms.join(', ');
    if (info.expressionSystems)
      map['Expression System'] = info.expressionSystems.join(', ');
    return new DG.Widget(ui.tableFromMap(map));
  } catch (e: any) {
    return new DG.Widget(ui.divText(`Failed to load info for ${pdbId}`));
  }
}
