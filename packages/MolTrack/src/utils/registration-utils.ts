import * as DG from 'datagrok-api/dg';

import {registerAssays, registerBulk, registerMolTrackProperties} from '../package';
import {Scope} from './constants';
import {fetchCsv, fetchSchema} from './fetch-utils';

export async function updateAllMolTrackSchemas(): Promise<void> {
  for (const scope of Object.values(Scope)) {
    try {
      const schemaPayload = await fetchSchema(`${scope}_schema.json`);
      await registerMolTrackProperties(schemaPayload);
    } catch (err) {
      console.error(`Error updating ${scope}:`, err);
    }
  }
}

export async function registerAllData(): Promise<void> {
  for (const scope of Object.values(Scope)) {
    if (scope !== Scope.ASSAYS) {
      const csvText = await fetchCsv(scope);
      const mappingJson = await fetchSchema(`${scope}_mapping.json`);
      const fileInfo = DG.FileInfo.fromString(`${scope}.csv`, csvText);
      await registerBulk(fileInfo, scope, mappingJson, 'reject_row');
    }
  }
}

export async function registerAssayData(): Promise<void> {
  const assayPayload = await fetchSchema(`${Scope.ASSAYS}.json`);
  await registerAssays(assayPayload);
}
