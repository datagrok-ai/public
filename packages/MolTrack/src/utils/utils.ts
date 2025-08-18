import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { registerAssays, registerBulk, updateMolTrackSchema } from '../package';
import { GITHUB_BASE_URL, Scope } from './constants';

export async function fetchSchema(fileName: string): Promise<any> {
  const url = `${GITHUB_BASE_URL}${fileName}`;

  const response = await fetch(url);
  if (!response.ok)
    throw new Error(`Failed to fetch schema for ${fileName}: ${response.statusText}`);
  return await response.json();
}

export async function fetchCsv(scope: Scope): Promise<string> {
  const url = `${GITHUB_BASE_URL}${scope}.csv`;
  const res = await fetch(url);
  if (!res.ok)
    throw new Error(`Failed to fetch CSV for ${scope}: ${res.statusText}`);
  return await res.text();
}

export async function updateAllMolTrackSchemas(): Promise<void> {
  for (const scope of Object.values(Scope)) {
    try {
      const schemaPayload = await fetchSchema(`${scope}_schema.json`);
      await updateMolTrackSchema(JSON.stringify(schemaPayload));
    } catch (err) {
      console.error(`Error updating ${scope}:`, err);
    }
  }
}

export async function registerAllData(): Promise<void> {
  for (const scope of Object.values(Scope)) {
    if (scope !== Scope.ASSAYS) {
      const csvText = await fetchCsv(scope);
      const fileInfo = DG.FileInfo.fromString(`${scope}.csv`, csvText);
      await registerBulk(fileInfo, scope, '', 'reject_row');
    }
  }
}

export async function registerAssayData(): Promise<void> {
  const assayPayload = await fetchSchema(`${Scope.ASSAYS}.json`);
  await registerAssays(JSON.stringify(assayPayload));
}

export function createPath(viewName: string) {
  let path = `${MOLTRACK_APP_PATH}/`;
  path += encodeURIComponent(viewName);
  return path;
}

const MOLTRACK_APP_PATH: string = 'apps/MolTrack';
