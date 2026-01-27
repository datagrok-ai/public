import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  DesirabilityProfile,
  PropertyDesirability,
  WeightedAggregation,
} from '@datagrok-libraries/statistics/src/mpo/mpo';

export type MpoProfileInfo = {
  fileName: string;
  name: string;
  description: string;
  properties: Record<string, PropertyDesirability>;
  file?: DG.FileInfo;
};

export enum MpoPathMode {
  List = 'list',
  Edit = 'edit',
  Create = 'create',
}

export const MPO_TEMPLATE_PATH = 'System:AppData/Chem/mpo';
export const MPO_PATH = 'Mpo';

export async function loadMpoProfiles(): Promise<MpoProfileInfo[]> {
  const files = await grok.dapi.files.list(MPO_TEMPLATE_PATH);
  const profiles: MpoProfileInfo[] = [];

  for (const file of files) {
    try {
      const text = await grok.dapi.files.readAsText(`${MPO_TEMPLATE_PATH}/${file.name}`);
      const content = JSON.parse(text) as DesirabilityProfile;

      profiles.push({
        file,
        fileName: file.name,
        name: content.name ?? file.name.replace(/\.json$/i, ''),
        description: content.description ?? '',
        properties: content.properties,
      });
    } catch (e) {
      grok.shell.warning(`Failed to load MPO profile "${file.name}": ${e instanceof Error ? e.message : e}`);
    }
  }

  return profiles;
}

export async function deleteMpoProfile(profile: MpoProfileInfo): Promise<void> {
  await grok.dapi.files.delete(`${MPO_TEMPLATE_PATH}/${profile.fileName}`);
}

export async function computeMpo(
  df: DG.DataFrame,
  profile: DesirabilityProfile,
  columnMapping: Record<string, string | null>,
  aggregation?: WeightedAggregation,
): Promise<string[]> {
  const mappedProperties: Record<string, PropertyDesirability> = {};
  for (const [propName, prop] of Object.entries(profile.properties)) {
    const columnName = columnMapping[propName] ?? propName;
    mappedProperties[columnName] = prop;
  }

  const [func] = await DG.Func.find({name: 'mpoTransformFunction'});
  const funcCall = await func.prepare({
    df,
    profileName: profile.name ?? 'MPO',
    currentProperties: mappedProperties,
    aggregation: aggregation ?? 'Average',
  }).call(undefined, undefined, {processed: false});

  return funcCall.getOutputParamValue() ?? [];
}

export function findSuitableProfiles(df: DG.DataFrame, profiles: MpoProfileInfo[]): MpoProfileInfo[] {
  const dfColumnNames = new Set(df.columns.names());

  return profiles.filter((p) => {
    return Object.keys(p.properties ?? {}).some((prop) =>
      dfColumnNames.has(prop),
    );
  });
}

export function updateMpoPath(
  view: DG.View,
  mode: MpoPathMode,
  profileName?: string,
): void {
  const url = new URL(window.location.href);
  const basePath = '/apps/Chem';

  if (!url.pathname.startsWith(basePath))
    url.pathname = basePath;

  url.pathname = url.pathname.replace(/\/$/, '');

  if (!url.pathname.includes(`/${MPO_PATH}`))
    url.pathname = `${basePath}/Mpo`;

  switch (mode) {
  case MpoPathMode.Edit: {
    const profileId = encodeURIComponent(profileName ?? '');
    url.searchParams.set('profileId', profileId);
    break;
  }

  case MpoPathMode.Create:
    url.pathname = `${basePath}/${MPO_PATH}/create-profile`;
    url.searchParams.delete('profileId');
    break;

  case MpoPathMode.List:
  default:
    url.pathname = `${basePath}/${MPO_PATH}`;
    url.searchParams.delete('profileId');
    break;
  }

  const newPath = url.pathname + (url.search ? url.search : '');
  if (newPath !== window.location.pathname + window.location.search)
    window.history.replaceState({}, '', newPath);

  view.path = newPath;
}
