import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  DesirabilityProfile,
  mpo,
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
export const MPO_PROFILE_CHANGED_EVENT = 'chem-mpo-profile-changed';

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

export function profileForEditing(profile: MpoProfileInfo): MpoProfileInfo {
  const copy = {...profile};
  delete copy.file;
  return copy;
}

export type MpoCalculationResult = {
  columnNames: string[];
  warnings: string[];
  error?: string;
};

export function calculateMpoCore(
  df: DG.DataFrame,
  profileName: string,
  currentProperties: Record<string, PropertyDesirability>,
  aggregation: WeightedAggregation,
): MpoCalculationResult {
  const columns: DG.Column[] = [];
  const warnings: string[] = [];

  for (const [propertyName, desirability] of Object.entries(currentProperties)) {
    const column = df.columns.byName(propertyName);
    if (!column) {
      warnings.push(`Column "${propertyName}" from profile not found in table. Skipping.`);
      continue;
    }
    column.setTag('desirabilityTemplate', JSON.stringify(desirability));
    columns.push(column);
  }

  if (columns.length === 0) {
    return {
      columnNames: [],
      warnings,
      error: 'No valid columns found matching the profile properties. Cannot calculate MPO score.',
    };
  }

  try {
    const resultCol = mpo(df, columns, profileName, aggregation);
    return {columnNames: resultCol ? [resultCol.name] : [], warnings};
  } catch (e) {
    console.error('MPO Calculation Error:', e);
    return {
      columnNames: [],
      warnings,
      error: `MPO calculation failed: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}

export async function computeMpo(
  df: DG.DataFrame,
  profile: DesirabilityProfile,
  columnMapping: Record<string, string | null>,
  aggregation?: WeightedAggregation,
  silent: boolean = false,
): Promise<string[]> {
  const mappedProperties: Record<string, PropertyDesirability> = {};
  for (const [propName, prop] of Object.entries(profile.properties)) {
    const columnName = columnMapping[propName] ?? propName;
    mappedProperties[columnName] = prop;
  }

  if (silent) {
    const result = calculateMpoCore(df, profile.name ?? 'MPO', mappedProperties, aggregation ?? 'Average');
    return result.columnNames;
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
  const mpoPath = `${basePath}/${MPO_PATH}`;

  if (!url.pathname.startsWith(basePath))
    url.pathname = basePath;

  url.pathname = url.pathname.replace(/\/$/, '');

  switch (mode) {
  case MpoPathMode.Edit:
    url.pathname = mpoPath;
    url.searchParams.set('profileId', encodeURIComponent(profileName ?? ''));
    break;

  case MpoPathMode.Create:
    url.pathname = `${mpoPath}/create-profile`;
    url.searchParams.delete('profileId');
    break;

  case MpoPathMode.List:
  default:
    url.pathname = mpoPath;
    url.searchParams.delete('profileId');
    break;
  }

  const newPath = url.pathname + (url.search ? url.search : '');
  if (newPath !== window.location.pathname + window.location.search)
    window.history.replaceState({}, '', newPath);

  view.path = newPath;
}

export function deepEqual<T>(current: T, original: T): boolean {
  if (current === original)
    return true;

  if (typeof original !== 'object')
    return false;

  if (Array.isArray(original)) {
    if (!Array.isArray(current) || current.length !== original.length)
      return false;
    for (let i = 0; i < original.length; i++) {
      if (!deepEqual(current[i], original[i]))
        return false;
    }
    return true;
  }

  for (const key in original) {
    if (!deepEqual(current[key], original[key]))
      return false;
  }
  return true;
}
