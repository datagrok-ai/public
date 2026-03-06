import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {calculateColumns} from '@datagrok-libraries/statistics/src/compute-functions/execution';
import {TemplateCompute, TemplateFunction} from '@datagrok-libraries/statistics/src/compute-functions/types';

import {
  DEFAULT_AGGREGATION,
  DESIRABILITY_PROFILE_TYPE,
  DesirabilityProfile,
  PropertyDesirability,
  WeightedAggregation,
  createDefaultNumerical,
} from '@datagrok-libraries/statistics/src/mpo/mpo';

export {MPO_PROFILE_CHANGED_EVENT, MPO_PROFILE_DELETED_EVENT} from '@datagrok-libraries/statistics/src/mpo/utils';

export type MpoProfileInfo = DesirabilityProfile & {
  fileName: string;
};

export enum MpoPathMode {
  List = 'list',
  Edit = 'edit',
  Create = 'create',
}

export const MPO_TEMPLATE_PATH = 'System:AppData/Chem/mpo';
export const MPO_PATH = 'MPOProfiles';
export const MAX_MPO_PROPERTIES = 20;

export async function loadMpoProfiles(): Promise<MpoProfileInfo[]> {
  const files = await grok.dapi.files.list(MPO_TEMPLATE_PATH);
  const profiles: MpoProfileInfo[] = [];

  for (const file of files) {
    try {
      const text = await grok.dapi.files.readAsText(`${MPO_TEMPLATE_PATH}/${file.name}`);
      const content = JSON.parse(text) as DesirabilityProfile;

      profiles.push({
        type: DESIRABILITY_PROFILE_TYPE,
        fileName: file.name,
        name: content.name ?? file.name.replace(/\.json$/i, ''),
        description: content.description ?? '',
        aggregation: content.aggregation,
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
  silent: boolean = false,
  processed: boolean = false,
): Promise<string[]> {
  // Execute per-property functions for missing columns
  await executePropertyFunctions(df, profile, columnMapping, silent);

  const mappedProperties: Record<string, PropertyDesirability> = {};
  for (const [propName, prop] of Object.entries(profile.properties)) {
    const columnName = columnMapping[propName] ?? propName;
    mappedProperties[columnName] = prop;
  }

  const profileName = profile.name || 'MPO';
  if (!processed) {
    const existingCol = df.col(profileName);
    if (existingCol)
      df.columns.remove(existingCol, false);
  }

  const resolvedAggregation = aggregation ?? profile.aggregation ?? DEFAULT_AGGREGATION;
  const call = await DG.Func.find({package: 'Chem', name: 'mpoTransformFunction'})[0].prepare({
    df,
    profileName,
    currentProperties: JSON.stringify(mappedProperties),
    aggregation: resolvedAggregation,
    silent,
  }).call(undefined, undefined, {processed});

  const result = call.getOutputParamValue() as DG.Column[];

  // Temporary fix until proper support for list<column> is implemented
  const colList = DG.DataFrame.fromColumns(result).columns;
  await DG.Func.find({package: 'Chem', name: 'mpoCalculate'})[0].prepare({
    df, columns: colList, profileName, aggregation: resolvedAggregation,
  }).call(undefined, undefined, {processed});
  return df.col(profileName) ? [profileName] : [];
}

export function parseCallString(callStr: string): TemplateFunction | null {
  const parenIdx = callStr.indexOf('(');
  const qualifiedName = parenIdx === -1 ? callStr : callStr.substring(0, parenIdx);
  const [pkg, name] = qualifiedName.split(':');
  if (!pkg || !name)
    return null;

  if (parenIdx === -1)
    return {package: pkg, name, args: {}};

  try {
    const argsJson = callStr.substring(parenIdx + 1, callStr.lastIndexOf(')'));
    const func = DG.Func.find({package: pkg, name})[0];
    if (!func)
      return {package: pkg, name, args: {}};
    const argValues = JSON.parse(`[${argsJson}]`);
    const args: Record<string, any> = {};
    for (let i = 0; i < argValues.length && i + 2 < func.inputs.length; i++)
      args[func.inputs[i + 2].name] = argValues[i];
    return {package: pkg, name, args};
  } catch (e) {
    console.warn(`Failed to parse function call string: ${callStr}`, e);
    return {package: pkg, name, args: {}};
  }
}

export function templateFromCallString(callStr: string): {compute: TemplateCompute} | undefined {
  const tf = parseCallString(callStr);
  if (!tf)
    return undefined;
  return {compute: {descriptors: {enabled: false, args: []}, functions: [tf]}};
}

export async function executePropertyFunctions(
  df: DG.DataFrame,
  profile: DesirabilityProfile,
  columnMapping: Record<string, string | null>,
  silent: boolean,
): Promise<void> {
  const molCol = [...df.columns].find((c) => c.semType === DG.SEMTYPE.MOLECULE);
  if (!molCol)
    return;

  // Collect unique functions for missing columns, parsing call strings into name + args
  const callStrings: string[] = [];
  for (const [propName, prop] of Object.entries(profile.properties)) {
    if (!prop.function)
      continue;
    const columnName = columnMapping[propName] ?? propName;
    if (df.col(columnName))
      continue;
    if (!callStrings.includes(prop.function))
      callStrings.push(prop.function);
  }

  if (callStrings.length === 0)
    return;

  const externals: Record<string, Record<string, any>> = {};
  for (const callStr of callStrings) {
    const tf = parseCallString(callStr);
    if (tf)
      externals[`${tf.package}:${tf.name}`] = tf.args;
  }

  const colsBefore = new Set(df.columns.names());
  await calculateColumns({descriptors: [], externals}, df, molCol.name);

  // Auto-map newly added columns to unmapped properties
  const newCols = df.columns.names().filter((n) => !colsBefore.has(n));
  if (newCols.length === 0)
    return;
  const newColsLower = new Map(newCols.map((n) => [n.toLowerCase(), n]));
  for (const [propName, prop] of Object.entries(profile.properties)) {
    if (!prop.function)
      continue;
    const mapped = columnMapping[propName];
    if (mapped != null && df.col(mapped))
      continue;
    const match = newColsLower.get(propName.toLowerCase());
    if (match)
      columnMapping[propName] = match;
  }
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

export function createDefaultProfile(): DesirabilityProfile {
  return {
    type: DESIRABILITY_PROFILE_TYPE,
    name: '',
    description: '',
    properties: {},
  };
}

export function createProfileForDf(df: DG.DataFrame): DesirabilityProfile {
  const props: {[key: string]: PropertyDesirability} = {};
  let count = 0;
  for (const col of df.columns.numerical) {
    if (count >= MAX_MPO_PROPERTIES)
      break;
    props[col.name] = createDefaultNumerical(1, col.min, col.max);
    count++;
  }
  return {type: DESIRABILITY_PROFILE_TYPE, name: '', description: '', properties: props};
}

export function mergeProfileWithDf(existing: DesirabilityProfile, df: DG.DataFrame): DesirabilityProfile {
  const merged: DesirabilityProfile = {
    ...existing,
    properties: structuredClone(existing.properties),
  };

  const existingLower = new Set(Object.keys(merged.properties).map((n) => n.toLowerCase()));
  for (const col of df.columns.numerical) {
    if (Object.keys(merged.properties).length >= MAX_MPO_PROPERTIES)
      break;
    if (!existingLower.has(col.name.toLowerCase()))
      merged.properties[col.name] = createDefaultNumerical(1, col.min, col.max);
  }

  return merged;
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
