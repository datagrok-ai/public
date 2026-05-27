import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {checkPackage} from '../utils/elemental-analysis-utils';

import {
  DEFAULT_AGGREGATION,
  DESIRABILITY_PROFILE_TYPE,
  DesirabilityProfile,
  MpoCalculator,
  MpoResult,
  PropertyDesirability,
  WeightedAggregation,
  createDefaultNumerical,
  migrateProfile,
} from '@datagrok-libraries/statistics/src/mpo/mpo';

export {MPO_PROFILE_CHANGED_EVENT, MPO_PROFILE_DELETED_EVENT} from '@datagrok-libraries/statistics/src/mpo/utils';

export const UNTITLED_PROFILE = 'Untitled Profile';

export enum MpoMethod {
  Manual = 'Manual',
  DataDriven = 'Data-driven',
}

export function isEdaPackageInstalled(): boolean {
  if (!checkPackage('EDA', 'getPmpoAppItems')) {
    grok.shell.warning('EDA package is not installed');
    return false;
  }
  return true;
}

export type MpoProfileInfo = DesirabilityProfile & {
  fileName: string;
};

export type MpoSaveResult = {saved: boolean; fileName: string};

export enum MpoPathMode {
  List = 'list',
  Edit = 'edit',
  Create = 'create',
}

export enum MpoUploadConflictAction {
  Replace = 'replace',
  KeepBoth = 'keep-both',
  Cancel = 'cancel',
}

export const MPO_TEMPLATE_PATH = 'System:AppData/Chem/mpo';
export const MPO_PATH = 'MPOProfiles';
export const MPO_PROFILES_NAME = 'MPO Profiles';
export const MAX_MPO_PROPERTIES = 20;

export async function loadMpoProfiles(): Promise<MpoProfileInfo[]> {
  const files = await grok.dapi.files.list(MPO_TEMPLATE_PATH);
  const profiles: MpoProfileInfo[] = [];

  for (const file of files) {
    try {
      const text = await grok.dapi.files.readAsText(`${MPO_TEMPLATE_PATH}/${file.name}`);
      const content = migrateProfile(JSON.parse(text) as DesirabilityProfile);

      profiles.push({
        type: DESIRABILITY_PROFILE_TYPE,
        version: content.version,
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
  createDesirabilityColumns: boolean = false,
  preview: boolean = false,
  calculator?: MpoCalculator,
): Promise<string[]> {
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

  // Preview path (interactive context panel): compute directly via the shared helpers, skipping the
  // two DG.Func.prepare().call() round-trips + the join(df) action. The OK path below keeps the
  // function-call machinery so the transform stays recorded/replayable. A passed-in calculator caches
  // per-column desirability across edits; one-shot callers get a transient one.
  if (preview) {
    const columns = applyDesirabilityTags(df, mappedProperties, true);
    if (columns.length > 0) {
      const isDifferent = df.rowCount !== columns[0].length;
      const result = (calculator ?? new MpoCalculator())
        .compute(df, columns, profileName, resolvedAggregation, isDifferent, createDesirabilityColumns);
      for (const c of collectMpoResultColumns(df, result, isDifferent)) {
        const existing = df.col(c.name);
        if (existing)
          df.columns.remove(existing.name);
        df.columns.add(c);
      }
    }
    return df.col(profileName) ? [profileName] : [];
  }

  const call = await DG.Func.find({package: 'Chem', name: 'mpoTransformFunction'})[0].prepare({
    df,
    profileName,
    currentProperties: JSON.stringify(mappedProperties),
    aggregation: resolvedAggregation,
    silent,
  }).call(undefined, undefined, {processed});

  const result = call.getOutputParamValue() as DG.Column[];

  await DG.Func.find({package: 'Chem', name: 'mpoCalculate'})[0].prepare({
    df, columns: result, profileName, aggregation: resolvedAggregation, createDesirabilityColumns,
  }).call(undefined, undefined, {processed});

  return df.col(profileName) ? [profileName] : [];
}

/// Sets the `desirabilityTemplate` tag on each column named in `properties` and returns the tagged
/// columns. Shared by Chem:mpoTransformFunction (the recorded transform) and the interactive preview.
export function applyDesirabilityTags(
  df: DG.DataFrame,
  properties: Record<string, PropertyDesirability>,
  silent: boolean = false,
): DG.Column[] {
  const columns: DG.Column[] = [];
  for (const [columnName, desirability] of Object.entries(properties)) {
    const column = df.col(columnName);
    if (!column) {
      if (!silent)
        grok.shell.warning(`Column "${columnName}" not found. Skipping.`);
      continue;
    }
    column.setTag('desirabilityTemplate', JSON.stringify(desirability));
    columns.push(column);
  }
  if (columns.length === 0 && !silent)
    grok.shell.error('No valid columns found matching the profile properties.');
  return columns;
}

/// Picks the columns to join back into `df` from an mpo() / MpoCalculator result: the score column when
/// it's new (or the row count differs), plus all desirability columns. Shared by Chem:mpoCalculate
/// (the recorded transform) and the interactive preview.
export function collectMpoResultColumns(df: DG.DataFrame, result: MpoResult, isDifferent: boolean): DG.Column[] {
  const columns: DG.Column[] = [];
  if (result.scoreColumn && (!df.col(result.scoreColumn.name) || isDifferent))
    columns.push(result.scoreColumn);
  if (result.desirabilityColumns)
    columns.push(...result.desirabilityColumns);
  return columns;
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
    name: UNTITLED_PROFILE,
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
  return {type: DESIRABILITY_PROFILE_TYPE, name: UNTITLED_PROFILE, description: '', properties: props};
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

  if (typeof original !== 'object' || original === null ||
      typeof current !== 'object' || current === null)
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

export function setupMpoBreadcrumbs(view: DG.ViewBase, lastSegment: string): void {
  const breadcrumbs = ui.breadcrumbs(['Home', MPO_PROFILES_NAME, lastSegment]);

  breadcrumbs.onPathClick.subscribe((path) => {
    const clicked = path[path.length - 1];
    if (clicked === lastSegment)
      return;
    if (clicked === MPO_PROFILES_NAME) {
      const listView = Array.from(grok.shell.views).find((v) => v.name === MPO_PROFILES_NAME);
      if (listView)
        grok.shell.v = listView;
    }
  });

  const homeEl = breadcrumbs.root.firstElementChild;
  if (homeEl) {
    const homeIcon = ui.iconFA('home', () => grok.shell.v = DG.View.createByType(DG.VIEW_TYPE.HOME), 'Home');
    homeEl.replaceWith(homeIcon);
  }

  const viewNameRoot = view.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
  if (viewNameRoot) {
    viewNameRoot.textContent = '';
    viewNameRoot.appendChild(breadcrumbs.root);
  }
}
