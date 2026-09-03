import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {CURRENT_MPO_VERSION, DESIRABILITY_PROFILE_TYPE, DesirabilityProfile, isDesirabilityProfile,
  lockProfileRanges, migrateProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';

import {_package} from '../package';
import {mpoDb, ProfileAggregation, ProfileInsert, ProfileRow} from '../generated/db';
import {MPO_PROFILE_CHANGED_EVENT, MPO_PROFILE_DELETED_EVENT, MpoProfileInfo, MpoProfileRef,
  UNTITLED_PROFILE} from './utils';

const MPO_PROFILE_LIMIT = 100;
const MPO_FOLDER = 'mpo';

export type MpoProfileParseResult = {profile: DesirabilityProfile} | {error: string};

export function parseMpoProfile(text: string, fallbackName?: string): MpoProfileParseResult {
  let parsed: unknown;
  try {
    parsed = JSON.parse(text);
  } catch (e) {
    return {error: `invalid JSON: ${e instanceof Error ? e.message : e}`};
  }
  if (!isDesirabilityProfile(parsed))
    return {error: 'not an MPO desirability profile'};
  if ((parsed.version ?? 0) > CURRENT_MPO_VERSION)
    return {error: 'created with a newer version of the application'};
  migrateProfile(parsed);
  if (!parsed.name && fallbackName)
    parsed.name = fallbackName;
  return {profile: parsed};
}

export function toRow(p: DesirabilityProfile): ProfileInsert {
  return {
    name: p.name?.trim() || UNTITLED_PROFILE,
    description: p.description,
    aggregation: (p.aggregation ?? null) as ProfileAggregation | undefined,
    format_version: CURRENT_MPO_VERSION,
    properties: JSON.stringify(p.properties),
  };
}

export function fromRow(r: ProfileRow): MpoProfileInfo {
  const version = r.format_version ?? 0;
  if (version > CURRENT_MPO_VERSION)
    throw new Error('created with a newer version of the application');

  const profile = migrateProfile({
    type: DESIRABILITY_PROFILE_TYPE,
    version,
    name: r.name,
    description: r.description ?? '',
    aggregation: r.aggregation,
    properties: JSON.parse(r.properties),
  });
  return {...profile, id: r.id, rowVersion: r.version};
}

export function toPlainProfile(profile: MpoProfileInfo): DesirabilityProfile {
  const {id: _id, rowVersion: _rowVersion, ...plain} = profile;
  return plain;
}

export class MpoProfileStore {
  private profiles: MpoProfileInfo[] = [];
  private loaded = false;

  get items(): MpoProfileInfo[] {
    return this.profiles;
  }

  async load(): Promise<MpoProfileInfo[]> {
    const rows = await mpoDb.profiles.query().orderBy('name').top(MPO_PROFILE_LIMIT);
    const profiles: MpoProfileInfo[] = [];
    for (const row of rows) {
      try {
        profiles.push(fromRow(row));
      } catch (e) {
        _package.logger.warning(`MPO profile "${row.name}" skipped: ${e instanceof Error ? e.message : e}`);
      }
    }
    this.profiles = profiles;
    this.loaded = true;
    return this.profiles;
  }

  async ensureLoaded(): Promise<MpoProfileInfo[]> {
    if (!this.loaded)
      await this.load();
    return this.profiles;
  }

  async findByName(name: string): Promise<MpoProfileRef | null> {
    const row = await mpoDb.profiles.getByKey({name});
    return row ? {id: row.id, rowVersion: row.version} : null;
  }

  /** Inserts, or updates the row `ref` identifies (a stale `ref.rowVersion` throws
   *  DG.DomainVersionConflictError). Normalizes the profile in place — migrated to the
   *  current format with ranges locked — so the caller's copy matches what was stored. */
  async save(profile: DesirabilityProfile, ref?: MpoProfileRef | null): Promise<MpoProfileRef> {
    migrateProfile(profile);
    lockProfileRanges(profile);
    let saved: MpoProfileRef;
    if (ref) {
      const r = await mpoDb.profiles.update(ref.id, toRow(profile), {version: ref.rowVersion});
      saved = {id: ref.id, rowVersion: r.version};
    } else {
      const [r] = await mpoDb.profiles.insert(toRow(profile), {errorOnDuplicate: true});
      saved = {id: r.id, rowVersion: r.version ?? 1};
    }
    await this.refresh();
    return saved;
  }

  async delete(id: string): Promise<void> {
    await mpoDb.profiles.delete(id);
    this.profiles = this.profiles.filter((p) => p.id !== id);
    grok.events.fireCustomEvent(MPO_PROFILE_CHANGED_EVENT, {});
    grok.events.fireCustomEvent(MPO_PROFILE_DELETED_EVENT, {id});
  }

  // Listeners re-render from `items`, so the cache must be refreshed before the event fires.
  private async refresh(): Promise<void> {
    await this.load();
    grok.events.fireCustomEvent(MPO_PROFILE_CHANGED_EVENT, {});
  }

  /** Grants all users access, then seeds the profiles from the package's `mpo` AppData folder;
   *  `batch` skips names that already exist, so re-running never touches saved profiles. */
  async seedDefaults(): Promise<string> {
    for (const permission of ['View', 'Edit', 'Delete'] as const)
      await mpoDb.profiles.grant(DG.Group.defaultGroupsIds['All users'], permission);

    const rows: ProfileInsert[] = [];
    for (const f of (await _package.files.list(MPO_FOLDER)).filter((f) => f.isFile)) {
      const result = parseMpoProfile(await _package.files.readAsText(`${MPO_FOLDER}/${f.name}`),
        f.name.replace(/\.json$/i, ''));
      if ('profile' in result)
        rows.push(toRow(result.profile));
      else
        _package.logger.warning(`MPO seed skipped ${f.name}: ${result.error}`);
    }

    const report = await mpoDb.profiles.batch(rows, {allOrNothing: false});
    await this.refresh();
    return `MPO profiles: ${report.inserted} added, ${report.skipped} already existed` +
      (report.errorCount ? `, ${report.errorCount} failed` : '');
  }
}

export const mpoProfileStore = new MpoProfileStore();
