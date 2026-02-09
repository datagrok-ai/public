import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {DesirabilityProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';

import {deleteMpoProfile, loadMpoProfiles, MPO_PROFILE_CHANGED_EVENT, MPO_PROFILE_DELETED_EVENT,
  MPO_TEMPLATE_PATH, MpoProfileInfo} from './utils';

class MpoProfileManagerImpl {
  private profiles: MpoProfileInfo[] = [];
  private loaded = false;

  get items(): MpoProfileInfo[] {
    return this.profiles;
  }

  get existingNames(): Set<string> {
    return new Set(this.profiles.map((p) => p.name));
  }

  get existingFileNames(): Set<string> {
    return new Set(this.profiles.map((p) => p.fileName));
  }

  async load(): Promise<MpoProfileInfo[]> {
    this.profiles = await loadMpoProfiles();
    this.loaded = true;
    return this.profiles;
  }

  async ensureLoaded(): Promise<MpoProfileInfo[]> {
    if (!this.loaded)
      await this.load();
    return this.profiles;
  }

  prepareClone(profile: MpoProfileInfo): {profile: MpoProfileInfo; fileName: string} {
    const baseName = this.getBaseName(profile.name);
    const baseFileName = this.getBaseFileName(profile.fileName);

    const clone = structuredClone(profile);
    clone.name = this.getNextAvailable(
      baseName,
      this.existingNames,
      (b, n) => n ? `${b} (Copy ${n})` : `${b} (Copy)`,
    );

    const cloneFileName = this.getNextAvailable(
      baseFileName,
      this.existingFileNames,
      (b, n) => n ? `${b}-copy-${n}.json` : `${b}-copy.json`,
    );

    return {profile: clone, fileName: cloneFileName};
  }

  confirmDelete(profile: MpoProfileInfo, onDeleted?: () => void): void {
    ui.dialog('Delete profile')
      .add(ui.divText(`Are you sure you want to delete profile "${profile.name}"?`))
      .onOK(async () => {
        try {
          await deleteMpoProfile(profile);
          this.profiles = this.profiles.filter((p) => p.fileName !== profile.fileName);
          this.fireChanged();
          grok.events.fireCustomEvent(MPO_PROFILE_DELETED_EVENT, {fileName: profile.fileName});
          onDeleted?.();
        } catch (e) {
          grok.shell.error(`Failed to delete profile "${profile.name}": ${e instanceof Error ? e.message : e}`);
        }
      })
      .show();
  }

  async save(profile: DesirabilityProfile, fileName: string): Promise<boolean> {
    try {
      await grok.dapi.files.writeAsText(`${MPO_TEMPLATE_PATH}/${fileName}`, JSON.stringify(profile));
      await this.load();
      this.fireChanged();
      grok.shell.info(`Profile "${profile.name}" saved.`);
      return true;
    } catch (e) {
      grok.shell.error(`Failed to save profile: ${e instanceof Error ? e.message : e}`);
      return false;
    }
  }

  fireChanged(): void {
    grok.events.fireCustomEvent(MPO_PROFILE_CHANGED_EVENT, {});
  }

  private getBaseName(name: string): string {
    return name.replace(/ \(Copy(?: \d+)?\)$/, '');
  }

  private getBaseFileName(fileName: string): string {
    return fileName.replace(/(-copy(?:-\d+)?)?\.json$/i, '');
  }

  private getNextAvailable(
    base: string,
    existing: Set<string>,
    format: (base: string, num?: number) => string,
  ): string {
    const first = format(base);
    if (!existing.has(first))
      return first;

    let num = 2;
    while (existing.has(format(base, num)))
      num++;

    return format(base, num);
  }
}

export const MpoProfileManager = new MpoProfileManagerImpl();
