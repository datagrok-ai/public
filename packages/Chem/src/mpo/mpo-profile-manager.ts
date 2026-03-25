import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {CURRENT_MPO_VERSION, DesirabilityProfile, isDesirabilityProfile, migrateProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {generateMpoFileName, getNextAvailable} from '@datagrok-libraries/statistics/src/mpo/utils';

import {deleteMpoProfile, loadMpoProfiles, MPO_PROFILE_CHANGED_EVENT, MPO_PROFILE_DELETED_EVENT,
  MPO_TEMPLATE_PATH, MpoProfileInfo, MpoSaveResult} from './utils';

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
    clone.name = getNextAvailable(
      baseName,
      this.existingNames,
      (b, n) => n ? `${b} (Copy ${n})` : `${b} (Copy)`,
    );

    const cloneFileName = getNextAvailable(
      baseFileName,
      this.existingFileNames,
      (b, n) => n ? `${b}-copy-${n}.json` : `${b}-copy.json`,
    );

    return {profile: clone, fileName: cloneFileName};
  }

  generateFileName(profileName: string): string {
    return generateMpoFileName(profileName, this.existingFileNames);
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

  async showSaveDialog(profile: DesirabilityProfile, existingFileName?: string | null): Promise<MpoSaveResult> {
    await this.ensureLoaded();
    const isExisting = !!existingFileName;
    return new Promise((resolve) => {
      const nameInput = ui.input.string('Name', {value: profile.name ?? '', nullable: false});
      const descInput = ui.input.textArea('Description', {value: profile.description ?? ''});
      const overrideInput = isExisting ? ui.input.bool('Override existing', {value: true}) : null;

      const inputs: HTMLElement[] = [nameInput.root, descInput.root];
      if (overrideInput)
        inputs.push(overrideInput.root);

      const dlg = ui.dialog({title: 'Save MPO Profile'})
        .add(ui.divV(inputs))
        .onOK(async () => {
          profile.name = nameInput.value!;
          profile.description = descInput.value || '';
          const fileName = isExisting && overrideInput?.value ? existingFileName! : this.generateFileName(profile.name);
          const saved = await this.save(profile, fileName);
          resolve({saved, fileName});
        })
        .onCancel(() => resolve({saved: false, fileName: ''}))
        .show();

      const okButton = dlg.getButton('OK');
      okButton.disabled = !nameInput.validate();
      nameInput.onInput.subscribe(() => okButton.disabled = !nameInput.validate());
    });
  }

  async save(profile: DesirabilityProfile, fileName: string): Promise<boolean> {
    try {
      await grok.dapi.files.writeAsText(`${MPO_TEMPLATE_PATH}/${fileName}`, JSON.stringify(profile));
      await this.load();
      this.fireChanged();
      grok.shell.info(`Profile "${profile.name}" saved as ${fileName}.`);
      return true;
    } catch (e) {
      grok.shell.error(`Failed to save profile: ${e instanceof Error ? e.message : e}`);
      return false;
    }
  }

  upload(): void {
    DG.Utils.openFile({accept: '.json', open: async (file) => {
      try {
        const text = await file.text();
        const parsed = JSON.parse(text);
        if (!isDesirabilityProfile(parsed)) {
          grok.shell.warning('Upload failed: not a valid MPO profile');
          return;
        }
        if ((parsed.version ?? 0) > CURRENT_MPO_VERSION) {
          grok.shell.warning('Upload failed: profile was created with a newer version of the application');
          return;
        }
        migrateProfile(parsed);
        if (!parsed.name)
          parsed.name = file.name.replace(/\.json$/i, '');
        await this.save(parsed, this.generateFileName(parsed.name));
      } catch (e) {
        grok.shell.warning('Upload failed: invalid JSON file');
      }
    }});
  }

  download(profile: MpoProfileInfo): void {
    const {fileName, ...data} = profile;
    DG.Utils.download(fileName, JSON.stringify(data, null, 2), 'application/json');
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
}

export const MpoProfileManager = new MpoProfileManagerImpl();
