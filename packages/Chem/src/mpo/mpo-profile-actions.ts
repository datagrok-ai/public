import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {DesirabilityProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {generateMpoFileName, getNextAvailable} from '@datagrok-libraries/statistics/src/mpo/utils';

import {mpoProfileStore, parseMpoProfile, toPlainProfile} from './mpo-profile-store';
import {MpoProfileInfo, MpoProfileRef} from './utils';

enum UploadConflictChoice {
  Replace = 'Replace existing profile',
  KeepBoth = 'Keep both profiles',
}

function takenNames(): Set<string> {
  return new Set(mpoProfileStore.items.map((p) => p.name));
}

function profileErrorMessage(e: unknown, action: string, profileName: string): string {
  return e instanceof DG.DomainValidationError && e.isDuplicate ? `A profile named "${profileName}" already exists.` :
    e instanceof DG.DomainVersionConflictError ? 'This profile was changed by someone else. Reload and try again.' :
      e instanceof DG.DomainForbiddenError ? 'You do not have permission to change MPO profiles.' :
        `Failed to ${action} profile "${profileName}": ${e instanceof Error ? e.message : e}`;
}

export async function saveProfileInteractive(profile: DesirabilityProfile, ref?: MpoProfileRef | null):
    Promise<MpoProfileRef | null> {
  try {
    const saved = await mpoProfileStore.save(profile, ref);
    grok.shell.info(`Profile "${profile.name}" saved.`);
    return saved;
  } catch (e) {
    grok.shell.error(profileErrorMessage(e, 'save', profile.name));
    return null;
  }
}

export function confirmDeleteProfile(profile: MpoProfileInfo, onDeleted?: () => void): void {
  ui.dialog('Delete profile')
    .add(ui.divText(`Are you sure you want to delete profile "${profile.name}"?`))
    .onOK(async () => {
      try {
        await mpoProfileStore.delete(profile.id);
        onDeleted?.();
      } catch (e) {
        grok.shell.error(profileErrorMessage(e, 'delete', profile.name));
      }
    })
    .show();
}

export function prepareProfileClone(profile: MpoProfileInfo): DesirabilityProfile {
  const clone = toPlainProfile(profile);
  clone.name = getNextAvailable(
    profile.name.replace(/ \(Copy(?: \d+)?\)$/, ''),
    takenNames(),
    (b, n) => n ? `${b} (Copy ${n})` : `${b} (Copy)`,
  );
  return clone;
}

export function uploadProfile(): void {
  DG.Utils.openFile({accept: '.json', open: async (file) => {
    await importProfile(await file.text(), file.name.replace(/\.json$/i, ''));
  }});
}

export async function importProfile(text: string, fallbackName: string): Promise<boolean> {
  const result = parseMpoProfile(text, fallbackName);
  if ('error' in result) {
    grok.shell.warning(`Upload failed: ${result.error}`);
    return false;
  }
  const profile = result.profile;

  const existing = await mpoProfileStore.findByName(profile.name);
  if (existing) {
    const choice = await promptUploadConflict(profile.name);
    if (choice === null)
      return false;
    if (choice === UploadConflictChoice.Replace)
      return await saveProfileInteractive(profile, existing) != null;
    await mpoProfileStore.ensureLoaded();
    profile.name = getNextAvailable(profile.name, takenNames(), (b, n) => n ? `${b} (${n})` : b);
  }
  return await saveProfileInteractive(profile) != null;
}

function promptUploadConflict(profileName: string): Promise<UploadConflictChoice | null> {
  return new Promise((resolve) => {
    const choice = ui.input.radio('', {value: UploadConflictChoice.KeepBoth,
      items: [UploadConflictChoice.Replace, UploadConflictChoice.KeepBoth], nullable: false});
    const dlg = ui.dialog('Upload options')
      .add(ui.divText(`A profile named "${profileName}" already exists.`))
      .add(choice)
      .addButton('Upload', () => {
        resolve(choice.value as UploadConflictChoice);
        dlg.close();
      })
      .onCancel(() => resolve(null))
      .show();
  });
}

export function downloadProfile(profile: MpoProfileInfo): void {
  DG.Utils.download(generateMpoFileName(profile.name, new Set()),
    JSON.stringify(toPlainProfile(profile), null, 2), 'application/json');
}
