import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {MpoProfileCreateView} from './mpo-create-profile';
import {deleteMpoProfile, MPO_PROFILE_CHANGED_EVENT, MpoProfileInfo} from './utils';

export function cloneMpoProfile(profile: MpoProfileInfo): void {
  const clone = structuredClone(profile);
  clone.name = `${profile.name} (Copy)`;

  const cloneFileName = profile.fileName.replace(/\.json$/i, '-copy.json');
  const view = new MpoProfileCreateView(clone, false, cloneFileName);
  grok.shell.v = grok.shell.addView(view.view);
}

export function confirmDeleteMpoProfile(profile: MpoProfileInfo, onDeleted?: () => void): void {
  ui.dialog('Delete profile')
    .add(ui.divText(`Are you sure you want to delete profile "${profile.name}"?`))
    .onOK(async () => {
      try {
        await deleteMpoProfile(profile);
        grok.events.fireCustomEvent(MPO_PROFILE_CHANGED_EVENT, {});
        onDeleted?.();
      } catch (e) {
        grok.shell.error(`Failed to delete profile "${profile.name}": ${e instanceof Error ? e.message : e}`);
      }
    })
    .show();
}
