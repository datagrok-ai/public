import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import type {ShareTarget} from './sharing';
import {SpacePicker} from './space-picker';

// Older platform builds lack the spaces client entirely
export function isWorkspaceSharingAvailable(): boolean {
  return (grok.dapi as any).spaces != null;
}

// Dialog that links a run into an existing space, saving it to history first when needed
export async function shareRunToWorkspace(target: ShareTarget): Promise<void> {
  const picker = await SpacePicker.create();
  const runName = target.defaultName?.() ?? 'run';
  const unsaved = target.savedCallId?.() == null;
  const selectedLabel = ui.divText('No space selected', {style: {color: 'var(--grey-4)'}});
  const dlg = ui.dialog({title: 'Share to workspace'})
    .add(ui.divText(`Link "${runName}" into a space${unsaved ? ' (the run will be saved to history first)' : ''}`))
    .add(picker.root)
    .add(selectedLabel)
    .onOK(() => {
      const space = picker.selected;
      if (!space) return;
      // Deferred past this dialog's teardown: opening a chained centered dialog (the save form)
      // while this one is closing leaks the platform's d4-dialog-center-wrapper, which then
      // swallows every click on the page.
      setTimeout(async () => {
        const id = target.savedCallId?.() ?? await target.saveRun?.();
        if (!id) return;
        try {
          await grok.dapi.spaces.id(space.id).addEntity(id, true);
          grok.shell.info(`Linked into "${space.friendlyName}"`);
        } catch (e: any) {
          grok.shell.error(`Could not link into "${space.friendlyName}": ${e?.message ?? e}`);
        }
      });
    });
  dlg.show({center: true, width: 400});
  const ok = dlg.getButton('OK') as HTMLButtonElement | null;
  if (ok) ok.disabled = true;
  picker.onChanged = (space) => {
    selectedLabel.textContent = space ? `Selected: ${space.friendlyName || space.name}` : 'No space selected';
    if (ok) ok.disabled = space == null;
  };
}
