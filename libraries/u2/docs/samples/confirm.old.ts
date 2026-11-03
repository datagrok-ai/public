import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

interface Campaign {
  code: string;
  name: string;
}

export function deleteCampaign(c: Campaign, remove: (code: string) => Promise<void>): void {
  const error = ui.divText('');
  error.style.color = 'var(--failure)';
  const dialog = ui.dialog('Delete campaign')
    .add(ui.divH([ui.iconFA('trash-alt'), ui.divV([
      ui.divText(`Are you sure you want to delete campaign ${c.name} (${c.code})?`),
      ui.divText('Its molecules and files go with it. This cannot be undone.')])]))
    .add(error);
  let busy = false;
  dialog.addButton('Delete', async () => {
    if (busy)
      return;
    busy = true;
    error.textContent = '';
    ui.setUpdateIndicator(dialog.root, true, 'Deleting...');
    try {
      await remove(c.code);
      dialog.close();
      grok.shell.info(`Campaign ${c.code} deleted`);
    }
    catch (e) {
      error.textContent = e instanceof Error ? e.message : String(e);
    }
    finally {
      busy = false;
      ui.setUpdateIndicator(dialog.root, false);
    }
  });
  dialog.show();
}
