import {Control, Toolbar, component, confirmDialog, notify} from '@datagrok-libraries/u2';

interface Campaign {
  code: string;
  name: string;
}

export async function deleteCampaign(c: Campaign, remove: (code: string) => Promise<void>): Promise<void> {
  if (await confirmDialog({title: 'Delete campaign', message: `Are you sure you want to delete campaign ${c.name} (${c.code})?`,
    detail: 'Its molecules and files go with it. This cannot be undone.', confirmText: 'Delete', destructive: true,
    onConfirm: () => remove(c.code)}))
    notify.info(`Campaign ${c.code} deleted`);
}

// Live demo
async function removeCampaign(code: string): Promise<void> {
  await new Promise((resolve) => setTimeout(resolve, 800));
  if (code === 'KHL-13')
    throw new Error('KHL-13 has submitted molecules: archive it instead');
}

export function render(host: HTMLElement): Control {
  const demo = component(() => new Toolbar({ariaLabel: 'Campaigns'})
    .addButton('Delete KHL-12...', () => deleteCampaign({code: 'KHL-12', name: 'PI3Kδ quick wins'}, removeCampaign),
      {icon: 'trash-alt'})
    .addButton('Delete KHL-13...', () => deleteCampaign({code: 'KHL-13', name: 'BTK follow-up'}, removeCampaign),
      {icon: 'trash-alt'}));
  host.append(demo.root);
  return demo;
}
