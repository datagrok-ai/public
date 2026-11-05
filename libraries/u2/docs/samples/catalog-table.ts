import {CatalogTable, Control, component, confirmDialog, divH, dot, h2, notify, signal, span, timestamp}
  from '@datagrok-libraries/u2';
import type {BadgeVariant} from '@datagrok-libraries/u2';

interface Campaign {
  code: string;
  name: string;
  template: string;
  status: 'In progress' | 'Submitted' | 'Archived';
  molecules: number;
  created: Date;
  author: string;
}

const STATUS: Record<Campaign['status'], BadgeVariant> = {'In progress': 'accent', 'Submitted': 'success', 'Archived': 'default'};

export function campaignCatalog(campaigns: Campaign[], open: (c: Campaign) => void): CatalogTable<Campaign> {
  const items = signal(campaigns);
  const remove = async (c: Campaign) => {
    await confirmDialog({title: 'Delete campaign', message: `Delete campaign ${c.name} (${c.code})?`,
      confirmText: 'Delete', destructive: true, onConfirm: () => {
        items.value = items.value.filter((x) => x !== c);
      }});
  };
  return new CatalogTable<Campaign>({
    name: 'campaigns', label: 'Campaigns', items, keyOf: (c) => c.code,
    columns: [
      {name: 'Name', value: (c) => c.name, pinned: true, link: 'Continue the campaign', maxWidth: '240px'},
      {name: 'Code', value: (c) => c.code},
      {name: 'Created', value: (c) => c.created, render: (c) => timestamp(c.created)},
      {name: 'Molecules', value: (c) => c.molecules, align: 'right'},
      {name: 'Status', value: (c) => c.status, render: (c) => divH([dot(STATUS[c.status]), span(c.status)])},
      {name: 'Author', value: (c) => c.author},
    ],
    groupings: [{name: 'Template', key: (c) => c.template}, {name: 'Status', key: (c) => c.status}],
    defaultColumns: ['Code', 'Created', 'Molecules', 'Status'],
    defaultOrder: (a, b) => b.created.getTime() - a.created.getTime(),
    actions: (c) => [{text: 'Delete campaign', icon: 'trash-alt', destructive: true, run: () => remove(c)}],
    onActivate: open,
    empty: 'No campaigns yet. Start one under New campaign.',
    footer: (shown) => `${shown.length} campaigns`,
    persist: 'build-with-u2.campaigns',
  });
}

// Live demo
const TEMPLATES = ['Kinase hit-to-lead', 'GPCR screen', 'Fragment follow-up'];
const STATUSES: Campaign['status'][] = ['In progress', 'Submitted', 'Archived'];
const AUTHORS = ['Ada Lovelace', 'Alan Turing', 'Grace Hopper'];

export function render(host: HTMLElement): Control {
  const demo = component(() => {
    const table = campaignCatalog(Array.from({length: 12}, (_, i) => ({code: `KHL-${i + 1}`,
      name: `${['PI3Kδ', 'BTK', 'JAK2', 'EGFR'][i % 4]} ${['quick wins', 'follow-up', 'selectivity', 'scale-up'][i % 3]}`,
      template: TEMPLATES[i % 3], status: STATUSES[i % 3 === 0 ? 0 : i % 4 === 1 ? 1 : 2], molecules: 40 + (i * 53) % 900,
      created: new Date(2026, i % 10, 1 + (i * 7) % 27), author: AUTHORS[i % 3]})), (c) => notify.info(`Opening ${c.name}`));
    return [divH([h2('Campaigns'), ...table.chrome()]), table];
  });
  host.append(demo.root);
  return demo;
}
