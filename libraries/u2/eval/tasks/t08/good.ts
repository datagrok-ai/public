import {CatalogTable} from '@datagrok-libraries/u2';

type Campaign = {code: string, name: string, status: string, created: Date};

export function campaignsTable(campaigns: Campaign[]): CatalogTable<Campaign> {
  return new CatalogTable<Campaign>({
    name: 'campaigns', label: 'Campaigns', items: campaigns, keyOf: (c) => c.code,
    columns: [
      {name: 'Name', value: (c) => c.name, pinned: true},
      {name: 'Code', value: (c) => c.code},
      {name: 'Status', value: (c) => c.status},
      {name: 'Created', value: (c) => c.created},
    ],
    groupings: [{name: 'Status', key: (c) => c.status}],
    state: {groupBy: 'Status'},
  });
}
