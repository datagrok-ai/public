import {CatalogTable} from '@datagrok-libraries/u2';

type Campaign = {code: string, name: string, kickoff: Date, owner: {login: string, friendlyName: string}};

const day = (d: Date) => `${d.getFullYear()}/${d.getMonth() + 1}/${d.getDate()}`;

export function campaignCatalog(campaigns: Campaign[]): CatalogTable<Campaign> {
  return new CatalogTable<Campaign>({
    name: 'campaigns', items: campaigns, keyOf: (c) => c.code,
    columns: [
      {name: 'Name', value: (c) => c.name, pinned: true},
      {name: 'Code', value: (c) => c.code},
      {name: 'Kickoff', value: (c) => c.kickoff, render: (c) => day(c.kickoff)},
      {name: 'Owner', value: (c) => c.owner.friendlyName},
    ],
    groupings: [
      {name: 'Kickoff', key: (c) => day(c.kickoff)},
      {name: 'Owner', key: (c) => c.owner.login, groupLabel: (_, first) => first.owner.friendlyName},
    ],
  });
}
