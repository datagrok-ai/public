import * as grok from 'datagrok-api/grok';
import {AsyncValue, CatalogTable, Control, Form, TextInput, button, div, each, signal} from '@datagrok-libraries/u2';
import './good.css';

type Campaign = {code: string, name: string, status: string, kickoff: Date,
  owner: {login: string, friendlyName: string}, reference: {id: string, smiles: string}};

const day = (d: Date) => `${d.getFullYear()}/${d.getMonth() + 1}/${d.getDate()}`;

function details(campaign: Campaign, open: (campaign: Campaign) => unknown): HTMLElement {
  const structure = new TextInput({label: 'Structure', value: campaign.reference.smiles, readOnly: true});
  structure.addOptions(grok.chem.drawMolecule(campaign.reference.smiles, 160, 100));
  return new Form().addAll([
    new TextInput({label: 'Code', value: campaign.code, readOnly: true}),
    new TextInput({label: 'Name', value: campaign.name, readOnly: true}),
    new TextInput({label: 'Owner', value: campaign.owner.friendlyName, readOnly: true}),
    structure,
  ]).addButtons((row) => row.append(button('OPEN', () => open(campaign), {primary: true}))).root;
}

export function campaignsHome(load: () => Promise<Campaign[]>, open: (campaign: Campaign) => unknown): Control {
  const home = new Control();
  home.root.classList.add('campaigns-home');
  home.build(() => {
    const chosen = signal<Campaign | null>(null);
    const table = new CatalogTable<Campaign>({
      name: 'campaigns', label: 'Campaigns', items: new AsyncValue(() => load()), keyOf: (c) => c.code,
      columns: [
        {name: 'Name', value: (c) => c.name, pinned: true, link: 'Show the details'},
        {name: 'Code', value: (c) => c.code},
        {name: 'Status', value: (c) => c.status},
        {name: 'Kickoff', value: (c) => c.kickoff, render: (c) => day(c.kickoff)},
        {name: 'Owner', value: (c) => c.owner.friendlyName},
      ],
      groupings: [
        {name: 'Status', key: (c) => c.status},
        {name: 'Kickoff', key: (c) => day(c.kickoff)},
        {name: 'Owner', key: (c) => c.owner.login, groupLabel: (_, first) => first.owner.friendlyName},
      ],
      state: {groupBy: 'Status'},
      onActivate: (c) => chosen.value = c,
    });
    const pane = div([], 'campaigns-home-details');
    each(pane, () => chosen.value === null ? [] : [chosen.value], (c) => c.code, (c) => details(c, open));
    home.root.append(table.root, pane);
  });
  return home;
}
