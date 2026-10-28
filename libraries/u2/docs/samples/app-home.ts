import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {AppHeader, AsyncValue, CatalogTable, divH, h1} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

interface Campaign {
  code: string;
  name: string;
  template: string;
  status: string;
  created: string;
}

async function loadCampaigns(): Promise<Campaign[]> {
  return JSON.parse(await grok.dapi.files.readAsText('System:AppData/MyApp/campaigns.json'));
}

function openCampaign(c: Campaign): void {
  grok.shell.info(`Opening ${c.code}`);
}

//name: My App
//tags: app
//output: view result
export function myApp(): DG.ViewBase {
  return ViewHost.create({name: 'My App', content: () => {
    const campaigns = new AsyncValue(loadCampaigns);
    const table = new CatalogTable<Campaign>({
      name: 'campaigns', items: campaigns, keyOf: (c) => c.code,
      columns: [
        {name: 'Name', value: (c) => c.name, pinned: true, link: 'Continue the campaign'},
        {name: 'Code', value: (c) => c.code},
        {name: 'Created', value: (c) => c.created},
        {name: 'Status', value: (c) => c.status},
      ],
      groupings: [{name: 'Template', key: (c) => c.template}, {name: 'Status', key: (c) => c.status}],
      onActivate: openCampaign,
      empty: 'No campaigns yet.',
      persist: 'my-app.campaigns',
    });
    return [
      new AppHeader({icon: 'flask', items: ['Design molecules with your team', 'Score them as you go'],
        learnMore: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README.md'}),
      divH([h1('Continue campaign'), ...table.chrome()]),
      table,
    ];
  }}).view;
}
