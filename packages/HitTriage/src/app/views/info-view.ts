import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {HitTriageBaseView} from './base-view';
import {_package} from '../../package';
import {ITemplate} from '../types';
import {createTemplateDialog} from '../dialogs/templateDialog';
import {CampaignIdKey, CampaignJsonName} from '../consts';
import {ICampaign} from '../types';

export class InfoView extends HitTriageBaseView {
  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Hit Triage';
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
    this.init();
    this.checkCampaign();
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
  }

  async init() {
    const wikiLink = ui.link('wiki', _package.webRoot + 'README.md');
    const textLink = ui.inlineText(['For more details, see our ', wikiLink, '.']);

    const appDescription = ui.divV([
      ui.h1('Process, analyse and filter molecules for your needs using Hit Triage:'),
      ui.list([
        '-  Configure your own workflow using the template editor.',
        '-  Calculate differnet molecular properties.',
        '-  Filter molecules using different criteria.',
        '-  Submit processed dataframe to the function of your choice.',
        '-  Save campaigns and continue any time from where you left off.',
      ]),
    ],
    );

    const templates = (await _package.files.list('templates')).map((file) => file.name.slice(0, -5));
    const templatesInput = ui.choiceInput('Select template', templates[0], templates, null);
    const useTemplateButton = ui.button('Use template', () => {
      _package.files.readAsText('templates/' + templatesInput.value + '.json').then((template) => {
        const templateJson: ITemplate = JSON.parse(template);
        this.app.setTemplate(templateJson);
      });
    }, '');

    const campaignFolders = await _package.files.list('campaigns');
    const campaignNamesMap: {[name: string]: string} = {};
    for (const folder of campaignFolders) {
      const campaignJson = JSON.parse(await _package.files
        .readAsText(`campaigns/${folder.name}/${CampaignJsonName}`));
      campaignNamesMap[campaignJson.name] = folder.name;
    }
    const campaignsInput = ui.choiceInput('Select campaign', Object.keys(campaignNamesMap)[0],
      Object.keys(campaignNamesMap), null);
    const useCampaignButton = ui.button('Continue campaign', () => {
      const campaignId = campaignNamesMap[campaignsInput.value!];
      const url = location.href.split('?')[0] + '?' + CampaignIdKey + '=' + campaignId;
      location.href = url;
    });
    this.root.appendChild(ui.divV([
      appDescription,
      ui.info([textLink]),
      ui.divH([
        templatesInput.root,
        useTemplateButton,
        ui.button('Create new template', () => createTemplateDialog().then((t) => this.app.setTemplate(t)), ''),
      ]),
      ui.divH([
        campaignsInput.root,
        useCampaignButton,
      ]),
      await this.getCampaignsTable(),
    ]));
  }

  private async checkCampaign() {
    const url = location.search;
    const urlParams = new URLSearchParams(url);
    if (!urlParams.has(CampaignIdKey))
      return;
    const campaignId = urlParams.get(CampaignIdKey);
    // check if such campaign exists
    if (!await _package.files.exists(`campaigns/${campaignId}/${CampaignJsonName}`))
      return;
    const campaign: ICampaign =
      JSON.parse(await _package.files.readAsText(`campaigns/${campaignId}/${CampaignJsonName}`));
    // Load the template and modify it
    const template: ITemplate = JSON.parse(await _package.files.readAsText(`templates/${campaign.templateName}.json`));
    // modify the template with path to the campaign's precalculated table
    this.app.setTemplate(template, campaign.filters, campaignId!, campaign.ingest);
  }

  private async getCampaignsTable() {
    const campaignFolders = await _package.files.list('campaigns');
    const campaignNamesMap: {[name: string]: ICampaign} = {};
    for (const folder of campaignFolders) {
      const campaignJson: ICampaign = JSON.parse(await _package.files
        .readAsText(`campaigns/${folder.name}/${CampaignJsonName}`));
      campaignNamesMap[campaignJson.name] = campaignJson;
    }

    const campaignsInfo = Object.values(campaignNamesMap).map((campaign) =>
      ({name: campaign.name, template: campaign.templateName}));
    const table = ui.table(campaignsInfo, (info) => ([info.name, info.template]), ['Campaign', 'Template']);
    return table;
  }
};
