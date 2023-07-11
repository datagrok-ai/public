import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {HitTriageBaseView} from './base-view';
import {_package} from '../../package';
import {ITemplate} from '../types';
import {createTemplateDialog} from '../dialogs/templateDialog';
import {CampaignIdKey, CampaignJsonName, CampaignTableName} from '../consts';
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

    const appDescription = ui.info([
      ui.list([
        '-  ipsum lorem dolor sit amet, consectetur adipiscing elit, ut labore et dolore magna aliqua.',
        '-  suscipit urna quis, placerat dui. Aliquam erat volutpat.',
        '-  Mauris sit amet orci eleifend nunc viverra varius',
        '-  Donec auctor, nunc vel tempor aliquam, nisl nunc ultricies nunc, quis aliquam nunc nunc nec nunc.',
        '-  Maecenas vehicula nunc vel augue lobortis, ut ultrices nibh suscipit.',
        '-  roin nec lectus tempus, ultrices eros auctor',
      ]),
    ], 'Process, analyse and filter molecules for your needs using Hit Triage:',
    );

    const templates = (await _package.files.list('templates')).map((file) => file.name.slice(0, -5));
    const templatesInput = ui.choiceInput('Select template', templates[0], templates, null);
    const useTemplateButton = ui.button('Use template', () => {
      _package.files.readAsText('templates/' + templatesInput.value + '.json').then((template) => {
        const templateJson: ITemplate = JSON.parse(template);
        this.app.setTemplate(templateJson);
      });
    }, '');
    this.root.appendChild(ui.divV([
      appDescription,
      ui.info([textLink]),
      ui.divH([
        templatesInput.root,
        useTemplateButton,
        ui.button('Create new template', () => createTemplateDialog().then((t) => this.app.setTemplate(t)), ''),
      ]),
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
    template.ingest.query = `System:AppData/HitTriage/campaigns/${campaignId}/${CampaignTableName}`;
    this.app.setTemplate(template, campaign.filters, campaignId!);
  }
};
