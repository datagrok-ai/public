import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {HitTriageBaseView} from './base-view';
import {_package} from '../../package';
import {ITemplate} from '../types';
import {CampaignIdKey, CampaignJsonName} from '../consts';
import {ICampaign} from '../types';
import '../../../css/hit-triage.css';
import {hideComponents, modifyUrl} from '../utils';
import {newCampaignAccordeon} from '../accordeons/new-campaign-accordeon';
import $ from 'cash-dom';
import {createTemplateAccordeon} from '../accordeons/new-template-accordeon';

export class InfoView extends HitTriageBaseView {
  newItemHeader: HTMLElement = ui.h1('Create New Campaign');
  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Hit Triage';
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
    this.init();
    this.checkCampaign().then((c) => this.app.campaign = c);
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
  }

  async init() {
    $(this.root).empty();
    this.newItemHeader.style.display = 'none';
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

    const campaignAccordionDiv = ui.div();
    const templatesDiv = ui.divH([], {classes: 'hit-triage-templates-input-div'});
    const startNewCampaignButton = ui.button('Start new campaign', () => {
      this.startNewCampaign(campaignAccordionDiv, templatesDiv, [startNewCampaignButton.style, campaignsTable.style]);
    });
    startNewCampaignButton.classList.add('hit-triage-start-new-campaign-button');
    const campaignsTable = await this.getCampaignsTable();
    this.root.appendChild(ui.divV([
      appDescription,
      ui.info([textLink]),
      campaignsTable,
      startNewCampaignButton,
      this.newItemHeader,
      templatesDiv,
      campaignAccordionDiv,
    ]));
  }

  private async startNewCampaign(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, toRemove: CSSStyleDeclaration[],
    presetTemplate?: ITemplate) {
    hideComponents(toRemove);
    this.newItemHeader.style.display = 'block';
    this.newItemHeader.innerText = 'Create New Campaign';
    const templates = (await _package.files.list('templates')).map((file) => file.name.slice(0, -5));
    // if the template is just created and saved, it may not be in the list of templates
    if (presetTemplate && !templates.includes(presetTemplate.name))
      templates.push(presetTemplate.name);

    const onTemmplateChange = async () => {
      const templateName = templatesInput.value;
      const template: ITemplate = presetTemplate && presetTemplate.name === templateName ? presetTemplate :
        JSON.parse(await _package.files.readAsText('templates/' + templateName + '.json'));
      const newCampaignAccordeon = await this.getNewCampaignAccordeon(template);
      $(containerDiv).empty();
      containerDiv.appendChild(newCampaignAccordeon);
    };

    const templatesInput = ui.choiceInput('Select template', presetTemplate?.name ?? templates[0], templates,
      async () => {
        await onTemmplateChange();
      });
    const createNewtemplateButton = ui.link('Create new template', () => {
      this.createNewTemplate(containerDiv, templateInputDiv, toRemove);
    }, undefined, {style: {marginLeft: '15px'}});
    templateInputDiv.style.paddingBottom = '10px';
    templateInputDiv.style.borderBottom = '1px solid gray';
    await onTemmplateChange();
    $(templateInputDiv).empty();
    templateInputDiv.appendChild(templatesInput.root);
    templateInputDiv.appendChild(createNewtemplateButton);
  }
  private async checkCampaign(campId?: string) {
    const url = location.search;
    const urlParams = new URLSearchParams(url);
    if (!urlParams.has(CampaignIdKey) && !campId)
      return;
    const campaignId = campId ?? urlParams.get(CampaignIdKey);
    // check if such campaign exists
    if (!await _package.files.exists(`campaigns/${campaignId}/${CampaignJsonName}`))
      return;
    const campaign: ICampaign =
      JSON.parse(await _package.files.readAsText(`campaigns/${campaignId}/${CampaignJsonName}`));
    // Load the template and modify it
    const template: ITemplate = JSON.parse(await _package.files.readAsText(`templates/${campaign.templateName}.json`));
    // modify the template with path to the campaign's precalculated table
    this.app.setTemplate(template, campaign.filters, campaignId!, campaign.ingest);

    if (campId)
      modifyUrl(CampaignIdKey, campId);

    return campaign;
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
      ({name: campaign.name, createDate: campaign.createDate, status: campaign.status}));
    const table = ui.table(campaignsInfo, (info) => ([ui.link(info.name, () => this.setCampaign(info.name)),
      info.createDate, info.status]), ['Campaign', 'Create date', 'Status']);
    table.classList.add('hit-triage-table');
    return ui.div(table, {classes: 'hit-triage-table-container'});
  }

  private async setCampaign(campaignName: string) {
    const campaign = await this.checkCampaign(campaignName);
    this.app.campaign = campaign;
  }

  private async getNewCampaignAccordeon(template: ITemplate) {
    const {root, promise, cancelPromise} = newCampaignAccordeon(template);
    promise.then((camp) => {
      this.app.dataFrame = camp.df;
      this.app._fileInputType = camp.type;
      this.app.setTemplate(template);
      this.app.campaignProps = camp.campaignProps;
      this.newItemHeader.style.display = 'none';
    });

    cancelPromise.then(() => {
      this.init();
    });
    return root;
  }

  private async createNewTemplate(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, toRemove: CSSStyleDeclaration[]) {
    hideComponents(toRemove);
    this.newItemHeader.style.display = 'block';
    this.newItemHeader.innerText = 'Create New Template';
    const newTemplateAccordeon = await createTemplateAccordeon();
    $(containerDiv).empty();
    $(templateInputDiv).empty();
    containerDiv.appendChild(newTemplateAccordeon.root);
    newTemplateAccordeon.template.then((t) => this.startNewCampaign(containerDiv, templateInputDiv, toRemove, t));
    newTemplateAccordeon.cancelPromise.then(() => this.init());
  }
};
