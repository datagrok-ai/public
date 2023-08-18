import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {HitTriageBaseView} from './base-view';
import {_package} from '../../package';
import {HitTriageTemplate} from '../types';
import {CampaignIdKey, CampaignJsonName, i18n} from '../consts';
import {HitTriageCampaign} from '../types';
import '../../../css/hit-triage.css';
import {addBreadCrumbsToRibbons, hideComponents, modifyUrl} from '../utils';
import {newCampaignAccordeon} from '../accordeons/new-campaign-accordeon';
import $ from 'cash-dom';
import {createTemplateAccordeon} from '../accordeons/new-template-accordeon';

export class InfoView extends HitTriageBaseView {
  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Hit Triage';
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
    this.checkCampaign().then((c) => {this.app.campaign = c; this.init();});
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
  }

  async init(presetTemplate?: HitTriageTemplate): Promise<void> {
    $(this.root).empty();
    this.root.style.flexDirection = 'column';
    const wikiLink = ui.link('Read more', _package.webRoot + 'README.md');
    const textLink = ui.inlineText([wikiLink, '.']);
    const continueCampaignsHeader = ui.h1(i18n.continueCampaigns);
    const createNewCampaignHeader = ui.h1(i18n.createNewCampaignHeader, {style: {marginLeft: '10px'}});
    const appDescription = ui.divV([
      ui.h1('Process, analyse and filter molecules for your needs using Hit Triage:'),
      ui.list([
        '-  Configure your own workflow using the template editor.',
        '-  Calculate different molecular properties.',
        '-  Filter molecules using different criteria.',
        '-  Submit processed dataframe to the function of your choice.',
        '-  Save campaigns and continue any time from where you left off.',
      ]), textLink, continueCampaignsHeader,
    ], {style: {marginLeft: '10px'}});

    const campaignAccordionDiv = ui.div();
    const templatesDiv = ui.divH([], {classes: 'hit-triage-templates-input-div ui-form'});

    const campaignsTable = await this.getCampaignsTable();
    this.root.appendChild(ui.divV([
      appDescription,
      campaignsTable,
      createNewCampaignHeader,
      templatesDiv,
      campaignAccordionDiv,
    ]));
    this.startNewCampaign(campaignAccordionDiv, templatesDiv,
      [campaignsTable.style, continueCampaignsHeader.style, createNewCampaignHeader.style, appDescription.style],
      presetTemplate);
  }

  private async startNewCampaign(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, toRemove: CSSStyleDeclaration[],
    presetTemplate?: HitTriageTemplate) {
    //hideComponents(toRemove);
    const templates = (await _package.files.list('Hit Triage/templates')).map((file) => file.name.slice(0, -5));
    // if the template is just created and saved, it may not be in the list of templates
    if (presetTemplate && !templates.includes(presetTemplate.name))
      templates.push(presetTemplate.name);

    const onTemmplateChange = async () => {
      const templateName = templatesInput.value;
      const template: HitTriageTemplate = presetTemplate && presetTemplate.name === templateName ? presetTemplate :
        JSON.parse(await _package.files.readAsText('Hit Triage/templates/' + templateName + '.json'));
      const newCampaignAccordeon = await this.getNewCampaignAccordeon(template);
      $(containerDiv).empty();
      containerDiv.appendChild(newCampaignAccordeon);
    };

    const templatesInput = ui.choiceInput(i18n.selectTemplate, presetTemplate?.name ?? templates[0], templates,
      async () => {
        await onTemmplateChange();
      });
    const createNewtemplateButton = ui.icons.add(() => {
      this.createNewTemplate(containerDiv, templateInputDiv, toRemove);
    }, i18n.createNewTemplate);
    createNewtemplateButton.style.marginLeft = '15px';
    createNewtemplateButton.style.color = '#2083d5';
    templatesInput.root.style.width = '100%';
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
    if (!await _package.files.exists(`Hit Triage/campaigns/${campaignId}/${CampaignJsonName}`))
      return;
    const campaign: HitTriageCampaign =
      JSON.parse(await _package.files.readAsText(`Hit Triage/campaigns/${campaignId}/${CampaignJsonName}`));
    // Load the template and modify it
    const template: HitTriageTemplate = JSON.parse(
      await _package.files.readAsText(`Hit Triage/templates/${campaign.templateName}.json`),
    );
    // modify the template with path to the campaign's precalculated table
    this.app.setTemplate(template, campaign.filters, campaignId!, campaign.ingest);

    if (campId)
      modifyUrl(CampaignIdKey, campId);

    return campaign;
  }

  private async getCampaignsTable() {
    const campaignFolders = await _package.files.list('Hit Triage/campaigns');
    const campaignNamesMap: {[name: string]: HitTriageCampaign} = {};
    for (const folder of campaignFolders) {
      const campaignJson: HitTriageCampaign = JSON.parse(await _package.files
        .readAsText(`Hit Triage/campaigns/${folder.name}/${CampaignJsonName}`));
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

  private async getNewCampaignAccordeon(template: HitTriageTemplate) {
    const {root, promise, cancelPromise} = newCampaignAccordeon(template);
    promise.then((camp) => {
      this.app.dataFrame = camp.df;
      this.app._fileInputType = camp.type;
      this.app.setTemplate(template);
      this.app.campaignProps = camp.campaignProps;
    });

    cancelPromise.then(() => {
      this.init();
    });
    return root;
  }

  private async createNewTemplate(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, toRemove: CSSStyleDeclaration[]) {
    hideComponents(toRemove);
    const newTemplateAccordeon = await createTemplateAccordeon();
    $(containerDiv).empty();
    $(templateInputDiv).empty();
    const {breadcrumbs, sub} = addBreadCrumbsToRibbons(grok.shell.v, i18n.createNewTemplate, () => {
      this.init();
    });
    //this.root.prepend(breadcrumbs.root);
    containerDiv.appendChild(newTemplateAccordeon.root);
    newTemplateAccordeon.template.then((t) => {
      sub.unsubscribe();
      this.init(t);
      $(breadcrumbs.root).empty();
      $(breadcrumbs.root).remove();
    });
    newTemplateAccordeon.cancelPromise.then(() => this.init());
  }
};
