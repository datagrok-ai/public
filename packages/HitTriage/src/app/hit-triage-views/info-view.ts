import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {_package} from '../../package';
import {HitTriageTemplate} from '../types';
import {CampaignIdKey, CampaignJsonName, i18n} from '../consts';
import {HitTriageCampaign} from '../types';
import '../../../css/hit-triage.css';
import {addBreadCrumbsToRibbons, modifyUrl, popRibbonPannels} from '../utils';
import {newCampaignAccordeon} from '../accordeons/new-campaign-accordeon';
import $ from 'cash-dom';
import {createTemplateAccordeon} from '../accordeons/new-template-accordeon';
import {HitBaseView} from '../base-view';

export class InfoView extends HitBaseView<HitTriageTemplate, HitTriageApp> {
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
    const continueCampaignsHeader = ui.h1(i18n.continueCampaigns);
    const createNewCampaignHeader = ui.h1(i18n.createNewCampaignHeader, {style: {marginLeft: '10px'}});
    const description = '- Configure your own workflow using the template editor.\n'+
    '- Calculate different molecular properties.\n'+
    '- Filter molecules using different criteria.\n'+
    '- Submit processed dataframe to the function of your choice.\n'+
    '- Save campaigns and continue any time from where you left off.\n ';
    const appDescription = ui.div([
      ui.h1('Process, analyse and filter molecules for your needs using Hit Triage:'),
      ui.div(ui.markdown(description), {style: {color: 'var(--grey-5)'}, classes: 'mb-small'}),
      ui.link('Read more', 'https://github.com/datagrok-ai/public/tree/master/packages/HitTriage'),
    ]);

    const campaignAccordionDiv = ui.div();
    const templatesDiv = ui.div([], {classes: 'ui-form'});

    const campaignsTable = await this.getCampaignsTable();
    $(this.root).empty();
    this.root.appendChild(ui.div([
      appDescription,
      continueCampaignsHeader,
      campaignsTable,
      createNewCampaignHeader,
      templatesDiv,
      campaignAccordionDiv,
    ]));
    this.startNewCampaign(campaignAccordionDiv, templatesDiv, presetTemplate).then(() => this.app.resetBaseUrl());
  }

  private async startNewCampaign(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, presetTemplate?: HitTriageTemplate,
  ) {
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
      containerDiv.className = 'ui-form';
      containerDiv.appendChild(newCampaignAccordeon);
    };

    const templatesInput = ui.choiceInput(i18n.selectTemplate, presetTemplate?.name ?? templates[0], templates,
      async () => {
        await onTemmplateChange();
      });
    const createNewtemplateButton = ui.icons.add(() => {
      this.createNewTemplate();
    }, i18n.createNewTemplate);
    templatesInput.addOptions(createNewtemplateButton);

    await onTemmplateChange();
    $(templateInputDiv).empty();
    templateInputDiv.appendChild(templatesInput.root);
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
      ({name: campaign.name, createDate: campaign.createDate,
        rowCount: campaign.rowCount, filtered: campaign.filteredRowCount, status: campaign.status}));
    const table = ui.table(campaignsInfo, (info) =>
      ([ui.link(info.name, () => this.setCampaign(info.name), '', ''),
        info.createDate, info.rowCount, info.filtered, info.status]),
    ['Name', 'Created', 'Total', 'Selected', 'Status']);
    table.style.color = 'var(--grey-5)';
    return table;
  }

  private async setCampaign(campaignName: string) {
    const campaign = await this.checkCampaign(campaignName);
    this.app.campaign = campaign;
  }

  private async getNewCampaignAccordeon(template: HitTriageTemplate) {
    const {root, promise, cancelPromise} = newCampaignAccordeon(template);
    promise.then(async (camp) => {
      this.app.dataFrame = camp.df;
      this.app._fileInputType = camp.type;
      await this.app.setTemplate(template);
      this.app.campaignProps = camp.campaignProps;
      this.app.saveCampaign(undefined, false);
    });

    cancelPromise.then(() => {
      this.init();
    });
    return root;
  }

  private async createNewTemplate() {
    const newTemplateAccordeon = await createTemplateAccordeon();
    // hideComponents(toRemove);
    // $(containerDiv).empty();
    // $(templateInputDiv).empty();

    const newView = new DG.ViewBase();
    const curView = grok.shell.v;
    newView.name = 'New Template';
    newView.root.appendChild(newTemplateAccordeon.root);
    grok.shell.addView(newView);
    newView.path = new URL(this.app.baseUrl).pathname + '/new-template';
    newView.parentView = curView;
    const {sub} = addBreadCrumbsToRibbons(grok.shell.v, 'Hit Triage', i18n.createNewTemplate, () => {
      grok.shell.v = curView;
      newView.close();
    });
    //this.root.prepend(breadcrumbs.root);
    //containerDiv.appendChild(newTemplateAccordeon.root);

    newTemplateAccordeon.template.then(async (t) => {
      await this.init(t);
      newView.close();
      popRibbonPannels(newView);
      grok.shell.v = curView;
      sub.unsubscribe();
    });
    newTemplateAccordeon.cancelPromise.then(() => {
      sub.unsubscribe();
      newView.close();
      grok.shell.v = curView;
    });
  }
};
