import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import $ from 'cash-dom';
import {CampaignJsonName, HitDesignCampaignIdKey, i18n} from '../consts';
import {HitDesignCampaign, HitDesignTemplate} from '../types';
import {addBreadCrumbsToRibbons, hideComponents, modifyUrl} from '../utils';
import {newHitDesignCampaignAccordeon} from '../accordeons/new-hit-design-campaign-accordeon';
import {newHitDesignTemplateAccordeon} from '../accordeons/new-hit-design-template-accordeon';
import {HitBaseView} from '../base-view';

export class HitDesignInfoView extends HitBaseView<HitDesignTemplate, HitDesignApp> {
  newItemHeader: HTMLElement = ui.h1(i18n.startNewCampaign);
  constructor(app: HitDesignApp) {
    super(app);
    this.name = 'Hit Design';
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README.md'); // TODO: Separate readme for Hit Design
    this.checkCampaign().then((c) => {this.app.campaign = c; this.init();});
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
  }

  async init(presetTemplate?: HitDesignTemplate) {
    $(this.root).empty();
    this.newItemHeader.style.display = 'none';
    const wikiLink = ui.link('Read more', _package.webRoot + 'README.md'); // TODO: Separate readme for Hit Design
    const textLink = ui.inlineText([wikiLink, '.']);
    const continueCampaignsHeader = ui.h1(i18n.continueCampaigns);
    const createNewCampaignHeader = ui.h1(i18n.createNewCampaignHeader, {style: {marginLeft: '10px'}});
    const appDescription = ui.divV([
      ui.h1('Create, Process, analyse and filter molecules for your needs using Hit Design:'),
      ui.list([
        '-  Configure your own workflow using the template editor.',
        '-  Calculate differnet molecular properties.',
        '-  Add molecule rows using sketcher.',
        '-  Move molecules between stages using drag and drop in tile viewer.',
        '-  Submit processed dataframe to the function of your choice.',
        '-  Save campaigns and continue any time from where you left off.',
      ]),
    ]);
    const campaignAccordionDiv = ui.div();
    const templatesDiv = ui.divH([], {classes: 'hit-triage-templates-input-div ui-form'});

    const campaignsTable = await this.getCampaignsTable();
    this.root.appendChild(ui.divV([
      ui.divV([appDescription, textLink, continueCampaignsHeader], {style: {marginLeft: '10px'}}),
      campaignsTable,
      createNewCampaignHeader,
      this.newItemHeader,
      templatesDiv,
      campaignAccordionDiv,
    ]));
    this.startNewCampaign(campaignAccordionDiv, templatesDiv,
      [campaignsTable.style, continueCampaignsHeader.style, createNewCampaignHeader.style], presetTemplate);
  }

  private async startNewCampaign(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, toRemove: CSSStyleDeclaration[],
    presetTemplate?: HitDesignTemplate) {
    // hideComponents(toRemove);
    // this.newItemHeader.style.display = 'block';
    this.newItemHeader.innerText = i18n.startNewCampaign;
    const templates = (await _package.files.list('Hit Design/templates')).map((file) => file.name.slice(0, -5));
    // if the template is just created and saved, it may not be in the list of templates
    if (presetTemplate && !templates.includes(presetTemplate.name))
      templates.push(presetTemplate.name);

    const onTemmplateChange = async () => {
      const templateName = templatesInput.value;
      const template: HitDesignTemplate = presetTemplate && presetTemplate.name === templateName ? presetTemplate :
        JSON.parse(await _package.files.readAsText('Hit Design/templates/' + templateName + '.json'));
      const newCampaignAccordeon = await this.getNewCampaignAccordeon(template);
      $(containerDiv).empty();
      containerDiv.appendChild(newCampaignAccordeon);
    };

    const templatesInput = ui.choiceInput('Template', presetTemplate?.name ?? templates[0], templates,
      async () => {
        await onTemmplateChange();
      });
    templatesInput.root.style.width = '100%';
    const createNewtemplateButton = ui.icons.add(() => {
      this.createNewTemplate(containerDiv, templateInputDiv, toRemove);
    }, i18n.createNewTemplate);
    createNewtemplateButton.style.marginLeft = '15px';
    createNewtemplateButton.style.color = '#2083d5';

    await onTemmplateChange();
    $(templateInputDiv).empty();
    templateInputDiv.appendChild(templatesInput.root);
    templateInputDiv.appendChild(createNewtemplateButton);
  }


  private async checkCampaign(campId?: string) {
    const url = location.search;
    const urlParams = new URLSearchParams(url);
    if (!urlParams.has(HitDesignCampaignIdKey) && !campId)
      return;
    const campaignId = campId ?? urlParams.get(HitDesignCampaignIdKey);
    // check if such campaign exists
    if (!await _package.files.exists(`Hit Design/campaigns/${campaignId}/${CampaignJsonName}`))
      return;
    const campaign: HitDesignCampaign =
      JSON.parse(await _package.files.readAsText(`Hit Design/campaigns/${campaignId}/${CampaignJsonName}`));
    // Load the template and modify it
    const template: HitDesignTemplate = JSON.parse(
      await _package.files.readAsText(`Hit Design/templates/${campaign.templateName}.json`),
    );
    // modify the template with path to the campaign's precalculated table
    this.app.setTemplate(template, campaignId!);

    if (campId)
      modifyUrl(HitDesignCampaignIdKey, campId);

    return campaign;
  }

  private async getCampaignsTable() {
    const campaignFolders = await _package.files.list('Hit Design/campaigns');
    const campaignNamesMap: {[name: string]: HitDesignCampaign} = {};
    for (const folder of campaignFolders) {
      const campaignJson: HitDesignCampaign = JSON.parse(await _package.files
        .readAsText(`Hit Design/campaigns/${folder.name}/${CampaignJsonName}`));
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
  private async getNewCampaignAccordeon(template: HitDesignTemplate) {
    const {root, promise, cancelPromise} = newHitDesignCampaignAccordeon(template);
    promise.then((camp) => {
      this.app.dataFrame = camp.df;
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
    this.newItemHeader.innerText = i18n.createNewTemplate;
    const newTemplateAccordeon = await newHitDesignTemplateAccordeon();
    const {breadcrumbs, sub} = addBreadCrumbsToRibbons(grok.shell.v, 'Hit design', i18n.createNewTemplate, () => {
      this.init();
    });
    $(containerDiv).empty();
    $(templateInputDiv).empty();
    containerDiv.appendChild(newTemplateAccordeon.root);
    newTemplateAccordeon.template.then((t) => {
      $(breadcrumbs.root).empty();
      $(breadcrumbs.root).remove();
      sub.unsubscribe();
      this.init(t);
    });
  }
}
