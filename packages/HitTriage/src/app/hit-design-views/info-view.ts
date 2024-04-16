import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import $ from 'cash-dom';
import {CampaignJsonName, HitDesignCampaignIdKey, i18n} from '../consts';
import {HitDesignCampaign, HitDesignTemplate} from '../types';
import {addBreadCrumbsToRibbons, loadCampaigns, modifyUrl, popRibbonPannels} from '../utils';
import {newHitDesignCampaignAccordeon} from '../accordeons/new-hit-design-campaign-accordeon';
import {newHitDesignTemplateAccordeon} from '../accordeons/new-hit-design-template-accordeon';
import {HitBaseView} from '../base-view';

export class HitDesignInfoView extends HitBaseView<HitDesignTemplate, HitDesignApp> {
  constructor(app: HitDesignApp) {
    super(app);
    this.name = 'Hit Design';
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README_HD.md'); // TODO: Separate readme for Hit Design
    this.checkCampaign().then((c) => {this.app.campaign = c; this.init();});
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README_HD.md');
  }

  async init(presetTemplate?: HitDesignTemplate) {
    ui.setUpdateIndicator(this.root, true);
    try {
      const continueCampaignsHeader = ui.h1(i18n.continueCampaigns);
      const createNewCampaignHeader = ui.h1(i18n.createNewCampaignHeader, {style: {marginLeft: '10px'}});
      const appHeader = u2.appHeader({
        iconPath: _package.webRoot + '/images/icons/hit-design-icon.png',
        learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HD.md',
        description:
        '-  Configure your own workflow using the template editor\n' +
        '-  Sketch molecules in the molecular spreadsheet\n' +
        '-  Annotate and share ideas with the team\n' +
        '-  Calculate different molecular properties\n' +
        '-  Save campaigns and continue from where you left off\n' +
        '-  Submit final selection to the function of your choice',
      });

      const campaignAccordionDiv = ui.div();
      const templatesDiv = ui.divH([]);
      const contentDiv = ui.div([templatesDiv, campaignAccordionDiv], 'ui-form');

      const campaignsTable = await this.getCampaignsTable();
      $(this.root).empty();
      this.root.appendChild(ui.div([
        ui.divV([appHeader, continueCampaignsHeader], {style: {marginLeft: '10px'}}),
        campaignsTable,
        createNewCampaignHeader,
        contentDiv,
      ]));
      await this.startNewCampaign(campaignAccordionDiv, templatesDiv, presetTemplate);
      this.app.resetBaseUrl();
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  private async startNewCampaign(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, presetTemplate?: HitDesignTemplate,
  ) {
    const templates = (await _package.files.list('Hit Design/templates'))
      .filter((file) => file.name.endsWith('.json'))
      .map((file) => file.name.slice(0, -5));
    // if the template is just created and saved, it may not be in the list of templates
    if (presetTemplate && !templates.includes(presetTemplate.name))
      templates.push(presetTemplate.name);

    let selectedTemplate: HitDesignTemplate | null = null;
    const onTemmplateChange = async () => {
      const templateName = templatesInput.value;
      const template: HitDesignTemplate = presetTemplate && presetTemplate.name === templateName ? presetTemplate :
        JSON.parse(await _package.files.readAsText('Hit Design/templates/' + templateName + '.json'));
      selectedTemplate = template;
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
      this.createNewTemplate();
    }, i18n.createNewTemplate);

    const cloneTemplateButton = ui.icons.copy(() => {
      if (selectedTemplate)
        this.createNewTemplate(selectedTemplate);
    }, 'Clone template');
    createNewtemplateButton.style.color = '#2083d5';
    templatesInput.addOptions(cloneTemplateButton);
    templatesInput.addOptions(createNewtemplateButton);
    await onTemmplateChange();
    $(templateInputDiv).empty();
    templateInputDiv.appendChild(templatesInput.root);
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
    if (campaign)
      this.app.campaign = campaign;
    // Load the template and modify it
    const template: HitDesignTemplate = campaign.template ?? JSON.parse(
      await _package.files.readAsText(`Hit Design/templates/${campaign.templateName}.json`),
    );
    // modify the template with path to the campaign's precalculated table
    await this.app.setTemplate(template, campaignId!);

    if (campId)
      modifyUrl(HitDesignCampaignIdKey, campId);

    campaign.template = template;
    return campaign;
  }

  private async getCampaignsTable() {
    const campaignNamesMap = await loadCampaigns('Hit Design', this.deletedCampaigns);

    const table = ui.table(Object.values(campaignNamesMap), (info) =>
      ([ui.link(info.name, () => this.setCampaign(info.name), '', ''),
        info.createDate,
        info.rowCount,
        //info.filteredRowCount,
        info.status,
        // ui.icons.copy(async () => {
        //   const template = info.template ?? JSON.parse(
        //     await _package.files.readAsText(`Hit Design/templates/${info.templateName}.json`),
        //   );
        //   const df = await _package.files.readCsv(`Hit Design/campaigns/${info.name}/${CampaignTableName}`);
        //   this.app.dataFrame = df;
        //   await this.app.setTemplate(template);
        //   this.app.campaignProps = info.campaignFields;
        //   await this.app.saveCampaign(undefined, false);
        // }, 'Clone campaign'),
        ui.icons.delete(async () => {
          ui.dialog('Delete campaign')
            .add(ui.divText(`Are you sure you want to delete campaign ${info.name}?`))
            .onOK(async () => {
              await this.deleteCampaign('Hit Design', info.name);
              this.deletedCampaigns.push(info.name);
              await this.init();
            })
            .show();
        }, 'Delete campaign')]),
    ['Name', 'Created', 'Molecules', 'Status', '']);
    table.style.color = 'var(--grey-5)';
    table.style.marginLeft = '24px';
    return table;
  }

  public async setCampaign(campaignName: string) {
    const campaign = await this.checkCampaign(campaignName);
    this.app.campaign = campaign;
  }
  private async getNewCampaignAccordeon(template: HitDesignTemplate) {
    const {root, promise, cancelPromise} = newHitDesignCampaignAccordeon(template);
    promise.then(async (camp) => {
      this.app.dataFrame = camp.df;
      await this.app.setTemplate(template);
      this.app.campaignProps = camp.campaignProps;
      this.app.saveCampaign(undefined, false);
    });

    cancelPromise.then(() => {
      this.init();
    });
    return root;
  }

  private async createNewTemplate(preset?: HitDesignTemplate) {
    ui.setUpdateIndicator(this.root, true);
    try {
      const newTemplateAccordeon = await newHitDesignTemplateAccordeon(this.app, preset);
      // hideComponents(toRemove);
      // $(containerDiv).empty();
      // $(templateInputDiv).empty();

      const newView = new DG.ViewBase();
      const curView = grok.shell.v;
      newView.name = 'New Template';
      newView.root.appendChild(newTemplateAccordeon.root);
      newView.parentCall = this.app.parentCall;
      grok.shell.addView(newView);
      newView.path = new URL(this.app.baseUrl).pathname + '/new-template';
      const {sub} = addBreadCrumbsToRibbons(newView, 'Hit Design', i18n.createNewTemplate, async () => {
        grok.shell.v = curView;
        newView.close();
      });
      //containerDiv.appendChild(newTemplateAccordeon.root);
      newTemplateAccordeon.template.then(async (t) => {
        await this.init(t);
        newView.close();
        popRibbonPannels(newView);
        grok.shell.v = curView;
        sub.unsubscribe();
      });
      newTemplateAccordeon.cancelPromise.then(async () => {
        sub.unsubscribe();
        newView.close();
        grok.shell.v = curView;
      });
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }
}
