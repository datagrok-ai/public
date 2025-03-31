import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {_package} from '../../package';
import {HitTriageTemplate} from '../types';
import {CampaignIdKey, CampaignJsonName, i18n} from '../consts';
import {HitTriageCampaign} from '../types';
import '../../../css/hit-triage.css';
import {addBreadCrumbsToRibbons, checkEditPermissions, checkViewPermissions, modifyUrl, popRibbonPannels} from '../utils';
import {newCampaignAccordeon} from '../accordeons/new-campaign-accordeon';
import $ from 'cash-dom';
import {createTemplateAccordeon} from '../accordeons/new-template-accordeon';
import {HitBaseView} from '../base-view';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {defaultPermissions} from '../dialogs/permissions-dialog';

export class InfoView extends HitBaseView<HitTriageTemplate, HitTriageApp> {
  public readmePath = _package.webRoot + 'README_HT.md';
  private dataSourceFunctionsMap: {[key: string]: DG.Func | DG.DataQuery} = {};
  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Hit Triage';
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(this.readmePath);
    this.checkCampaign().then((c) => {this.app.campaign = c; this.init();});
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(this.readmePath);
  }

  async init(presetTemplate?: HitTriageTemplate): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      const continueCampaignsHeader = ui.h1(i18n.continueCampaigns);
      const createNewCampaignHeader = ui.h1(i18n.createNewCampaignHeader, {style: {marginLeft: '10px'}});
      const appHeader = u2.appHeader({
        iconPath: _package.webRoot + '/images/icons/hit-triage-icon.png',
        learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HT.md',
        description: '- Configure your own workflow using the template editor.\n'+
        '- Calculate different molecular properties.\n'+
        '- Filter molecules using different criteria.\n'+
        '- Submit processed dataframe to the function of your choice.\n'+
        '- Save campaigns and continue any time from where you left off.\n ',
      });

      const campaignAccordionDiv = ui.div();
      const templatesDiv = ui.div([], {classes: 'ui-form'});

      const campaignsTable = await this.getCampaignsTable();

      await this.startNewCampaign(campaignAccordionDiv, templatesDiv, presetTemplate);
      $(this.root).empty();
      this.root.appendChild(ui.div([
        ui.divV([appHeader, continueCampaignsHeader], {style: {marginLeft: '10px'}}),
        campaignsTable,
        createNewCampaignHeader,
        templatesDiv,
        campaignAccordionDiv,
      ]));
    } catch (e) {
      ui.setUpdateIndicator(this.root, false);
      throw e;
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
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
      const newCampaignAccordeon = await this.getNewCampaignAccordeon(template, containerDiv);
      $(containerDiv).empty();
      containerDiv.className = 'ui-form';
      containerDiv.appendChild(newCampaignAccordeon);
    };

    const templatesInput =
      ui.input.choice(i18n.selectTemplate, {value: presetTemplate?.name ?? templates[0], items: templates,
        onValueChanged: async () => {
          await onTemmplateChange();
        }});
    const createNewtemplateButton = ui.icons.add(async () => {
      ui.setUpdateIndicator(this.root, true);
      await this.createNewTemplate();
      ui.setUpdateIndicator(this.root, false);
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

    if (campaign) {
      // in case if the link was opened and user has no permissions to view the campaign
      if (campaign.authorUserId && campaign.permissions &&
          !await checkViewPermissions(campaign.authorUserId, campaign.permissions)
      ) {
        this.app.setBaseUrl();
        return;
      }
      this.app.campaign = campaign;
    }
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
    const campaignNamesMap = await _package.loadCampaigns('Hit Triage', this.deletedCampaigns);

    const deleteCampaignIcon = (info: HitTriageCampaign) => {
      const icon = ui.icons.delete(async () => {
        ui.dialog('Delete campaign')
          .add(ui.divText(`Are you sure you want to delete campaign ${info.name}?`))
          .onOK(async () => {
            await this.deleteCampaign('Hit Triage', info.name);
            this.deletedCampaigns.push(info.name);
            await this.init();
          })
          .show();
      }, 'Delete campaign');
      icon.style.display = 'none';
      const authorId = info.authorUserId ?? DG.User.current().id;
      const perms = info.permissions ?? defaultPermissions;
      checkEditPermissions(authorId, perms).then((canEdit) => {
        if (canEdit)
          icon.style.display = 'inline-block';
      });
      return icon;
    };

    const campaignsInfo = Object.values(campaignNamesMap);
    const table = ui.table(campaignsInfo, (info) =>
      ([ui.link(info.name, () => this.setCampaign(info.name), '', ''),
        info.createDate, info.rowCount, info.filteredRowCount, info.status,
        deleteCampaignIcon(info),
      ]),
    ['Name', 'Created', 'Total', 'Selected', 'Status', '']);
    table.style.color = 'var(--grey-5)';
    return table;
  }

  public async setCampaign(campaignName: string) {
    const campaign = await this.checkCampaign(campaignName);
    this.app.campaign = campaign;
  }

  private async getNewCampaignAccordeon(template: HitTriageTemplate, campaignDetailsDiv: HTMLElement) {
    ui.setUpdateIndicator(campaignDetailsDiv, true);
    const {root, promise, cancelPromise} = await newCampaignAccordeon(template, this.dataSourceFunctionsMap);
    ui.setUpdateIndicator(campaignDetailsDiv, false);
    promise.then(async (camp) => {
      this.app.dataFrame = camp.df;
      this.app._fileInputType = camp.type;
      await this.app.setTemplate(template);
      this.app.campaignProps = camp.campaignProps;
      await this.app.saveCampaign(undefined, false);
      if (template.layoutViewState && this.app.campaign)
        this.app.campaign.layout = template.layoutViewState;
    });

    cancelPromise.then(() => {
      this.init();
    });
    return root;
  }

  private async createNewTemplate() {
    const newTemplateAccordeon = await createTemplateAccordeon(this.app, this.dataSourceFunctionsMap);

    const newView = new DG.ViewBase();
    const curView = grok.shell.v;
    newView.name = 'New Template';
    newView.root.appendChild(newTemplateAccordeon.root);
    newView.parentCall = this.app.parentCall;
    grok.shell.addView(newView);
    newView.path = new URL(this.app.baseUrl).pathname + '/new-template';
    newView.parentView = curView;
    const {sub} = addBreadCrumbsToRibbons(grok.shell.v, 'Hit Triage', i18n.createNewTemplate, () => {
      grok.shell.v = curView;
      newView.close();
    });

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
