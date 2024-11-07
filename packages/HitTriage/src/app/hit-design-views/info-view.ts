/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import $ from 'cash-dom';
import {CampaignGroupingType, CampaignJsonName, HitDesignCampaignIdKey, i18n} from '../consts';
import {HitDesignCampaign, HitDesignTemplate} from '../types';
import {addBreadCrumbsToRibbons, checkEditPermissions,
  checkViewPermissions, getGroupedCampaigns, getSavedCampaignsGrouping,
  loadCampaigns, modifyUrl, popRibbonPannels,
  processGroupingTable,
  setSavedCampaignsGrouping} from '../utils';
import {newHitDesignCampaignAccordeon} from '../accordeons/new-hit-design-campaign-accordeon';
import {newHitDesignTemplateAccordeon} from '../accordeons/new-hit-design-template-accordeon';
import {HitBaseView} from '../base-view';
import {defaultPermissions, PermissionsDialog} from '../dialogs/permissions-dialog';

export class HitDesignInfoView
  <T extends HitDesignTemplate = HitDesignTemplate, K extends HitDesignApp = HitDesignApp>
  extends HitBaseView<T, K> {
  currentSorting: string = 'None';
  constructor(app: K) {
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

  getAppHeader() {
    return u2.appHeader({
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
  }

  async init(presetTemplate?: T) {
    ui.setUpdateIndicator(this.root, true);
    try {
      const continueCampaignsHeader = ui.h1(i18n.continueCampaigns);

      const createNewCampaignHeader = ui.h1(i18n.createNewCampaignHeader, {style: {marginLeft: '10px'}});
      const appHeader = this.getAppHeader();

      const campaignAccordionDiv = ui.div();
      const templatesDiv = ui.divH([]);
      const contentDiv = ui.div([templatesDiv, campaignAccordionDiv], 'ui-form');

      const campaignsTable = await this.getCampaignsTable();
      const tableRoot = ui.div([campaignsTable], {style: {position: 'relative'}});

      const sortIcon = ui.iconFA('sort', () => {
        const menu = DG.Menu.popup();
        Object.values(CampaignGroupingType).forEach((i) => {
          menu.item(i, async () => {
            setSavedCampaignsGrouping(i as CampaignGroupingType);
            ui.setUpdateIndicator(tableRoot, true);
            try {
              const t = await this.getCampaignsTable();
              ui.setUpdateIndicator(tableRoot, false);
              ui.empty(tableRoot);
              tableRoot.appendChild(t);
            } catch (e) {
              grok.shell.error('Failed to update campaigns table');
              console.error(e);
            } finally {
              ui.setUpdateIndicator(tableRoot, false);
            }
          });
          menu.show({element: sortingHeader, x: 100, y: sortingHeader.offsetTop + 30});
        });
      });
      sortIcon.style.marginBottom = '12px';
      sortIcon.style.marginLeft = '5px';
      sortIcon.style.fontSize = '15px';
      ui.tooltip.bind(sortIcon, () => `Group Campaigns. Current: ${this.currentSorting}`);
      const sortingHeader = ui.divH([continueCampaignsHeader, sortIcon], {style: {alignItems: 'center'}});
      $(this.root).empty();
      this.root.appendChild(ui.div([
        ui.divV([appHeader, sortingHeader], {style: {marginLeft: '10px'}}),
        tableRoot,
        createNewCampaignHeader,
        contentDiv,
      ]));
      await this.startNewCampaign(campaignAccordionDiv, templatesDiv, presetTemplate);
      this.app.resetBaseUrl();
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  protected async startNewCampaign(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, presetTemplate?: T,
  ) {
    const templates = (await _package.files.list(`${this.app.appName}/templates`))
      .filter((file) => file.name.endsWith('.json'))
      .map((file) => file.name.slice(0, -5));
    // if the template is just created and saved, it may not be in the list of templates
    if (presetTemplate && !templates.includes(presetTemplate.name))
      templates.push(presetTemplate.name);

    let selectedTemplate: T | null = null;
    let templateChangeFlag = true;
    const onTemmplateChange = async () => {
      if (!templateChangeFlag)
        return;
      const templateName = templatesInput.value;
      const template: T = presetTemplate && presetTemplate.name === templateName ? presetTemplate :
        JSON.parse(await _package.files.readAsText(`${this.app.appName}/templates/${templateName}.json`));
      selectedTemplate = template;
      const newCampaignAccordeon = await this.getNewCampaignAccordeon(template);
      $(containerDiv).empty();
      containerDiv.appendChild(newCampaignAccordeon);
    };

    const templatesInput = ui.input.choice('Template', {value: presetTemplate?.name ?? templates[0], items: templates,
      onValueChanged: async () => {
        await onTemmplateChange();
      }});
    templatesInput.root.style.width = '100%';
    const createNewtemplateButton = ui.icons.add(() => {
      this.createNewTemplate();
    }, i18n.createNewTemplate);

    const cloneTemplateButton = ui.icons.copy(() => {
      if (selectedTemplate)
        this.createNewTemplate(selectedTemplate);
    }, 'Clone template');

    const deleteTempleteButton = ui.icons.delete(() => {
      if (!templatesInput.value || templatesInput.items.length < 2)
        return;
      ui.dialog('Delete template')
        .add(ui.divText(`Are you sure you want to delete template ${templatesInput.value}?`))
        .onOK(async () => {
          try {
            const prevValue = templatesInput.value;
            await _package.files.delete(`${this.app.appName}/templates/${prevValue}.json`);
            templateChangeFlag = false;
            templatesInput.items = templatesInput.items.filter((item) => item !== prevValue);
            templatesInput.value = templatesInput.items[0];
            templateChangeFlag = true;
            await onTemmplateChange();
            grok.shell.info('Template ' + prevValue + ' was deleted successfully');
          } catch (e) {
            templateChangeFlag = true;
            grok.shell.error('Failed to delete template ' + templatesInput.value);
            console.error(e);
          }
        })
        .show();
    }, 'Delete template');
    createNewtemplateButton.style.color = '#2083d5';
    templatesInput.addOptions(deleteTempleteButton);
    templatesInput.addOptions(cloneTemplateButton);
    templatesInput.addOptions(createNewtemplateButton);
    await onTemmplateChange();
    $(templateInputDiv).empty();
    templateInputDiv.appendChild(templatesInput.root);
  }


  protected async checkCampaign(campId?: string) {
    const url = location.search;
    const urlParams = new URLSearchParams(url);
    if (!urlParams.has(HitDesignCampaignIdKey) && !campId)
      return;
    const campaignId = campId ?? urlParams.get(HitDesignCampaignIdKey);
    // check if such campaign exists
    if (!await _package.files.exists(`${this.app.appName}/campaigns/${campaignId}/${CampaignJsonName}`))
      return;
    const campaign: HitDesignCampaign =
      JSON.parse(await _package.files.readAsText(`${this.app.appName}/campaigns/${campaignId}/${CampaignJsonName}`));
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
    const template: T = campaign.template ?? JSON.parse(
      await _package.files.readAsText(`${this.app.appName}/templates/${campaign.templateName}.json`),
    );
    // modify the template with path to the campaign's precalculated table
    await this.app.setTemplate(template, campaignId!);

    if (campId)
      modifyUrl(HitDesignCampaignIdKey, campId);

    campaign.template = template;
    return campaign;
  }

  private async getCampaignsTable() {
    const campaignNamesMap = await loadCampaigns(this.app.appName, this.deletedCampaigns);
    const grouppingMode = getSavedCampaignsGrouping();
    const grouppedCampaigns = getGroupedCampaigns<HitDesignCampaign>(Object.values(campaignNamesMap), grouppingMode);
    this.currentSorting = grouppingMode;
    this.app.existingStatuses = Array.from(new Set(Object.values(campaignNamesMap).map((c) => c.status).filter((s) => !!s)));
    const deleteAndShareCampaignIcons = (info: HitDesignCampaign) => {
      const deleteIcon = ui.icons.delete(async () => {
        ui.dialog('Delete campaign')
          .add(ui.divText(`Are you sure you want to delete campaign ${info.name}?`))
          .onOK(async () => {
            await this.deleteCampaign(this.app.appName, info.name);
            this.deletedCampaigns.push(info.name);
            await this.init();
          })
          .show();
      }, 'Delete campaign');
      const shareIcon = ui.iconFA('share', async () => {
        await (new PermissionsDialog(info.permissions)).show(async (res) => {
          try {
            info.permissions = res;
            info.authorUserId ??= grok.shell.user.id;
            await _package.files.writeAsText(
              `${this.app.appName}/campaigns/${info.name}/${CampaignJsonName}`, JSON.stringify(info));
            grok.shell.info('Permissions updated for campaign ' + info.name);
          } catch (e) {
            grok.shell.error('Failed to update permissions for campaign ' + info.name);
            console.error(e);
          }
        });
      }, 'Manage campaign permissions');
      deleteIcon.style.display = 'none';
      shareIcon.style.display = 'none';
      const authorId = info.authorUserId ?? DG.User.current().id;
      const perms = info.permissions ?? defaultPermissions;
      checkEditPermissions(authorId, perms).then((canEdit) => {
        if (canEdit) {
          deleteIcon.style.display = 'inline-block';
          shareIcon.style.display = 'inline-block';
        }
      });
      return [shareIcon, deleteIcon];
    };
    const table = ui.table(Object.values(campaignNamesMap), (info) =>
      ([ui.link(info.name, () => this.setCampaign(info.name), '', ''),
        info.createDate,
        info.rowCount,
        info.status,
        ...(deleteAndShareCampaignIcons(info)),
      ]),
    ['Name', 'Created', 'Molecules', 'Status', '']);
    table.style.color = 'var(--grey-5)';
    table.style.marginLeft = '24px';
    processGroupingTable(table, grouppedCampaigns);
    return table;
  }

  public async setCampaign(campaignName: string) {
    const campaign = await this.checkCampaign(campaignName);
    this.app.campaign = campaign;
  }
  async getNewCampaignAccordeon(template: T) {
    const {root, promise, cancelPromise} = newHitDesignCampaignAccordeon(template);
    promise.then(async (camp) => {
      this.app.dataFrame = camp.df;
      await this.app.setTemplate(template);
      this.app.campaignProps = camp.campaignProps;
      await this.app.saveCampaign(false);
      if (template.layoutViewState && this.app.campaign)
        this.app.campaign.layout = template.layoutViewState;
    });

    cancelPromise.then(() => {
      this.init();
    });
    return root;
  }

  private async createNewTemplate(preset?: T) {
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
      const {sub} = addBreadCrumbsToRibbons(newView, this.app.appName, i18n.createNewTemplate, async () => {
        grok.shell.v = curView;
        newView.close();
      });
      //containerDiv.appendChild(newTemplateAccordeon.root);
      newTemplateAccordeon.template.then(async (t) => {
        await this.init(t as any);
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
