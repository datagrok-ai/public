/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import $ from 'cash-dom';
import {CampaignGrouping, CampaignGroupingType, CampaignJsonName, CampaignTableColumns, DefaultCampaignTableInfoGetters, HitDesignCampaignIdKey, i18n} from '../consts';
import {HitDesignCampaign, HitDesignTemplate} from '../types';
import {addBreadCrumbsToRibbons, checkEditPermissions,
  checkViewPermissions, getGroupedCampaigns, getSavedCampaignsGrouping, getSavedCampaignsSorting, getSavedCampaignTableColumns, modifyUrl, popRibbonPannels,
  processGroupingTable,
  SavedCampaignsTableSorting,
  setSavedCampaignsGrouping,
  setSavedCampaignsSorting,
  setSavedCampaignTableColumns,
  sortCampaigns} from '../utils';
import {newHitDesignCampaignAccordeon} from '../accordeons/new-hit-design-campaign-accordeon';
import {newHitDesignTemplateAccordeon} from '../accordeons/new-hit-design-template-accordeon';
import {HitBaseView} from '../base-view';
import {defaultPermissions, PermissionsDialog} from '../dialogs/permissions-dialog';

export class HitDesignInfoView
  <T extends HitDesignTemplate = HitDesignTemplate, K extends HitDesignApp = HitDesignApp>
  extends HitBaseView<T, K> {
  currentGroupping: string = 'None';
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

  protected campaignsTableRoot: HTMLElement | null = null;

  protected async refreshCampaignsTable() {
    if (this.campaignsTableRoot == null)
      return;
    ui.setUpdateIndicator(this.campaignsTableRoot, true);
    try {
      const t = await this.getCampaignsTable();
      ui.setUpdateIndicator(this.campaignsTableRoot, false);
      ui.empty(this.campaignsTableRoot);
      this.campaignsTableRoot.appendChild(t);
    } catch (e) {
      grok.shell.error('Failed to update campaigns table');
      console.error(e);
    } finally {
      ui.setUpdateIndicator(this.campaignsTableRoot, false);
    }
  };
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
      this.campaignsTableRoot = ui.div([campaignsTable], {style: {position: 'relative'}});

      // grouping via different campaign properties
      const groupIcon = ui.iconFA('layer-group', async () => {
        const menu = DG.Menu.popup();
        Object.values(CampaignGrouping).forEach((i) => {
          menu.item(i, async () => {
            setSavedCampaignsGrouping(i as CampaignGroupingType);
            await this.refreshCampaignsTable();
          });
        });
        const campaignFieldsGroup = menu.group('Campaign Fields');
        const campaignNamesMap = await _package.loadCampaigns(this.app.appName, this.deletedCampaigns);
        const customFields = new Set<string>();
        Object.values(campaignNamesMap).forEach((c) => {
          if (c.campaignFields)
            Object.keys(c.campaignFields).forEach((field) => customFields.add(field));
        });
        Array.from(customFields).forEach((field) => {
          campaignFieldsGroup.item(field, async () => {
            setSavedCampaignsGrouping(`campaignFields.${field}`);
            await this.refreshCampaignsTable();
          });
        });

        menu.show({element: sortingHeader, x: 120, y: sortingHeader.offsetTop + 30});
      });
      groupIcon.style.marginBottom = '9px';
      groupIcon.style.marginLeft = '8px';
      groupIcon.style.fontSize = '15px';
      groupIcon.style.color = 'var(--blue-1)';
      ui.tooltip.bind(groupIcon, () => `Group Campaigns. Current: ${this.currentGroupping}`);

      const editColumnsIcon = ui.iconFA('eye', async () => {
        const menu = DG.Menu.popup();
        const campaignNamesMap = await _package.loadCampaigns(this.app.appName, this.deletedCampaigns);
        const getters: {[key in CampaignTableColumns]: (a: HitDesignCampaign) => string} = {...DefaultCampaignTableInfoGetters};
        // remove Name as it is always present
        // @ts-ignore
        delete getters['Name'];
        Object.values(campaignNamesMap).forEach((c) => {
          if (c.campaignFields) {
            Object.keys(c.campaignFields).forEach((field) => {
              getters[`campaignFields.${field}`] = (a) => a.campaignFields?.[field] ?? '';
            });
          }
        });
        const savedColumns = getSavedCampaignTableColumns().reduce((acc: {[code in CampaignTableColumns]: boolean}, col) => {acc[col] = true; return acc;}, {} as any);
        // first add the main columns
        Object.keys(getters).filter((g) => !g.startsWith('campaignFields.')).forEach((col) => {
          menu.item(col, () => {
            savedColumns[col as CampaignTableColumns] = !savedColumns[col as CampaignTableColumns];
            setSavedCampaignTableColumns(Object.keys(savedColumns).filter((c) => savedColumns[c as CampaignTableColumns]) as CampaignTableColumns[]);
            this.refreshCampaignsTable();
          }, null, {check: savedColumns[col as CampaignTableColumns] ?? false});
        });
        // then add a new group
        const campaignsPropsGroup = menu.group('Campaign Fields');

        Object.keys(getters).filter((g) => g.startsWith('campaignFields.')).forEach((col) => {
          campaignsPropsGroup.item(col.replace('campaignFields.', ''), () => {
            savedColumns[col as CampaignTableColumns] = !savedColumns[col as CampaignTableColumns];
            setSavedCampaignTableColumns(Object.keys(savedColumns).filter((c) => savedColumns[c as CampaignTableColumns]) as CampaignTableColumns[]);
            this.refreshCampaignsTable();
          }, null, {check: savedColumns[col as CampaignTableColumns] ?? false});
        });

        // then add the custom fields
        menu.show({element: sortingHeader, x: 100, y: sortingHeader.offsetTop + 30});
      });
      ui.tooltip.bind(editColumnsIcon, 'Edit visible columns in campaigns table');
      editColumnsIcon.style.marginBottom = '9px';
      editColumnsIcon.style.marginLeft = '8px';
      editColumnsIcon.style.fontSize = '15px';
      editColumnsIcon.style.color = 'var(--blue-1)';

      const refreshIcon = ui.iconFA('sync', async () => {
        await this.refreshCampaignsTable();
      });
      refreshIcon.style.marginBottom = '9px';
      refreshIcon.style.marginLeft = '8px';
      refreshIcon.style.color = 'var(--blue-1)';
      ui.tooltip.bind(refreshIcon, () => 'Refresh campaigns table');
      const sortingHeader = ui.divH([continueCampaignsHeader, editColumnsIcon, groupIcon, refreshIcon], {style: {alignItems: 'center'}});
      $(this.root).empty();
      this.root.appendChild(ui.div([
        ui.divV([appHeader, sortingHeader], {style: {marginLeft: '10px'}}),
        this.campaignsTableRoot,
        createNewCampaignHeader,
        contentDiv,
      ], {classes: 'hit-triage-info-view-container'}));
      await this.startNewCampaign(campaignAccordionDiv, templatesDiv, presetTemplate);
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
    deleteTempleteButton.style.color = 'var(--blue-1)';
    createNewtemplateButton.style.color = '#2083d5';
    templatesInput.addOptions(createNewtemplateButton);
    templatesInput.addOptions(cloneTemplateButton);
    templatesInput.addOptions(deleteTempleteButton);
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
    const campaignNamesMap = await _package.loadCampaigns(this.app.appName, this.deletedCampaigns);
    const grouppingMode = getSavedCampaignsGrouping();
    const campaignSorting = getSavedCampaignsSorting();
    const allCampaigns = Object.values(campaignNamesMap);
    // sort all campaigns
    if (campaignSorting)
      sortCampaigns(allCampaigns, campaignSorting);
    const grouppedCampaigns = getGroupedCampaigns<HitDesignCampaign>(allCampaigns, grouppingMode);
    const shownColumns = getSavedCampaignTableColumns();

    this.currentGroupping = grouppingMode;
    this.app.existingStatuses = Array.from(new Set(allCampaigns.map((c) => c.status).filter((s) => !!s)));
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
            await _package.saveCampaignJson(this.app.appName, info);
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
      deleteIcon.style.color = 'var(--blue-1)';
      shareIcon.style.color = 'var(--blue-1)';
      return [shareIcon, deleteIcon];
    };

    const table = ui.table(Object.values(grouppedCampaigns).flat(), (info) =>
      ([ui.link(info.friendlyName ?? info.name, () => this.setCampaign(info.name), 'Continue Campaign', ''),
        ...(shownColumns.map((col) => col.startsWith('campaignFields.') ?
          info.campaignFields?.[col.replace('campaignFields.', '')] ?? '' :
          DefaultCampaignTableInfoGetters[col as unknown as keyof typeof DefaultCampaignTableInfoGetters](info) ?? '')),
        ...(deleteAndShareCampaignIcons(info)),
      ]),
    ['Name', ...(shownColumns.map((a) => a.replace('campaignFields.', ''))), '', '']);
    table.style.color = 'var(--grey-5)';
    table.style.marginLeft = '24px';
    // add sorting support
    table.querySelectorAll('tr.header > td').forEach((el, i) => {
      if (!el.textContent)
        return;
      const field = i === 0 ? 'Name' : shownColumns[i - 1];
      if (campaignSorting && field === campaignSorting.columnName) {
        // add icon
        const sortIcon = campaignSorting.ascending ? '↑' : '↓';
        el.appendChild(ui.span([sortIcon], {style: {fontSize: '12px', marginLeft: '2px'}}));
      }
      (el as HTMLElement).addEventListener('dblclick', () => {
        let nextSorting: SavedCampaignsTableSorting | null = {columnName: field, ascending: true};
        if (campaignSorting && field === campaignSorting.columnName)
          nextSorting = campaignSorting.ascending ? {columnName: field, ascending: false} : null;
        setSavedCampaignsSorting(nextSorting);
        this.refreshCampaignsTable();
      });
    });


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
      this.app.clearCampaign();// make sure the previous campaign is cleared
      this.app.dataFrame = camp.df;
      await this.app.setTemplate(template);
      this.app.campaignProps = camp.campaignProps;
      const campaignName = !!camp.name?.trim() ? camp.name?.trim() : undefined;
      await this.app.saveCampaign(false, true, {friendlyName: campaignName});
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
