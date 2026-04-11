/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from '../hit-triage-app';
import {_package} from '../../package';
import {HitTriageTemplate} from '../types';
import {CampaignIdKey, CampaignJsonName, CampaignGrouping, CampaignGroupingType,
  HTCampaignTableColumns, HTDefaultCampaignTableInfoGetters, HitTriageDataSourceTag, i18n} from '../consts';
import {HitTriageCampaign} from '../types';
import '../../../css/hit-triage.css';
import {addBreadCrumbsToRibbons, checkEditPermissions, checkViewPermissions,
  getGroupedCampaigns, getHTSavedCampaignsGrouping, getHTSavedCampaignsSorting,
  getHTSavedCampaignTableColumns, loadTemplate, modifyUrl, popRibbonPannels,
  processGroupingTable, setHTSavedCampaignsGrouping, setHTSavedCampaignsSorting,
  setHTSavedCampaignTableColumns, sortCampaigns,
  HTSavedCampaignsTableSorting} from '../utils';
import {newCampaignAccordeon} from '../accordeons/new-campaign-accordeon';
import $ from 'cash-dom';
import {createTemplateAccordeon} from '../accordeons/new-template-accordeon';
import {HitBaseView} from '../base-view';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {defaultPermissions, PermissionsDialog} from '../dialogs/permissions-dialog';

export class InfoView extends HitBaseView<HitTriageTemplate, HitTriageApp> {
  public readmePath = _package.webRoot + 'README_HT.md';
  private dataSourceFunctionsMap: {[key: string]: DG.Func | DG.DataQuery} = {};
  currentGroupping: string = 'None';
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

  protected getAppHeader() {
    return u2.appHeader({
      iconPath: _package.webRoot + '/images/icons/hit-triage-icon.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HT.md',
      description: '- Configure your own workflow using the template editor.\n'+
      '- Calculate different molecular properties.\n'+
      '- Filter molecules using different criteria.\n'+
      '- Submit processed dataframe to the function of your choice.\n'+
      '- Save campaigns and continue any time from where you left off.\n ',
    });
  }

  protected getTemplatesFolder(): string {
    return `${this.app.appName}/templates`;
  }

  protected getDataSourceTag(): string {
    return HitTriageDataSourceTag;
  }

  protected getDataSourceFunctionsMap(): {[key: string]: DG.Func | DG.DataQuery} {
    return this.dataSourceFunctionsMap;
  }

  async init(presetTemplate?: HitTriageTemplate): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      const continueCampaignsHeader = ui.h1(i18n.continueCampaigns);
      const createNewCampaignHeader = ui.h1(i18n.createNewCampaignHeader, {style: {marginLeft: '10px'}});
      const appHeader = this.getAppHeader();

      const campaignAccordionDiv = ui.div();
      const templatesDiv = ui.div([], {classes: 'ui-form'});

      const campaignsTable = await this.getCampaignsTable();
      this.campaignsTableRoot = ui.div([campaignsTable], {style: {position: 'relative'}});

      // grouping via different campaign properties
      const groupIcon = ui.iconFA('layer-group', async () => {
        const menu = DG.Menu.popup();
        Object.values(CampaignGrouping).forEach((i) => {
          menu.item(i, async () => {
            setHTSavedCampaignsGrouping(i as CampaignGroupingType);
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
            setHTSavedCampaignsGrouping(`campaignFields.${field}`);
            await this.refreshCampaignsTable();
          });
        });

        menu.show({element: sortingHeader, x: 120, y: sortingHeader.offsetTop + 30});
      }, 'Group by');
      groupIcon.style.marginBottom = '9px';
      groupIcon.style.marginLeft = '8px';
      groupIcon.style.fontSize = '15px';
      groupIcon.style.color = 'var(--blue-1)';
      ui.tooltip.bind(groupIcon, () => `Group Campaigns. Current: ${this.currentGroupping}`);

      const editColumnsIcon = ui.iconFA('eye', async () => {
        const menu = DG.Menu.popup();
        const campaignNamesMap = await _package.loadCampaigns(this.app.appName, this.deletedCampaigns);
        const getters: {[key in HTCampaignTableColumns]: (a: HitTriageCampaign) => string} = {...HTDefaultCampaignTableInfoGetters};
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
        const savedColumns = getHTSavedCampaignTableColumns().reduce((acc: {[code in HTCampaignTableColumns]: boolean}, col) => {acc[col] = true; return acc;}, {} as any);
        // first add the main columns
        Object.keys(getters).filter((g) => !g.startsWith('campaignFields.')).forEach((col) => {
          menu.item(col, () => {
            savedColumns[col as HTCampaignTableColumns] = !savedColumns[col as HTCampaignTableColumns];
            setHTSavedCampaignTableColumns(Object.keys(savedColumns).filter((c) => savedColumns[c as HTCampaignTableColumns]) as HTCampaignTableColumns[]);
            this.refreshCampaignsTable();
          }, null, {check: savedColumns[col as HTCampaignTableColumns] ?? false});
        });
        // then add a new group for custom campaign fields
        const campaignsPropsGroup = menu.group('Campaign Fields');

        Object.keys(getters).filter((g) => g.startsWith('campaignFields.')).forEach((col) => {
          campaignsPropsGroup.item(col.replace('campaignFields.', ''), () => {
            savedColumns[col as HTCampaignTableColumns] = !savedColumns[col as HTCampaignTableColumns];
            setHTSavedCampaignTableColumns(Object.keys(savedColumns).filter((c) => savedColumns[c as HTCampaignTableColumns]) as HTCampaignTableColumns[]);
            this.refreshCampaignsTable();
          }, null, {check: savedColumns[col as HTCampaignTableColumns] ?? false});
        });

        menu.show({element: sortingHeader, x: 100, y: sortingHeader.offsetTop + 30});
      }, 'View');
      ui.tooltip.bind(editColumnsIcon, 'Edit visible columns in campaigns table');
      editColumnsIcon.style.marginBottom = '9px';
      editColumnsIcon.style.marginLeft = '8px';
      editColumnsIcon.style.fontSize = '15px';
      editColumnsIcon.style.color = 'var(--blue-1)';

      const refreshIcon = ui.iconFA('sync', async () => {
        await this.refreshCampaignsTable();
      }, 'Refresh');
      refreshIcon.style.marginBottom = '9px';
      refreshIcon.style.marginLeft = '8px';
      refreshIcon.style.color = 'var(--blue-1)';
      ui.tooltip.bind(refreshIcon, () => 'Refresh campaigns table');
      const sortingHeader = ui.divH([continueCampaignsHeader, editColumnsIcon, groupIcon, refreshIcon], {style: {alignItems: 'center'}});

      await this.startNewCampaign(campaignAccordionDiv, templatesDiv, presetTemplate);
      $(this.root).empty();
      this.root.appendChild(ui.div([
        ui.divV([appHeader, sortingHeader], {style: {marginLeft: '10px'}}),
        this.campaignsTableRoot,
        createNewCampaignHeader,
        templatesDiv,
        campaignAccordionDiv,
      ], {classes: 'hit-triage-info-view-container'}));
    } catch (e) {
      ui.setUpdateIndicator(this.root, false);
      throw e;
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  protected async startNewCampaign(
    containerDiv: HTMLElement, templateInputDiv: HTMLElement, presetTemplate?: HitTriageTemplate,
  ) {
    const templates = (await _package.files.list(this.getTemplatesFolder()))
      .filter((file) => file.name.endsWith('.json'))
      .map((file) => file.name.slice(0, -5));
    // if the template is just created and saved, it may not be in the list of templates
    if (presetTemplate && !templates.includes(presetTemplate.name))
      templates.push(presetTemplate.name);

    let selectedTemplate: HitTriageTemplate | null = null;
    let templateChangeFlag = true;
    const onTemmplateChange = async () => {
      if (!templateChangeFlag)
        return;
      const templateName = templatesInput.value;
      const template: HitTriageTemplate = presetTemplate && presetTemplate.name === templateName ? presetTemplate :
        await loadTemplate<HitTriageTemplate>(this.getTemplatesFolder() + '/' + templateName + '.json');
      selectedTemplate = template;
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
    templatesInput.root.style.width = '100%';
    const createNewtemplateButton = ui.icons.add(async () => {
      ui.setUpdateIndicator(this.root, true);
      await this.createNewTemplate();
      ui.setUpdateIndicator(this.root, false);
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
            await _package.files.delete(`${this.getTemplatesFolder()}/${prevValue}.json`);
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
    if (!urlParams.has(CampaignIdKey) && !campId)
      return;
    const campaignId = campId ?? urlParams.get(CampaignIdKey);
    // check if such campaign exists
    if (!await _package.files.exists(`${this.app.appName}/campaigns/${campaignId}/${CampaignJsonName}`))
      return;
    const campaign: HitTriageCampaign =
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
    const template: HitTriageTemplate = await loadTemplate<HitTriageTemplate>(
      `${this.getTemplatesFolder()}/${campaign.templateName}.json`,
    );
    // modify the template with path to the campaign's precalculated table
    this.app.setTemplate(template, campaign.filters, campaignId!, campaign.ingest);

    if (campId)
      modifyUrl(CampaignIdKey, campId);

    return campaign;
  }

  private async getCampaignsTable() {
    const campaignNamesMap = await _package.loadCampaigns(this.app.appName, this.deletedCampaigns) as {[name: string]: HitTriageCampaign};
    const grouppingMode = getHTSavedCampaignsGrouping();
    const campaignSorting = getHTSavedCampaignsSorting();
    const allCampaigns = Object.values(campaignNamesMap);
    // sort all campaigns
    if (campaignSorting)
      sortCampaigns(allCampaigns, campaignSorting, HTDefaultCampaignTableInfoGetters as any);
    const grouppedCampaigns = getGroupedCampaigns<HitTriageCampaign>(allCampaigns, grouppingMode);
    const shownColumns = getHTSavedCampaignTableColumns();

    this.currentGroupping = grouppingMode;

    const deleteAndShareCampaignIcons = (info: HitTriageCampaign) => {
      const deleteIcon = ui.icons.delete(async () => {
        ui.dialog('Delete campaign')
          .add(ui.divText(`Are you sure you want to delete campaign ${info.friendlyName ?? info.name}?`))
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
            grok.shell.info('Permissions updated for campaign ' + (info.friendlyName ?? info.name));
          } catch (e) {
            grok.shell.error('Failed to update permissions for campaign ' + (info.friendlyName ?? info.name));
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

    const campaignsInfo = Object.values(grouppedCampaigns).flat();
    const table = ui.table(campaignsInfo, (info) =>
      ([ui.link(info.friendlyName ?? info.name, () => this.setCampaign(info.name), 'Continue Campaign', ''),
        ...(shownColumns.map((col) => col.startsWith('campaignFields.') ?
          info.campaignFields?.[col.replace('campaignFields.', '')] ?? '' :
          HTDefaultCampaignTableInfoGetters[col as unknown as keyof typeof HTDefaultCampaignTableInfoGetters]?.(info) ?? '')),
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
        let nextSorting: HTSavedCampaignsTableSorting | null = {columnName: field, ascending: true};
        if (campaignSorting && field === campaignSorting.columnName)
          nextSorting = campaignSorting.ascending ? {columnName: field, ascending: false} : null;
        setHTSavedCampaignsSorting(nextSorting);
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

  protected async getNewCampaignAccordeon(template: HitTriageTemplate, campaignDetailsDiv: HTMLElement) {
    ui.setUpdateIndicator(campaignDetailsDiv, true);
    const {root, promise, cancelPromise} = await newCampaignAccordeon(template, this.getDataSourceFunctionsMap(), this.getDataSourceTag());
    ui.setUpdateIndicator(campaignDetailsDiv, false);
    promise.then(async (camp) => {
      this.app.dataFrame = camp.df;
      this.app._fileInputType = camp.type;
      await this.app.setTemplate(template, undefined, undefined, undefined, camp.friendlyName);
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

  protected async createNewTemplate(preset?: HitTriageTemplate) {
    const newTemplateAccordeon = await createTemplateAccordeon(this.app, this.getDataSourceFunctionsMap(), preset);

    const newView = new DG.ViewBase();
    const curView = grok.shell.v;
    newView.name = 'New Template';
    newView.root.appendChild(newTemplateAccordeon.root);
    newView.parentCall = this.app.parentCall;
    grok.shell.addView(newView);
    newView.path = new URL(this.app.baseUrl).pathname + '/new-template';
    newView.parentView = curView;
    const {sub} = addBreadCrumbsToRibbons(grok.shell.v, this.app.appName, i18n.createNewTemplate, () => {
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
