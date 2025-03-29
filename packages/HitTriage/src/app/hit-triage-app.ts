import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageCampaign, IFunctionArgs,
  HitTriageTemplate, HitTriageTemplateIngest, IngestType, HitTriageCampaignStatus,
  HitDesignTemplate} from './types';
import {InfoView} from './hit-triage-views/info-view';
import {SubmitView} from './hit-triage-views/submit-view';
import {CampaignIdKey, CampaignJsonName, CampaignTableName,
  HTQueryPrefix, HTScriptPrefix, HitSelectionColName, i18n} from './consts';
import {addBreadCrumbsToRibbons, checkRibbonsHaveSubmit, modifyUrl, toFormatedDateString} from './utils';
import {_package} from '../package';
import '../../css/hit-triage.css';
import {chemFunctionsDialog} from './dialogs/functions-dialog';
import {HitAppBase} from './hit-app-base';
import {HitBaseView} from './base-view';
import {saveCampaignDialog} from './dialogs/save-campaign-dialog';
import {calculateColumns} from './utils/calculate-single-cell';
import {defaultPermissions, PermissionsDialog} from './dialogs/permissions-dialog';

export class HitTriageApp extends HitAppBase<HitTriageTemplate> {
  multiView: DG.MultiView;

  private _infoView: InfoView;
  get infoView(): InfoView {return this._infoView;}
  private _pickView?: DG.TableView;
  private _submitView?: SubmitView;

  private _filterViewName = 'Hit triage | Pick';
  private _campaignFilters?: {[key: string]: any}[];
  private _campaignId?: string;
  private _dfName?: string;
  private _molColName?: string;
  public _fileInputType?: IngestType;

  private _campaign?: HitTriageCampaign;
  protected _filterDescriptions: string[] = [];
  public campaignProps: {[key: string]: any} = {};
  private currentPickViewId?: string;
  private _pickViewPromise?: Promise<void> | null = null;

  constructor(c: DG.FuncCall) {
    super(c, 'Hit Triage');
    this._infoView = new InfoView(this);
    this.multiView = new DG.MultiView({viewFactories: {[this._infoView.name]: () => this._infoView}});
    this.multiView.tabs.onTabChanged.subscribe((_) => {
      if (this.multiView.currentView instanceof HitBaseView)
        (this.multiView.currentView as HitBaseView<HitTriageTemplate, HitTriageApp>).onActivated();
    });
    this.multiView.parentCall = c;
    //grok.shell.addView(this.multiView);

    grok.events.onCurrentViewChanged.subscribe(() => {
      try {
        if (grok.shell.v?.name === this.currentPickViewId) {
          grok.shell.windows.showHelp = false;
          this.setBaseUrl();
          modifyUrl(CampaignIdKey, this._campaignId ?? this._campaign?.name ?? '');
        }
      } catch (e) {
        console.error(e);
      }
    });
  }

  public async setTemplate(template: HitTriageTemplate, presetFilters?: {[key: string]: any}[],
    campaignId?: string, ingestProps?: HitTriageTemplateIngest) {
    this._pickView?.dataFrame && grok.shell.closeTable(this._pickView?.dataFrame);
    this._pickView = undefined;
    if (!campaignId) {
      campaignId = await this.getNewCampaignName('Hit Triage/campaigns', template.key);
      modifyUrl(CampaignIdKey, campaignId);
    } else if (ingestProps) {
      this._fileInputType = ingestProps.type;
      if (ingestProps.type === 'File')
        this.dataFrame = await grok.dapi.files.readCsv(ingestProps.query);
      else
        this.dataFrame = await grok.functions.call(ingestProps.query);
    }
    if (!this.dataFrame) {
      console.error('Dataframe is empty.');
      return;
    }
    if (this._campaign?.columnSemTypes) {
      Object.entries(this._campaign.columnSemTypes).forEach(([colName, semtype]) => {
        const col = this.dataFrame!.columns.byName(colName);
        if (col && col.semType !== semtype)
          col.semType = semtype;
      });
    };
    if (this._campaign?.columnTypes) {
      Object.entries(this._campaign.columnTypes).forEach(([colName, type]) => {
        const col = this.dataFrame!.columns.byName(colName);
        if (col && col.type !== type)
          try {col.convertTo(type);} catch (e) {console.error(e);}
      });
    }
    await this.dataFrame.meta.detectSemanticTypes();
    this._molColName = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? undefined;
    this._dfName = this.dataFrame.name ?? ingestProps?.query;
    this.dataFrame.name = this._dfName ?? this.dataFrame.name ?? 'Molecules';
    this._campaignId = campaignId;
    this.template = this.campaign?.template ?
      {...this.campaign.template, dataSourceType: template.dataSourceType,
        queryFunctionName: template.queryFunctionName} : template;
    this._campaignFilters = presetFilters;

    if (!presetFilters) {
      const funcs: {[_: string]: IFunctionArgs} = {};
      const scripts: {[_: string]: IFunctionArgs} = {};
      const queries: {[_: string]: IFunctionArgs} = {};
      template.compute.functions.forEach((func) => {
        const fName = `${func.package}:${func.name}`;
        funcs[fName] = func.args;
      });
      template.compute.scripts?.forEach((script) => {
        const fName = `${HTScriptPrefix}:${script.name}:${script.id}`;
        scripts[fName] = script.args;
      });
      template.compute.queries?.forEach((query) => {
        const fName = `${HTQueryPrefix}:${query.name}:${query.id}`;
        queries[fName] = query.args;
      });

      await calculateColumns({
        descriptors: template.compute.descriptors.enabled ? template.compute.descriptors.args : [],
        externals: funcs, scripts, queries,
      }, this.dataFrame!, this._molColName!);
    };

    if (this.campaign)
      await this.setCanEdit(this.campaign);
    else
      this.hasEditPermission = true; // if the campaign is new, obviously the user can edit it

    const curView = grok.shell.v;
    const pickV = grok.shell.addView(this.pickView);
    this.currentPickViewId = pickV.name;
    this._submitView = new SubmitView(this);
    this.setBaseUrl();
    modifyUrl(CampaignIdKey, this._campaignId ?? this._campaign?.name ?? '');

    const newView = pickV;

    setTimeout(() => {
      this._pickViewPromise && this._pickViewPromise.then(() => {
        const {sub} = addBreadCrumbsToRibbons(newView, 'Hit Triage', 'Pick', () => {
          grok.shell.v = curView;
          newView.close();
          sub.unsubscribe();
        });
      });
    }, 300);
    // this.multiView.addView(this._submitView.name, () => this._submitView!, false);
    grok.shell.windows.showHelp = false;
  }

  get filterSettings(): {[key: string]: any}[] | undefined {return this._campaignFilters;}

  get campaignId(): string | undefined {return this._campaignId;}

  get pickView(): DG.TableView {return this._pickView = this.getFilterView();}

  get molColName() {return this._molColName ??= this.dataFrame?.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name;}

  get fileInputType(): IngestType | undefined {return this._fileInputType;}

  get campaign(): HitTriageCampaign | undefined {return this._campaign;}

  set campaign(campaign: HitTriageCampaign | undefined) {this._campaign = campaign;}

  /**
   * A view that lets you filter the molecules using either molecules, or
   * their properties derived at the enrichment step.
   * @return {DG.TableView}
   * */
  getFilterView(): DG.TableView {
    let resolvePromise: Function | null = null;
    this._pickViewPromise = new Promise<void>((res) => {resolvePromise = res;});
    const getComputeDialog = async () => {
      chemFunctionsDialog(this, async (resultMap) => {
        const oldDescriptors = this.template!.compute.descriptors.args;
        const oldFunctions = this.template!.compute.functions;
        const oldScripts = this.template!.compute.scripts ?? [];
        const oldQueries = this.template!.compute.queries ?? [];
        const newDescriptors = resultMap.descriptors;
        const newComputeObj = {
          descriptors: {
            enabled: !!resultMap?.descriptors?.length,
            args: resultMap?.descriptors ?? [],
          },
          functions: Object.entries(resultMap?.externals ?? {}).map(([funcName, args]) => {
            const splitFunc = funcName.split(':');
            return ({
              name: splitFunc[1],
              package: splitFunc[0],
              args: args,
            });
          }),
          scripts: Object.entries(resultMap?.scripts ?? {})
            .filter(([name, _]) => name.startsWith(HTScriptPrefix) && name.split(':').length === 3)
            .map(([scriptId, args]) => {
              const scriptNameParts = scriptId.split(':');
              return ({
                name: scriptNameParts[1] ?? '',
                id: scriptNameParts[2] ?? '',
                args: args,
              });
            }),
          queries: Object.entries(resultMap?.queries ?? {})
            .filter(([name, _]) => name.startsWith(HTQueryPrefix) && name.split(':').length === 3)
            .map(([scriptId, args]) => {
              const scriptNameParts = scriptId.split(':');
              return ({
                name: scriptNameParts[1] ?? '',
                id: scriptNameParts[2] ?? '',
                args: args,
              });
            }),
        };
        this.template!.compute = newComputeObj;
        const uncalculatedDescriptors = newDescriptors.filter((d) => !oldDescriptors.includes(d));
        const uncalculatedFunctions = newComputeObj.functions.filter((func) => {
          if (!oldFunctions.some((f) => f.name === func.name && f.package === func.package))
            return true;
          const oldFunc = oldFunctions.find((f) => f.name === func.name && f.package === func.package)!;
          return !Object.entries(func.args).every(([key, value]) => oldFunc.args[key] === value);
        });
        const uncalculatedScripts = newComputeObj.scripts.filter((func) => {
          if (!oldScripts.some((f) => f.id === func.id ))
            return true;
          const oldScript = oldScripts.find((f) => f.id === func.id)!;
          return !Object.entries(func.args).every(([key, value]) => oldScript.args[key] === value);
        });
        const uncalculatedQueries = newComputeObj.queries.filter((query) => {
          if (!oldQueries.some((f) => f.id === query.id))
            return true;
          const oldQuery = oldQueries.find((f) => f.id === query.id)!;
          return !Object.entries(query.args).every(([key, value]) => oldQuery.args[key] === value);
        });

        const externalFuncs: {[_: string]: IFunctionArgs} = {};
        uncalculatedFunctions.forEach((func) => {
          const fName = `${func.package}:${func.name}`;
          externalFuncs[fName] = func.args;
        });
        const externalScripts: {[_: string]: IFunctionArgs} = {};
        uncalculatedScripts.forEach((script) => {
          const fName = `${HTScriptPrefix}:${script.name}:${script.id}`;
          externalScripts[fName] = script.args;
        });
        const externalQueries: {[_: string]: IFunctionArgs} = {};
        uncalculatedQueries.forEach((query) => {
          const fName = `${HTQueryPrefix}:${query.name}:${query.id}`;
          externalQueries[fName] = query.args;
        });
        await calculateColumns(
          {descriptors: uncalculatedDescriptors, externals: externalFuncs,
            scripts: externalScripts, queries: externalQueries},
          this.dataFrame!, this._molColName!);
        this.saveCampaign('In Progress', false);
      }, () => null, this.template!, true);
    };
    if (!this.dataFrame!.col(HitSelectionColName))
      this.dataFrame!.columns.addNewBool(HitSelectionColName).init(false);

    const view = DG.TableView.create(this.dataFrame!, false);
    const ribbons = view.getRibbonPanels();
    const calculateRibbon = ui.iconFA('wrench', getComputeDialog, 'Calculate additional properties');
    const submitButton = ui.bigButton('Submit', () => {
      const dialogContent = this._submitView?.render();
      if (dialogContent) {
        const dlg = ui.dialog('Submit');
        dlg.add(dialogContent);
        dlg.addButton('Save', ()=>{this.saveCampaign(); dlg.close();});
        dlg.addButton('Submit', ()=>{this._submitView?.submit(); dlg.close();});
        dlg.show();
      }
    });
    const permissionsButton = ui.iconFA('share', async () => {
      await (new PermissionsDialog(this.campaign?.permissions)).show((res) => {
        this.campaign!.permissions = res;
        this.saveCampaign(undefined, true);
      });
    }, 'Edit permissions');
    submitButton.classList.add('hit-design-submit-button');
    const hasSubmit = checkRibbonsHaveSubmit(ribbons);
    ribbons.push([
      calculateRibbon,
      ...(this.hasEditPermission ? [permissionsButton] : []),
      ...(hasSubmit ? [] : [submitButton])]);
    view.setRibbonPanels(ribbons);
    view.name = this._filterViewName;
    setTimeout(async () => {
      view._onAdded();
      await new Promise((r) => setTimeout(r, 1000));

      const f = view.filters(this._campaignFilters ? {filters: this._campaignFilters} : undefined);
      const layoutViewState = this._campaign?.layout ?? this.template?.layoutViewState;
      if (layoutViewState) {
        try {
          const layout = DG.ViewLayout.fromViewState(layoutViewState);
          view.loadLayout(layout);
        } catch (e) {
          console.error(e);
        }
      }
      view.dataFrame.onFilterChanged
        .subscribe((_) => {
          this._filterDescriptions = Array.from(view.dataFrame.rows.filters);
          this._campaignFilters = f.getOptions().look.filters;
        });
      setTimeout(() => {
        this._filterDescriptions = Array.from(view.dataFrame.rows.filters);
        this._campaignFilters = f.getOptions().look.filters;
      }, 300);
      resolvePromise && resolvePromise();
    }, 300);
    view.parentCall = this.parentCall;
    return view;
  }

  getSummary(): {[_: string]: any} {
    const campaignProps = this.campaign?.campaignFields ?? this.campaignProps;
    return {
      'Template': this.template?.name ?? 'Molecules',
      'File path': ui.divH([ui.link(this._dfName ?? '-', () => this.download(this.dataFrame!, 'molecules'),
        i18n.download)],
      {style: {alignItems: 'center'}}),
      ...campaignProps,
      'Number of molecules': this.dataFrame!.rowCount.toString(),
      'Enrichment methods': [this.template!.compute.descriptors.enabled ? 'descriptors' : '',
        ...this.template!.compute.functions.map((func) => func.name)].filter((f) => f && f.trim() !== '').join(', '),
      'Filters': this._filterDescriptions.join(', '),
      'Result Molecules': ui.divH([ui.link(this.dataFrame?.filter.trueCount.toString() ?? '0', () => {
        this.download(this.dataFrame!, 'hits', true);
      }, i18n.download)], {style: {alignItems: 'center'}}),
    };
  }

  async saveCampaign(status?: HitTriageCampaignStatus, notify = true): Promise<any> {
    const campaignId = this.campaignId!;
    const filters = this.filterSettings!;
    const templateName = this.template!.name;
    const enrichedDf = this.dataFrame!;
    const campaignPrefix = `System:AppData/HitTriage/Hit Triage/campaigns/${campaignId}/`;
    const campaignName = campaignId ?? await saveCampaignDialog(campaignId);
    const columnSemTypes: {[_: string]: string} = {};
    enrichedDf.columns.toList().forEach((col) => columnSemTypes[col.name] = col.semType);

    const colTypeMap: {[_: string]: string} = {};
    enrichedDf.columns.toList().forEach((col) => colTypeMap[col.name] = col.type);

    // if its first time save author as current user, else keep the same
    const authorUserId = this.campaign?.authorUserId ?? grok.shell.user.id;
    const permissions = this.campaign?.permissions ?? defaultPermissions;


    const campaign: HitTriageCampaign = {
      name: campaignName,
      templateName,
      filters: filters ?? {},
      ingest: {
        type: 'File',
        query: `${campaignPrefix}${CampaignTableName}`,
        molColName: this.molColName!,
      },
      status: status ?? this.campaign?.status ?? 'In Progress',
      createDate: this.campaign?.createDate ?? toFormatedDateString(new Date()),
      campaignFields: this.campaign?.campaignFields ?? this.campaignProps,
      columnSemTypes,
      rowCount: enrichedDf.rowCount,
      filteredRowCount: enrichedDf.filter.trueCount,
      template: this.template as HitDesignTemplate | undefined,
      columnTypes: colTypeMap,
      authorUserId,
      permissions,
    };

    this.campaign = campaign;
    if (!this.hasEditPermission) {
      grok.shell.error('You do not have permission to modify this campaign');
      return campaign;
    }
    const csvDf = DG.DataFrame.fromColumns(
      enrichedDf.columns.toList().filter((col) => !col.name.startsWith('~')),
    ).toCsv();
    await _package.files.writeAsText(`Hit Triage/campaigns/${campaignId}/${CampaignTableName}`, csvDf);
    const newLayout = this._pickView!.saveLayout();
    if (!newLayout)
      grok.shell.warning('Layout cound not be saved');
    else
      campaign.layout = newLayout.viewState;
    await _package.saveCampaignJson(this.appName, campaign);
    notify && grok.shell.info('Campaign saved successfully.');
    return campaign;
  }
}
