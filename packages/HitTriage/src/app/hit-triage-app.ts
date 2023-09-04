import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageCampaign, IComputeDialogResult, IFunctionArgs,
  HitTriageTemplate, HitTriageTemplateIngest, IngestType, HitTriageCampaignStatus} from './types';
import {InfoView} from './hit-triage-views/info-view';
import {SubmitView} from './hit-triage-views/submit-view';
import {CampaignIdKey, CampaignJsonName, CampaignTableName, HitSelectionColName, i18n} from './consts';
import {modifyUrl, toFormatedDateString} from './utils';
import {_package} from '../package';
import '../../css/hit-triage.css';
import {chemFunctionsDialog} from './dialogs/functions-dialog';
import {HitAppBase} from './hit-app-base';
import {HitBaseView} from './base-view';
import {saveCampaignDialog} from './dialogs/save-campaign-dialog';
export class HitTriageApp extends HitAppBase<HitTriageTemplate> {
  multiView: DG.MultiView;

  private _infoView: InfoView;
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
  constructor() {
    super();
    this._infoView = new InfoView(this);
    this.multiView = new DG.MultiView({viewFactories: {[this._infoView.name]: () => this._infoView}});
    this.multiView.tabs.onTabChanged.subscribe((_) => {
      if (this.multiView.currentView instanceof HitBaseView)
        (this.multiView.currentView as HitBaseView<HitTriageTemplate, HitTriageApp>).onActivated();
    });
    grok.shell.addView(this.multiView);
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
        if (col)
          col.semType = semtype;
      });
    }
    await this.dataFrame.meta.detectSemanticTypes();
    this._molColName = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? undefined;
    this._dfName = this.dataFrame.name ?? ingestProps?.query;
    this.dataFrame.name = this._dfName ?? this.dataFrame.name ?? 'Molecules';
    this._campaignId = campaignId;
    this.template = template;
    this._campaignFilters = presetFilters;

    if (!presetFilters) {
      const funcs: {[_: string]: IFunctionArgs} = {};
      template.compute.functions.forEach((func) => {
        const fName = `${func.package}:${func.name}`;
        funcs[fName] = func.args;
      } );
      await this.calculateColumns({
        descriptors: template.compute.descriptors.enabled ? template.compute.descriptors.args : [],
        externals: funcs,
      });
    };

    const pickV = grok.shell.addView(this.pickView);
    this.currentPickViewId = pickV.name;
    this._submitView ??= new SubmitView(this);
    this.setBaseUrl();
    modifyUrl(CampaignIdKey, this._campaignId ?? this._campaign?.name ?? '');
    // this.multiView.addView(this._submitView.name, () => this._submitView!, false);
    grok.shell.windows.showHelp = false;
  }

  get filterSettings(): {[key: string]: any}[] | undefined {return this._campaignFilters;}

  get campaignId(): string | undefined {return this._campaignId;}

  get pickView(): DG.TableView {return this._pickView ??= this.getFilterView();}

  get molColName() {return this._molColName ??= this.dataFrame?.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name;}

  get fileInputType(): IngestType | undefined {return this._fileInputType;}

  get campaign(): HitTriageCampaign | undefined {return this._campaign;}

  set campaign(campaign: HitTriageCampaign | undefined) {this._campaign = campaign;}

  public async calculateColumns(resultMap: IComputeDialogResult, view?: DG.TableView) {
    const previousColumns = this.dataFrame!.columns.names();
    if (resultMap.descriptors && resultMap.descriptors.length > 0)
      await grok.chem.descriptors(this.dataFrame!, this.molColName!, resultMap.descriptors);

    for (const funcName of Object.keys(resultMap.externals)) {
      const props = resultMap.externals[funcName];
      props['table'] = this.dataFrame!;
      props['molecules'] = this.molColName!;
      if (props)
        await grok.functions.call(funcName, props);
    };

    if (view) {
      const filterGroup = view.getFiltersGroup();
    this.dataFrame!.columns.names().filter((colName) => !previousColumns.includes(colName))
      .forEach((colName) => {filterGroup.add({type: this.getFilterType(colName), column: colName});});
    }
  }

  /**
   * A view that lets you filter the molecules using either molecules, or
   * their properties derived at the enrichment step.
   * @return {DG.TableView}
   * */
  getFilterView(): DG.TableView {
    const getComputeDialog = async () => {
      chemFunctionsDialog(async (resultMap) => {
        this.calculateColumns(resultMap, view);
      }, () => null, this.template!, true);
    };
    if (!this.dataFrame!.col(HitSelectionColName))
      this.dataFrame!.columns.addNewBool(HitSelectionColName).init(false);

    const view = DG.TableView.create(this.dataFrame!, false);
    const ribbons = view.getRibbonPanels();
    const calculateRibbon = ui.icons.add(getComputeDialog, 'Calculate additional properties');
    const submitButton = ui.div(ui.bigButton('Submit', () => {
      const dialogContent = this._submitView?.render();
      if (dialogContent)
        ui.dialog('Submit').add(dialogContent).show();
    }));

    ribbons.push([calculateRibbon, submitButton]);
    view.setRibbonPanels(ribbons);
    view.name = this._filterViewName;
    setTimeout(async () => {
      view._onAdded();
      await new Promise((r) => setTimeout(r, 1000));
      // const ribbons = view.getRibbonPanels();
      // ribbons[0].push(ui.div(ui.icons.add(() => null), {classes: 'd4-ribbon-item'}));
      // console.log(ribbons);
      // we need to wait for chem package to be initialized first to be able to use chem filters

      //const f = view.filters();
      const f = view.filters(this._campaignFilters ? {filters: this._campaignFilters} : undefined);

      view.dataFrame.onFilterChanged
        .subscribe((_) => {
          this._filterDescriptions = Array.from(view.dataFrame.rows.filters);
          this._campaignFilters = f.getOptions().look.filters;
        });
      setTimeout(() => {
        this._filterDescriptions = Array.from(view.dataFrame.rows.filters);
        this._campaignFilters = f.getOptions().look.filters;
      }, 300);
    }, 300);
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
    };
    await _package.files.writeAsText(`Hit Triage/campaigns/${campaignId}/${CampaignJsonName}`,
      JSON.stringify(campaign));

    const csvDf = DG.DataFrame.fromColumns(
      enrichedDf.columns.toList().filter((col) => !col.name.startsWith('~')),
    ).toCsv();
    await _package.files.writeAsText(`Hit Triage/campaigns/${campaignId}/${CampaignTableName}`, csvDf);
    notify && grok.shell.info('Campaign saved successfully.');
  }
}
