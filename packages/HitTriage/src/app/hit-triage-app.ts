import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageCampaign, IComputeDialogResult, IFunctionArgs,
  HitTriageTemplate, HitTriageTemplateIngest, IngestType} from './types';
import {InfoView} from './hit-triage-views/info-view';
import {HitTriageBaseView} from './hit-triage-views/base-view';
import {SubmitView} from './hit-triage-views/submit-view';
import {CampaignIdKey, HitSelectionColName, i18n} from './consts';
import {modifyUrl} from './utils';
import {_package} from '../package';
import '../../css/hit-triage.css';
import {chemFunctionsDialog} from './dialogs/functions-dialog';
export class HitTriageApp {
  template?: HitTriageTemplate;
  dataFrame?: DG.DataFrame;
  multiView: DG.MultiView;

  private _infoView: InfoView;
  private _pickView?: DG.TableView;
  private _submitView?: SubmitView;

  private _filterViewName = 'Pick';
  private _campaignFilters?: {[key: string]: any}[];
  private _campaignId?: string;
  private _dfName?: string;
  private _molColName?: string;
  public _fileInputType?: IngestType;

  private _campaign?: HitTriageCampaign;
  protected _filterDescriptions: string[] = [];
  public campaignProps: {[key: string]: any} = {};
  constructor() {
    this._infoView = new InfoView(this);
    this.multiView = new DG.MultiView({viewFactories: {[this._infoView.name]: () => this._infoView}});
    this.multiView.tabs.onTabChanged.subscribe((_) => {
      if (this.multiView.currentView instanceof HitTriageBaseView)
        (this.multiView.currentView as HitTriageBaseView).onActivated();
    });
    grok.shell.addView(this.multiView);
  }

  public async setTemplate(template: HitTriageTemplate, presetFilters?: {[key: string]: any}[],
    campaignId?: string, ingestProps?: HitTriageTemplateIngest) {
    if (!campaignId) {
      campaignId = await this.getNewCampaignName(template.key);
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

    //this.currentView.root.appendChild(this.pickView.root);
    this.multiView.addView(this._filterViewName, () => this.pickView, true);
    this._submitView ??= new SubmitView(this);
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

  private getFilterType(colName: string): DG.FILTER_TYPE {
    const col = this.dataFrame!.col(colName);
    if (col?.semType === DG.SEMTYPE.MOLECULE)
      return DG.FILTER_TYPE.SUBSTRUCTURE;
    if (col?.type === DG.COLUMN_TYPE.BOOL)
      return DG.FILTER_TYPE.BOOL_COLUMNS;
    if (col?.type === DG.COLUMN_TYPE.STRING)
      return DG.FILTER_TYPE.CATEGORICAL;
    return DG.FILTER_TYPE.HISTOGRAM;
  }


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
      // view.getFiltersGroup().add()
      // const group = view.getFiltersGroup();
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

  download(df: DG.DataFrame, name: string, onlyFiltered = false): void {
    const element = document.createElement('a');
    const result = DG.DataFrame.fromColumns(df.columns.toList().filter((c) => !c.name.startsWith('~')))
      .toCsv({filteredRowsOnly: onlyFiltered});
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
    element.setAttribute('download', name + '.csv');
    element.click();
  }

  async getNewCampaignName(templateKey: string) {
    const templateCampaigns = (await _package.files.list('Hit Triage/campaigns'))
      .map((file) => file.name)
      .filter((name) => name.startsWith(templateKey));
    if (templateCampaigns.length === 0)
      return templateKey + '-1';
    const postFixes = templateCampaigns.map((c) => c.split('-')[1]).filter(Boolean).map((c) => parseInt(c, 10)).sort();
    return templateKey + '-' + ((postFixes[postFixes.length - 1] + 1).toString());
  }
}
