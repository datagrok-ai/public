import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ICampaign, ITemplate, ITemplateIngest, IngestType} from './types';
import {InfoView} from './views/info-view';
import {HitTriageBaseView} from './views/base-view';
import {ComputeView} from './views/compute-view';
import {SubmitView} from './views/submit-view';
import {CampaignIdKey, HitSelectionColName} from './consts';
import {modifyUrl} from './utils';
import {_package} from '../package';
import '../../css/hit-triage.css';

export class HitTriageApp {
  template?: ITemplate;
  dataFrame?: DG.DataFrame;
  multiView: DG.MultiView;

  private _infoView: InfoView;
  private _computeView?: ComputeView;
  private _pickView?: DG.TableView;
  private _submitView?: SubmitView;

  private _filterViewName = 'Pick';
  private _campaignFilters?: {[key: string]: any}[];
  private _campaignId?: string;
  private _dfName?: string;
  private _molColName?: string;
  public _fileInputType?: IngestType;

  private _campaign?: ICampaign;
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

  public async setTemplate(template: ITemplate, presetFilters?: {[key: string]: any}[],
    campaignId?: string, ingestProps?: ITemplateIngest) {
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
    this._campaignId = campaignId;
    this.template = template;
    this._campaignFilters = presetFilters;
    await this.dataFrame.meta.detectSemanticTypes();
    this._computeView ??= new ComputeView(this);
    this.multiView.addView(this._computeView.name, () => this._computeView!, !presetFilters);
    this.multiView.addView(this._filterViewName, () => this.pickView, !!presetFilters);
    this._submitView ??= new SubmitView(this);
    this.multiView.addView(this._submitView.name, () => this._submitView!, false);
    grok.shell.windows.showHelp = false;
  }

  get filterSettings(): {[key: string]: any}[] | undefined {return this._campaignFilters;}

  get campaignId(): string | undefined {return this._campaignId;}

  get pickView(): DG.TableView {return this._pickView ??= this.getFilterView();}

  get molColName(): string | undefined {return this._molColName;}

  get fileInputType(): IngestType | undefined {return this._fileInputType;}

  get campaign(): ICampaign | undefined {return this._campaign;}

  set campaign(campaign: ICampaign | undefined) {this._campaign = campaign;}
  /**
   * A view that lets you filter the molecules using either molecules, or
   * their properties derived at the enrichment step.
   * @return {DG.TableView}
   * */
  getFilterView(): DG.TableView {
    if (!this.dataFrame!.col(HitSelectionColName))
      this.dataFrame!.columns.addNewBool(HitSelectionColName).init(false);

    const view = DG.TableView.create(this.dataFrame!, false);
    view.name = this._filterViewName;
    setTimeout(async () => {
      console.log(view.getRibbonPanels());
      view._onAdded();
      // we need to wait for chem package to be initialized first to be able to use chem filters
      // const someChemFunction = DG.Func.find({package: 'Chem', name: 'substructureFilter'})[0];
      // await someChemFunction.package.init();
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
    }, 100);
    return view;
  }

  getSummary(): {[_: string]: any} {
    const getStyledDownloadButton = (callBackFn: () => void): HTMLElement => {
      const db = ui.iconFA('arrow-to-bottom', callBackFn);
      db.classList.add('hit-triage-download-button');
      return db;
    };
    const campaignProps = this.campaign?.campaignFields ?? this.campaignProps;
    return {
      'Template': this.template?.name ?? 'Molecules',
      'File path': ui.divH([ui.divText(this._dfName ?? ''), getStyledDownloadButton(
        () => this.download(this.dataFrame!, 'molecules'))], {style: {alignItems: 'center'}}),
      ...campaignProps,
      'Number of molecules': this.dataFrame!.rowCount.toString(),
      'Enrichment methods': [this.template!.compute.descriptors.enabled ? 'descriptors' : '',
        ...this.template!.compute.functions.map((func) => func.name)].filter((f) => f && f.trim() !== '').join(', '),
      'Filters': this._filterDescriptions.join(', '),
      'Result Molecules': ui.divH([ui.divText(this.dataFrame?.filter.trueCount.toString() ?? '0'),
        getStyledDownloadButton(
          () => this.download(this.dataFrame!, 'hits', true))], {style: {alignItems: 'center'}}),
    };
  }

  download(df: DG.DataFrame, name: string, onlyFiltered = false): void {
    const element = document.createElement('a');
    const result = df.toCsv({filteredRowsOnly: onlyFiltered});
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
    element.setAttribute('download', name + '.csv');
    element.click();
  }

  async getNewCampaignName(templateKey: string) {
    const templateCampaigns = (await _package.files.list('campaigns'))
      .map((file) => file.name)
      .filter((name) => name.startsWith(templateKey));
    if (templateCampaigns.length === 0)
      return templateKey + '-1';
    const postFixes = templateCampaigns.map((c) => c.split('-')[1]).filter(Boolean).map((c) => parseInt(c, 10)).sort();
    return templateKey + '-' + ((postFixes[postFixes.length - 1] + 1).toString());
  }
}
