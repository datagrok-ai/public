import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitDesignCampaign, HitDesignTemplate} from './types';
import {HitDesignInfoView} from './hit-design-views/info-view';
import {CampaignTableName, EmptyStageCellValue, HitDesignCampaignIdKey,
  HitDesignMolColName, TileCategoriesColName, ViDColName, i18n} from './consts';
import {HitDesignBaseView} from './hit-design-views/base-view';
import {calculateSingleCellValues, getNewVid} from './utils/calculate-single-cell';
import '../../css/hit-triage.css';
import {_package} from '../package';
import {modifyUrl} from './utils';
import {HitDesignSubmitView} from './hit-design-views/submit-view';
import {HitDesignTilesView} from './hit-design-views/tiles-view';

export class HitDesignApp {
  template?: HitDesignTemplate;
  dataFrame?: DG.DataFrame;
  multiView: DG.MultiView;

  private _infoView: HitDesignInfoView;
  private _designView?: DG.TableView;
  public _submitView?: HitDesignSubmitView;
  private _tilesView?: HitDesignTilesView;
  private _designViewName = 'Design';
  private _campaignId?: string;
  private _dfName?: string;
  private _molColName: string = HitDesignMolColName;
  private _campaign?: HitDesignCampaign;
  public campaignProps: {[key: string]: any} = {};
  private processedValues: string[] = [];
  private _extraStageColsCount = 0;
  constructor() {
    this._infoView = new HitDesignInfoView(this);
    this.multiView = new DG.MultiView({viewFactories: {[this._infoView.name]: () => this._infoView}});
    this.multiView.tabs.onTabChanged.subscribe((_) => {
      if (this.multiView.currentView instanceof HitDesignBaseView)
        (this.multiView.currentView as HitDesignBaseView).onActivated();
    });
    grok.shell.addView(this.multiView);
  }

  public async setTemplate(template: HitDesignTemplate, campaignId?: string) {
    if (!campaignId) {
      campaignId = await this.getNewCampaignName(template.key);
      modifyUrl(HitDesignCampaignIdKey, campaignId);
    } else {
      this.dataFrame = await _package.files.readCsv(`Hit Design/campaigns/${campaignId}/${CampaignTableName}`);
      await this.dataFrame.meta.detectSemanticTypes();
    }

    if (!this.dataFrame) {
      console.error('DataFrame is empty');
      return;
    }
    this._molColName = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? HitDesignMolColName;
    this._dfName = this.dataFrame.name ??= 'Molecules';
    this._campaignId = campaignId;
    this.template = template;
    this.multiView.addView(this._designViewName, () => this.designView, true);
    this._tilesView ??= new HitDesignTilesView(this);
    this.multiView.addView(this._tilesView.name, () => this._tilesView!, false);
    this._submitView ??= new HitDesignSubmitView(this);
    // this.multiView.addView(this._submitView.name, () => this._submitView!, false);
    grok.shell.windows.showHelp = false;
    //add empty rows to define stages, used for tile categories;
    const stagesRow = this.dataFrame.getCol(TileCategoriesColName);
    if (stagesRow) {
      const categories = stagesRow.categories;
      if (categories && categories.length) {
        template.stages.forEach((s) => {
          if (!categories.includes(s)) {
            const newRow = this.dataFrame!.rows.addNew();
            const idx = newRow.idx;
            this.dataFrame!.set(TileCategoriesColName, idx, s);
            this.dataFrame!.set(ViDColName, idx, EmptyStageCellValue);
          }
        });
      }
    }
    this.dataFrame.rows.filter((r) => r[ViDColName] !== EmptyStageCellValue);
    this._extraStageColsCount = this.dataFrame!.rowCount - this.dataFrame.filter.trueCount;
  }

  get campaignId(): string | undefined {return this._campaignId;}

  get designView(): DG.TableView {return this._designView ??= this.getDesignView();}

  get molColName() {
    return this._molColName ??= this.dataFrame?.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? HitDesignMolColName;
  }

  get campaign(): HitDesignCampaign | undefined {return this._campaign;}

  set campaign(campaign: HitDesignCampaign | undefined) {this._campaign = campaign;}

  private getDesignView(): DG.TableView {
    const isNew = this.dataFrame!.rowCount === 1;
    const view = DG.TableView.create(this.dataFrame!, false);
    view.name = this._designViewName;
    this.processedValues = this.dataFrame!.getCol(this.molColName).toList();
    setTimeout(async () => {
      view._onAdded();
      await new Promise((r) => setTimeout(r, 1000)); // needed for substruct filter
      if (isNew)
        grok.functions.call('Chem:editMoleculeCell', {cell: view.grid.cell(this._molColName, 0)});
      //const _f = view.filters();

      this.dataFrame!.onRowsAdded.subscribe(() => { // TODO, insertion of rows in the middle
        const newRowsNum = this.dataFrame!.rowCount - this.processedValues.length;
        this.processedValues.push(...new Array(newRowsNum).fill(''));
        //const lastCell = view.grid.cell(this.molColName, this.dataFrame!.rowCount - 1);
        //view.grid.onCellValueEdited
      });
      this.dataFrame?.onRowsRemoved.subscribe(() => {
        this.processedValues = this.dataFrame!.getCol(this.molColName).toList();
      });

      this.dataFrame!.onValuesChanged.subscribe(async () => {
        let newValueIdx: number | null = null;
        for (let i = 0; i < this.processedValues.length; i++) {
          if (this.processedValues[i] !== this.dataFrame!.get(this.molColName, i)) {
            newValueIdx = i;
            break;
          }
        }
        if (newValueIdx == null) return;
        this.processedValues[newValueIdx] = this.dataFrame!.get(this.molColName, newValueIdx);
        const newCellValue: string = this.processedValues[newValueIdx];
        this.dataFrame!.col(TileCategoriesColName)!.set(newValueIdx, this.template!.stages[0], false);
        if (!this.dataFrame!.col(ViDColName)?.get(newValueIdx))
          this.dataFrame!.col(ViDColName)!.set(newValueIdx, getNewVid(this.dataFrame!.col(ViDColName)!), false);

        const computeObj = this.template!.compute;
        const calcDf =
          await calculateSingleCellValues(newCellValue, computeObj.descriptors.args, computeObj.functions);
        for (const col of calcDf.columns.toList()) {
          if (col.name === HitDesignMolColName) continue;
          if (!this.dataFrame!.columns.contains(col.name))
            this.dataFrame!.columns.addNew(col.name, col.type);
          this.dataFrame!.col(col.name)!.set(newValueIdx, col.get(0), false);
        }
        this.dataFrame!.fireValuesChanged();
      });
    }, 300);
    const ribbons = view?.getRibbonPanels();
    if (ribbons) {
      const submitButton = ui.div(ui.bigButton('Submit', () => {
        const dialogContent = this._submitView?.render();
        if (dialogContent)
          ui.dialog('Submit').add(dialogContent).show();
      }));
      ribbons.push([submitButton]);
      view.setRibbonPanels(ribbons);
    }
    return view;
  }
  async getNewCampaignName(templateKey: string) {
    const templateCampaigns = (await _package.files.list('Hit Design/campaigns'))
      .map((file) => file.name)
      .filter((name) => name.startsWith(templateKey));
    if (templateCampaigns.length === 0)
      return templateKey + '-1';
    const postFixes = templateCampaigns.map((c) => c.split('-')[1]).filter(Boolean).map((c) => parseInt(c, 10)).sort();
    return templateKey + '-' + ((postFixes[postFixes.length - 1] + 1).toString());
  }


  getSummary(): {[_: string]: any} {
    const campaignProps = this.campaign?.campaignFields ?? this.campaignProps;
    return {
      'Template': this.template?.name ?? 'Molecules',
      'File path': ui.divH([ui.link(this._dfName ?? '-',
        () => this.download(this.dataFrame!, 'Molecules'), i18n.download)],
      {style: {alignItems: 'center'}}),
      ...campaignProps,
      'Number of molecules': (this.dataFrame!.rowCount - this._extraStageColsCount).toString(),
      'Enrichment methods': [this.template!.compute.descriptors.enabled ? 'descriptors' : '',
        ...this.template!.compute.functions.map((func) => func.name)].filter((f) => f && f.trim() !== '').join(', '),
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
}
