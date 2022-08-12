/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageBaseView} from "./views/hit-triage-base-view";
import {InfoView} from "./views/0-info-view";
import {GetMoleculesView} from "./views/1-ingest-view";
import {ComputeView} from "./views/2-compute-view";
import {SubmitView} from "./views/4-submit-view";

export class HitTriageApp {
  template: HitTriageTemplate;
  multiView: DG.MultiView;

  private _infoView?: InfoView;
  private _getMoleculesView?: GetMoleculesView;
  private _computeView?: ComputeView;
  private _filterView?: DG.TableView;
  private _submitView?: SubmitView;

  /** Creates and starts the application. */
  constructor(session: HitTriageTemplate) {
    this.template = session;

    this.multiView = new DG.MultiView({ viewFactories: this.viewFactories });

    this.multiView.tabs.onTabChanged.subscribe((_) => {
      if (this.multiView.currentView instanceof HitTriageBaseView)
        (this.multiView.currentView as HitTriageBaseView).onActivated();
    });

    grok.shell.addView(this.multiView);
  }

  get infoView(): InfoView { return this._infoView ??= new InfoView(this); }
  get getMoleculesView(): GetMoleculesView { return this._getMoleculesView ??= new GetMoleculesView(this); }
  get computeView(): ComputeView { return this._computeView ??= new ComputeView(this); }
  get filterView(): DG.TableView { return this._filterView ??= this.getFilterView(); }
  get submitView(): SubmitView { return this._submitView ??= new SubmitView(this); }

  get viewFactories() {
    return {
      '0. Info': () => this.infoView,
      '1. Ingest': () => this.getMoleculesView,
      '2. Compute': () => this.computeView,
      '3. Pick': () => this.filterView,
      '4. Submit': () => this.submitView,
    }
  }

  /**
   * A view that lets you filter the molecules using either molecules, or
   * their properties derived at the enrichment step.
   * */
  getFilterView(): DG.TableView {
    const template = this.template;
    const view = DG.TableView.create(template.hitsTable!, false);
    view.name = 'Hit Triage | Pick'
    setTimeout(function () {
      view._onAdded();
      view.scatterPlot();
      view.filters();
      view.dataFrame.onFilterChanged.subscribe((_) => template.filterDescriptions = Array.from(view.dataFrame.rows.filters))
    }, 100);
    return view;
  }
}



export class HitTriageTemplate {
  project: string = 'New project';

  hitsQuery: string = '';
  hitsTable?: DG.DataFrame;

  /** Name of the column that contains molecules */
  hitsMolColumn: string = '';
  hitsIdColumn: string = '';
  sourceType: string = 'file';
  sourceDescription: string = '';

  hitsTargetsQuery: string = '';

  /** Format: Hit, Target */
  hitsTargetsTable?: DG.DataFrame;

  enrichedTable?: DG.DataFrame;
  enrichmentSteps: string[] = [];
  enrichmentDescriptions: string[] = [];

  filterDescriptions: string[] = [];

  getSummary() {
    return {
      'From': this.sourceType,
      'Path': this.sourceDescription,
      'Total Molecules': this.hitsTable!.rowCount,
      'Enrich': this.enrichmentDescriptions.join('\n'),
      'Filter': this.filterDescriptions,
      'Result Molecules': this.hitsTable!.filter.trueCount
    }
  }

  async loadData(): Promise<void> {
    this.hitsTable = await grok.dapi.files.readCsv(this.hitsQuery);
    this.hitsTargetsTable = await grok.dapi.files.readCsv(this.hitsTargetsQuery);
    await this.hitsTable.meta.detectSemanticTypes();
    await this.hitsTargetsTable.meta.detectSemanticTypes();
  }

  static demo(): HitTriageTemplate {
    const session = new HitTriageTemplate();
    session.project = 'Demo project';
    session.hitsQuery = 'System:AppData/HitTriage/campaigns/chembl-101/hits.csv';
    session.hitsTargetsQuery = 'System:AppData/HitTriage/campaigns/chembl-101/hits-targets.csv';
    session.hitsTable = grok.data.demo.molecules(20000);
    session.hitsMolColumn = 'smiles';
    session.hitsTable.meta.detectSemanticTypes().then((_) => {});
    session.sourceType = 'file';
    session.sourceDescription = 'AppData:/HitTriage/campaigns/chembl-101';
    return session;
  }
}


export class HitTriageSession {

}