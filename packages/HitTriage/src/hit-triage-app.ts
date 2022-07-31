/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageBaseView} from "./views/hit-triage-base-view";
import {InfoView} from "./views/0-info-view";
import {GetMoleculesView} from "./views/1-get-molecules-view";
import {EnrichView} from "./views/2-enrich-view";
import {SubmitView} from "./views/4-submit-view";

export class HitTriageApp {
  template: HitTriageTemplate;
  multiView: DG.MultiView;

  private _infoView?: InfoView;
  private _getMoleculesView?: GetMoleculesView;
  private _enrichView?: EnrichView;
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
  get enrichView(): EnrichView { return this._enrichView ??= new EnrichView(this); }
  get filterView(): DG.TableView { return this._filterView ??= this.getFilterView(); }
  get submitView(): SubmitView { return this._submitView ??= new SubmitView(this); }

  get viewFactories() {
    return {
      '0. Info': () => this.infoView,
      '1. Get molecules': () => this.getMoleculesView,
      '2. Enrich': () => this.enrichView,
      '3. Filter': () => this.filterView,
      '4. Submit': () => this.submitView,
    }
  }

  /**
   * A view that lets you filter the molecules using either molecules, or
   * their properties derived at the enrichment step.
   * */
  getFilterView(): DG.TableView {
    const template = this.template;
    const view = DG.TableView.create(template.sourceDataFrame!, false);
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

  sourceQuery?: DG.FuncCall;
  sourceDataFrame?: DG.DataFrame;
  sourceMoleculeColumn: string = '';
  sourceType: string = 'file';
  sourceDescription: string = '';

  enrichedDataFrame?: DG.DataFrame;
  enrichmentSteps: string[] = [];
  enrichmentDescriptions: string[] = [];

  filterDescriptions: string[] = [];

  getSummary() {
    return {
      'From': this.sourceType,
      'Path': this.sourceDescription,
      'Total Molecules': this.sourceDataFrame!.rowCount,
      'Enrich': this.enrichmentDescriptions.join('\n'),
      'Filter': this.filterDescriptions,
      'Result Molecules': this.sourceDataFrame!.filter.trueCount
    }
  }

  static demo(): HitTriageTemplate {
    const session = new HitTriageTemplate();
    session.project = 'Demo project';
    session.sourceDataFrame = grok.data.demo.molecules(20000);
    session.sourceMoleculeColumn = 'smiles';
    session.sourceDataFrame.meta.detectSemanticTypes().then((_) => {});
    session.sourceType = 'file';
    session.sourceDescription = 'AppData:/HitTriage/campaigns/bfg9000.csv';
    return session;
  }
}
