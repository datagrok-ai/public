/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from "./package";
import {HitTriageSession} from "./hit-triage-session";
import {divText} from "datagrok-api/ui";

const session = HitTriageSession.demo();

export function hitTriageView(): DG.MultiView {
  const multiView = new DG.MultiView({
    viewFactories: {
      '0. Info': () => new InfoView(),
      '1. Get molecules': () => new GetMoleculesView(),
      '2. Enrich': () => new EnrichView(),
      '3. Filter': getFilterView,
      '4. Submit': () => new SubmitView(),
    }
  });

  multiView.tabs.onTabChanged.subscribe((_) => {
    if (multiView.currentView instanceof HitTriageBaseView)
      (multiView.currentView as HitTriageBaseView).onActivated();
  });

  return multiView;
}


class HitTriageBaseView extends DG.ViewBase {
  constructor() {
    super();
    this.root.classList.add('grok-hit-triage-view');
    this.root.style.display = 'flex';
    this.statusBarPanels = [divText('Hit Triage')]
  }

  /** Override to initialize the view based on the session. */
  onActivated(): void {}

  async process(): Promise<any> {}
}


export class InfoView extends HitTriageBaseView {

  constructor() {
    super();
    _package.files.readAsText('README.md').then((md) => {
      this.root.appendChild(ui.markdown(md));
    });
  }
}

export class GetMoleculesView extends HitTriageBaseView {

  constructor() {
    super();
    this.name = 'Source Molecules';

    const from = ui.choiceInput('From', 'file', ['file', 'database', 'webservice']);
    const content = ui.divV([
      ui.divH([
        ui.divText('Ingest', {style: {'font-weight': 'bold'}}),
        from.root,
        ui.divText(session.sourceDescription)],
        {style: {'display': 'flex', 'align-items': 'center', 'gap': '12px'}}
      ),
      session.sourceDataFrame!.plot.grid()
    ])

    this.root.appendChild(content);
  }
}


/**
 * Enrichment of the molecular dataset.
 **/
export class EnrichView extends HitTriageBaseView {

  grid: DG.Grid;

  constructor() {
    super();

    this.grid = session.sourceDataFrame!.plot.grid()
    const content = ui.divV([
      ui.h1('This is where we enrich our data'),
      this.grid
    ])

    this.root.appendChild(content);
  }

  onActivated() {
    super.onActivated();
    this.process().then(() => {});
  }

  async process(): Promise<any> {
    session.enrichedDataFrame = session.sourceDataFrame;
    await session.enrichedDataFrame!.columns.addNewCalculated("length", "length(${smiles})");
    this.grid.dataFrame = session.enrichedDataFrame!;
  }
}


/**
 * A view that lets you filter the molecules using either molecules, or
 * their properties derived at the enrichment step.
 * */
function getFilterView() {
  const view = DG.TableView.create(session.sourceDataFrame!, false);
  setTimeout(function () {
    view._onAdded();
    view.scatterPlot();
    view.filters();
    view.dataFrame.onFilterChanged.subscribe((_) => session.filterDescriptions = Array.from(view.dataFrame.rows.filters))
  }, 100);
  return view;
}


export class SubmitView extends HitTriageBaseView {

  constructor() {
    super();
    this.render();
  }

  render(): void {
    ui.empty(this.root);
    const content = ui.divV([
      ui.h1('Summary'),
      ui.div([ui.tableFromMap(session.getSummary())]),
      ui.divH([ui.bigButton('SUBMIT', () => grok.shell.info('Good job!'))])
    ])
    this.root.appendChild(content);
  }
}