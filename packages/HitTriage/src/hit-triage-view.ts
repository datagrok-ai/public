/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from "./package";
import {HitTriageSession} from "./hit-triage-session";

const session = HitTriageSession.demo();

export function hitTriageView(): DG.MultiView {
  return new DG.MultiView({
    viewFactories: {
      '0. Info': () => new InfoView(),
      '1. Get molecules': () => new GetMoleculesView(),
      '2. Enrich': () => new EnrichView(),
      '3. Filter': getFilterView,
      '4. Submit': () => new SubmitView(),
    }
  });
}


class HitTriageBaseView extends DG.ViewBase {

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

    const content = ui.divV([
      ui.h1('Getting molecules'),
      session.sourceDataFrame!.plot.grid()
    ])

    this.root.appendChild(content);
  }
}


/**
 * Enrichment of the molecular dataset.
 **/
export class EnrichView extends HitTriageBaseView {

  constructor() {
    super();

    const content = ui.divV([
      ui.h1('This is where we enrich our data'),
      session.sourceDataFrame!.plot.grid()
    ])

    this.root.appendChild(content);
  }
}


/**
 * A view that lets you filter the molecules using either molecules, or
 * their properties derived at the enrichment step.
 * */
function getFilterView() {
  const view = DG.TableView.create(session.sourceDataFrame!, false);
  view.scatterPlot();
  view.filters();
  return view;
}


export class SubmitView extends HitTriageBaseView {

  constructor() {
    super();

    const content = ui.divV([
      ui.h1('This is it!'),
      ui.bigButton('SUBMIT', () => grok.shell.info('Good job!'))
    ])

    this.root.appendChild(content);
  }
}


