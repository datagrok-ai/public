import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilter} from "../filter2";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {UaQueryViewer} from "../viewers/abstract/ua-query-viewer";
import {TopErrorsViewer} from "../drilldown_viewers/errors/top-errors-viewer";
import {TopErrorSourcesViewer} from "../drilldown_viewers/errors/top-error-sources-viewer";

export class ErrorsView extends UaView {

  constructor(uaToolbox: UaToolbox) {
    super('Errors', uaToolbox);
  }

  async initViewers() : Promise<void> {
    let errorsViewer = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Errors',
        'Errors1',
        (t: DG.DataFrame) => DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions).root
    );
    this.viewers.push(errorsViewer);

    let topErrorsViewer = new TopErrorsViewer('Errors', 'TopErrors', this.uaToolbox.filterStream);
    this.viewers.push(topErrorsViewer);

    let topNotErrorsViewer = new TopErrorsViewer('Disabled Errors', 'TopDisabledErrors', this.uaToolbox.filterStream);
    this.viewers.push(topNotErrorsViewer);

    let topErrorSourcesViewer = new TopErrorSourcesViewer(this.uaToolbox.filterStream);
    this.viewers.push(topErrorSourcesViewer);

    this.root.append(ui.divV([
      ui.div([ui.block([errorsViewer.root])]),
      ui.divH([ui.block50([topErrorsViewer.root]), ui.block50([topNotErrorsViewer.root])]),
      ui.div([ui.block50([topErrorSourcesViewer.root])])
    ]));

  }
}