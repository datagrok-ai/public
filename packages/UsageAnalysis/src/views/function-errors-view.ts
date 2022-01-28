import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";
import {UaQueryViewer} from "../viewers/ua-query-viewer";
import {TopFunctionErrorsViewer} from "../drilldown_viewers/function_errors/top-function-errors-viewer";

export class FunctionErrorsView extends UaView {

  constructor(uaToolbox: UaToolbox) {
    super('Function Errors', uaToolbox);
  }

  async initViewers() : Promise<void> {
    let functionErrorsViewer = new UaFilterableViewer(
        this.uaToolbox.filterStream,
        'Function Errors',
        'FunctionErrors',
        (t: DG.DataFrame) => DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions).root
    );
    this.viewers.push(functionErrorsViewer);

    let topFunctionErrorsViewer = new TopFunctionErrorsViewer('Function Errors', 'TopFunctionErrors', this.uaToolbox.filterStream);
    this.viewers.push(topFunctionErrorsViewer);

    let topFunctionNotErrorsViewer = new TopFunctionErrorsViewer('Function Not Errors', 'TopFunctionNotErrors', this.uaToolbox.filterStream);
    this.viewers.push(topFunctionNotErrorsViewer);

    this.root.append(ui.divV([
      ui.div([ui.block([functionErrorsViewer.root])]),
      ui.div([ui.block50([topFunctionErrorsViewer.root]), ui.block50([topFunctionNotErrorsViewer.root])])
    ]));

  }
}