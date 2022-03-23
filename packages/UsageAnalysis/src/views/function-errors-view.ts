import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {UaQueryViewer} from "../viewers/abstract/ua-query-viewer";
import {TopFunctionErrorsViewer} from "../drilldown_viewers/function_errors/top-function-errors-viewer";
import {TopPackagesViewer} from "../drilldown_viewers/events/top-packages-viewer";

export class FunctionErrorsView extends UaView {

  constructor(uaToolbox: UaToolbox) {
    super('Function Errors', uaToolbox);
  }

  async initViewers() : Promise<void> {
    let functionErrorsViewer = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Function Errors',
        'FunctionErrors',
        (t: DG.DataFrame) => DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions).root
    );
    this.viewers.push(functionErrorsViewer);

    let topFunctionErrorsViewer = new TopFunctionErrorsViewer('Function Errors', 'TopFunctionErrors', this.uaToolbox.filterStream);
    this.viewers.push(topFunctionErrorsViewer);

    let topFunctionNotErrorsViewer = new TopFunctionErrorsViewer('Function Disabled Errors', 'TopFunctionDisabledErrors', this.uaToolbox.filterStream);
    this.viewers.push(topFunctionNotErrorsViewer);

    let topPackagesByErrorsViewer = new TopPackagesViewer('Packages By Errors', 'TopPackagesByError', this.uaToolbox.filterStream);
    this.viewers.push(topPackagesByErrorsViewer);

    this.root.append(ui.divV([
      ui.div([ui.block([functionErrorsViewer.root])]),
      ui.div([ui.block50([topFunctionErrorsViewer.root]), ui.block50([topFunctionNotErrorsViewer.root])]),
      ui.div([ui.block50([topPackagesByErrorsViewer.root])])
    ]));

  }
}