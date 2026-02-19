/*
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';

export class FunctionsView extends UaView {
  static viewName = 'Functions';

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
  }

  async initViewers() : Promise<void> {
    const functionErrorsViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Function Errors',
      queryName: 'FunctionErrors',
      viewerFunction: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions);
        viewer.root.style.maxHeight = '150px';
        return viewer;
      },
    });
    this.viewers.push(functionErrorsViewer);

    const topFunctionErrorsViewer = new TopFunctionErrorsViewer(
      'Function Errors', 'TopFunctionErrors', this.uaToolbox.filterStream);
    this.viewers.push(topFunctionErrorsViewer);

    const topFunctionNotErrorsViewer = new TopFunctionErrorsViewer(
      'Function Disabled Errors', 'TopFunctionDisabledErrors', this.uaToolbox.filterStream);
    this.viewers.push(topFunctionNotErrorsViewer);

    const topPackagesByErrorsViewer = new TopPackagesViewer(
      'Packages By Errors', 'TopPackagesByError', this.uaToolbox.filterStream);
    this.viewers.push(topPackagesByErrorsViewer);

    this.root.append(ui.divV([
      ui.div([ui.block([functionErrorsViewer.root])]),
      ui.div([ui.block50([topFunctionErrorsViewer.root]), ui.block50([topFunctionNotErrorsViewer.root])]),
      ui.div([ui.block50([topPackagesByErrorsViewer.root])]),
    ]));
  }
}
*/
