import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class DistributionProfilerViewer extends DG.JsViewer {
  constructor() {
    super();
    this.init();
  }

  init() {
    let df = grok.data.demo.demog();

    var header = ui.inputs([
      ui.columnInput('Split by', df, df.col('sex'))
    ]);

    let host = ui.divV([
      header,
      DG.Viewer.grid(df).root,
      DG.Viewer.scatterPlot(df).root,
    ]);

    this.root.appendChild(host);
  }
}
