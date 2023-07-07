import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class DistributionProfilerViewer extends DG.JsViewer {
  constructor() {
    super();
    this.init();
  }

  init() {
    const df = grok.data.demo.demog();

    const header = ui.inputs([
      ui.columnInput('Split by', df, df.col('sex')),
    ]);

    const host = ui.divV([
      header,
      DG.Viewer.grid(df).root,
      DG.Viewer.scatterPlot(df).root,
    ]);

    this.root.appendChild(host);
  }
}
