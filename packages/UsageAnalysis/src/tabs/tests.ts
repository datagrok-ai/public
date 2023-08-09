import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';

export class TestsView extends UaView {
  loader = ui.div([ui.loader()], 'grok-wait');

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Tests';
  }

  async initViewers(): Promise<void> {
    this.root.appendChild(this.loader);
    const df: DG.DataFrame = await grok.data.query('UsageAnalysis:TestsToday');
    df.getCol('id').name = '~id';
    // const tv = DG.TableView.create(df, false);
    df.getCol('status').colors.setCategorical({passed: '#3CB173', failed: '#EB6767', skipped: '#FFA24A'});
    // const sl = DG.Viewer.fromType('Line chart', df);
    this.root.removeChild(this.loader);
    this.root.append(ui.splitV([
      // ui.box(sl.root, {style: {maxHeight: '20%'}}),
      df.plot.grid().root,
    ]));
    // tv.getFiltersGroup();
    df.onCurrentRowChanged.subscribe(async () => {
      const acc = DG.Accordion.create();
      const history: DG.DataFrame = await grok.data.query('UsageAnalysis:ScenarioHistory', {id: df.currentRow.get('~id')});
      acc.addPane('History', () => ui.waitBox(async () => {
        history.getCol('status').colors.setCategorical({passed: '#3CB173', failed: '#EB6767', skipped: '#FFA24A'});
        const grid = history.plot.grid();
        grid.col('date')!.format = 'MM/dd/yyyy HH:mm:ss';
        return grid.root;
      }), true);
      acc.addPane('Execution time', () => ui.waitBox(async () => {
        const lct = DG.Viewer.fromType('Line chart', history, execTimeStyle);
        return lct.root;
      }), true);
      grok.shell.o = acc.root;
    });
  }
}

const execTimeStyle = {
  x: 'date',
  y: 'ms',
  // showMouseOverRowLine: false,
  showXSelector: false,
  showYSelectors: false,
  showAggrSelectors: false,
  showSplitSelector: false,
  showXAxis: true,
  showYAxis: true,
  showMarkers: 'Never',
  lineWidth: 2,
  autoLayout: false,
};
