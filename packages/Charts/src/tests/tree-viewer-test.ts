import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {getOptions} from './utils';


category('TreeViewer', () => {
  const TYPE = 'TreeViewer';
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog(20);
    tv = grok.shell.addTableView(df);
  });

  test('Creation', async () => {
    const viewer = DG.Viewer.fromType(TYPE, df);
    expect(viewer instanceof DG.JsViewer, true);
    expect(viewer.type, TYPE);
    expect(viewer.table.id, df.id);
  });

  test('Properties', async () => {
    const viewer = DG.Viewer.fromType(TYPE, df);
    const standardOptions = {
      top: '5px',
      left: '5px',
      bottom: '5px',
      right: '5px',
      animationDuration: 750,
      animationDurationUpdate: 750,
      animation: true,
      layout: 'orthogonal',
      orient: 'LR',
      expandAndCollapse: true,
      initialTreeDepth: 3,
      edgeShape: 'curve',
      symbol: 'emptyCircle',
      symbolSize: 7,
      sizeColumnName: null,
      sizeAggrType: DG.AGG.AVG,
      colorColumnName: null,
      colorAggrType: DG.AGG.AVG,
      hierarchyColumnNames: ['sex', 'disease', 'site'],
    };

    expect(JSON.stringify(standardOptions), JSON.stringify(await getOptions(viewer)));

    viewer.setOptions({
      symbolSize: 5,
      sizeColumnName: 'age',
      colorColumnName: 'weight',
      hierarchyColumnNames: ['disease', 'sex'],
    });

    const options = await getOptions(viewer);
    expect(options.symbolSize, 5);
    expect(options.sizeColumnName, 'age');
    expect(options.colorColumnName, 'weight');
    expectArray(options.hierarchyColumnNames, ['disease', 'sex']);
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
