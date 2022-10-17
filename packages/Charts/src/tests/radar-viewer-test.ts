import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {getOptions} from './utils';


category('Radar', () => {
  const TYPE = 'RadarViewer';
  let df: DG.DataFrame;
  let tv: DG.TableView;
  const subs: Subscription[] = [];
  
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
      min: '5',
      max: '95',
      showCurrentRow: false,
      showTooltip: true,
      showAllRows: false,
      backgroundMinColor: 0xFFB0D7FF,
      backgroundMaxColor: 0xFFBCE2F5,
      showMin: true,
      showMax: true,
      showValues: true,
      columnNames: null,
    };
    expect(JSON.stringify(standardOptions), JSON.stringify(await getOptions(viewer)));
  });

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    tv.close();
    grok.shell.closeTable(df);
  });
});
