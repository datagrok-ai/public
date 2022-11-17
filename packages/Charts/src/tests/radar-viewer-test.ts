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
    df = DG.DataFrame.fromCsv(
        `USUBJID, SEX, AGE, COUNTRYID
      s1, F, 48, 12
      s2, M, 51, 3
      s3, F, 39, 5
      s4, M, 43, 89`);
    tv = grok.shell.addTableView(df);
  });

  test('Creation', async () => {
    const viewer = DG.Viewer.fromType(TYPE, df);
    expect(viewer instanceof DG.JsViewer, true);
    expect(viewer.type, TYPE);
    expect(viewer.table.id, df.id);
  });

  test('Standard Properties', async () => {
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
      valuesColumnNames: null,
    };
    expect(JSON.stringify(standardOptions), JSON.stringify(await getOptions(viewer)));
  });

  test('Changed properties', async () => {
    const changedViewer = DG.Viewer.fromType(TYPE, df, {
      min: '10',
      showTooltip: false,
      backgroundMaxColor: 0xFF4F616F,
      showValues: false,
    });

    const options = await getOptions(changedViewer);
    expect(options.min, '10');
    expect(options.showTooltip, false);
    expect(options.backgroundMaxColor, 0xFF4F616F);
    expect(options.showValues, false);
  });

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    tv.close();
    grok.shell.closeTable(df);
  });
});
