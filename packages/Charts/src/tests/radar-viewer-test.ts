import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {Subscription} from 'rxjs';
import {getOptions} from './utils';


category('Radar', () => {
  const TYPE = 'Radar';
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
      rowSource: 'Filtered',
      filter: '',
      title: 'Radar',
      min: '5',
      max: '95',
      showCurrentRow: true,
      showMouseOverRow: true,
      showMouseOverRowGroup: true, 
      showTooltip: true,
      colorColumnName: null,
      backgroundMinColor: 0xFFBB845D,
      backgroundMaxColor: 0xFFE7CDCD,
      currentRowColor: 0xFF00FF00,
      mouseOverRowColor: 0xAAAAAA,
      lineColor: 0xADD8E6,
      showMin: false,
      showMax: false,
      showValues: false,
      valuesColumnNames: ['AGE', 'COUNTRYID'],
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
