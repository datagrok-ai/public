import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
//import $ from 'cash-dom';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {Subscription} from 'rxjs';
import {getOptions} from './utils';


category('Timelines', () => {
  const TYPE = 'Timelines';
  let df: DG.DataFrame;
  let tv: DG.TableView;
  const subs: Subscription[] = [];

  before(async () => {
    df = DG.DataFrame.fromCsv(
      `USUBJID, AESTDY, AEENDY, SEX, AGE
      s1, 10/02/2018, 10/09/2018, F, 48
      s2, 10/04/2018, 10/07/2018, M, 51
      s3, 10/02/2018, 10/05/2018, F, 39
      s4, 10/07/2018, 10/08/2018, M, 43`);
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
      rowSource: 'Filtered',
      filter: '',
      splitByColumnName: 'USUBJID',
      startColumnName: 'AESTDY',
      endColumnName: 'AEENDY',
      colorColumnName: 'SEX',
      showOpenIntervals: false,
      eventColumnName: 'SEX',
      eventsColumnNames: null,
      showEventInTooltip: true,
      marker: 'circle',
      markerSize: 6,
      markerPosition: 'main line',
      lineWidth: 3,
      dateFormat: '{MMM} {d}',
      axisPointer: 'shadow',
      showZoomSliders: true,
      legendVisibility: 'Auto',
    };

    expect(JSON.stringify(standardOptions), JSON.stringify(await getOptions(viewer)));

    viewer.setOptions({
      colorColumnName: 'USUBJID',
      markerSize: 4,
      dateFormat: '{yyyy}-{MM}-{dd}',
      showZoomSliders: false,
    });

    const options = await getOptions(viewer);
    expect(options.colorColumnName, 'USUBJID');
    expect(options.markerSize, 4);
    expect(options.dateFormat, '{yyyy}-{MM}-{dd}');
    expect(options.showZoomSliders, false);
  });

  // test('Context menu', () => new Promise((resolve, reject) => {
  //   const viewer = DG.Viewer.fromType(TYPE, df);
  //   subs.push(viewer.onContextMenu.subscribe((menu) => {
  //     const customCommand = menu.find('Reset View');
  //     if (customCommand == null)
  //       reject('Command not found');
  //     else
  //       resolve('OK');
  //   }));

  //   setTimeout(() => reject('Timeout'), 50);
  //   $(viewer.root).trigger('contextmenu');
  // }));

  after(async () => {
    subs.forEach((sub) => sub.unsubscribe());
    tv.close();
    grok.shell.closeTable(df);
  });
});
