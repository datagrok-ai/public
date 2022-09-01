import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';


category('Timelines', () => {
  const TYPE = 'TimelinesViewer';
  let df: DG.DataFrame;
  let tv: DG.TableView;
  
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
      splitByColumnName: 'USUBJID',
      startColumnName: 'AESTDY',
      endColumnName: 'AEENDY',
      colorByColumnName: 'SEX',
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
      colorByColumnName: 'USUBJID',
      markerSize: 4,
      dateFormat: '{yyyy}-{MM}-{dd}',
      showZoomSliders: false,
    });

    const options = await getOptions(viewer);
    expect(options.colorByColumnName, 'USUBJID');
    expect(options.markerSize, 4);
    expect(options.dateFormat, '{yyyy}-{MM}-{dd}');
    expect(options.showZoomSliders, false);
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});

async function getOptions(viewer: DG.Viewer): Promise<{ [key: string]: any; }> {
  await delay(100);
  const options = (viewer.getOptions(true) as {id: string, type: string, look: {[key: string]: any}}).look;
  delete options['#type'];
  delete options.table;
  return options;
}
