import {after, before, category, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('Benchmarks', () => {
  let tv: DG.TableView;
  const df = grok.data.demo.randomWalk(1000000, 1);
  const dfSmall = grok.data.demo.randomWalk(100, 1);
  let detectors: DG.Func[];

  before(async () => {
    detectors = await grok.dapi.functions.filter('#semTypeDetector').list();
    tv = grok.shell.addTableView(dfSmall);
    const colSmall = dfSmall.columns.byIndex(0);
    for (const d of detectors) d.apply({col: colSmall});
    tv.close();
    grok.shell.closeTable(dfSmall);
  });

  test('detectors', async () => {
    const violatingDetectors: string[] = [];
    const res: any = [];
    let start: number;
    const col = df.columns.byIndex(0);
    for (const d of detectors) {
      start = Date.now();
      d.apply({col: col}).then(res.push([d.name, Date.now() - start]));
    }
    for (const i of res) if (i[1] > 10) violatingDetectors.push(`${i[0]}: ${i[1]}`);
    if (violatingDetectors.length > 0) throw new Error(violatingDetectors.join(', '));
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
