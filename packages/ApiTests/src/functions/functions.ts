import './date-functions';
import './math-functions';
import './text-functions';
import './logical-functions';
import './conversion-functions';
import './stats-functions';
import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('Functions: General', () => {
  test('Eval', async () => {
    let dfList: DG.DataFrame[] = await grok.functions.eval('OpenServerFile("System:AppData/ApiTests/datasets/demog.csv")');
    expect(dfList[0].columns instanceof DG.ColumnList, true);
  });
  test('Call', async () => {
    let dfList: DG.DataFrame[] = await grok.functions.call('OpenServerFile', {'fullPath': 'System:AppData/ApiTests/datasets/demog.csv'});
    expect(dfList[0].columns instanceof DG.ColumnList, true);
  });
  test('Def param', async () => {
      await grok.functions.call('AddNewColumn', {table: grok.data.demo.demog(), expression: 'test', name: 'test'});
  });
});
