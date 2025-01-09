import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, expectTable } from '@datagrok-libraries/utils/src/test';

const GDF = grok.dapi.functions;

category('Dapi: functions annotations', async () => {

  test('Join Df', async () => {
    await testAnnotation('ApiTests:testOutputAnnotationJoinDf', 'joined');
  });

  test('Join Column', async () => {
    await testAnnotation('ApiTests:testOutputAnnotationJoinCol', 'joined');
  });

  test('Join Column list', async () => {
    await testAnnotation('ApiTests:testOutputAnnotationJoinColList', 'joined');
  });

  test('Replace Df', async () => {
    await testAnnotation('ApiTests:testOutputAnnotationReplaceDf', 'val');
  });

  test('Replace Column', async () => {
    await testAnnotation('ApiTests:testOutputAnnotationReplaceCol', 'val');
  });

  test('Replace Column list', async () => {
    await testAnnotation('ApiTests:testOutputAnnotationReplaceColList', 'val');
  });

  //to check that in case we return column, parent dataframe is not created for it and stays null
  test('Return column without action', async () => {
    await testAnnotation('ApiTests:testOutputWithoutAction');
  });

}, { owner: 'aparamonov@datagrok.ai' });

async function testAnnotation(functionName: string, colName?: string): Promise<void> {
    var df = DG.DataFrame.fromCsv(`x, y, val
        1, 2, a
        4, 5, b
        7, 8, c`);
  
  // Run script and validate output
    const res = await grok.functions.call(functionName, {data: df, col: df.col('val')});
    var len = colName == 'joined' ? 4 : 3;
    expect(df.columns.length, len, "Incorrect number of columns in dataframe");
    if (!colName) {
      expect(res.dataFrame == null, true, "Parent dataframe shoul be null");
      return;
    }
    expect(df.col(colName!)!.get(0), 'a_abc', "Incorrect data in joined column");
    expect(df.col(colName!)!.get(1), 'b_abc', "Incorrect data in joined column");
    expect(df.col(colName!)!.get(2), 'c_abc', "Incorrect data in joined column");
  }