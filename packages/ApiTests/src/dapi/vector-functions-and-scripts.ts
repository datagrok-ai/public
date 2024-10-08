import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, expectTable, awaitCheck } from '@datagrok-libraries/utils/src/test';

const GDF = grok.dapi.functions;

category('Dapi: vector scripts and functions', async () => {

  test('Vector script', async () => {
    await testVectorFunc('ApiTests:vectorScript(${x}, ${y})', ['1_2', '4_5', '7_8'], 'x', [3, 6, 9], ['3_2', '6_5', '9_8']);
  });

  test('Vector script with non-vectorizable param', async () => {
    await testVectorFunc('ApiTests:vectorScriptNonVectorizableParam(${x}, ${y}, ${val})',
        ['3_a', '9_a', '15_a'], 'x', [3, 6, 9], ['5_a', '11_a', '17_a']);
  }, {skipReason: 'GROK-16773'});

  test('Vector function', async () => {
    await testVectorFunc('ApiTests:testVectorFunc(${x}, \'abc\')',
        ['abc_1', 'abc_4', 'abc_7'], 'x', [3, 6, 9], ['abc_3', 'abc_6', 'abc_9']);
  });

  test('Vector function with non-vectorizable param', async () => {
    await testVectorFunc('ApiTests:testVectorFuncNonVectorizableParam(${x}, \'abc\', ${val})',
        ['abc_1_a', 'abc_4_b', 'abc_7_c'], 'x', [3, 6, 9], ['abc_3_a', 'abc_6_b', 'abc_9_c']);
  });

});

async function testVectorFunc(formula: string, results: string[], colToChange: string, newVals: number[], resultsAfterValChange: string[]): Promise<void> {
    var df = DG.DataFrame.fromCsv(`x, y, val
        1, 2, a
        4, 5, b
        7, 8, c`);
  
    const newColName = 'new';
    await df.columns.addNewCalculated(newColName, formula);
    //check the results in calculated column
    for (let i = 0; i < df.rowCount; i++) {
        await awaitCheck(() => df.col(newColName)?.get(i) === results[i],
            `incorrect value in ${i} row, actual: ${df.col(newColName)?.get(i)}, expected: ${results[i]}`, 1000);
    }
    //check updates after value in initial col changed
    for (let i = 0; i < df.rowCount; i++) {
        df.col(colToChange)!.set(i, newVals[i], true);
        await awaitCheck(() => df.col(newColName)?.get(i) === resultsAfterValChange[i],
            `incorrect value in ${i} row after initial value changed, actual: ${df.col(newColName)?.get(i)}, expected: ${resultsAfterValChange[i]}`, 1000);
    }
  }