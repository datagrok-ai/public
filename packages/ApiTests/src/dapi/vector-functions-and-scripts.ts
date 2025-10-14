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

  test('Vector func: cascade columns modification', async () => {
    const df = DG.DataFrame.fromCsv(`x, y, val
      1, 2, a
      4, 5, b
      7, 8, c`);

     
    await df.columns.addNewCalculated('newColScalar', '${val} + \'_scalar\'');
    await df.columns.addNewCalculated('newColVector', 'ApiTests:testVectorFunc(${newColScalar}, \'vector\')');
    await df.columns.addNewCalculated('newColScalar2', '${newColVector} + \'_scalar2\'');

    df.col('val')!.set(0, 'b', true);

    await awaitCheck(() => df.col('newColScalar')?.get(0) === 'b_scalar',
        `incorrect value in newColScalar after initial value changed, actual: ${df.col('newColScalar')?.get(0)}, expected: b_scalar`, 1000);
    await awaitCheck(() => df.col('newColVector')?.get(0) === 'vector_b_scalar',
        `incorrect value in newColVector after initial value changed, actual: ${df.col('newColVector')?.get(0)}, expected: vector_b_scalar`, 1000);
    await awaitCheck(() => df.col('newColScalar2')?.get(0) === 'vector_b_scalar_scalar2',
        `incorrect value in newColScalar2 after initial value changed, actual: ${df.col('newColScalar2')?.get(0)}, expected: vector_b_scalar_scalar2`, 1000);
  })

  test('Vector script: cascade columns modification', async () => {
    const df = DG.DataFrame.fromCsv(`x, y, val
      1, 2, a
      4, 5, b
      7, 8, c`);

     
    await df.columns.addNewCalculated('newColScalar', '${x} + 1');
    await df.columns.addNewCalculated('newColVector', 'ApiTests:vectorScript(${newColScalar}, 1)');
    await df.columns.addNewCalculated('newColScalar2', '${newColVector} + \'1\'');

    df.col('x')!.set(0, 2, true);

    await awaitCheck(() => df.col('newColScalar')?.get(0) === 3,
        `incorrect value in newColScalar after initial value changed, actual: ${df.col('newColScalar')?.get(0)}, expected: 3`, 1000);
    await awaitCheck(() => df.col('newColVector')?.get(0) === '3_1',
        `incorrect value in newColVector after initial value changed, actual: ${df.col('newColVector')?.get(0)}, expected: 3_1`, 1000);
    await awaitCheck(() => df.col('newColScalar2')?.get(0) === '3_11',
        `incorrect value in newColScalar2 after initial value changed, actual: ${df.col('newColScalar2')?.get(0)}, expected: 3_11`, 1000);
  })

}, {owner: 'mdolotova@datagrok.ai'});

async function testVectorFunc(formula: string, results: string[], colToChange: string, newVals: number[], resultsAfterValChange: string[]): Promise<void> {
    const df = DG.DataFrame.fromCsv(`x, y, val
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