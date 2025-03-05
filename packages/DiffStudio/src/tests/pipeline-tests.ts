// Tests for computational pipelines

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test, timeout} from '@datagrok-libraries/utils/src/test';

import * as DSL from '@datagrok/diff-grok';

import {TEMPLATES, ENERGY_N_CONTROL} from '../templates';
import {USE_CASES} from '../use-cases';
import {DF_NAME} from '../constants';

const TIMEOUT = 10000;
const MIN_ROWS = 1;
const TINY = 0.1;

const getDeviation = (solutionTable: DG.DataFrame, solutionArrs: Float64Array[], hasStageCol: boolean): number => {
  let deviation = 0;
  const rowCount = solutionTable.rowCount;
  const cols = solutionTable.columns;
  const colsCount = cols.length - (hasStageCol ? 1 : 0);

  if (colsCount !== solutionArrs.length)
    throw new Error('Non-equal solution columns count');

  for (let i = 0; i < colsCount; ++i) {
    const raw = cols.byIndex(i).getRawData();
    const arr = solutionArrs[i];

    if (arr.length !== rowCount)
      throw new Error('Non-equal solution rows count');

    for (let j = 0; j < rowCount; ++j)
      deviation = Math.max(deviation, Math.abs(raw[j] - arr[j]));
  }

  return deviation;
};

/** Template for testing solving ODEs */
const testPipelineTemplate = (modelName: string, problem: string) => {
  test(modelName, async () => {
    // Main thread features
    let ivp: DSL.IVP | null = null;
    let code: string | null = null;
    let script: DG.Script | null = null;
    let params: Record<string, number> | null = null;
    let call: DG.FuncCall | null = null;
    let solutionDf: DG.DataFrame | null = null;

    // Worker features
    let ivpWW: DSL.IVP2WebWorker | null = null;
    let inputVector: Float64Array | null = null;
    let creator: DSL.PipelineCreator | null = null;
    let pipeline: DSL.Pipeline | null = null;
    let solutionArrs: Float64Array[] | null = null;

    // Comparison
    let sameColsCount = false;
    let sameRowsCount = false;
    let arrsSameLength = true;
    let deviation = TINY;

    try {
      // Apply computations via the package funcs
      ivp = DSL.getIVP(problem);
      code = DSL.getScriptLines(ivp).join('\n');
      script = DG.Script.create(code);
      params = DSL.getScriptParams(ivp);
      call = script.prepare(params);
      await call.call();
      solutionDf = call.outputs[DF_NAME];
      const solutionCols = solutionDf.columns;

      // Remove segments column
      const hasStageCol = (ivp.updates !== null);

      // Apply computations via the pipeline features
      ivpWW = DSL.getIvp2WebWorker(ivp);
      inputVector = DSL.getInputVector(params, ivp);
      creator = DSL.getPipelineCreator(ivp);
      pipeline = creator.getPipeline(inputVector);
      solutionArrs = DSL.applyPipeline(pipeline, ivpWW, inputVector);

      // Check length
      const arrRowCount = solutionArrs[0].length;
      const arrColCount = solutionArrs.length;

      for (let i = 1; i < arrColCount; ++i)
        arrsSameLength = arrsSameLength && (solutionArrs[i].length === arrRowCount);

      // Compare results
      sameColsCount = ((solutionCols.length - (hasStageCol ? 1 : 0)) === arrColCount);
      sameRowsCount = (solutionDf.rowCount === arrRowCount);

      if (arrsSameLength && sameColsCount && sameRowsCount)
        deviation = getDeviation(solutionDf, solutionArrs, hasStageCol);
    } catch (e) {}

    expect(ivp !== null, true, `Failed to parse equations: ${modelName}`);
    expect(code !== null, true, `Failed to create JS-code: ${modelName}`);
    expect(script !== null, true, `Failed to create the platform script: ${modelName}`);
    expect(params !== null, true, `Failed to get script params: ${modelName}`);
    expect(call !== null, true, `Failed to get funccall: ${modelName}`);
    expect(solutionDf !== null, true, `Failed to apply numerical method: ${modelName}`);
    expect(solutionDf.rowCount >= MIN_ROWS, true, `Solving failed: an empty dataframe: ${modelName}`);

    expect(ivpWW !== null, true, `Failed to create IVP-object for in-webworkers run: ${modelName}`);
    expect(inputVector !== null, true, `Failed to get input vector for in-webworkers run: ${modelName}`);
    expect(creator !== null, true, `Failed to build pipeline creator: ${modelName}`);
    expect(pipeline!== null, true, `Failed to create a pipeline: ${modelName}`);
    expect(solutionArrs !== null, true, `Failed to apply a pipeline: ${modelName}`);

    expect(arrsSameLength, true, `Non-equal lengths of solution arrays: ${modelName}`);
    expect(sameColsCount, true, `Wrong solution arrays count: ${modelName}.`);
    expect(sameRowsCount, true, `Wrong solution arrays length: ${modelName}`);
    expect(deviation < TINY, true, `Too large deviation between solution via package and pipelines: ${modelName}.
    Obtained deviation: ${deviation}. Expected: < ${TINY}`);
  }, {timeout: TIMEOUT, benchmark: true});
}; // testPipelineTemplate

// Correctness tests
category('Pipelines', () => {
  testPipelineTemplate('Basic template', TEMPLATES.BASIC);
  testPipelineTemplate('Advanced template', TEMPLATES.ADVANCED);
  testPipelineTemplate('Extended template', TEMPLATES.EXTENDED);
  testPipelineTemplate('Chem react', USE_CASES.CHEM_REACT);
  testPipelineTemplate('Robertson', USE_CASES.ROBERTSON);
  testPipelineTemplate('Fermentation', USE_CASES.FERMENTATION);
  testPipelineTemplate('Multistage model', USE_CASES.ACID_PROD);
  testPipelineTemplate('Nimotuzumab', USE_CASES.NIMOTUZUMAB);
  testPipelineTemplate('Bioreactor', USE_CASES.BIOREACTOR);
  testPipelineTemplate('Pollution', USE_CASES.POLLUTION);
  testPipelineTemplate('Energy-n-control', ENERGY_N_CONTROL);
}); // Correctness
