/* eslint-disable valid-jsdoc */
// sobol-sensitivity-analysis.ts

/* Tools that perform variance-based sensitivity analysis (VSA).

   WARNING. When computing Sobol' indices, negative values may occur.
            The reason is given in the discussion [3].

   References
     [1] A. Saltelli, et. al., "Variance based sensitivity analysis of model output. Design and estimator
         for the total sensitivity index", DOI: https://doi.org/10.1016/j.cpc.2009.09.018

     [2] Variance-based sensitivity analysis,
         LINK: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis

     [3] Negative Sobol indices, https://github.com/SALib/SALib/issues/102

     [4] Unexpected Sobol indices,
         www.researchgate.net/post/Is-it-ever-possible-to-have-the-sum-of-first-order-Sobol-indices-greater-than-one
*/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {VariedNumericalInputInfo, FixedInputItem, getVariedNumericalInputColumnsForSobolAnalysis} from './input-tools';
import {checkSize, getCalledFuncCalls} from './utils';
import {getOutput, OutputDataFromUI,
  SensitivityAnalysisResult, getDataFrameFromInputsOutputs} from './sa-outputs-routine';

import {DiffGrok} from '../fitting-view';
import {ModelEvaluator} from './diff-studio/model-evaluator';
import {DIFF_GROK_OUT_IDX} from './constants';

type VariedNumericalInputValues = VariedNumericalInputInfo & {column: DG.Column};

type SobolIndices = {
  firstOrder: DG.Column,
  totalOrder: DG.Column
};

const DEFAULT_VALUE_OF_SOBOL_INDEX = 0;

export type ResultOfSobolAnalysis = {
  firstOrderSobolIndices: DG.DataFrame,
  totalOrderSobolIndices: DG.DataFrame
} & SensitivityAnalysisResult;

/** Sobol sensitivity analysis */
export class SobolAnalysis {
  private samplesCount: number;
  private dimension: number;
  private runsCount: number;
  private fixedInputs: FixedInputItem[];

  private variedInputs: VariedNumericalInputValues[];
  private func: DG.Func;
  private funcCalls: DG.FuncCall[] | null = null;
  private funcInputs: any[];
  private diffGrok: DiffGrok | undefined;

  private outputsOfInterest: OutputDataFromUI[];

  constructor(
    func: DG.Func,
    fixedInputs: FixedInputItem[],
    variedInputs: VariedNumericalInputInfo[],
    outputsOfInterest: OutputDataFromUI[],
    samplesCount: number,
    diffGrok: DiffGrok | undefined,
  ) {
    // check size
    checkSize(samplesCount);

    this.samplesCount = samplesCount;
    this.fixedInputs = fixedInputs;

    this.func = func;

    /* Get values for the varied inputs:
         1) random data with the uniform on [0, 1] distribution is generated
            - the approach given in [1] is implemented;
         2) scale [0, 1]-to-[min, max] is applied.

       Further, these values are put to the function analyzed. */
    const numericalColumns = getVariedNumericalInputColumnsForSobolAnalysis(samplesCount, variedInputs);

    this.dimension = variedInputs.length;

    this.variedInputs = [...Array(this.dimension).keys()].map((i) =>
      ({prop: variedInputs[i].prop,
        min: variedInputs[i].min,
        max: variedInputs[i].max,
        column: numericalColumns[i]}),
    );

    // number of the function runs (this formula is given in [1])
    this.runsCount = samplesCount * (this.dimension + 2);

    this.funcCalls = (diffGrok === undefined) ? [] : null;
    this.funcInputs = [];

    // create an array of funccalls
    for (let i = 0; i < this.runsCount; ++i) {
      const inputs: any = {};

      for (const input of fixedInputs)
        inputs[input.name] = input.value;

      for (const input of this.variedInputs)
        inputs[input.prop.name] = input.column.get(i);

      if (this.funcCalls !== null)
        this.funcCalls.push(func.prepare(inputs));
      this.funcInputs.push(inputs);
    }

    this.diffGrok = diffGrok;

    this.outputsOfInterest = outputsOfInterest;
  } // constructor

  /** Returns 1-st and total order Sobol' indices. */
  private getSobolIndeces(outputColumn: DG.Column): SobolIndices {
    /* 1-st order and total order Sobol' indices are defined by
       the formulas (2) and (4) respectively [1]. Computations requires:
         - the variance V(Y);
         - the quantities V_i and E_i.
       The formulas (b),(f) (see Table 2 of [1]) are implemented in order
       to obtain V_i & E_i.

       Further, we use notations from the paper [1]. */

    const d = this.dimension;
    const N = this.samplesCount;

    const firstOrderIndeces = new Float32Array(d);
    const totalOrderIndeces = new Float32Array(d);

    // Compute variance, i.e. V(Y).
    const arr = outputColumn.getRawData();
    const len = 2 * N;

    let sum = 0;
    let sumOfSquares = 0;

    for (let i = 0; i < len; ++i) {
      sum += arr[i];
      sumOfSquares += arr[i] * arr[i];
    }

    const mean = sum / len;
    const variance = sumOfSquares / len - mean * mean;

    // Compute Sobol' indices

    // check variance: default values are used if variance is zero
    if (variance === 0) {
      firstOrderIndeces.fill(DEFAULT_VALUE_OF_SOBOL_INDEX);
      totalOrderIndeces.fill(DEFAULT_VALUE_OF_SOBOL_INDEX);
    } else {
      for (let i = 0; i < d; ++i) {
        let sumForFisrtOrderIndex = 0;
        let sumForTotalOrderIndex = 0;

        let buf = 0;

        for (let j = 0; j < N; ++j) {
          buf = arr[(2 + i) * N + j] - arr[j];
          sumForFisrtOrderIndex += arr[N + j] * buf;
          sumForTotalOrderIndex += buf * buf;
        }

        firstOrderIndeces[i] = sumForFisrtOrderIndex / (N * variance);
        totalOrderIndeces[i] = sumForTotalOrderIndex / (2 * N * variance);
      }
    }

    return {
      firstOrder: DG.Column.fromFloat32Array(outputColumn.name, firstOrderIndeces),
      totalOrder: DG.Column.fromFloat32Array(outputColumn.name, totalOrderIndeces),
    };
  }

  // Performs variance-based sensitivity analysis
  async perform(): Promise<ResultOfSobolAnalysis> {
    // columns with the varied inputs values
    const inputCols = this.variedInputs.map((varInput) => varInput.column as DG.Column);

    // columns with outputs
    let outputCols: DG.Column[];

    if (this.diffGrok !== undefined) {
      const evaluator = new ModelEvaluator(
        this.diffGrok,
        this.funcInputs,
        this.outputsOfInterest[DIFF_GROK_OUT_IDX],
      );

      outputCols = await evaluator.getResults();
    } else {
      if (this.funcCalls ===null)
        throw new Error('Failed Sobol analysis: empty list of funccalls');

      this.funcCalls = await getCalledFuncCalls(this.funcCalls);
      outputCols = getOutput(this.funcCalls, this.outputsOfInterest).columns.toList();
    }

    // create table with the varied inputs
    const funcEvalResults = getDataFrameFromInputsOutputs(inputCols, outputCols);
    funcEvalResults.name = `Sensitivity Analysis of ${this.func.friendlyName}`;

    // compute 1-st & total order Sobol' indices
    const sobolIndeces: SobolIndices[] = outputCols.map((col) => this.getSobolIndeces(col));

    // create dataframes with 1-st & total order Sobol' indices

    const inputNames = DG.Column.fromStrings('input', inputCols.map((col) => (col.name)));

    const firstOrderSobolIndecesCols: DG.Column[] = [inputNames];
    const totalOrderSobolIndecesCols: DG.Column[] = [inputNames];

    for (const item of sobolIndeces) {
      firstOrderSobolIndecesCols.push(item.firstOrder);
      totalOrderSobolIndecesCols.push(item.totalOrder);
    }

    const firstOrderSobolIndicesDF = DG.DataFrame.fromColumns(firstOrderSobolIndecesCols);
    firstOrderSobolIndicesDF.name = 'First order Sobol\' indices';

    const totalOrderSobolIndicesDF = DG.DataFrame.fromColumns(totalOrderSobolIndecesCols);
    totalOrderSobolIndicesDF.name = 'Total order Sobol\' indices';

    return {
      funcEvalResults: funcEvalResults,
      funcInputs: this.funcInputs,
      firstOrderSobolIndices: firstOrderSobolIndicesDF,
      totalOrderSobolIndices: totalOrderSobolIndicesDF,
    };
  }
};
