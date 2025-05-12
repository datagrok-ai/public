// random-sensitivity-analysis.ts

/* Random or Monte Carlo sensitivity analysis:
   the function is evaluated with respect to random variation
   of the selected inputs within the specified ranges.
*/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {VariedNumericalInputInfo, FixedInputItem, getVariedNumericalInputColumnsForRandomAnalysis} from './input-tools';
import {checkSize, getCalledFuncCalls} from './utils';
import {OutputDataFromUI, getOutput,
  SensitivityAnalysisResult, getDataFrameFromInputsOutputs} from './sa-outputs-routine';

import {DiffGrok} from '../fitting-view';
import {ModelEvaluator} from './diff-studio/model-evaluator';
import {DIFF_GROK_OUT_IDX} from './constants';

type VariedNumericalInputValues = VariedNumericalInputInfo & {column: DG.Column};

/** Monte-Carlo analyzer */
export class RandomAnalysis {
  private dimension: number;

  private variedInputs: VariedNumericalInputValues[];
  private func: DG.Func;
  private funcCalls: DG.FuncCall[] | null = null;
  private outputsOfInterest: OutputDataFromUI[];
  private funcInputs: any[];
  private diffGrok: DiffGrok | undefined;

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
    this.func = func;
    const numericalColumns = getVariedNumericalInputColumnsForRandomAnalysis(samplesCount, variedInputs);
    this.dimension = variedInputs.length;

    this.variedInputs = [...Array(this.dimension).keys()].map((i) =>
      ({prop: variedInputs[i].prop,
        min: variedInputs[i].min,
        max: variedInputs[i].max,
        column: numericalColumns[i]}),
    );

    this.funcInputs = [];
    this.diffGrok = diffGrok;
    this.funcCalls = (diffGrok === undefined) ? [] : null;

    // create an array of funccalls
    for (let i = 0; i < samplesCount; ++i) {
      const inputs: any = {};

      for (const input of fixedInputs)
        inputs[input.name] = input.value;

      for (const input of this.variedInputs)
        inputs[input.prop.name] = input.column.get(i);

      if (this.funcCalls !== null)
        this.funcCalls.push(func.prepare(inputs));

      this.funcInputs.push(inputs);
    }

    this.outputsOfInterest = outputsOfInterest;
  } // constructor

  /** Performs variance-based sensitivity analysis */
  async perform(): Promise<SensitivityAnalysisResult> {
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
      if (this.funcCalls === null)
        throw new Error('Failed Monte-Carlo analysis: empty list of funccalls');

      this.funcCalls = await getCalledFuncCalls(this.funcCalls);
      outputCols = getOutput(this.funcCalls, this.outputsOfInterest).columns.toList();
    }

    // create table with the varied inputs
    const funcEvalResults = getDataFrameFromInputsOutputs(inputCols, outputCols);

    funcEvalResults.name = `Sensitivity Analysis of ${this.func.friendlyName}`;

    return {funcEvalResults: funcEvalResults, funcInputs: this.funcInputs};
  }
};
