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
import {OutputInfo, getOutput, SensitivityAnalysisResult, getDataFrameFromInputsOutputs} from './sa-outputs-routine';

type VariedNumericalInputValues = VariedNumericalInputInfo & {column: DG.Column};

type OutputDataFromUI = {
  prop: DG.Property,
  value: {
    row: number,
    columns: string | null
  },
};

const DEFAULT_VALUE_OF_SOBOL_INDEX = 0;

export class RandomAnalysis {
  private dimension: number;

  private variedInputs: VariedNumericalInputValues[];
  private func: DG.Func;
  private funcCalls: DG.FuncCall[];

  private outputInfo: OutputInfo[];

  constructor(
    func: DG.Func,
    fixedInputs: FixedInputItem[],
    variedInputs: VariedNumericalInputInfo[],
    outputsOfInterest: OutputDataFromUI[],
    samplesCount: number,
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

    this.funcCalls = [];

    // create an array of funccalls
    for (let i = 0; i < samplesCount; ++i) {
      const inputs: any = {};

      for (const input of fixedInputs)
        inputs[input.name] = input.value;

      for (const input of this.variedInputs)
        inputs[input.prop.name] = input.column.get(i);

      this.funcCalls.push(func.prepare(inputs));
    }

    this.outputInfo = outputsOfInterest.map((output) => ({
      prop: output.prop,
      elements: [],
      row: output.value.row,
    }
    ));
  }

  // Runs the function with each inputs set
  private async run(): Promise<void> {
    await Promise.all(this.funcCalls.map((call) => call.call()));
  }

  // Performs variance-based sensitivity analysis
  async perform(): Promise<SensitivityAnalysisResult> {
    //await this.run();

    this.funcCalls = await getCalledFuncCalls(this.funcCalls);

    // columns with the varied inputs values
    const inputCols = this.variedInputs.map((varInput) => varInput.column as DG.Column);

    // columns with outputs
    const outputCols = getOutput(this.funcCalls, this.outputInfo).columns.toList();    

    // create table with the varied inputs
    const funcEvalResults = getDataFrameFromInputsOutputs(inputCols, outputCols);

    funcEvalResults.name = `Sensitivity Analysis of ${this.func.friendlyName}`; 

    return {funcEvalResults: funcEvalResults, funcCalls: this.funcCalls};
  }
};
