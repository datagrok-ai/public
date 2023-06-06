// sensitivityAnalysis.ts

/* Tools that perform variance-based sensitivity analysis (VSA).

   References
     [1] A. Saltelli, et. al., "Variance based sensitivity analysis of model output. Design and estimator
         for the total sensitivity index", DOI: https://doi.org/10.1016/j.cpc.2009.09.018

     [2] Variance-based sensitivity analysis,
         LINK: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis 
*/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {VariedNumericalInputInfo, FixedInputItem, getVariedNumericalInputColumns} from './inputTools';
import {getOutputColumns} from './outputTools';
import {checkSize} from './utils';

type VariedNumericalInputValues = VariedNumericalInputInfo & {column: DG.Column};

export class VarianceBasedSenstivityAnalysis {
  private samplesCount: number;
  private runsCount: number;
  private fixedInputs: FixedInputItem[];

  private variedInputs: VariedNumericalInputValues[];
  private func: DG.Func;
  private funcCalls: DG.FuncCall[];

  constructor(
    func: DG.Func,
    fixedInputs: FixedInputItem[],
    variedInputs: VariedNumericalInputInfo[],
    samplesCount: number,
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
    const numericalColumns = getVariedNumericalInputColumns(samplesCount, variedInputs);

    const dimension = variedInputs.length;

    this.variedInputs = [...Array(dimension).keys()].map((i) =>
      ({prop: variedInputs[i].prop,
        min: variedInputs[i].min,
        max: variedInputs[i].max,
        column: numericalColumns[i]}),
    );

    // number of the function runs (this formula is given in [1])
    this.runsCount = samplesCount * (dimension + 2);

    this.funcCalls = [];

    // create an array of funccalls
    for (let i = 0; i < this.runsCount; ++i) {
      const inputs: any = {};

      for (const input of fixedInputs)
        inputs[input.name] = input.value;

      for (const input of this.variedInputs)
        inputs[input.prop.name] = input.column.get(i);

      this.funcCalls.push(func.prepare(inputs));
    }
  }

  // Runs the function with each inputs set
  private async run(): Promise<void> {
    await Promise.all(this.funcCalls.map((call) => call.call()));
  }

  // Performs variance-based sensitivity analysis
  async perform(): Promise<DG.DataFrame> {
    await this.run();

    // columns with the varied inputs values
    const columns = this.variedInputs.map((varInput) => varInput.column as DG.Column);

    // create table with the varied inputs
    const table = DG.DataFrame.fromColumns(columns);
    table.name = `Sensitivity Analysis of ${this.func.friendlyName}`;

    // add columns with outputs
    for (const col of getOutputColumns(this.funcCalls))
      table.columns.add(col);

    // TODO (Viktor Makarichev): add Sobol' indeces computation.

    return table;
  }
};
