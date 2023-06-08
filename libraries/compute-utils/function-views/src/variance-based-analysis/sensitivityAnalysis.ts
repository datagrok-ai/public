// sensitivityAnalysis.ts

/* Tools that perform variance-based sensitivity analysis (VSA).

   WARNING. When computing Sobol' indeces, negative values may occur. 
            The reason is given in the discussion [3].

   References
     [1] A. Saltelli, et. al., "Variance based sensitivity analysis of model output. Design and estimator
         for the total sensitivity index", DOI: https://doi.org/10.1016/j.cpc.2009.09.018

     [2] Variance-based sensitivity analysis,
         LINK: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis 

     [3] Negative Sobol indices, https://github.com/SALib/SALib/issues/102
*/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {VariedNumericalInputInfo, FixedInputItem, getVariedNumericalInputColumns} from './inputTools';
import {getOutputColumns} from './outputTools';
import {checkSize} from './utils';

type VariedNumericalInputValues = VariedNumericalInputInfo & {column: DG.Column};

type SobolIndeces = {
  firstOrder: DG.Column,
  totalOrder: DG.Column
};

export type ResultOfVarianceBasedSenstivityAnalysis = {
  funcEvalResults: DG.DataFrame,
  firstOrderSobolIndeces: DG.DataFrame,
  totalOrderSobolIndeces: DG.DataFrame
};

export class VarianceBasedSenstivityAnalysis {
  private samplesCount: number;
  private dimension: number;
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

    this.dimension = variedInputs.length;

    this.variedInputs = [...Array(this.dimension).keys()].map((i) =>
      ({prop: variedInputs[i].prop,
        min: variedInputs[i].min,
        max: variedInputs[i].max,
        column: numericalColumns[i]}),
    );

    // number of the function runs (this formula is given in [1])
    this.runsCount = samplesCount * (this.dimension + 2);

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

  // Returns 1-st and totoal order Sobol' indeces.
  private getSobolIndeces(outputColumn: DG.Column): SobolIndeces {

    /* 1-st order and total order Sobol' indeces are defined by
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

    // TODO (Viktor Makarichev): the following expressions can be reduced.

    const mean = sum / len;
    const variance = sumOfSquares / len - mean * mean;   
    
    // TODO (Viktor Makarichev): add processing the case variance = 0 or its small values
   
    // Compute Sobol' indeces 
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

    return {
      firstOrder: DG.Column.fromFloat32Array(outputColumn.name, firstOrderIndeces), 
      totalOrder: DG.Column.fromFloat32Array(outputColumn.name, totalOrderIndeces)
    };
  }

  // Performs variance-based sensitivity analysis
  async perform(): Promise<ResultOfVarianceBasedSenstivityAnalysis> {
    await this.run();

    // columns with the varied inputs values
    const inputColumns = this.variedInputs.map((varInput) => varInput.column as DG.Column);

    // create table with the varied inputs
    const funcEvalResults = DG.DataFrame.fromColumns(inputColumns);
    funcEvalResults.name = `Sensitivity Analysis of ${this.func.friendlyName}`;

    const outputColumns = getOutputColumns(this.funcCalls);

    // add columns with outputs
    for (const col of outputColumns)
      funcEvalResults.columns.add(col);

    // compute 1-st & total order Sobol' indeces
    const sobolIndeces: SobolIndeces[] = outputColumns.map((col) => this.getSobolIndeces(col));

    /*for (const el of sobolIndeces) {
      console.log(el.firstOrder);
      console.log(el.totalOrder);
    }*/

    // create dataframes with 1-st & total order Sobol' indeces

    const inputNames = DG.Column.fromStrings('input', inputColumns.map((col) => (col.name)));

    const firstOrderSobolIndecesCols: DG.Column[] = [inputNames];
    const totalOrderSobolIndecesCols: DG.Column[] = [inputNames];

    for (const item of sobolIndeces) {
      firstOrderSobolIndecesCols.push(item.firstOrder);
      totalOrderSobolIndecesCols.push(item.totalOrder);
    }

    const firstOrderSobolIndecesDF = DG.DataFrame.fromColumns(firstOrderSobolIndecesCols);
    firstOrderSobolIndecesDF.name = "First order Sobol' indeces";    

    const totalOrderSobolIndecesDF = DG.DataFrame.fromColumns(totalOrderSobolIndecesCols);
    totalOrderSobolIndecesDF.name = "Total order Sobol' indeces";

    return {
      funcEvalResults: funcEvalResults,
      firstOrderSobolIndeces: firstOrderSobolIndecesDF,
      totalOrderSobolIndeces: totalOrderSobolIndecesDF
    };
  }
};
