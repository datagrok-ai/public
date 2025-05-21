// random-sensitivity-analysis.ts

/* Random or Monte Carlo sensitivity analysis:
   the function is evaluated with respect to random variation
   of the selected inputs within the specified ranges.
*/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getCalledFuncCalls} from './utils';
import {OutputDataFromUI, getOutput, getDataFrameFromInputsOutputs} from './sa-outputs-routine';

import {DiffGrok} from '../fitting-view';
import {ModelEvaluator} from './diff-studio/model-evaluator';
import {DIFF_GROK_OUT_IDX} from './constants';

/** Grid analyzer */
export class GridAnalysis {
  private func: DG.Func;
  private funcCalls: DG.FuncCall[] | null = null;
  private outputsOfInterest: OutputDataFromUI[];
  private diffGrok: DiffGrok | undefined;
  private funcInputs: Record<string, any>[];
  private inputCols: DG.Column[]= [];

  constructor(
    func: DG.Func,
    funcInputs: Record<string, any>[],
    variedInputsInf: {name: string, caption: string, type: DG.TYPE}[],
    outputsOfInterest: OutputDataFromUI[],
    diffGrok: DiffGrok | undefined,
  ) {
    this.func = func;
    this.diffGrok = diffGrok;
    this.funcInputs = funcInputs;
    this.funcCalls = (diffGrok === undefined) ?
      funcInputs.map((inp) => func.prepare(inp)) :
      null;

    variedInputsInf.forEach((info) => {
      const col = DG.Column.fromList(
        info.type as unknown as DG.COLUMN_TYPE,
        info.caption,
        funcInputs.map((inputs) => inputs[info.name]),
      );
      this.inputCols.push(col);
    });

    this.outputsOfInterest = outputsOfInterest;
  } // constructor

  /** Performs non-random sensitivity analysis */
  async perform(): Promise<DG.DataFrame> {
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
        throw new Error('Failed Grid analysis: empty list of funccalls');

      this.funcCalls = await getCalledFuncCalls(this.funcCalls);
      outputCols = getOutput(this.funcCalls, this.outputsOfInterest).columns.toList();
    }

    // create table with the varied inputs
    const funcEvalResults = getDataFrameFromInputsOutputs(this.inputCols, outputCols);

    funcEvalResults.name = `Sensitivity Analysis of ${this.func.friendlyName}`;

    return funcEvalResults;
  }
};
