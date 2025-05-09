import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DiffGrok} from '../../fitting-view';
import {OutputDataFromUI} from '../sa-outputs-routine';

import * as DGL from '@datagrok/diff-grok';
import {DEFAULT_NUM} from './defs';
import {getAnalysis} from './sens-analysis-utils';

export class ModelEvaluator {
  private diffGrok: DiffGrok;
  private funcInputs: Record<string, number>[];
  private rowIdx: number = DEFAULT_NUM;
  private argColIdx: number = DEFAULT_NUM;
  private argVal: number = DEFAULT_NUM;
  private outputNames: string[];

  constructor(
    diffGrok: DiffGrok,
    funcInputs: Record<string, number>[],
    outputsOfInterest: OutputDataFromUI) {
    this.diffGrok = diffGrok;
    this.funcInputs = funcInputs;
    this.outputNames = DGL.getOutputNames(this.diffGrok.ivp);

    if (outputsOfInterest.value.row !== null)
      this.rowIdx = outputsOfInterest.value.row;
    else {
      const name = outputsOfInterest.value.colName;
      this.argColIdx = this.outputNames.indexOf(name);

      if (this.argColIdx < 0)
        throw new Error(`Incorrect arg column name: ${name}`);

      this.argVal = outputsOfInterest.value.colValue;
    }
  }

  public async getResults(): Promise<DG.Column[]> {
    const samplesCount = this.funcInputs.length;

    const rawArrs: Float64Array[] = [];

    for (let i = 0; i < this.outputNames.length; ++i)
      rawArrs.push(new Float64Array(samplesCount));

    for (let rowIdx = 0; rowIdx < samplesCount; ++rowIdx) {
      const inputVec = DGL.getInputVector(
        this.funcInputs[rowIdx],
        this.diffGrok.ivp,
      );

      const rowVals = getAnalysis({
        ivp2ww: this.diffGrok.ivpWW,
        pipeline: this.diffGrok.pipelineCreator.getPipeline(inputVec),
        inputVec: inputVec,
        rowIdx: this.rowIdx,
        argColIdx: this.argColIdx,
        argVal: this.argVal,
      });

      for (let colIdx = 0; colIdx < rowVals.length; ++colIdx)
        rawArrs[colIdx][rowIdx] = rowVals[colIdx];
    }

    return this.outputNames.map((name, idx) => DG.Column.fromFloat64Array(
      name,
      rawArrs[idx],
    ));
  }
};
