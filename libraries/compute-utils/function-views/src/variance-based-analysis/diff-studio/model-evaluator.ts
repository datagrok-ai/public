import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DiffGrok} from '../../fitting-view';
import {OutputDataFromUI, SensitivityAnalysisResult} from '../sa-outputs-routine';

export class ModelEvaluator {
  private diffGrok: DiffGrok;
  private funcInputs: Record<string, number>[];
  private outputsOfInterest: OutputDataFromUI;

  constructor(
    diffGrok: DiffGrok,
    funcInputs: Record<string, number>[],
    outputsOfInterest: OutputDataFromUI) {
    this.diffGrok = diffGrok;
    this.funcInputs = funcInputs;
    this.outputsOfInterest = outputsOfInterest;
  }

  public async getResults(): Promise<SensitivityAnalysisResult> {
    return {
      funcEvalResults: grok.data.demo.randomWalk(this.funcInputs.length),
      funcInputs: this.funcInputs,
    };
  }
};
