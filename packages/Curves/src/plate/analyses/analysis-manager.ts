import {DoseRatioAnalysis} from './dose-ratio/dose-ratio-analysis';
import {DrcAnalysis} from './drc/drc-analysis';
import {IPlateAnalysis} from './base-analysis';
import {QpcrAnalysis} from './qpcr/qpcr-analysis';

export class AnalysisManager {
  private static _instance: AnalysisManager;
  private _analyses: IPlateAnalysis[] = [];

  private constructor() {
    this._analyses.push(new DrcAnalysis());
    this._analyses.push(new DoseRatioAnalysis());
    this._analyses.push(new QpcrAnalysis());
    //and you add new analyses here the same way
  }

  public static get instance(): AnalysisManager {
    if (!AnalysisManager._instance)
      AnalysisManager._instance = new AnalysisManager();
    return AnalysisManager._instance;
  }

  public async init(): Promise<void> {
    for (const analysis of this._analyses) {
      if ('registerProperties' in analysis)
        await (analysis as any).registerProperties();
    }
  }

  public get analyses(): IPlateAnalysis[] {
    return this._analyses;
  }
  public byFriendlyName(name: string): IPlateAnalysis | undefined {
    return this._analyses.find((a) => a.friendlyName === name);
  }
}
