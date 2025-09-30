import {DrcAnalysis} from './drc/drc-analysis';
// import { DoseRatioAnalysis } from './dose-ratio/dose-ratio-analysis';
// import { QpcrAnalysis } from './qpcr/qpcr-analysis';
import {IPlateAnalysis} from './plate-analysis';

/**
 * A singleton manager to register and provide access to all available plate analyses.
 */
export class AnalysisManager {
  private static _instance: AnalysisManager;
  private _analyses: IPlateAnalysis[] = [];

  private constructor() {
    // --- Register all your analyses here ---
    this._analyses.push(new DrcAnalysis());
    // this._analyses.push(new DoseRatioAnalysis()); // Uncomment when refactored
    // this._analyses.push(new QpcrAnalysis());    // Uncomment when refactored
  }

  /** Gets the singleton instance of the manager. */
  public static get instance(): AnalysisManager {
    if (!AnalysisManager._instance)
      AnalysisManager._instance = new AnalysisManager();
    return AnalysisManager._instance;
  }

  /**
     * Initializes all registered analyses by registering their properties with the database.
     * This should be called once when your package starts.
     */
  public async init(): Promise<void> {
    for (const analysis of this._analyses) {
      if ('registerProperties' in analysis)
        await (analysis as any).registerProperties();
    }
  }

  /** Returns the list of all registered analysis instances. */
  public get analyses(): IPlateAnalysis[] {
    return this._analyses;
  }

  /** Finds a specific analysis by its friendly name (e.g., 'Dose Response'). */
  public byFriendlyName(name: string): IPlateAnalysis | undefined {
    return this._analyses.find((a) => a.friendlyName === name);
  }
}
