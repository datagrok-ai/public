// Bookkeeping for early-stopping logic shared by `MainExecutor` (sync loop)
// and `WorkerExecutor` (async pool dispatch). Owns the valid/above-threshold
// buckets and the count vs. `stopAfter`; each executor keeps its own loop
// control and uses `accept()` as the stop signal.
//
// `accept()` returns `true` once the caller should stop dispatching seeds.
// It stays correct under the worker-arm race where multiple in-flight seeds
// complete after the stop has been signaled — the `validPointsCount`
// guard short-circuits before pushing.

import {Extremum} from './optimizer-misc';
import {EarlyStoppingSettings} from './constants';

export class EarlyStopTracker {
  private readonly validExtremums: Extremum[] = [];
  private readonly allExtremums: Extremum[] = [];
  private validPointsCount = 0;
  private readonly useEarlyStopping: boolean;
  private readonly useAboveThresholdPoints: boolean;
  private readonly threshold: number | undefined;
  private readonly maxValidPoints: number;

  constructor(settings: EarlyStoppingSettings) {
    this.useEarlyStopping = settings.useEarlyStopping;
    this.useAboveThresholdPoints = settings.useAboveThresholdPoints;
    this.threshold = this.useEarlyStopping ? settings.costFuncThreshold : undefined;
    this.maxValidPoints = settings.stopAfter;
  }

  /** Bucket the extremum and return `true` once the caller should stop. */
  accept(extremum: Extremum): boolean {
    if (!this.useEarlyStopping) {
      this.validExtremums.push(extremum);
      return false;
    }
    if (this.validPointsCount >= this.maxValidPoints) return true;
    if (extremum.cost <= this.threshold!) {
      this.validExtremums.push(extremum);
      ++this.validPointsCount;
      return this.validPointsCount >= this.maxValidPoints;
    }
    this.allExtremums.push(extremum);
    return false;
  }

  /** Final extremum list — applies the above-threshold spillover rule. */
  finalize(): Extremum[] {
    if (this.useEarlyStopping && this.useAboveThresholdPoints &&
        this.maxValidPoints > this.validExtremums.length)
      return [...this.validExtremums, ...this.allExtremums];
    return this.validExtremums;
  }
}
