import {LoopProgress} from '../engines/engine-calls';

export type QueuedTraining<T> = {outcome: 'completed'; result: T} | {outcome: 'failed' | 'cancelled'; error: unknown} |
  {outcome: 'superseded'};

interface Request { isSuperseded: boolean }

/** Runs one training at a time. A newer request supersedes the one running (it stops before its next fit) and
 * any still waiting (they never start); a superseded training reports `superseded`, whatever it returned. */
export class TrainingQueue {
  private latest: Request | null = null;
  private previous: Promise<unknown> = Promise.resolve();

  /** Runs [train] once the trainings before it have settled. [indicator] is created when it starts, so a request
   * superseded while waiting has none; its progress is cancelled when the request is superseded or the indicator is
   * cancelled, and `update` goes to the indicator. */
  run<T>(train: (progress: LoopProgress) => Promise<T>, indicator?: () => LoopProgress): Promise<QueuedTraining<T>> {
    this.supersede();
    const request: Request = {isSuperseded: false};
    this.latest = request;
    const previous = this.previous;
    const settled = (async () => {
      await previous;
      try {
        return await this.execute(request, train, indicator);
      } finally {
        if (this.latest === request)
          this.latest = null;
      }
    })();
    this.previous = settled;
    return settled;
  }

  /** The data changed: the latest training ends as `superseded` before its next fit, without a newer request. */
  supersede(): void {
    if (this.latest !== null)
      this.latest.isSuperseded = true;
  }

  private async execute<T>(request: Request, train: (progress: LoopProgress) => Promise<T>,
    indicatorOf?: () => LoopProgress): Promise<QueuedTraining<T>> {
    if (request.isSuperseded)
      return {outcome: 'superseded'};
    const indicator = indicatorOf?.();
    const progress: LoopProgress = {
      get canceled() {
        return request.isSuperseded || (indicator?.canceled ?? false);
      },
      update: (percent, description) => indicator?.update(percent, description),
    };
    try {
      const result = await train(progress);
      return request.isSuperseded ? {outcome: 'superseded'} : {outcome: 'completed', result};
    } catch (error) {
      // The outcome is the queue's answer: the error is returned, not thrown.
      if (request.isSuperseded)
        return {outcome: 'superseded'};
      return {outcome: progress.canceled ? 'cancelled' : 'failed', error};
    }
  }
}
