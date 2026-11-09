import * as DG from 'datagrok-api/dg';
import {Hyperparameters} from '../engines/engine';
import {EnginePrediction, LoopProgress, yieldToEventLoop} from '../engines/engine-calls';
import {errorMessage, ForgeError} from '../forge-error';
import {classesOf, ClassProbabilities, METRIC_LABELS, MetricId, metricsOf, taskOf} from '../metrics/metrics';
import {SCORE_METRIC_IDS} from '../metrics/report';
import {replayPostprocessing} from '../preparation/pipeline';
import {validationFolds} from './k-fold';
import {MAX_DUPLICATES, MAX_FAILED_WITHOUT_COMPLETE, PRUNE_AFTER_FOLDS, PRUNE_MIN_TRIALS, SearchExecutor,
  SearchOutcome, SearchPlan, SearchStop, Trial, trialPoints} from './search';
import {foldMasks, foldPrediction, outOfFoldColumn, outOfFoldProbabilities, positionsOf} from './train-model';

/** At most [slots] holders at a time; the others wait in turn. */
class Semaphore {
  private waiters: (() => void)[] = [];

  constructor(private slots: number) {}

  async acquire(): Promise<void> {
    if (this.slots > 0)
      this.slots--;
    else
      await new Promise<void>((resolve) => this.waiters.push(resolve));
  }

  release(): void {
    const next = this.waiters.shift();
    if (next === undefined)
      this.slots++;
    else
      next();
  }
}

const CANCELLED = new ForgeError('The search was cancelled.');

/** The search in the browser: the trials' fold fits through the method's `train` and `apply`, up to the plan's
 * parallelism at a time, each on its own clones of the fold rows. */
export class BrowserSearchExecutor implements SearchExecutor {
  async run(plan: SearchPlan, progress: LoopProgress): Promise<SearchOutcome> {
    const started = performance.now();
    const elapsed = () => performance.now() - started;
    const {request, fold, metric, direction} = plan;
    const {columns, target, options} = request.prepared;
    const task = taskOf(target);
    const classes = task === 'classification' ? classesOf(request.target) : [];
    const positiveClass = classes.length === 2 ? classes[0] : undefined;
    const label = METRIC_LABELS[metric];
    const folds = validationFolds(request.validation);
    const masks = Array.from({length: folds}, (_, i) => foldMasks(fold, i));
    const heldOut = masks.map(({isHeldOut}) => target.clone(isHeldOut));
    const position = positionsOf(fold);
    const semaphore = new Semaphore(plan.parallelism);
    const trials: Trial[] = [];
    let completed = 0;
    let settled = 0;
    let best: Trial | undefined;
    // Each trial's first fold scores once it passes them unpruned, else undefined.
    const early: Promise<number[] | undefined>[] = [];
    let launched = 0;
    let failures = 0;

    // Tune: the one candidate's points in turn; a point already tried is skipped.
    const candidate = plan.candidates[0];
    const points = trialPoints(candidate, plan.seed, task);
    const tried = new Set<string>();
    let isExhausted = false;
    const nextTrial = (): Hyperparameters | null => {
      for (let duplicates = 0; !isExhausted && duplicates < MAX_DUPLICATES; duplicates++) {
        const next = points.next();
        if (next.done)
          break;
        const key = JSON.stringify((candidate.space ?? []).map((d) => next.value[d.name]));
        if (!tried.has(key)) {
          tried.add(key);
          return next.value;
        }
      }
      isExhausted = true;
      return null;
    };

    const needsProbabilities = `The metric ${label} needs class probabilities; ${candidate.engine.name} gives none.`;
    const isScoreMetric = (SCORE_METRIC_IDS as readonly MetricId[]).includes(metric);
    const scoreOf = (actual: DG.Column, labels: DG.Column, probabilities?: ClassProbabilities) =>
      metricsOf(task, actual, labels, positiveClass, isScoreMetric ? probabilities : undefined)[metric];
    const foldScore = (prediction: EnginePrediction, i: number): number => {
      const score = scoreOf(heldOut[i], replayPostprocessing(prediction.prediction, options, prediction.probabilities),
        prediction.probabilities);
      if (score === undefined)
        throw new ForgeError(needsProbabilities);
      return score;
    };

    // Only the lower-numbered trials that passed their first folds count, so the finish order changes nothing.
    const isPruned = (passed: number[][], scores: number[]): boolean => {
      if (passed.length < PRUNE_MIN_TRIALS)
        return false;
      const meanOf = (s: number[]) => direction * s.slice(0, PRUNE_AFTER_FOLDS).reduce((a, b) => a + b, 0) /
        PRUNE_AFTER_FOLDS;
      const means = passed.map(meanOf).sort((a, b) => a - b);
      const middle = means.length >> 1;
      const median = means.length % 2 === 1 ? means[middle] : (means[middle - 1] + means[middle]) / 2;
      return meanOf(scores) < median;
    };

    const runTrial = async (index: number, hyperparameters: Hyperparameters): Promise<void> => {
      const {engine} = candidate;
      const scores: number[] = [];
      const results: EnginePrediction[] = [];
      let fitMs = 0;
      // A failed fold fails the trial: its folds still waiting for a slot skip their fits.
      let isFailed = false;
      const fitWave = async (folds: number[]) => {
        const outcomes = await Promise.allSettled(folds.map(fitFold));
        const failed = outcomes.find((o): o is PromiseRejectedResult => o.status === 'rejected');
        if (failed !== undefined)
          throw failed.reason;
      };
      const fitFold = async (i: number) => {
        await semaphore.acquire();
        try {
          await yieldToEventLoop();
          if (progress.canceled)
            throw CANCELLED;
          if (isFailed)
            return;
          const t0 = performance.now();
          results[i] = await foldPrediction(engine, columns, target, hyperparameters, masks[i]);
          fitMs += performance.now() - t0;
          scores[i] = foldScore(results[i], i);
        } catch (e) {
          isFailed = true;
          throw e;
        } finally {
          semaphore.release();
        }
      };
      const settle = (trial: Omit<Trial, 'candidate' | 'number' | 'hyperparameters' | 'fitMs'>) => {
        trials[index] = {candidate: 0, number: index, hyperparameters, fitMs, ...trial};
        settled++;
        if (trial.state === 'complete') {
          completed++;
          if (best === undefined || direction * (trial.score! - best.score!) > 0)
            best = trials[index];
        } else if (trial.state === 'failed')
          failures++;
      };
      let pass: (scores?: number[]) => void = () => {};
      early[index] = new Promise((resolve) => pass = resolve);
      const first = Array.from({length: Math.min(PRUNE_AFTER_FOLDS, folds)}, (_, i) => i);
      const rest = Array.from({length: folds - first.length}, (_, i) => first.length + i);
      // A trial is judged only when folds are left after the first ones: never on a hold out.
      const canPrune = rest.length > 0 && index >= PRUNE_MIN_TRIALS;
      try {
        await fitWave(first);
        const passed = canPrune ?
          (await Promise.all(early.slice(0, index))).filter((s): s is number[] => s !== undefined) : [];
        const forecast = rest.length === 0 ? 0 :
          fitMs / first.length * rest.length / Math.min(plan.parallelism, rest.length);
        if (completed > 0 && elapsed() + forecast > plan.budgetMs)
          settle({state: 'abandoned', foldScores: [...scores]});
        else if (isPruned(passed, scores))
          settle({state: 'pruned', foldScores: [...scores]});
        else {
          pass([...scores]);
          await fitWave(rest);
          const outOfFold = outOfFoldColumn(results, fold, position);
          const probabilities = results[0].probabilities === undefined ? undefined :
            outOfFoldProbabilities(results, results[0].probabilities, fold, position);
          const score = scoreOf(target, replayPostprocessing(outOfFold, options, probabilities), probabilities);
          if (score === undefined)
            throw new ForgeError(needsProbabilities);
          settle({state: 'complete', score, foldScores: [...scores]});
        }
      } catch (e) {
        const foldScores = scores.filter(() => true);
        settle(e === CANCELLED ? {state: 'abandoned', foldScores} : {state: 'failed', foldScores,
          error: errorMessage(e)});
      } finally {
        pass();
      }
      progress.update(Math.min(99, 100 * elapsed() / plan.budgetMs), best === undefined ? `${settled} trials` :
        `${settled} trials, best ${label} ${best.score!.toFixed(3)} (${engine.name})`);
    };

    const running = new Set<Promise<void>>();
    let stopped: SearchStop;
    while (true) {
      if (progress.canceled) {
        stopped = 'cancelled';
        break;
      }
      if (elapsed() >= plan.budgetMs) {
        stopped = 'budget';
        break;
      }
      if (plan.maxTrials !== undefined && completed >= plan.maxTrials) {
        stopped = 'trials';
        break;
      }
      if (completed === 0 && failures >= MAX_FAILED_WITHOUT_COMPLETE) {
        stopped = 'failed';
        break;
      }
      const isFull = running.size >= plan.parallelism ||
        (plan.maxTrials !== undefined && completed + running.size >= plan.maxTrials) ||
        (completed === 0 && failures + running.size >= MAX_FAILED_WITHOUT_COMPLETE);
      const hyperparameters = isFull ? null : nextTrial();
      if (hyperparameters === null) {
        if (running.size === 0) {
          stopped = 'exhausted';
          break;
        }
        await Promise.race(running);
        continue;
      }
      const p: Promise<void> = runTrial(launched++, hyperparameters).finally(() => running.delete(p));
      running.add(p);
    }
    await Promise.all(running);
    return {trials, elapsedMs: elapsed(), stopped};
  }
}

export const browserSearch = new BrowserSearchExecutor();
