import * as DG from 'datagrok-api/dg';
import {Engine, Hyperparameters} from '../engines/engine';
import {EnginePrediction, LoopProgress, yieldToEventLoop} from '../engines/engine-calls';
import {errorMessage, ForgeError} from '../forge-error';
import {classesOf, ClassProbabilities, METRIC_LABELS, metricsOf, taskOf} from '../metrics/metrics';
import {isScoreMetric} from '../metrics/report';
import {replayPostprocessing} from '../preparation/pipeline';
import {validationFolds} from './k-fold';
import {ExecutorState, Halving, MAX_DUPLICATES, MAX_FAILED_WITHOUT_COMPLETE, PRUNE_AFTER_FOLDS, PRUNE_MIN_TRIALS,
  SearchExecutor, SearchOutcome, SearchPlan, SearchSpace, SearchStop, searchTitle, Trial, trialPoints} from './search';
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
/** The shortest time between two progress updates. */
const REPORT_MS = 250;

/** The search in the browser: the trials' fold fits through the method's `train` and `apply`, up to the plan's
 * parallelism at a time, each on its own clones of the fold rows. */
export class BrowserSearchExecutor implements SearchExecutor {
  async run(plan: SearchPlan, progress: LoopProgress, previous?: SearchOutcome): Promise<SearchOutcome> {
    const started = performance.now();
    const elapsed = () => performance.now() - started;
    const resume = previous?.resume;
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
    const trials: Trial[] = [...previous?.trials ?? []];
    const before = trials.length;
    const isComplete = (t: Trial) => t.state === 'complete';
    let completed = trials.filter(isComplete).length;
    let settled = before;
    // By candidate and its own trial number: the first fold scores once they pass unpruned, else null; a trial in
    // flight also has the promise of its first folds' outcome.
    const passedScores = plan.candidates.map((_, k) => [...resume?.passed[k] ?? []]);
    const early = plan.candidates.map((): Promise<number[] | undefined>[] => []);
    const inFlight = plan.candidates.map(() => new Set<number>());
    let launched = 0;
    let failures = 0;

    // Each candidate's points in turn, a continued one's from where it stopped; a point already tried is skipped.
    const states = plan.candidates.map((c, k) => {
      const space = c.space ?? [];
      const points = trialPoints(c, plan.seed, task);
      const drawn = resume?.drawn[k] ?? 0;
      for (let i = 0; i < drawn; i++)
        points.next();
      const own = trials.filter((t) => t.candidate === k);
      return {points, space, drawn, tried: new Set(own.map((t) => keyOf(space, t.hyperparameters))), next: own.length,
        isExhausted: resume?.exhausted[k] ?? false};
    });
    const drawOf = (k: number): Hyperparameters | null => {
      const state = states[k];
      for (let duplicates = 0; !state.isExhausted && duplicates < MAX_DUPLICATES; duplicates++) {
        const next = state.points.next();
        if (next.done)
          break;
        state.drawn++;
        const key = keyOf(state.space, next.value);
        if (!state.tried.has(key)) {
          state.tried.add(key);
          return next.value;
        }
      }
      state.isExhausted = true;
      return null;
    };

    // Per candidate: its complete trials over all runs, its failed ones in this run.
    const ownCompleted = states.map((_, k) => trials.filter((t) => t.candidate === k && isComplete(t)).length);
    const ownFailures = states.map(() => 0);
    // An equal split (Find best model): the candidates not exhausted, in plan (cost) order, one after another; each
    // gets its share of the trials or the time left when it starts (rounded up: the remainder to the first ones), so
    // what a candidate leaves unused goes to the ones after it. A share ends at its count or time, at exhaustion, or
    // after MAX_FAILED_WITHOUT_COMPLETE failures without a complete trial; the next starts once its trials settled.
    const isEqual = plan.strategy === 'equal';
    const order = isEqual ? states.map((_, k) => k).filter((k) => !states[k].isExhausted) : [];
    const {budgetMs, maxTrials} = plan;
    const quotaOf = (at: number) => maxTrials === undefined || at >= order.length ? Infinity :
      Math.ceil((maxTrials - launched) / (order.length - at));
    const shareEnd = (at: number) => budgetMs === undefined || at >= order.length ? Infinity :
      elapsed() + (budgetMs - elapsed()) / (order.length - at);
    const equal = {at: 0, given: states.map(() => 0), running: states.map(() => 0), endsMs: shareEnd(0),
      quota: quotaOf(0)};
    const isShareOver = (k: number) => equal.given[k] >= equal.quota || elapsed() >= equal.endsMs ||
      (ownCompleted[k] === 0 && ownFailures[k] >= MAX_FAILED_WITHOUT_COMPLETE);

    // One candidate (Tune) draws in turn. Several (Find best model): an equal split, or successive halving, then the
    // finalists in turn.
    const halving: Halving | undefined = states.length === 1 || isEqual ? undefined : resume?.halving === undefined ?
      {round: 0, pool: states.map((_, k) => k), given: states.map(() => 0), eliminated: states.map(() => null),
        next: 0} :
      {...resume.halving, pool: [...resume.halving.pool], given: [...resume.halving.given],
        eliminated: [...resume.halving.eliminated], finalists: resume.halving.finalists?.slice()};
    let pending = 0;
    const bestScoreOf = (k: number): number | undefined => {
      let score: number | undefined;
      for (const t of trials) {
        if (t !== undefined && t.candidate === k && isComplete(t) &&
          (score === undefined || direction * (t.score! - score) > 0))
          score = t.score;
      }
      return score;
    };
    // The pool with a space, not exhausted and with a complete trial, best first; the worse half is eliminated, but
    // never down to fewer than topK; topK or fewer become the finalists.
    const closeRound = (h: Halving) => {
      const best = new Map(h.pool.map((k) => [k, bestScoreOf(k)]));
      const ranked = h.pool.filter((k) => states[k].space.length > 0 && !states[k].isExhausted &&
        best.get(k) !== undefined).sort((a, b) => direction * (best.get(b)! - best.get(a)!) || a - b);
      const keep = Math.min(ranked.length, Math.max(plan.topK, Math.ceil(ranked.length / 2)));
      for (const k of ranked.slice(keep))
        h.eliminated[k] = h.round;
      h.pool = ranked.slice(0, keep);
      h.given = states.map(() => 0);
      h.round++;
      if (h.pool.length <= plan.topK)
        h.finalists = [...h.pool];
    };
    const nextTrial = (): {candidate: number; hyperparameters: Hyperparameters} | null => {
      while (isEqual && equal.at < order.length) {
        const k = order[equal.at];
        if (!isShareOver(k)) {
          if (ownCompleted[k] === 0 && ownFailures[k] + equal.running[k] >= MAX_FAILED_WITHOUT_COMPLETE)
            return null;
          const hyperparameters = drawOf(k);
          if (hyperparameters !== null) {
            equal.given[k]++;
            equal.running[k]++;
            return {candidate: k, hyperparameters};
          }
        }
        if (equal.running[k] > 0)
          return null;
        equal.at++;
        equal.endsMs = shareEnd(equal.at);
        equal.quota = quotaOf(equal.at);
      }
      if (isEqual)
        return null;
      if (halving === undefined) {
        const hyperparameters = drawOf(0);
        return hyperparameters === null ? null : {candidate: 0, hyperparameters};
      }
      while (halving.finalists === undefined) {
        const quota = 2 ** halving.round;
        let k: number | undefined;
        while ((k = halving.pool.filter((c) => !states[c].isExhausted && halving.given[c] < quota)
          .sort((a, b) => halving.given[a] - halving.given[b])[0]) !== undefined) {
          const hyperparameters = drawOf(k);
          if (hyperparameters !== null) {
            halving.given[k]++;
            pending++;
            return {candidate: k, hyperparameters};
          }
        }
        // A trial of the round in flight: wait for it, never exhaustion.
        if (pending > 0)
          return null;
        closeRound(halving);
      }
      const cycle = halving.finalists;
      while (cycle.length > 0) {
        halving.next %= cycle.length;
        const k = cycle[halving.next];
        const hyperparameters = drawOf(k);
        if (hyperparameters !== null) {
          halving.next++;
          return {candidate: k, hyperparameters};
        }
        cycle.splice(halving.next, 1);
      }
      return null;
    };

    const action = searchTitle(plan.mode, plan.candidates[0].engine.name);
    // At most one update per REPORT_MS: a report within it is sent when it ends; [isFinal] sends at once.
    let reportedAt = -Infinity;
    let trailing: ReturnType<typeof setTimeout> | undefined;
    const report = (isFinal = false) => {
      clearTimeout(trailing);
      trailing = undefined;
      const wait = reportedAt + REPORT_MS - performance.now();
      if (!isFinal && wait > 0) {
        trailing = setTimeout(() => report(), wait);
        return;
      }
      reportedAt = performance.now();
      const finished = settled - before;
      const running = launched - finished;
      const seconds = budgetMs === undefined ? 0 : Math.max(0, Math.ceil((budgetMs - elapsed()) / 1000));
      // The time is up and the last trial settled: `finishing 1 trial` stays until the search ends.
      if (budgetMs !== undefined && seconds === 0 && running === 0)
        return;
      const limit = budgetMs === undefined ? `${finished}/${maxTrials}` :
        seconds === 0 && running > 0 ? `finishing ${running} ${running === 1 ? 'trial' : 'trials'}` :
          seconds >= 60 ? `${Math.ceil(seconds / 60)} min left` : `${seconds} s left`;
      const share = budgetMs === undefined ? finished / maxTrials! : elapsed() / budgetMs;
      progress.update(Math.min(99, 100 * share), `${action} (${limit})`);
    };

    const needsProbabilities = (engine: Engine) =>
      `The metric ${label} needs class probabilities; ${engine.name} gives none.`;
    const isScore = isScoreMetric(metric);
    const scoreOf = (actual: DG.Column, labels: DG.Column, probabilities?: ClassProbabilities) =>
      metricsOf(task, actual, labels, positiveClass, isScore ? probabilities : undefined)[metric];
    const foldScore = (prediction: EnginePrediction, i: number, engine: Engine): number => {
      const score = scoreOf(heldOut[i], replayPostprocessing(prediction.prediction, options, prediction.probabilities),
        prediction.probabilities);
      if (score === undefined)
        throw new ForgeError(needsProbabilities(engine));
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

    const runTrial = async (index: number, candidate: number, hyperparameters: Hyperparameters): Promise<void> => {
      const {engine} = plan.candidates[candidate];
      const number = states[candidate].next++;
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
          scores[i] = foldScore(results[i], i, engine);
        } catch (e) {
          isFailed = true;
          throw e;
        } finally {
          semaphore.release();
        }
      };
      const settle = (trial: Omit<Trial, 'candidate' | 'number' | 'hyperparameters' | 'fitMs'>) => {
        trials[index] = {candidate, number, hyperparameters, fitMs, ...trial};
        settled++;
        if (trial.state === 'complete') {
          completed++;
          ownCompleted[candidate]++;
        } else if (trial.state === 'failed') {
          failures++;
          ownFailures[candidate]++;
        }
      };
      let pass: (scores?: number[]) => void = () => {};
      early[candidate][number] = new Promise((resolve) => pass = resolve);
      passedScores[candidate][number] = null;
      inFlight[candidate].add(number);
      const first = Array.from({length: Math.min(PRUNE_AFTER_FOLDS, folds)}, (_, i) => i);
      const rest = Array.from({length: folds - first.length}, (_, i) => first.length + i);
      // A trial is judged only when folds are left after the first ones (never on a hold out), against its own
      // candidate's earlier trials; a candidate's first trials never are.
      const canPrune = rest.length > 0 && number >= PRUNE_MIN_TRIALS;
      try {
        await fitWave(first);
        // The earlier trials still in flight pass or not first; the settled ones' scores are kept.
        if (canPrune)
          await Promise.all([...inFlight[candidate]].filter((j) => j < number).map((j) => early[candidate][j]));
        const passed = canPrune ?
          passedScores[candidate].slice(0, number).filter((s): s is number[] => Array.isArray(s)) : [];
        const forecast = rest.length === 0 ? 0 :
          fitMs / first.length * rest.length / Math.min(plan.parallelism, rest.length);
        // Abandoned when it would end after the time: an equal split's after its share, its candidate's first
        // complete trial never.
        const deadline = isEqual ? equal.endsMs : budgetMs ?? Infinity;
        if ((isEqual ? ownCompleted[candidate] : completed) > 0 && elapsed() + forecast > deadline)
          settle({state: 'abandoned', foldScores: [...scores]});
        else if (isPruned(passed, scores))
          settle({state: 'pruned', foldScores: [...scores]});
        else {
          passedScores[candidate][number] = [...scores];
          pass([...scores]);
          await fitWave(rest);
          const outOfFold = outOfFoldColumn(results, fold, position);
          const probabilities = results[0].probabilities === undefined ? undefined :
            outOfFoldProbabilities(results, results[0].probabilities, fold, position);
          const score = scoreOf(target, replayPostprocessing(outOfFold, options, probabilities), probabilities);
          if (score === undefined)
            throw new ForgeError(needsProbabilities(engine));
          settle({state: 'complete', score, foldScores: [...scores]});
        }
      } catch (e) {
        const foldScores = scores.filter(() => true);
        settle(e === CANCELLED ? {state: 'abandoned', foldScores} : {state: 'failed', foldScores,
          error: errorMessage(e)});
      } finally {
        pass();
        inFlight[candidate].delete(number);
      }
      report();
    };

    const running = new Set<Promise<void>>();
    let stopped: SearchStop;
    while (true) {
      if (progress.canceled) {
        stopped = 'cancelled';
        break;
      }
      if (budgetMs !== undefined && elapsed() >= budgetMs) {
        stopped = 'budget';
        break;
      }
      if (maxTrials !== undefined && launched >= maxTrials) {
        stopped = 'trials';
        break;
      }
      // An equal split ends a failing candidate's share instead.
      const isFailing = !isEqual && completed === 0;
      if (isFailing && failures >= MAX_FAILED_WITHOUT_COMPLETE) {
        stopped = 'failed';
        break;
      }
      const isFull = running.size >= plan.parallelism ||
        (isFailing && failures + running.size >= MAX_FAILED_WITHOUT_COMPLETE);
      // A null while a round or a share is unsettled always has a trial of it in flight: a wait, not exhaustion.
      const next = isFull ? null : nextTrial();
      if (next === null) {
        if (running.size === 0) {
          stopped = isEqual && completed === 0 && failures > 0 ? 'failed' : 'exhausted';
          break;
        }
        await Promise.race(running);
        continue;
      }
      const index = before + launched++;
      const round = halving?.finalists === undefined ? halving?.round : undefined;
      const p: Promise<void> = runTrial(index, next.candidate, next.hyperparameters).finally(() => {
        running.delete(p);
        if (isEqual)
          equal.running[next.candidate]--;
        if (round === undefined || round !== halving!.round)
          return;
        pending--;
        if (trials[index].state === 'abandoned')
          halving!.given[next.candidate]--;
      });
      running.add(p);
      report();
    }
    await Promise.all(running);
    report(true);
    const resumeState: ExecutorState = {drawn: states.map((s) => s.drawn), exhausted: states.map((s) => s.isExhausted),
      passed: passedScores, ...(halving === undefined ? {} : {halving})};
    return {trials, elapsedMs: (previous?.elapsedMs ?? 0) + elapsed(), stopped, resume: resumeState};
  }
}

/** The values of [point] in the dimensions of [space], as text: points with the same key are one. */
function keyOf(space: SearchSpace, point: Hyperparameters): string {
  return JSON.stringify(space.map((d) => point[d.name]));
}

export const browserSearch = new BrowserSearchExecutor();
