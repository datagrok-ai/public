import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
  OptimizationResult,
  CommonSettings,
  IterationCallback,
  IterationState,
} from './types';
import {applyPenalty, applyPenaltyAsync} from './penalty';

/**
 * Abstract base for every optimisation algorithm.
 *
 * Subclasses implement:
 *   - `runInternal`       — synchronous hot path
 *   - `runInternalAsync`  — asynchronous path (API calls, simulations, etc.)
 *   - `withDefaults`      — fill in algorithm-specific defaults
 *
 * The base class owns input validation, the iteration-callback
 * mechanism, and constraint handling so that every optimizer
 * behaves consistently.
 */
export abstract class Optimizer<S extends CommonSettings = CommonSettings> {
  readonly name: string;

  protected constructor(name: string) {
    this.name = name;
  }

  /* ------------------------------------------------------------------ */
  /*  Public API — synchronous                                           */
  /* ------------------------------------------------------------------ */

  /** Find the point that minimises `fn`. */
  minimize(
    fn: ObjectiveFunction,
    x0: Float64Array,
    settings: S,
  ): OptimizationResult {
    this.validate(x0);
    const effective = this.applyConstraints(fn, settings);
    return this.runInternal(effective, x0, this.withDefaults(settings));
  }

  /** Find the point that maximises `fn`. */
  maximize(
    fn: ObjectiveFunction,
    x0: Float64Array,
    settings: S,
  ): OptimizationResult {
    this.validate(x0);
    const negated: ObjectiveFunction = (x) => -fn(x);
    const effective = this.applyConstraints(negated, settings);
    const result = this.runInternal(effective, x0, this.withDefaults(settings));

    const ch = result.costHistory;
    for (let i = 0; i < ch.length; i++) ch[i] = -ch[i];

    return {
      ...result,
      value: -result.value,
      costHistory: ch,
    };
  }

  /* ------------------------------------------------------------------ */
  /*  Public API — asynchronous                                          */
  /* ------------------------------------------------------------------ */

  /** Find the point that minimises async `fn`. */
  async minimizeAsync(
    fn: AsyncObjectiveFunction,
    x0: Float64Array,
    settings: S,
  ): Promise<OptimizationResult> {
    this.validate(x0);
    const effective = this.applyConstraintsAsync(fn, settings);
    return this.runInternalAsync(effective, x0, this.withDefaults(settings));
  }

  /** Find the point that maximises async `fn`. */
  async maximizeAsync(
    fn: AsyncObjectiveFunction,
    x0: Float64Array,
    settings: S,
  ): Promise<OptimizationResult> {
    this.validate(x0);
    const negated: AsyncObjectiveFunction = async (x) => -(await fn(x));
    const effective = this.applyConstraintsAsync(negated, settings);
    const result = await this.runInternalAsync(effective, x0, this.withDefaults(settings));

    const ch = result.costHistory;
    for (let i = 0; i < ch.length; i++) ch[i] = -ch[i];

    return {
      ...result,
      value: -result.value,
      costHistory: ch,
    };
  }

  /* ------------------------------------------------------------------ */
  /*  To be implemented by each algorithm                                */
  /* ------------------------------------------------------------------ */

  protected abstract runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    settings: S,
  ): OptimizationResult;

  protected abstract runInternalAsync(
    fn: AsyncObjectiveFunction,
    x0: Float64Array,
    settings: S,
  ): Promise<OptimizationResult>;

  /** Return settings with algorithm-specific defaults filled in. */
  protected abstract withDefaults(settings: S): S;

  /* ------------------------------------------------------------------ */
  /*  Shared helpers available to all optimizers                         */
  /* ------------------------------------------------------------------ */

  /** Fire the callback; returns `true` when the optimizer should stop. */
  protected notify(cb: IterationCallback | undefined, state: IterationState): boolean {
    return cb?.(state) === true;
  }

  /* ------------------------------------------------------------------ */
  /*  Validation & constraints                                           */
  /* ------------------------------------------------------------------ */

  private validate(x0: Float64Array): void {
    if (x0.length === 0)
      throw new Error(`${this.name}: x0 must have at least one element`);
  }

  /** Wrap sync `fn` with penalty if constraints are provided in settings. */
  private applyConstraints(fn: ObjectiveFunction, settings: S): ObjectiveFunction {
    if (!settings.constraints || settings.constraints.length === 0)
      return fn;
    return applyPenalty(fn, settings.constraints, settings.penaltyOptions);
  }

  /** Wrap async `fn` with penalty if constraints are provided in settings. */
  private applyConstraintsAsync(fn: AsyncObjectiveFunction, settings: S): AsyncObjectiveFunction {
    if (!settings.constraints || settings.constraints.length === 0)
      return fn;
    return applyPenaltyAsync(fn, settings.constraints, settings.penaltyOptions);
  }
}
