import type {
  ObjectiveFunction,
  OptimizationResult,
  CommonSettings,
  IterationCallback,
  IterationState,
} from './types';
import {applyPenalty} from './penalty';

/**
 * Abstract base for every optimisation algorithm.
 *
 * Subclasses only need to implement `runInternal` and `withDefaults`.
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
  /*  Public API                                                         */
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
    // Negate the objective first, then apply penalty — penalty always increases
    // the value, so the solver is correctly pushed away from infeasible regions.
    const negated: ObjectiveFunction = (x) => -fn(x);
    const effective = this.applyConstraints(negated, settings);
    const result = this.runInternal(effective, x0, this.withDefaults(settings));

    return {
      ...result,
      value: -result.value,
      costHistory: result.costHistory.map((v) => -v),
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

  /** Wrap `fn` with penalty if constraints are provided in settings. */
  private applyConstraints(fn: ObjectiveFunction, settings: S): ObjectiveFunction {
    if (!settings.constraints || settings.constraints.length === 0)
      return fn;
    return applyPenalty(fn, settings.constraints, settings.penaltyOptions);
  }
}
