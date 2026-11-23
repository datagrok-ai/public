import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DEFAULT_BUDGET_SECONDS} from '../training/search';

export const TIME = 'Time';
export const TRIALS = 'Trials';
export type LimitKind = typeof TIME | typeof TRIALS;
const MINUTES_RANGE = [1, 60];
const TRIALS_RANGE = [1, 10000];
const LIMIT_TOOLTIP = 'What ends the search: a time or a trial limit';
const MINUTES_TOOLTIP = 'Minutes the search may run';
const TRIALS_TOOLTIP = 'Maximum number of trials, pruned and failed ones included';

/** The one limit of a search as remembered per user. */
export interface StoredLimit { limit: LimitKind; minutes: number; trials: number }

/** The limit rows of a Search options grid (Tune's and Find best model's): **Limit** (Time / Trials) and under it the
 * one input it shows, **Max time** in whole minutes with a clicker or **Max trials**; each a platform int input with
 * its range (`min`, `max`, not nullable) in a grid cell, its caption in the grid's label column. The platform checks
 * a value and shows a problem (a red underline, the message in the input's tooltip under its description); only the
 * shown input counts. */
export class LimitInputs {
  readonly kind: DG.InputBase<string | null>;
  readonly minutes: DG.InputBase<number | null>;
  readonly trials: DG.InputBase<number | null>;
  /** The grid's label and cell of each row. */
  private readonly rows: Record<'limit' | LimitKind, HTMLElement[]>;

  /** [stored] gives the values, [defaults] what an invalid one is remembered as; [onChange] follows every change;
   * [trialsProblem] is the caller's own problem of a trial count in range, a validator of Max trials. */
  constructor(stored: StoredLimit, private readonly defaults: StoredLimit, onChange: () => void,
    trialsProblem: (trials: number) => string | null = () => null) {
    this.kind = ui.input.radio('Limit', {items: [TIME, TRIALS], value: stored.limit, nullable: false,
      onValueChanged: () => {
        this.showKind();
        onChange();
      }});
    this.kind.root.classList.add('forge-inline-radio');
    this.minutes = ui.input.int('Max time', {value: stored.minutes, min: 1, max: 60, nullable: false, step: 1,
      showPlusMinus: true, tooltipText: MINUTES_TOOLTIP, onValueChanged: () => onChange()});
    this.minutes.addPostfix('min');
    this.trials = ui.input.int('Max trials', {value: stored.trials, min: 1, max: 10000, nullable: false,
      showPlusMinus: false, tooltipText: TRIALS_TOOLTIP, onValueChanged: () => onChange()});
    this.trials.addValidator(() => {
      const trials = LimitInputs.valid(this.trials.value, TRIALS_RANGE);
      return trials === undefined ? null : trialsProblem(trials);
    });
    LimitInputs.bindOptionTooltips(this.kind, () => LIMIT_TOOLTIP);
    const rowOf = (input: DG.InputBase, tooltip: string) => {
      const {label, cell} = LimitInputs.row(input, tooltip);
      return [label, cell];
    };
    this.rows = {limit: rowOf(this.kind, LIMIT_TOOLTIP), [TIME]: rowOf(this.minutes, MINUTES_TOOLTIP),
      [TRIALS]: rowOf(this.trials, TRIALS_TOOLTIP)};
    this.showKind();
  }

  /** The grid's label and cell of each row, in order: Limit, Max time, Max trials. */
  get elements(): HTMLElement[] {
    return [...this.rows.limit, ...this.rows[TIME], ...this.rows[TRIALS]];
  }

  /** The platform inputs. */
  get inputs(): DG.InputBase[] {
    return [this.kind, this.minutes, this.trials];
  }

  get isTime(): boolean {
    return this.kind.value !== TRIALS;
  }

  /** The minutes in seconds; the default while they are invalid. */
  get seconds(): number {
    const minutes = LimitInputs.valid(this.minutes.value, MINUTES_RANGE);
    return minutes === undefined ? DEFAULT_BUDGET_SECONDS : minutes * 60;
  }

  /** The search's one limit, as the shown input has it: the time in seconds, or the trials. Read while it is valid. */
  limit(): {budgetSeconds: number} | {maxTrials: number} {
    return this.isTime ? {budgetSeconds: this.seconds} : {maxTrials: this.trials.value!};
  }

  /** The limit to remember: an invalid value gives way to the default. */
  stored(): StoredLimit {
    return {limit: this.isTime ? TIME : TRIALS,
      minutes: LimitInputs.valid(this.minutes.value, MINUTES_RANGE) ?? this.defaults.minutes,
      trials: LimitInputs.valid(this.trials.value, TRIALS_RANGE) ?? this.defaults.trials};
  }

  /** The shown input's first validation message, captioned (`Max time: ...`); null while it is valid. */
  problem(): string | null {
    const input = this.shown;
    input.validate();
    const [message] = input.validationResults;
    return message === undefined ? null : `${input.caption}: ${message}`;
  }

  /** Validates both inputs, so each shows its state; true while the shown one is valid. */
  validate(): boolean {
    (this.isTime ? this.trials : this.minutes).validate();
    return this.shown.validate();
  }

  setLocked(isLocked: boolean): void {
    for (const input of this.inputs)
      input.enabled = !isLocked;
  }

  private get shown(): DG.InputBase<number | null> {
    return this.isTime ? this.minutes : this.trials;
  }

  /** Shows the row of the chosen kind, hides the other; the host validates it with the change. */
  private showKind(): void {
    for (const kind of [TIME, TRIALS] as const) {
      for (const element of this.rows[kind])
        ui.setDisplay(element, (kind === TIME) === this.isTime);
    }
  }

  /** The stored limit of a `grok.userSettings` record [value]: its kind, minutes and trials, a malformed one falling
   * back to [defaults]; a record of the old time select (`seconds`) is a time limit of its whole minutes, rounded up
   * and clamped to 1..60. */
  static load(value: {[key: string]: unknown}, defaults: StoredLimit): StoredLimit {
    const {limit, minutes, trials, seconds} = value;
    const [low, high] = MINUTES_RANGE;
    const isOldTime = typeof seconds === 'number' && seconds > 0;
    const migrated = isOldTime ? Math.min(high, Math.max(low, Math.ceil(seconds / 60))) : defaults.minutes;
    return {limit: limit === TIME || limit === TRIALS ? limit : isOldTime ? TIME : defaults.limit,
      minutes: LimitInputs.valid(minutes, MINUTES_RANGE) ?? migrated,
      trials: LimitInputs.valid(trials, TRIALS_RANGE) ?? defaults.trials};
  }

  /** The grid's label and cell of [input], whose own label is hidden, the label with [tooltip]. */
  static row(input: DG.InputBase, tooltip: string): {label: HTMLElement; cell: HTMLElement} {
    input.root.classList.add('forge-tune-input');
    const label = ui.div([input.caption], 'forge-tune-label');
    ui.tooltip.bind(label, tooltip);
    return {label, cell: ui.div([input.root], 'forge-tune-cell')};
  }

  /** Binds to each option of the radio [input] the tooltip [tooltipOf] gives it; the options. */
  static bindOptionTooltips(input: DG.InputBase, tooltipOf: (option: HTMLElement) => string | null): HTMLElement[] {
    const options = Array.from(input.root.querySelectorAll<HTMLElement>('.ui-radio-button'));
    for (const option of options)
      ui.tooltip.bind(option, tooltipOf(option));
    return options;
  }

  /** [n] when it is a whole number in [low, high]; undefined otherwise. */
  private static valid(n: unknown, [low, high]: number[]): number | undefined {
    return typeof n === 'number' && Number.isInteger(n) && n >= low && n <= high ? n : undefined;
  }
}
