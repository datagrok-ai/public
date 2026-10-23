import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {missingColumnsOf, MissingValuesSettings} from '../preparation/missing-values';
import {cutoffProblem, DEFAULT_CUTOFF, DEFAULT_STEPS, defaultOneHot, hasUniqueCategories, PreparationSteps,
  twoClasses} from '../preparation/pipeline';
import {CollapsibleGroup} from './collapsible-group';
import {MissingValuesInputs} from './missing-values-inputs';

/** The **Preparation** group of the Train view: Missing values, One-hot encoding, Skip unique categories, Predict
 * probability and its Positive class cutoff, each shown only while it applies to the selection; the group is hidden
 * while none does. An input that hides goes back to its default, so it reappears as it first appeared; One-hot
 * encoding gets its default when it appears and follows it until the user toggles it. */
export class PreparationInputs {
  readonly group: CollapsibleGroup;
  readonly missingValues: MissingValuesInputs;
  readonly oneHot: DG.InputBase<boolean>;
  readonly skipUniqueCategories: DG.InputBase<boolean>;
  readonly predictProbability: DG.InputBase<boolean>;
  readonly cutoff: DG.InputBase<number | null>;
  private features: DG.Column[] = [];
  private hasTwoClasses = false;
  private oneHotDefault = false;
  private isOneHotChosen = false;

  /** [onChanged] follows every input but the cutoff, which calls [onCutoffChanged]; [missingValuesProblem] is the
   * validator of **Missing values**. */
  constructor(onChanged: () => void, onCutoffChanged: () => void, missingValuesProblem: () => string | null) {
    this.missingValues = new MissingValuesInputs(onChanged, missingValuesProblem);
    this.oneHot = ui.input.bool('One-hot encoding', {value: false,
      tooltipText: 'Lets the model use categories as numbers without implying an order: each category gets its own ' +
        '0/1 column.',
      onValueChanged: (value) => {
        if (value !== this.oneHotDefault)
          this.isOneHotChosen = true;
        onChanged();
      }});
    const onStepChanged = () => {
      this.updateVisibility();
      onChanged();
    };
    this.skipUniqueCategories = ui.input.bool('Skip unique categories', {value: DEFAULT_STEPS.skipUniqueCategories,
      tooltipText: 'Leaves out text features whose values are all different, such as ids.',
      onValueChanged: onStepChanged});
    this.predictProbability = ui.input.bool('Predict probability', {value: DEFAULT_STEPS.predictProbability,
      onValueChanged: onStepChanged});
    this.cutoff = ui.input.float('Positive class cutoff', {value: DEFAULT_CUTOFF, min: 0, max: 1, step: 0.01,
      showSlider: true, nullable: false, onValueChanged: onCutoffChanged});
    this.cutoff.addValidator(() => this.cutoffProblem);
    for (const input of this.stepInputs)
      ui.setDisplay(input.root, false);
    this.group = new CollapsibleGroup('Preparation', this.inputs.map((input) => input.root));
    this.updateVisibility();
  }

  get inputs(): DG.InputBase[] {
    return [...this.missingValues.inputs, ...this.stepInputs];
  }

  get cutoffValue(): number {
    return this.cutoff.value ?? DEFAULT_CUTOFF;
  }

  /** Why the cutoff cannot be used; null while it is a number from 0 to 1. */
  get cutoffProblem(): string | null {
    return cutoffProblem(this.cutoff.value);
  }

  /** The steps as chosen; a hidden input's step is off. */
  steps(): PreparationSteps {
    return {oneHot: this.oneHot.value, skipUniqueCategories: this.hasUniqueText && this.skipUniqueCategories.value,
      predictProbability: this.predictProbability.value, cutoff: this.cutoffValue};
  }

  /** Shows the inputs that apply to the checked [features] and the [target]. */
  update(features: DG.Column[], target: DG.Column | null): void {
    this.missingValues.update(missingColumnsOf(features));
    this.features = features;
    const classes = target === null ? null : twoClasses(target);
    this.hasTwoClasses = classes !== null;
    if (classes !== null) {
      this.predictProbability.setTooltip('Trains a regression model on the two classes and predicts the probability ' +
        `of ${classes[0]}. The score is the regression output, not a calibrated probability.`);
      this.cutoff.setTooltip(`Probability from which a row is predicted as ${classes[0]}.`);
    }
    this.updateVisibility();
  }

  /** Hides **Missing values** while the method takes the rows with missing values as they are ([keeps]). */
  setKeepsMissingValues(keeps: boolean): void {
    this.missingValues.setKeeps(keeps);
    this.updateVisibility();
  }

  /** Sets the inputs shown for the current selection to [steps] and [missingValues], as a user would; a hidden input
   * keeps its default. Skip unique categories goes first, as it changes the One-hot default. */
  restore(steps: PreparationSteps, missingValues: MissingValuesSettings): void {
    this.missingValues.restore(missingValues);
    for (const [input, value] of [[this.skipUniqueCategories, steps.skipUniqueCategories],
      [this.oneHot, steps.oneHot], [this.predictProbability, steps.predictProbability]] as const) {
      if (PreparationInputs.isShown(input))
        input.value = value;
    }
    if (PreparationInputs.isShown(this.cutoff))
      this.cutoff.value = steps.cutoff;
  }

  private get hasText(): boolean {
    return this.features.some((c) => c.isCategorical);
  }

  private get hasUniqueText(): boolean {
    return this.features.some(hasUniqueCategories);
  }

  private get stepInputs(): DG.InputBase[] {
    return [this.oneHot, this.skipUniqueCategories, this.predictProbability, this.cutoff];
  }

  private updateVisibility(): void {
    this.updateOneHot();
    PreparationInputs.show(this.skipUniqueCategories, this.hasUniqueText, true);
    PreparationInputs.show(this.predictProbability, this.hasTwoClasses, false);
    PreparationInputs.show(this.cutoff, this.hasTwoClasses && this.predictProbability.value, DEFAULT_CUTOFF);
    const isApplicable = this.missingValues.visibleInputs.length > 0 || this.hasText || this.hasTwoClasses;
    ui.setDisplay(this.group.root, isApplicable);
  }

  /** Shows One-hot encoding while a feature is text or yes/no, checked by its {@link defaultOneHot}; it follows that
   * until the user toggles it, and is unchecked while hidden. */
  private updateOneHot(): void {
    if (!this.hasText) {
      this.oneHotDefault = false;
      PreparationInputs.show(this.oneHot, false, false);
      this.isOneHotChosen = false;
      return;
    }
    ui.setDisplay(this.oneHot.root, true);
    this.oneHotDefault = defaultOneHot(this.features, this.skipUniqueCategories.value);
    if (!this.isOneHotChosen)
      this.oneHot.value = this.oneHotDefault;
  }

  /** Shows or hides [input]; a hidden input gets [defaultValue] back. */
  private static show<T>(input: DG.InputBase<T>, isShown: boolean, defaultValue: T): void {
    const wasShown = PreparationInputs.isShown(input);
    ui.setDisplay(input.root, isShown);
    if (wasShown && !isShown)
      input.value = defaultValue;
  }

  private static isShown(input: DG.InputBase): boolean {
    return input.root.style.display !== 'none';
  }
}
