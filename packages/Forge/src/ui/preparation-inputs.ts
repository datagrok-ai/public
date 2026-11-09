import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BEST_THRESHOLD_METRICS, BestThresholdMetric, classesOf, METRIC_DESCRIPTIONS, METRIC_LABELS}
  from '../metrics/metrics';
import {missingColumnsOf, MissingValuesSettings} from '../preparation/missing-values';
import {DEFAULT_STEPS, DEFAULT_THRESHOLD, DEFAULT_THRESHOLD_METRIC, defaultOneHot, hasUniqueCategories,
  PreparationSteps, thresholdProblem} from '../preparation/pipeline';
import {CollapsibleGroup} from './collapsible-group';
import {MissingValuesInputs} from './missing-values-inputs';
import {CUSTOMIZED, TuneSettings} from './tune-settings';

const BALANCED_ACCURACY_TOOLTIP = `${METRIC_DESCRIPTIONS.balancedAccuracy} For two classes it is ` +
  '(1 + Youden\'s J) / 2, so it picks the threshold Youden\'s J picks.';

/** The **Preparation** group of the Train view: Missing values, One-hot encoding, Skip unique categories, Predict
 * probability, its Classification threshold, and its **Best threshold by** metric with **Set**, each shown only while
 * it applies to the selection; the group is hidden while none does. An input that hides goes back to its default, so it
 * reappears as it first appeared; One-hot encoding gets its default when it appears and follows it until the user
 * toggles it. The threshold follows the chosen metric's best (the view sets it on every result) until it is edited by
 * hand; **Set** sets the best again. An italic status under the threshold says which metric's best it is, or
 * `Customized`. The metric, the threshold and its hand edit outlive their hiding; other classes of the target reset
 * them. */
export class PreparationInputs {
  readonly group: CollapsibleGroup;
  readonly missingValues: MissingValuesInputs;
  readonly oneHot: DG.InputBase<boolean>;
  readonly skipUniqueCategories: DG.InputBase<boolean>;
  readonly predictProbability: DG.InputBase<boolean>;
  /** **Best threshold by**: the metric's label. */
  readonly thresholdMetric: DG.ChoiceInput<string | null>;
  /** **Set**, at the right of the metric's row: sets the threshold to the chosen metric's best; shown while the
   * threshold differs from it. */
  readonly setThresholdButton: HTMLButtonElement;
  readonly threshold: DG.InputBase<number | null>;
  /** The status under the threshold: `Best for <metric>` or `Customized`; hidden while no best exists. */
  readonly thresholdStatus: HTMLSpanElement;
  private readonly thresholdStatusRow: HTMLDivElement;
  /** The metric the threshold was last set to the best of. */
  private optimizedFor: BestThresholdMetric = DEFAULT_THRESHOLD_METRIC;
  private features: DG.Column[] = [];
  private classes: string[] = [];
  /** Why the method cannot predict probabilities; null while it can. */
  private probabilityReason: string | null = null;
  private oneHotDefault = false;
  private isOneHotChosen = false;
  /** The threshold was edited by hand: it stays until **Set**, and the steps carry no metric. */
  private isThresholdEdited = false;

  /** [onChanged] follows every input but Predict probability, which calls [onProbabilityChanged], and the threshold,
   * which calls [onThresholdChanged] (**Set** too); [bestThreshold] gives a metric's best threshold on the last
   * result, undefined while there is none; [missingValuesProblem] is the validator of **Missing values**. */
  constructor(onChanged: () => void, onProbabilityChanged: () => void, onThresholdChanged: () => void,
    private readonly bestThreshold: (metric: BestThresholdMetric) => number | undefined,
    missingValuesProblem: () => string | null) {
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
      onValueChanged: () => {
        this.updateVisibility();
        onProbabilityChanged();
      }});
    this.thresholdMetric = ui.input.choice<string | null>('Best threshold by', {
      items: BEST_THRESHOLD_METRICS.map((m) => METRIC_LABELS[m]), value: METRIC_LABELS[DEFAULT_THRESHOLD_METRIC],
      nullable: false, tooltipText: 'The metric the classification threshold maximizes on the validation rows.',
      onValueChanged: () => this.updateSet()});
    for (const option of this.thresholdMetric.input.querySelectorAll('option')) {
      const metric = BEST_THRESHOLD_METRICS.find((m) => METRIC_LABELS[m] === option.textContent);
      if (metric !== undefined)
        option.title = PreparationInputs.metricTooltip(metric);
    }
    this.setThresholdButton = ui.button('Set', () => {
      const best = this.bestThreshold(this.metric);
      if (best === undefined || !this.thresholdMetric.enabled)
        return;
      this.setThreshold(best, this.metric);
      onThresholdChanged();
    });
    ui.tooltip.bind(this.setThresholdButton,
      () => `Set the threshold that maximizes ${METRIC_LABELS[this.metric]} on the validation rows.`);
    this.thresholdMetric.root.classList.add('forge-threshold-metric');
    this.thresholdMetric.root.append(this.setThresholdButton);
    this.threshold = ui.input.float('Classification threshold', {value: DEFAULT_THRESHOLD, min: 0, max: 1,
      step: 0.01, showSlider: true, nullable: false, format: '3 significant digits', onValueChanged: () => {
        this.isThresholdEdited = true;
        onThresholdChanged();
      }});
    this.threshold.addValidator(() => this.thresholdProblem);
    const status = TuneSettings.statusLine();
    [this.thresholdStatus, this.thresholdStatusRow] = [status.text, status.row];
    for (const input of this.stepInputs)
      ui.setDisplay(input.root, false);
    this.group = new CollapsibleGroup('Preparation', this.inputs.map((input) => input.root));
    this.threshold.root.after(this.thresholdStatusRow);
    this.updateVisibility();
  }

  get inputs(): DG.InputBase[] {
    return [...this.missingValues.inputs, ...this.stepInputs];
  }

  get thresholdValue(): number {
    return this.threshold.value ?? DEFAULT_THRESHOLD;
  }

  /** Why the threshold cannot be used; null while it is a number from 0 to 1. */
  get thresholdProblem(): string | null {
    return thresholdProblem(this.threshold.value);
  }

  /** The metric chosen in **Best threshold by**. */
  get metric(): BestThresholdMetric {
    return BEST_THRESHOLD_METRICS.find((m) => METRIC_LABELS[m] === this.thresholdMetric.value) ??
      DEFAULT_THRESHOLD_METRIC;
  }

  /** The steps as chosen; a hidden input's step is off (a locked form keeps them). The threshold is optimized for the
   * chosen metric unless it was edited by hand. */
  steps(): PreparationSteps {
    return {oneHot: this.oneHot.value, skipUniqueCategories: this.hasUniqueText && this.skipUniqueCategories.value,
      predictProbability: this.predictProbability.value,
      threshold: this.thresholdValue, ...(this.isThresholdEdited ? {} : {thresholdMetric: this.metric})};
  }

  /** Sets the threshold to [value], optimized for [metric]; without one, [value] counts as a hand edit. The view is
   * not told. */
  setThreshold(value: number, metric?: BestThresholdMetric): void {
    this.threshold.notify = false;
    try {
      this.threshold.value = value;
    } finally {
      this.threshold.notify = true;
    }
    if (metric !== undefined)
      this.thresholdMetric.value = METRIC_LABELS[metric];
    this.isThresholdEdited = metric === undefined;
    if (metric !== undefined)
      this.optimizedFor = metric;
    this.updateSet();
  }

  /** **Set** and the status follow the chosen metric's best on the last result: without one both hide, Set shows
   * while the threshold differs from it. */
  updateSet(): void {
    const best = PreparationInputs.isShown(this.thresholdMetric) ? this.bestThreshold(this.metric) : undefined;
    ui.setDisplay(this.setThresholdButton, best !== undefined && this.threshold.value !== best);
    ui.setDisplay(this.thresholdStatusRow, best !== undefined);
    this.thresholdStatus.textContent = best === undefined ? '' : this.statusOf(best);
  }

  /** Shows the inputs that apply to the checked [features] and the [target]. */
  update(features: DG.Column[], target: DG.Column | null): void {
    this.missingValues.update(missingColumnsOf(features));
    this.features = features;
    const classes = target !== null && target.isCategorical ? classesOf(target) : [];
    if (classes.length !== this.classes.length || classes.some((c, i) => c !== this.classes[i]))
      this.setThreshold(DEFAULT_THRESHOLD, DEFAULT_THRESHOLD_METRIC);
    this.classes = classes;
    this.updateProbabilityTooltip();
    this.updateThresholdTooltip();
    this.updateVisibility();
  }

  /** Enables Predict probability for a method that gives class probabilities ([gives]); otherwise unchecks and
   * disables it, its tooltip saying why. */
  setMethod(name: string, gives: boolean): void {
    this.probabilityReason = gives ? null : `${name} gives no class probabilities`;
    if (!gives && this.predictProbability.value)
      this.predictProbability.value = false;
    this.predictProbability.enabled = gives;
    this.updateProbabilityTooltip();
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
    this.restoreThreshold(steps);
  }

  /** Sets a shown threshold to the one of [steps] with its metric (a hand edit without one). */
  restoreThreshold({threshold, thresholdMetric}: PreparationSteps): void {
    if (PreparationInputs.isShown(this.threshold))
      this.setThreshold(threshold, thresholdMetric);
  }

  private updateProbabilityTooltip(): void {
    const positive = this.classes.length === 2 ? this.classes[0] : undefined;
    this.predictProbability.setTooltip(this.probabilityReason ?? (positive === undefined ?
      'Adds the probability of each class when the model is applied.' :
      'Adds the probability of each class when the model is applied; for two classes, lets you choose the ' +
      `classification threshold of ${positive}.`));
  }

  /** The input's own tooltip (its caption's and editor's: a tooltip bound to the root never shows over them), the
   * classes in bold. */
  private updateThresholdTooltip(): void {
    const [positive, negative] = this.classes.map((c) => {
      const b = ui.element('b');
      b.textContent = c;
      return b.outerHTML;
    });
    this.threshold.setTooltip(`Rows whose probability of ${positive} is at least this value are predicted as ` +
      `${positive}, the others as ${negative}.`);
  }

  private get hasText(): boolean {
    return this.features.some((c) => c.isCategorical);
  }

  private get hasUniqueText(): boolean {
    return this.features.some(hasUniqueCategories);
  }

  private get stepInputs(): DG.InputBase[] {
    return [this.oneHot, this.skipUniqueCategories, this.predictProbability, this.threshold, this.thresholdMetric];
  }

  /** `Best for <metric>`: the chosen one's while the threshold is its [best], else the one it was set to the best of;
   * `Customized` after a hand edit. */
  private statusOf(best: number): string {
    if (this.threshold.value === best)
      return `Best for ${METRIC_LABELS[this.metric]}`;
    return this.isThresholdEdited ? CUSTOMIZED : `Best for ${METRIC_LABELS[this.optimizedFor]}`;
  }

  private updateVisibility(): void {
    this.updateOneHot();
    PreparationInputs.show(this.skipUniqueCategories, this.hasUniqueText, true);
    const hasClasses = this.classes.length >= 2;
    PreparationInputs.show(this.predictProbability, hasClasses, false);
    const hasThreshold = this.classes.length === 2 && this.predictProbability.value;
    for (const input of [this.thresholdMetric, this.threshold])
      ui.setDisplay(input.root, hasThreshold);
    this.updateSet();
    const isApplicable = this.missingValues.visibleInputs.length > 0 || this.hasText || hasClasses;
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

  private static metricTooltip(metric: BestThresholdMetric): string {
    return metric === 'balancedAccuracy' ? BALANCED_ACCURACY_TOOLTIP : METRIC_DESCRIPTIONS[metric];
  }
}
