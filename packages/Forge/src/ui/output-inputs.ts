import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BEST_THRESHOLD_METRICS, BestThresholdMetric, classesOf, METRIC_DESCRIPTIONS, METRIC_LABELS, MetricId,
  metricOfLabel} from '../metrics/metrics';
import {DEFAULT_STEPS, DEFAULT_THRESHOLD, DEFAULT_THRESHOLD_METRIC, PreparationSteps, thresholdProblem}
  from '../preparation/pipeline';
import {CollapsibleGroup} from './collapsible-group';
import {CUSTOMIZED, TuneSettings} from './tune-settings';

const BALANCED_ACCURACY_TOOLTIP = `${METRIC_DESCRIPTIONS.balancedAccuracy} For two classes it is ` +
  '(1 + Youden\'s J) / 2, so it picks the threshold Youden\'s J picks.';

/** The **Output** group of the Train view; the threshold follows the chosen metric's best until edited by hand. */
export class OutputInputs {
  readonly group: CollapsibleGroup;
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
  private classes: string[] = [];
  /** Why the method cannot predict probabilities; null while it can. */
  private probabilityReason: string | null = null;
  /** The threshold was edited by hand: it stays until **Set**, and the steps carry no metric. */
  private isThresholdEdited = false;

  /** Predict probability calls [onProbabilityChanged], the threshold [onThresholdChanged] (**Set** too);
   * [bestThreshold] gives a metric's best threshold on the last result, undefined while there is none. */
  constructor(onProbabilityChanged: () => void, onThresholdChanged: () => void,
    private readonly bestThreshold: (metric: BestThresholdMetric) => number | undefined) {
    this.predictProbability = ui.input.bool('Predict probability', {value: DEFAULT_STEPS.predictProbability,
      onValueChanged: () => {
        this.updateVisibility();
        onProbabilityChanged();
      }});
    this.thresholdMetric = ui.input.choice<string | null>('Best threshold by', {
      items: BEST_THRESHOLD_METRICS.map((m) => METRIC_LABELS[m]), value: METRIC_LABELS[DEFAULT_THRESHOLD_METRIC],
      nullable: false, tooltipText: 'The metric the classification threshold maximizes on the validation rows.',
      onValueChanged: () => this.updateSet()});
    TuneSettings.describeMetricOptions(this.thresholdMetric.input, OutputInputs.metricTooltip);
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
    // Hidden first: the first update must not reset Predict probability to its default.
    for (const input of this.inputs)
      input.visible = false;
    this.group = new CollapsibleGroup('Output', this.inputs.map((input) => input.root));
    this.threshold.root.after(this.thresholdStatusRow);
    this.updateVisibility();
  }

  get inputs(): DG.InputBase[] {
    return [this.predictProbability, this.threshold, this.thresholdMetric];
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
    const id = metricOfLabel(this.thresholdMetric.value);
    return BEST_THRESHOLD_METRICS.find((m) => m === id) ?? DEFAULT_THRESHOLD_METRIC;
  }

  /** The output steps as chosen. The threshold is optimized for the chosen metric unless it was edited by hand. */
  steps(): Pick<PreparationSteps, 'predictProbability' | 'threshold' | 'thresholdMetric'> {
    return {predictProbability: this.predictProbability.value, threshold: this.thresholdValue,
      ...(this.isThresholdEdited ? {} : {thresholdMetric: this.metric})};
  }

  /** Sets the threshold to [value], optimized for [metric]; without one, [value] counts as a hand edit. The view is
   * not told. */
  setThreshold(value: number, metric?: BestThresholdMetric): void {
    TuneSettings.silently(this.threshold, () => this.threshold.value = value);
    if (metric !== undefined) {
      this.thresholdMetric.value = METRIC_LABELS[metric];
      this.optimizedFor = metric;
    }
    this.isThresholdEdited = metric === undefined;
    this.updateSet();
  }

  /** **Set** and the status follow the chosen metric's best on the last result: without one both hide, Set shows
   * while the threshold differs from it. */
  updateSet(): void {
    const best = this.thresholdMetric.visible ?this.bestThreshold(this.metric) : undefined;
    ui.setDisplay(this.setThresholdButton, best !== undefined && this.threshold.value !== best);
    ui.setDisplay(this.thresholdStatusRow, best !== undefined);
    this.thresholdStatus.textContent = best === undefined ? '' : this.statusOf(best);
  }

  /** Shows the inputs that apply to the [target]. */
  update(target: DG.Column | null): void {
    const classes = target !== null && target.isCategorical ? classesOf(target) : [];
    if (classes.length !== this.classes.length || classes.some((c, i) => c !== this.classes[i]))
      this.setThreshold(DEFAULT_THRESHOLD, DEFAULT_THRESHOLD_METRIC);
    this.classes = classes;
    this.updateProbabilityTooltip();
    this.updateThresholdTooltip();
    this.updateVisibility();
  }

  /** Disables and unchecks Predict probability with [reason] as its tooltip; null enables it. */
  setProbabilityReason(reason: string | null): void {
    this.probabilityReason = reason;
    if (reason !== null && this.predictProbability.value)
      this.predictProbability.value = false;
    this.predictProbability.enabled = reason === null;
    this.updateProbabilityTooltip();
  }

  /** Sets the inputs shown for the current target to [steps], as a user would; a hidden input keeps its default. */
  restore(steps: PreparationSteps): void {
    if (this.predictProbability.visible)
      this.predictProbability.value = steps.predictProbability;
    this.restoreThreshold(steps);
  }

  /** Sets a shown threshold to the one of [steps] with its metric (a hand edit without one). */
  restoreThreshold({threshold, thresholdMetric}: PreparationSteps): void {
    if (this.threshold.visible)
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

  /** `Best for <metric>`: the chosen one's while the threshold is its [best], else the one it was set to the best of;
   * `Customized` after a hand edit. */
  private statusOf(best: number): string {
    if (this.threshold.value === best)
      return `Best for ${METRIC_LABELS[this.metric]}`;
    return this.isThresholdEdited ? CUSTOMIZED : `Best for ${METRIC_LABELS[this.optimizedFor]}`;
  }

  private updateVisibility(): void {
    const hasClasses = this.classes.length >= 2;
    CollapsibleGroup.showInput(this.predictProbability, hasClasses, false);
    const hasThreshold = this.classes.length === 2 && this.predictProbability.value;
    for (const input of [this.thresholdMetric, this.threshold])
      input.visible = hasThreshold;
    this.updateSet();
    ui.setDisplay(this.group.root, hasClasses);
  }

  private static metricTooltip(metric: MetricId): string {
    return metric === 'balancedAccuracy' ? BALANCED_ACCURACY_TOOLTIP : METRIC_DESCRIPTIONS[metric];
  }
}
