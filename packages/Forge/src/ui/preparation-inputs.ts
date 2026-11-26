import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {missingColumnsOf, MissingValuesSettings} from '../preparation/missing-values';
import {DEFAULT_STEPS, defaultOneHot, hasUniqueCategories, isEncodable, PreparationSteps}
  from '../preparation/pipeline';
import {CollapsibleGroup} from './collapsible-group';
import {MissingValuesInputs} from './missing-values-inputs';

/** The **Preparation** group of the Train view: Missing values, One-hot encoding and Skip unique categories, each
 * shown only while it applies to the selection; the group is hidden while none does. An input that hides goes back to
 * its default, so it reappears as it first appeared; One-hot encoding gets its default when it appears and follows it
 * until the user toggles it. */
export class PreparationInputs {
  readonly group: CollapsibleGroup;
  readonly missingValues: MissingValuesInputs;
  readonly oneHot: DG.InputBase<boolean>;
  readonly skipUniqueCategories: DG.InputBase<boolean>;
  private features: DG.Column[] = [];
  private oneHotDefault = false;
  private isOneHotChosen = false;

  /** [onChanged] follows every input; [missingValuesProblem] is the validator of **Missing values**; [semTypes]: the
   * semantic types the view's methods take as they are; such features are not text here. */
  constructor(onChanged: () => void, missingValuesProblem: () => string | null,
    private readonly semTypes: readonly string[]) {
    this.missingValues = new MissingValuesInputs(onChanged, missingValuesProblem);
    this.oneHot = ui.input.bool('One-hot encoding', {value: false,
      tooltipText: 'Lets the model use categories as numbers without implying an order: each category gets its own ' +
        '0/1 column.',
      onValueChanged: (value) => {
        if (value !== this.oneHotDefault)
          this.isOneHotChosen = true;
        onChanged();
      }});
    this.skipUniqueCategories = ui.input.bool('Skip unique categories', {value: DEFAULT_STEPS.skipUniqueCategories,
      tooltipText: 'Leaves out text features whose values are all different, such as ids.',
      onValueChanged: () => {
        this.updateVisibility();
        onChanged();
      }});
    // Hidden first: the first update must not reset the steps to their defaults.
    for (const input of this.stepInputs)
      input.visible = false;
    this.group = new CollapsibleGroup('Preparation', this.inputs.map((input) => input.root));
    this.updateVisibility();
  }

  get inputs(): DG.InputBase[] {
    return [...this.missingValues.inputs, ...this.stepInputs];
  }

  /** The steps as chosen; a hidden input's step is off (a locked form keeps them). */
  steps(): Pick<PreparationSteps, 'oneHot' | 'skipUniqueCategories'> {
    return {oneHot: this.oneHot.value, skipUniqueCategories: this.hasUniqueText && this.skipUniqueCategories.value};
  }

  /** Shows the inputs that apply to the checked [features]. */
  update(features: DG.Column[]): void {
    this.missingValues.update(missingColumnsOf(features));
    this.features = features;
    this.updateVisibility();
  }

  /** Hides **Missing values** while the method takes the rows with missing values as they are ([keeps]). */
  setKeepsMissingValues(keeps: boolean): void {
    this.missingValues.setKeeps(keeps);
    this.updateVisibility();
  }

  /** Sets **Missing values** to [settings] without telling the owner and shows the group as it then applies. */
  restoreMissingValues(settings: MissingValuesSettings): void {
    this.missingValues.restore(settings);
    this.updateVisibility();
  }

  /** Sets the inputs shown for the current selection to [steps] and [missingValues] ({@link restoreMissingValues}), the
   * steps as a user would; a hidden input keeps its default, and an input that already holds its value is not set, so
   * nothing is checked for it. Skip unique categories goes first, as it changes the One-hot default. */
  restore(steps: PreparationSteps, missingValues: MissingValuesSettings): void {
    this.restoreMissingValues(missingValues);
    for (const [input, value] of [[this.skipUniqueCategories, steps.skipUniqueCategories],
      [this.oneHot, steps.oneHot]] as const) {
      if (input.visible && input.value !== value)
        input.value = value;
    }
  }

  private get hasText(): boolean {
    return this.features.some((c) => isEncodable(c, this.semTypes));
  }

  private get hasUniqueText(): boolean {
    return this.features.some((c) => isEncodable(c, this.semTypes) && hasUniqueCategories(c));
  }

  private get stepInputs(): DG.InputBase[] {
    return [this.oneHot, this.skipUniqueCategories];
  }

  private updateVisibility(): void {
    this.updateOneHot();
    CollapsibleGroup.showInput(this.skipUniqueCategories, this.hasUniqueText, true);
    const isApplicable = this.missingValues.visibleInputs.length > 0 || this.hasText;
    ui.setDisplay(this.group.root, isApplicable);
  }

  /** Shows One-hot encoding while a feature is text or yes/no, checked by its {@link defaultOneHot}; it follows that
   * until the user toggles it, and is unchecked while hidden. */
  private updateOneHot(): void {
    if (!this.hasText) {
      this.oneHotDefault = false;
      CollapsibleGroup.showInput(this.oneHot, false, false);
      this.isOneHotChosen = false;
      return;
    }
    this.oneHot.visible = true;
    this.oneHotDefault = defaultOneHot(this.features, this.skipUniqueCategories.value, this.semTypes);
    if (!this.isOneHotChosen && this.oneHot.value !== this.oneHotDefault)
      this.oneHot.value = this.oneHotDefault;
  }
}
