import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {taskOf} from '../metrics/metrics';
import {DEFAULT_SEED, FOLDS, FOLDS_RANGE, TEST_FRACTION, TEST_FRACTION_RANGE, ValidationSettings}
  from '../training/k-fold';
import {MAX_SEED} from '../training/train-model';
import {CollapsibleGroup} from './collapsible-group';

const STRATIFIED = 'K-fold, stratified';
const K_FOLD = 'K-fold';
const HOLD_OUT = 'Hold-out';

/** The **Validation** group of the Train view, collapsed at first: Scheme, Folds (k-fold) or Test fraction (hold-out),
 * Shuffle and Split seed (enabled while Shuffle is on). K-fold, stratified is listed for a classification target only,
 * and is its default. */
export class ValidationInputs {
  readonly group: CollapsibleGroup;
  readonly scheme: DG.ChoiceInput<string | null>;
  readonly folds: DG.InputBase<number | null>;
  readonly testFraction: DG.InputBase<number | null>;
  readonly shuffle: DG.InputBase<boolean>;
  readonly seed: DG.InputBase<number | null>;
  private isClassification = true;
  /** The user picked the scheme; a target change keeps it while it is listed. */
  private isSchemeChosen = false;
  private isUpdating = false;

  constructor(onChanged: () => void) {
    this.scheme = ui.input.choice<string | null>('Scheme', {items: [STRATIFIED, K_FOLD, HOLD_OUT], value: STRATIFIED,
      nullable: false, tooltipText: 'How rows are split',
      onValueChanged: () => {
        if (!this.isUpdating)
          this.isSchemeChosen = true;
        this.updateVisibility();
        onChanged();
      }});
    const [fewestFolds, mostFolds] = FOLDS_RANGE;
    this.folds = ui.input.int('Folds', {value: FOLDS, min: fewestFolds, max: mostFolds, nullable: false,
      tooltipText: 'Number of cross-validation folds', onValueChanged: onChanged});
    const [fewest, most] = TEST_FRACTION_RANGE;
    this.testFraction = ui.input.float('Test fraction', {value: TEST_FRACTION, min: fewest, max: most, step: 0.05,
      nullable: false, tooltipText: 'Share of rows held out', onValueChanged: onChanged});
    this.shuffle = ui.input.bool('Shuffle', {value: true, tooltipText: 'Shuffle rows before splitting',
      onValueChanged: () => {
        this.updateVisibility();
        onChanged();
      }});
    this.seed = ui.input.int('Split seed', {value: DEFAULT_SEED, min: 0, max: MAX_SEED, nullable: false,
      showPlusMinus: false, tooltipText: 'Seed for random splitting', onValueChanged: onChanged});
    this.group = new CollapsibleGroup('Validation', this.inputs.map((input) => input.root), false);
    this.updateVisibility();
  }

  get inputs(): DG.InputBase[] {
    return [this.scheme, this.folds, this.testFraction, this.shuffle, this.seed];
  }

  /** The inputs whose values the validation uses now: Folds or Test fraction for the scheme, Split seed when
   * shuffled. */
  get usedInputs(): DG.InputBase[] {
    const isHoldOut = this.scheme.value === HOLD_OUT;
    return [this.scheme, isHoldOut ? this.testFraction : this.folds, this.shuffle,
      ...(this.shuffle.value ? [this.seed] : [])];
  }

  get seedValue(): number {
    return this.seed.value ?? DEFAULT_SEED;
  }

  /** The validation as chosen; the hold out of a classification is stratified. */
  settings(): ValidationSettings {
    const scheme = this.scheme.value;
    return {scheme: scheme === HOLD_OUT ? 'holdout' : 'kfold', folds: this.folds.value ?? FOLDS,
      testFraction: this.testFraction.value ?? TEST_FRACTION, stratified: scheme !== K_FOLD,
      shuffle: this.shuffle.value};
  }

  /** Lists the schemes of the [target]'s task: a regression has no stratified one, so K-fold replaces it; a
   * classification gets its default back unless the user chose K-fold. */
  update(target: DG.Column | null): void {
    const isClassification = target === null ? this.isClassification : taskOf(target) === 'classification';
    if (isClassification === this.isClassification)
      return;
    this.isClassification = isClassification;
    const value = this.scheme.value;
    this.isUpdating = true;
    try {
      this.scheme.items = isClassification ? [STRATIFIED, K_FOLD, HOLD_OUT] : [K_FOLD, HOLD_OUT];
      if (!isClassification && value === STRATIFIED)
        this.isSchemeChosen = false;
      this.scheme.value = value === HOLD_OUT ? HOLD_OUT :
        isClassification && !this.isSchemeChosen ? STRATIFIED : K_FOLD;
    } finally {
      this.isUpdating = false;
    }
  }

  /** Sets the inputs to [validation] and [seed], as a user would. */
  restore(validation: ValidationSettings, seed: number): void {
    this.scheme.value = validation.scheme === 'holdout' ? HOLD_OUT :
      validation.stratified && this.isClassification ? STRATIFIED : K_FOLD;
    this.folds.value = validation.folds;
    this.testFraction.value = validation.testFraction;
    this.shuffle.value = validation.shuffle;
    this.seed.value = seed;
  }

  private updateVisibility(): void {
    const isHoldOut = this.scheme.value === HOLD_OUT;
    ui.setDisplay(this.folds.root, !isHoldOut);
    ui.setDisplay(this.testFraction.root, isHoldOut);
    this.seed.enabled = this.shuffle.value;
  }
}
