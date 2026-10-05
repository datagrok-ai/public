import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {defaultValuesOf, Hyperparameters} from '../engines/engine';
import {imputeFunction, imputeSettingsOf, MissingColumn, MissingValuesSettings} from '../preparation/missing-values';

const SKIP_ROWS = 'Skip rows';
const IMPUTE = 'Impute';

/** The **Missing values** radio (Skip rows / Impute) and the imputation inputs under it, shown for gapped columns. */
export class MissingValuesInputs {
  readonly choice: DG.InputBase<string | null>;
  private readonly imputeInputs = new Map<string, DG.InputBase>();
  private readonly defaults: Hyperparameters;
  private hasGaps = false;
  private isImputeShown = false;

  constructor(onChanged: () => void, problem: () => string | null) {
    const func = imputeFunction();
    const properties = func === undefined ? [] : imputeSettingsOf(func);
    this.defaults = defaultValuesOf(properties);
    this.choice = ui.input.radio('Missing values', {items: func === undefined ? [SKIP_ROWS] : [SKIP_ROWS, IMPUTE],
      value: SKIP_ROWS, nullable: false, onValueChanged: () => {
        this.updateVisibility();
        onChanged();
      }});
    this.choice.addValidator(problem);
    this.choice.root.classList.add('forge-inline-radio');
    // The default goes in with the options, before the change handler: the owner is not built yet.
    for (const p of properties) {
      this.imputeInputs.set(p.name, ui.input.forProperty(p, null, {value: this.defaults[p.name],
        onValueChanged: () => onChanged()}));
    }
    this.updateVisibility();
  }

  get inputs(): DG.InputBase[] {
    return [this.choice, ...this.imputeInputs.values()];
  }

  get visibleInputs(): DG.InputBase[] {
    return [...(this.hasGaps ? [this.choice] : []), ...(this.isImpute ? this.imputeInputs.values() : [])];
  }

  get isImpute(): boolean {
    return this.hasGaps && this.choice.value === IMPUTE;
  }

  /** Shows the inputs when some columns have missing values ([missing], from `missingColumnsOf`); the tooltip lists
   * them with their counts. */
  update(missing: MissingColumn[]): void {
    this.hasGaps = missing.length > 0;
    this.choice.setTooltip(`Rows with missing values in: ${missing.map((m) => `${m.name}: ${m.count}`).join(', ')}.`);
    this.updateVisibility();
  }

  settings(): MissingValuesSettings {
    if (!this.isImpute)
      return {mode: 'skip'};
    const neighbors: unknown = this.imputeInputs.get('neighbors')?.value;
    const distance: unknown = this.imputeInputs.get('distance')?.value;
    return typeof neighbors === 'number' && typeof distance === 'string' ?
      {mode: 'impute', impute: {neighbors, distance}} : {mode: 'impute'};
  }

  private updateVisibility(): void {
    const isImputeShown = this.isImpute;
    ui.setDisplay(this.choice.root, this.hasGaps);
    for (const input of this.imputeInputs.values())
      ui.setDisplay(input.root, isImputeShown);
    const isHiding = this.isImputeShown && !isImputeShown;
    this.isImputeShown = isImputeShown;
    // A hidden setting goes back to its default, so a value left invalid cannot block the form unseen.
    if (isHiding) {
      for (const [name, input] of this.imputeInputs) {
        const value = this.defaults[name];
        if (value !== undefined)
          input.value = value;
      }
    }
  }
}
