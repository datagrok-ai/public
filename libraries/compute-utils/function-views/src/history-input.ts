import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {openRunDialog} from './custom-function-view';
import {getRunCreated, getRunName, getRunTitle} from '../../shared-utils/utils';

export interface HistoryInputOptions extends ui.input.IInputInitOptions<DG.FuncCall | null> {
  /** nqName of the model whose saved runs the dialog lists */
  model: string;
}

/** A saved run of a model, picked in Compute2's Open dialog; its value is loaded with its inputs and outputs. */
export class HistoryInput extends DG.JsInputBase<DG.FuncCall | null> {
  private run: DG.FuncCall | null = null;
  private readonly field = ui.element('input') as HTMLInputElement;
  private readonly clearIcon = ui.iconFA('times', () => this.pick(null), 'Clear');

  get inputType(): string {return 'History';}
  get dataType(): string {return DG.TYPE.FUNC_CALL;}

  constructor(name: string, private readonly options: HistoryInputOptions) {
    super();
    this.field.classList.add('ui-input-editor');
    this.field.readOnly = true;
    this.field.placeholder = 'No run chosen';
    this.clearIcon.style.display = 'none';
    this.field.addEventListener('click', () => this.openDialog());
    this.addCaption(name);
    this.addOptions(ui.divH([this.clearIcon, ui.iconFA('search', () => this.openDialog(), 'Choose a saved run')]));
    if (options.nullable !== undefined)
      this.nullable = options.nullable;
    if (options.tooltipText)
      this.setTooltip(options.tooltipText);
    this.setValue(options.value ?? null);
    if (options.onValueChanged)
      this.onChanged.subscribe(() => options.onValueChanged!(this.value, this));
    options.onCreated?.(this);
  }

  getInput(): HTMLElement {return this.field;}

  getValue(): DG.FuncCall | null {return this.run;}

  setValue(value: DG.FuncCall | null): void {
    if (value === this.run)
      return;
    this.run = value;
    this.field.value = this.getStringValue();
    this.clearIcon.style.display = value && this.nullable ? '' : 'none';
    this.fireChanged();
  }

  getStringValue(): string {
    if (!this.run)
      return '';
    const created = getRunCreated(this.run)?.format('YYYY-MM-DD HH:mm');
    return getRunName(this.run) && created ? `${getRunName(this.run)} · ${created}` : getRunTitle(this.run);
  }

  // a run can't be found by its label
  setStringValue(_value: string): void {}

  private pick(run: DG.FuncCall | null) {
    this.setValue(run);
    this.fireInput();
  }

  private async openDialog() {
    try {
      const run = await openRunDialog(DG.Func.byName(this.options.model));
      if (run)
        this.pick(run);
    } catch (e: any) {
      grok.shell.error(`Couldn't open the run: ${e?.message ?? e}`);
    }
  }
}

/** An input for a saved run of `options.model`, picked in Compute2's Open dialog. */
export function historyInput(name: string, options: HistoryInputOptions): HistoryInput {
  return new HistoryInput(name, options);
}
