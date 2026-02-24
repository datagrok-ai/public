import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {
  CronFields, CRON_PRESETS, cronToDescription, cronToString, dowTokenToNames,
  findMatchingPreset, formatNextRun, getNextRun, parseCron, validateCron, validateCronPerField,
} from './cron-utils';

const CUSTOM_VALUE = '__custom__';

export class CronInput extends DG.JsInputBase<string> {
  private _value: string = '* * * * *';

  private readonly presetSelect: HTMLSelectElement;
  private readonly toggleIcon: HTMLElement;
  private readonly infoLabel: HTMLElement;
  private readonly advancedPanel: HTMLDivElement;

  private fieldInputs: {[key: string]: HTMLInputElement} = {};
  private readonly rawInput: HTMLInputElement;
  private readonly descriptionLine: HTMLElement;
  private readonly dowHint: HTMLElement;
  private readonly editorRoot: HTMLDivElement;

  private _syncing = false;
  private _outsideHandler: ((e: MouseEvent) => void) | null = null;

  get inputType(): string { return 'Cron'; }
  get dataType(): string { return DG.TYPE.STRING; }

  constructor() {
    super();

    this.presetSelect = this.createPresetSelect();
    this.toggleIcon = this.createToggleIcon();
    this.infoLabel = ui.div([], 'pp-cron-info');
    this.advancedPanel = ui.div([], 'pp-cron-advanced') as HTMLDivElement;
    this.rawInput = ui.element('input') as HTMLInputElement;
    this.descriptionLine = ui.div([], 'pp-cron-description');
    this.dowHint = ui.div([], 'pp-cron-dow-hint');
    this.editorRoot = this.buildUI();

    this.addValidator((v) => validateCron(v));
    this.syncUI();
  }

  getInput(): HTMLElement { return this.editorRoot; }

  getValue(): string { return this._value; }

  setValue(value: string): void {
    if (value === this._value)
      return;
    this._value = value;
    this.syncUI();
    this.fireChanged();
  }

  getStringValue(): string { return this._value; }

  setStringValue(value: string): void { this.setValue(value); }

  private createPresetSelect(): HTMLSelectElement {
    const select = ui.element('select') as HTMLSelectElement;
    select.className = 'pp-cron-preset-select';
    for (const preset of CRON_PRESETS) {
      const opt = document.createElement('option');
      opt.value = preset.value;
      opt.textContent = preset.label;
      select.appendChild(opt);
    }
    const customOpt = document.createElement('option');
    customOpt.value = CUSTOM_VALUE;
    customOpt.textContent = 'Custom...';
    select.appendChild(customOpt);

    select.addEventListener('change', () => {
      if (this._syncing)
        return;
      if (select.value === CUSTOM_VALUE) {
        this.showAdvanced(true);
        return;
      }
      this._value = select.value;
      this.syncUI();
      this.fireInput();
      this.fireChanged();
    });
    return select;
  }

  private createToggleIcon(): HTMLElement {
    const icon = ui.iconFA('cog', () => {
      const visible = !this.advancedPanel.classList.contains('pp-cron-collapsed');
      this.showAdvanced(!visible);
    }, 'Toggle advanced editor');
    icon.className = 'pp-cron-toggle grok-icon fal fa-cog';
    return icon;
  }

  private showAdvanced(show: boolean): void {
    this.advancedPanel.classList.toggle('pp-cron-collapsed', !show);
    this.toggleIcon.classList.toggle('pp-cron-toggle-active', show);

    if (show) {
      requestAnimationFrame(() => {
        this._outsideHandler = (e: MouseEvent) => {
          if (!this.editorRoot.contains(e.target as Node))
            this.showAdvanced(false);
        };
        document.addEventListener('mousedown', this._outsideHandler);
      });
    }
    else if (this._outsideHandler) {
      document.removeEventListener('mousedown', this._outsideHandler);
      this._outsideHandler = null;
    }
  }

  private buildUI(): HTMLDivElement {
    const mainRow = ui.div([this.presetSelect, this.toggleIcon], 'pp-cron-main-row');

    const fieldDefs: {key: string; label: string; placeholder: string}[] = [
      {key: 'minute', label: 'Minute', placeholder: '*'},
      {key: 'hour', label: 'Hour', placeholder: '*'},
      {key: 'dayOfMonth', label: 'Day of Mo.', placeholder: '*'},
      {key: 'month', label: 'Month', placeholder: '*'},
      {key: 'dayOfWeek', label: 'Weekday', placeholder: '*'},
    ];

    const fieldElements: HTMLElement[] = [];
    for (const f of fieldDefs) {
      const input = ui.element('input') as HTMLInputElement;
      input.className = 'pp-cron-field-input';
      input.placeholder = f.placeholder;
      input.addEventListener('input', () => this.onFieldInput());
      this.fieldInputs[f.key] = input;

      const label = ui.div([], 'pp-cron-field-label');
      label.textContent = f.label;
      const fieldEl = ui.div([label, input], 'pp-cron-field');
      if (f.key === 'dayOfWeek')
        fieldEl.appendChild(this.dowHint);
      fieldElements.push(fieldEl);
    }
    const fieldsRow = ui.div(fieldElements, 'pp-cron-fields-row');

    this.rawInput.className = 'pp-cron-raw-input';
    this.rawInput.placeholder = '* * * * *';
    this.rawInput.readOnly = true;
    this.rawInput.addEventListener('focus', () => { this.rawInput.readOnly = false; });
    this.rawInput.addEventListener('blur', () => { this.rawInput.readOnly = true; });
    this.rawInput.addEventListener('input', () => this.onRawInput());
    const rawRow = ui.div([this.rawInput], 'pp-cron-raw-row');

    this.advancedPanel.appendChild(fieldsRow);
    this.advancedPanel.appendChild(rawRow);
    this.advancedPanel.appendChild(this.descriptionLine);
    this.advancedPanel.classList.add('pp-cron-collapsed');

    return ui.div([mainRow, this.infoLabel, this.advancedPanel], 'pp-cron-root') as HTMLDivElement;
  }

  private onFieldInput(): void {
    if (this._syncing)
      return;
    const expr = cronToString(this.readFieldInputs());
    this._value = expr;
    this._syncing = true;
    try {
      this.rawInput.value = expr;
      this.syncPreset();
      this.updateValidationState();
    }
    finally {
      this._syncing = false;
    }
    this.fireInput();
    this.fireChanged();
  }

  private onRawInput(): void {
    if (this._syncing)
      return;
    const expr = this.rawInput.value.trim();
    this._value = expr;
    this._syncing = true;
    try {
      const fields = parseCron(expr);
      if (fields)
        this.setFieldInputs(fields);
      this.syncPreset();
      this.updateValidationState();
    }
    finally {
      this._syncing = false;
    }
    this.fireInput();
    this.fireChanged();
  }

  private syncUI(): void {
    this._syncing = true;
    try {
      const fields = parseCron(this._value);
      if (fields)
        this.setFieldInputs(fields);
      this.rawInput.value = this._value;
      this.syncPreset();
      this.updateValidationState();
    }
    finally {
      this._syncing = false;
    }
  }

  private setFieldInputs(fields: CronFields): void {
    for (const key of Object.keys(this.fieldInputs) as (keyof CronFields)[])
      this.fieldInputs[key].value = fields[key];
  }

  private readFieldInputs(): CronFields {
    return {
      minute: this.fieldInputs['minute'].value || '*',
      hour: this.fieldInputs['hour'].value || '*',
      dayOfMonth: this.fieldInputs['dayOfMonth'].value || '*',
      month: this.fieldInputs['month'].value || '*',
      dayOfWeek: this.fieldInputs['dayOfWeek'].value || '*',
    };
  }

  private syncPreset(): void {
    const preset = findMatchingPreset(this._value);
    if (preset)
      this.presetSelect.value = preset.value;
    else
      this.presetSelect.value = CUSTOM_VALUE;
  }

  private updateValidationState(): void {
    const fields = parseCron(this._value);
    let errors: {[key: string]: string} = {};
    let firstError: string | null = null;

    if (!fields) {
      firstError = 'Expression must have exactly 5 fields';
      for (const key of Object.keys(this.fieldInputs))
        this.fieldInputs[key].classList.add('pp-cron-invalid');
    }
    else {
      errors = validateCronPerField(fields);
      const errorKeys = Object.keys(errors);
      if (errorKeys.length > 0)
        firstError = errors[errorKeys[0]];
      for (const key of Object.keys(this.fieldInputs))
        this.fieldInputs[key].classList.toggle('pp-cron-invalid', key in errors);
    }

    const hasErrors = firstError !== null;
    this.rawInput.classList.toggle('pp-cron-invalid', hasErrors);
    this.presetSelect.classList.toggle('pp-cron-invalid', hasErrors);

    if (hasErrors) {
      this.infoLabel.textContent = firstError;
      this.infoLabel.classList.add('pp-cron-info-error');
      this.descriptionLine.textContent = '';
    }
    else {
      this.infoLabel.classList.remove('pp-cron-info-error');
      const next = getNextRun(this._value);
      this.infoLabel.textContent = next ? 'Next run: ' + formatNextRun(next) : '';
      this.descriptionLine.textContent = cronToDescription(this._value);
    }

    const parsed = parseCron(this._value);
    const dowNames = parsed ? dowTokenToNames(parsed.dayOfWeek) : null;
    this.dowHint.textContent = dowNames ?? '';
  }
}
