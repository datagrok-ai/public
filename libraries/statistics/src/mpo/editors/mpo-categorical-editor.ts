import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';
import {CategoricalDesirability} from '../mpo';

export class MpoCategoricalEditor {
  root = ui.divV([], {style: {width: '300px'}});
  onChanged = new Subject<CategoricalDesirability>();
  supportsModeDialog: boolean = false;

  private _prop: CategoricalDesirability;
  private desirabilityInputs: DG.InputBase[] = [];
  private form: HTMLElement | null = null;

  constructor(prop: CategoricalDesirability) {
    this._prop = prop;
    this.buildForm();
  }

  private buildForm(): void {
    const categories = this._prop.categories ?? [];

    if (this.form) {
      ui.empty(this.root);
      this.form = null;
      this.desirabilityInputs = [];
    }

    this.desirabilityInputs = categories.map((cat) => {
      return ui.input.float(cat.name, {
        value: cat.desirability ?? 0.5,
        min: 0,
        max: 1,
        format: '#0.000',
        onValueChanged: (value: number) => {
          cat.desirability = value;
          this.onChanged.next(this._prop);
        },
      });
    });

    this.form = ui.form(this.desirabilityInputs);
    this.root.append(this.form);
  }

  redrawAll(notify: boolean = true): void {
    const categories = this._prop.categories ?? [];

    categories.forEach((cat, idx) => {
      const desirabilityInput = this.desirabilityInputs[idx];
      if (!desirabilityInput)
        return;

      desirabilityInput.value = cat.desirability ?? 0.5;
    });

    if (notify)
      this.onChanged.next(this._prop);
  }

  setColumn(col: DG.Column | null): void {
    if (!col)
      return;

    if (!this._prop.categories) {
      this._prop.categories = col.categories.map((c) => ({
        name: c,
        desirability: 1,
      }));
    }
    this.buildForm();
  }
}
