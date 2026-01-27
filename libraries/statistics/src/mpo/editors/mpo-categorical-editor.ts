import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';
import {PropertyDesirability} from '../mpo';

export class MpoCategoricalEditor {
  root = ui.divV([], {style: {width: '300px'}});
  onChanged = new Subject<PropertyDesirability>();
  supportsModeDialog: boolean = false;

  private _prop: PropertyDesirability;
  private sliders: DG.InputBase[] = [];
  private form: HTMLElement | null = null;

  constructor(prop: PropertyDesirability) {
    this._prop = prop;
    this.buildForm();
  }

  private buildForm(): void {
    const categories = this._prop.categories ?? [];

    if (this.form) {
      ui.empty(this.root);
      this.form = null;
      this.sliders = [];
    }

    this.sliders = categories.map((cat, idx) => {
      return ui.input.slider(cat.name, {
        value: cat.desirability ?? 0.5,
        min: 0,
        max: 1,
        step: 0.01,
        onValueChanged: (value: number) => {
          this._prop.categories![idx].desirability = value;
          this.onChanged.next(this._prop);
        },
      });
    });

    this.form = ui.form(this.sliders);
    this.root.append(this.form);
  }

  redrawAll(notify: boolean = true): void {
    const categories = this._prop.categories ?? [];

    categories.forEach((cat, idx) => {
      const slider = this.sliders[idx] as any;
      if (!slider) return;

      slider.value = cat.desirability ?? 0.5;
    });

    if (notify)
      this.onChanged.next(this._prop);
  }

  setColumn(col: DG.Column | null): void {
    if (!col) return;

    if (!this._prop.categories) {
      this._prop.categories = col.categories.map((c) => ({
        name: c,
        desirability: 1,
      }));
    }
    this.buildForm();
  }
}
